use arrow_array::{
    Array, ArrayRef, Int32Array, Int64Array, LargeStringArray, RecordBatch, UInt8Array,
};
use arrow_schema::{DataType, Field, Schema};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use parquet::format::KeyValue;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{self, Read};
use rust_htslib::faidx;
use std::collections::HashMap;
use std::env;
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::Instant;

#[path = "../profiler_baq.rs"]
mod baq;
#[path = "../profiler_stats.rs"]
mod stats;

type AnyResult<T> = Result<T, Box<dyn Error + Send + Sync>>;

#[derive(Clone, Copy, Debug, PartialEq)]
enum Inclusion {
    AllMapped,
    Paired,
    ProperPairs,
}

#[derive(Clone)]
struct Args {
    bam: PathBuf,
    bed: PathBuf,
    stb: PathBuf,
    null_model: PathBuf,
    gene_range: Option<PathBuf>,
    reference: Option<PathBuf>,
    contract: Option<PathBuf>,
    output_dir: PathBuf,
    threads: usize,
    num_chunks: usize,
    min_mapq: u8,
    min_baseq: u8,
    min_read_ani: Option<f64>,
    read_inclusion: Inclusion,
    min_freq: f64,
    max_depth: u32,
}

fn parse_args() -> AnyResult<Args> {
    let mut a = HashMap::new();
    let known = [
        "--bam",
        "--bed",
        "--stb",
        "--null-model",
        "--gene-range",
        "--reference",
        "--profiling-contract",
        "--output-dir",
        "--threads",
        "--num-chunks",
        "--min-mapq",
        "--min-baseq",
        "--min-read-ani",
        "--read-inclusion",
        "--min-freq",
        "--max-depth",
    ];
    let mut it = env::args().skip(1);
    while let Some(key) = it.next() {
        if !key.starts_with("--") {
            return Err(format!("unexpected argument: {key}").into());
        }
        if !known.contains(&key.as_str()) {
            return Err(format!("unknown option: {key}").into());
        }
        let val = it
            .next()
            .ok_or_else(|| format!("missing value for {key}"))?;
        a.insert(key, val);
    }
    let required = |key: &str| -> AnyResult<PathBuf> {
        Ok(PathBuf::from(
            a.get(key).ok_or_else(|| format!("missing {key}"))?,
        ))
    };
    let parse = |key: &str, default: &str| -> AnyResult<String> {
        Ok(a.get(key).cloned().unwrap_or_else(|| default.to_owned()))
    };
    let read_inclusion = match parse("--read-inclusion", "all-mapped")?.as_str() {
        "all-mapped" => Inclusion::AllMapped,
        "paired" => Inclusion::Paired,
        "proper-pairs" => Inclusion::ProperPairs,
        v => return Err(format!("invalid --read-inclusion: {v}").into()),
    };
    let args = Args {
        bam: required("--bam")?,
        bed: required("--bed")?,
        stb: required("--stb")?,
        null_model: required("--null-model")?,
        gene_range: a.get("--gene-range").map(PathBuf::from),
        reference: a.get("--reference").map(PathBuf::from),
        contract: a.get("--profiling-contract").map(PathBuf::from),
        output_dir: required("--output-dir")?,
        threads: parse("--threads", "4")?.parse()?,
        num_chunks: parse("--num-chunks", "24")?.parse()?,
        min_mapq: parse("--min-mapq", "0")?.parse()?,
        min_baseq: parse("--min-baseq", "13")?.parse()?,
        min_read_ani: a.get("--min-read-ani").map(|v| v.parse()).transpose()?,
        read_inclusion,
        min_freq: parse("--min-freq", "0.01")?.parse()?,
        max_depth: parse("--max-depth", "8000")?.parse()?,
    };
    if args.threads == 0
        || args.num_chunks == 0
        || args.max_depth == 0
        || args.max_depth > i32::MAX as u32
        || !(0.0..=1.0).contains(&args.min_freq)
        || args.min_read_ani.is_some_and(|v| !(0.0..=1.0).contains(&v))
    {
        return Err("threads/chunks/depth must be positive; frequencies must be in [0,1]".into());
    }
    Ok(args)
}

#[derive(Clone)]
struct Region {
    chrom: String,
    start: i64,
    end: i64,
}
#[derive(Clone)]
struct Gene {
    name: String,
    start: i32,
    end: i32,
}
struct Inputs {
    tasks: Vec<Vec<Region>>,
    regions: HashMap<String, Vec<Region>>,
    genome: HashMap<String, String>,
    genes: HashMap<String, Vec<Gene>>,
    null: Vec<i64>,
    genome_lengths: HashMap<String, i64>,
}

fn read_bed(path: &Path) -> AnyResult<Vec<Region>> {
    let mut rows = Vec::new();
    for (line_no, line) in BufReader::new(File::open(path)?).lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<_> = line.split_whitespace().collect();
        if cols.len() < 3 {
            return Err(format!("BED line {} has fewer than 3 fields", line_no + 1).into());
        }
        let start: i64 = cols[1].parse()?;
        let end: i64 = cols[2].parse()?;
        if start < 0 || end <= start {
            return Err(format!("invalid BED interval at line {}", line_no + 1).into());
        }
        rows.push(Region {
            chrom: cols[0].to_owned(),
            start,
            end,
        });
    }
    rows.sort_unstable_by(|x, y| (&x.chrom, x.start, x.end).cmp(&(&y.chrom, y.start, y.end)));
    let mut merged: Vec<Region> = Vec::new();
    for r in rows {
        if let Some(last) = merged.last_mut() {
            if last.chrom == r.chrom && r.start <= last.end {
                last.end = last.end.max(r.end);
                continue;
            }
        }
        merged.push(r);
    }
    if merged.is_empty() {
        return Err("BED file has no intervals".into());
    }
    Ok(merged)
}

fn split_tasks(regions: &[Region], count: usize) -> Vec<Vec<Region>> {
    let total: i64 = regions.iter().map(|r| r.end - r.start).sum();
    let target = (total + count as i64 - 1) / count as i64;
    let mut tasks = Vec::new();
    let mut current = Vec::new();
    let mut used = 0;
    for r in regions {
        let mut at = r.start;
        while at < r.end {
            if used == target {
                tasks.push(std::mem::take(&mut current));
                used = 0;
            }
            let end = r.end.min(at + target - used);
            current.push(Region {
                chrom: r.chrom.clone(),
                start: at,
                end,
            });
            used += end - at;
            at = end;
        }
    }
    if !current.is_empty() {
        tasks.push(current);
    }
    tasks
}

fn load_stb(path: &Path) -> AnyResult<HashMap<String, String>> {
    let mut map = HashMap::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let cols: Vec<_> = line.split('\t').map(str::trim).collect();
        if cols.len() >= 2 && !cols[0].is_empty() {
            map.insert(cols[0].to_owned(), cols[1].to_owned());
        }
    }
    Ok(map)
}

fn load_genes(path: Option<&Path>) -> AnyResult<HashMap<String, Vec<Gene>>> {
    let mut map: HashMap<String, Vec<Gene>> = HashMap::new();
    let Some(path) = path else { return Ok(map) };
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let cols: Vec<_> = line.split('\t').map(str::trim).collect();
        if cols.len() < 4 {
            continue;
        }
        let (Ok(start), Ok(end)) = (cols[2].parse::<i32>(), cols[3].parse::<i32>()) else {
            continue;
        };
        map.entry(cols[1].to_owned()).or_default().push(Gene {
            name: cols[0].to_owned(),
            start,
            end,
        });
    }
    for rows in map.values_mut() {
        rows.sort_by_key(|g| g.start);
    }
    Ok(map)
}

fn int_at(array: &dyn Array, i: usize) -> AnyResult<i64> {
    if let Some(a) = array.as_any().downcast_ref::<Int64Array>() {
        return Ok(a.value(i));
    }
    if let Some(a) = array.as_any().downcast_ref::<Int32Array>() {
        return Ok(a.value(i) as i64);
    }
    Err(format!("unsupported integer Arrow type: {:?}", array.data_type()).into())
}

fn load_null(path: &Path) -> AnyResult<Vec<i64>> {
    let file = File::open(path)?;
    let mut reader = ParquetRecordBatchReaderBuilder::try_new(file)?.build()?;
    let mut entries = Vec::new();
    while let Some(batch) = reader.next() {
        let batch = batch?;
        let cov_idx = batch.schema().index_of("cov")?;
        let max_idx = batch.schema().index_of("max_error_count")?;
        for i in 0..batch.num_rows() {
            let cov = usize::try_from(int_at(batch.column(cov_idx).as_ref(), i)?)?;
            let max = int_at(batch.column(max_idx).as_ref(), i)?;
            entries.push((cov, max));
        }
    }
    let limit = entries
        .iter()
        .map(|(cov, _)| *cov)
        .max()
        .ok_or("null model is empty")?;
    let mut null = vec![i64::MIN; limit + 1];
    for (cov, max) in entries {
        null[cov] = max;
    }
    // Samtools may emit a site with zero passing bases; ZipStrain's generated
    // null model starts at coverage 1 and its left join retains that zero row.
    if null[0] == i64::MIN {
        null[0] = 0;
    }
    Ok(null)
}

fn load_inputs(args: &Args) -> AnyResult<Inputs> {
    let bed = read_bed(&args.bed)?;
    let genome = load_stb(&args.stb)?;
    let genes = load_genes(args.gene_range.as_deref())?;
    let null = load_null(&args.null_model)?;
    let tasks = split_tasks(&bed, args.num_chunks);
    let mut regions: HashMap<String, Vec<Region>> = HashMap::new();
    let mut genome_lengths = HashMap::new();
    for r in bed {
        let gen = genome.get(&r.chrom).map(String::as_str).unwrap_or("NA");
        *genome_lengths.entry(gen.to_owned()).or_insert(0) += r.end - r.start;
        regions.entry(r.chrom.clone()).or_default().push(r);
    }
    Ok(Inputs {
        tasks,
        regions,
        genome,
        genes,
        null,
        genome_lengths,
    })
}

fn allowed(record: &bam::Record, args: &Args) -> bool {
    if record.is_unmapped() || record.mapq() < args.min_mapq {
        return false;
    }
    match args.read_inclusion {
        Inclusion::AllMapped => (),
        Inclusion::Paired if record.is_paired() => {
            if record.is_mate_unmapped() || record.mtid() != record.tid() {
                return false;
            }
        }
        Inclusion::ProperPairs => {
            if !record.is_paired() || !record.is_proper_pair() || record.is_mate_unmapped() {
                return false;
            }
        }
        _ => (),
    }
    if let Some(min_ani) = args.min_read_ani {
        let nm = match record.aux(b"NM") {
            Ok(Aux::I8(v)) => Some(v as i64),
            Ok(Aux::U8(v)) => Some(v as i64),
            Ok(Aux::I16(v)) => Some(v as i64),
            Ok(Aux::U16(v)) => Some(v as i64),
            Ok(Aux::I32(v)) => Some(v as i64),
            Ok(Aux::U32(v)) => Some(v as i64),
            _ => None,
        };
        if let Some(nm) = nm {
            let aligned: u64 = record
                .cigar()
                .iter()
                .map(|c| match c {
                    Cigar::Match(n) | Cigar::Ins(n) | Cigar::Equal(n) | Cigar::Diff(n) => *n as u64,
                    _ => 0,
                })
                .sum();
            if aligned == 0 || 1.0 - nm as f64 / (aligned as f64) < min_ani {
                return false;
            }
        }
    }
    true
}

fn ref_mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        _ => 0,
    }
}

fn profile_schema(with_ref: bool) -> Arc<Schema> {
    let mut fields = vec![
        Field::new("chrom", DataType::LargeUtf8, true),
        Field::new("genome", DataType::LargeUtf8, true),
        Field::new("gene", DataType::LargeUtf8, true),
        Field::new("pos", DataType::Int32, true),
        Field::new("A", DataType::Int32, true),
        Field::new("C", DataType::Int32, true),
        Field::new("G", DataType::Int32, true),
        Field::new("T", DataType::Int32, true),
    ];
    if with_ref {
        fields.push(Field::new("ref_base_bitmask", DataType::UInt8, true));
    }
    Arc::new(Schema::new(fields))
}

#[derive(Default)]
struct ProfileBatch {
    chrom: Vec<String>,
    genome: Vec<String>,
    gene: Vec<String>,
    pos: Vec<i32>,
    a: Vec<i32>,
    c: Vec<i32>,
    g: Vec<i32>,
    t: Vec<i32>,
    reference: Vec<u8>,
}
impl ProfileBatch {
    fn len(&self) -> usize {
        self.pos.len()
    }
    fn flush(
        &mut self,
        writer: &mut ArrowWriter<File>,
        schema: &Arc<Schema>,
        with_ref: bool,
    ) -> AnyResult<()> {
        if self.len() == 0 {
            return Ok(());
        }
        let mut arrays: Vec<ArrayRef> = vec![
            Arc::new(LargeStringArray::from_iter_values(self.chrom.drain(..))),
            Arc::new(LargeStringArray::from_iter_values(self.genome.drain(..))),
            Arc::new(LargeStringArray::from_iter_values(self.gene.drain(..))),
            Arc::new(Int32Array::from(std::mem::take(&mut self.pos))),
            Arc::new(Int32Array::from(std::mem::take(&mut self.a))),
            Arc::new(Int32Array::from(std::mem::take(&mut self.c))),
            Arc::new(Int32Array::from(std::mem::take(&mut self.g))),
            Arc::new(Int32Array::from(std::mem::take(&mut self.t))),
        ];
        if with_ref {
            arrays.push(Arc::new(UInt8Array::from(std::mem::take(
                &mut self.reference,
            ))));
        }
        writer.write(&RecordBatch::try_new(schema.clone(), arrays)?)?;
        Ok(())
    }
}

fn writer(
    path: &Path,
    schema: Arc<Schema>,
    metadata: Vec<KeyValue>,
) -> AnyResult<ArrowWriter<File>> {
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .set_max_row_group_size(100_000)
        .set_key_value_metadata(Some(metadata))
        .build();
    Ok(ArrowWriter::try_new(
        File::create(path)?,
        schema,
        Some(props),
    )?)
}

struct ReferenceWindow {
    chrom: String,
    start: usize,
    bases: Vec<u8>,
}

fn reference_base(
    reader: &faidx::Reader,
    window: &mut Option<ReferenceWindow>,
    chrom: &str,
    pos: usize,
) -> AnyResult<u8> {
    if !window
        .as_ref()
        .is_some_and(|w| w.chrom == chrom && pos >= w.start && pos < w.start + w.bases.len())
    {
        let start = (pos / 65_536) * 65_536;
        let end = start + 65_535;
        let bases = reader.fetch_seq(chrom, start, end)?.to_vec();
        *window = Some(ReferenceWindow {
            chrom: chrom.to_owned(),
            start,
            bases,
        });
    }
    let w = window.as_ref().ok_or("reference cache missing")?;
    Ok(ref_mask(
        *w.bases
            .get(pos - w.start)
            .ok_or("reference position out of bounds")?,
    ))
}

fn run_task(
    index: usize,
    args: &Args,
    inputs: &Inputs,
    schema: &Arc<Schema>,
    tmp: &Path,
) -> AnyResult<PathBuf> {
    let path = tmp.join(format!("chunk_{index:05}.parquet"));
    let mut output = writer(&path, schema.clone(), vec![])?;
    let mut batch = ProfileBatch::default();
    let mut bam = baq::BaqReader::new(args)?;
    let reference = args
        .reference
        .as_deref()
        .map(faidx::Reader::from_path)
        .transpose()?;
    let mut ref_window = None;
    let mut last_key: Option<(String, i32)> = None;
    let header = bam.header().to_owned();
    for region in &inputs.tasks[index] {
        bam.fetch(&region.chrom, region.start, region.end)?;
        let mut pileups = bam.pileup_with_overlaps(args.max_depth)?;
        while let Some(pileup) = pileups.next_position()? {
            let p0 = i64::from(pileup.pos());
            if p0 < region.start || p0 >= region.end {
                continue;
            }
            let chrom = std::str::from_utf8(header.tid2name(pileup.tid()))?;
            if chrom != region.chrom {
                return Err("indexed BAM returned an unexpected scaffold".into());
            }
            let pos = i32::try_from(p0 + 1)?;
            let key = (chrom.to_owned(), pos);
            if last_key.as_ref().is_some_and(|prev| prev >= &key) {
                return Err(format!("non-increasing chunk coordinates at {chrom}:{pos}").into());
            }
            last_key = Some(key);
            let mut counts = [0_i32; 4];
            for alignment in pileup.alignments() {
                let Some(qpos) = alignment.qpos() else {
                    continue;
                };
                let record = alignment.record();
                if record.flags() & 0x704 != 0 {
                    continue;
                }
                if record
                    .qual()
                    .get(qpos)
                    .is_none_or(|quality| *quality < args.min_baseq)
                {
                    continue;
                }
                let index = match record.seq()[qpos].to_ascii_uppercase() {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => continue,
                };
                counts[index] += 1;
            }
            let depth: usize = counts.iter().map(|v| *v as usize).sum();
            let null_limit = *inputs.null.get(depth).ok_or_else(|| {
                format!(
                    "observed depth {depth} exceeds null-model limit {} at {chrom}:{pos}",
                    inputs.null.len() - 1
                )
            })?;
            if null_limit == i64::MIN {
                return Err(format!("null model lacks coverage {depth}").into());
            }
            for count in &mut counts {
                if i64::from(*count) <= null_limit
                    || f64::from(*count) < depth as f64 * args.min_freq
                {
                    *count = 0;
                }
            }
            let gene = inputs
                .genes
                .get(chrom)
                .and_then(|genes| {
                    let at = genes.partition_point(|g| g.start <= pos);
                    (at > 0 && pos <= genes[at - 1].end).then(|| genes[at - 1].name.as_str())
                })
                .unwrap_or("NA");
            batch.chrom.push(chrom.to_owned());
            batch.genome.push(
                inputs
                    .genome
                    .get(chrom)
                    .map(String::as_str)
                    .unwrap_or("NA")
                    .to_owned(),
            );
            batch.gene.push(gene.to_owned());
            batch.pos.push(pos);
            batch.a.push(counts[0]);
            batch.c.push(counts[1]);
            batch.g.push(counts[2]);
            batch.t.push(counts[3]);
            if let Some(reader) = &reference {
                batch
                    .reference
                    .push(reference_base(reader, &mut ref_window, chrom, p0 as usize)?);
            }
            if batch.len() >= 10_000 {
                batch.flush(&mut output, schema, reference.is_some())?;
            }
        }
    }
    batch.flush(&mut output, schema, reference.is_some())?;
    output.close()?;
    Ok(path)
}

fn run_workers(
    args: &Args,
    inputs: &Inputs,
    schema: &Arc<Schema>,
    tmp: &Path,
) -> AnyResult<Vec<PathBuf>> {
    let next = AtomicUsize::new(0);
    let mut paths: Vec<Option<PathBuf>> = vec![None; inputs.tasks.len()];
    thread::scope(|scope| -> AnyResult<()> {
        let mut jobs = Vec::new();
        for _ in 0..args.threads.min(inputs.tasks.len()) {
            jobs.push(scope.spawn(|| -> AnyResult<Vec<(usize, PathBuf)>> {
                let mut completed = Vec::new();
                loop {
                    let id = next.fetch_add(1, Ordering::Relaxed);
                    if id >= inputs.tasks.len() {
                        break;
                    }
                    completed.push((id, run_task(id, args, inputs, schema, tmp)?));
                }
                Ok(completed)
            }));
        }
        for job in jobs {
            for (id, path) in job.join().map_err(|_| "profile worker panicked")?? {
                paths[id] = Some(path);
            }
        }
        Ok(())
    })?;
    paths
        .into_iter()
        .map(|p| p.ok_or_else(|| "missing chunk output".into()))
        .collect()
}

fn metadata(args: &Args) -> AnyResult<Vec<KeyValue>> {
    let mut pairs = vec![KeyValue::new(
        "zipstrain_sorted_by".to_owned(),
        Some("chrom,pos".to_owned()),
    )];
    if let Some(contract) = &args.contract {
        let value: serde_json::Value = serde_json::from_reader(File::open(contract)?)?;
        for (json_key, parquet_key) in [
            ("reference_hash", "zipstrain_reference_hash"),
            ("gene_hash", "zipstrain_gene_hash"),
            ("null_model_hash", "zipstrain_null_model_hash"),
            ("stb_hash", "zipstrain_stb_hash"),
        ] {
            let hash = value.get(json_key).and_then(|v| v.as_str()).unwrap_or("NA");
            pairs.push(KeyValue::new(parquet_key.to_owned(), Some(hash.to_owned())));
        }
    } else {
        for key in [
            "zipstrain_reference_hash",
            "zipstrain_gene_hash",
            "zipstrain_null_model_hash",
            "zipstrain_stb_hash",
        ] {
            pairs.push(KeyValue::new(key.to_owned(), Some("NA".to_owned())));
        }
    }
    Ok(pairs)
}

fn merge_chunks(
    paths: &[PathBuf],
    args: &Args,
    schema: &Arc<Schema>,
    output: &Path,
) -> AnyResult<stats::Accumulator> {
    let temp = output.with_extension("parquet.tmp");
    let mut writer = writer(&temp, schema.clone(), metadata(args)?)?;
    let mut stats = stats::Accumulator::default();
    let mut last: Option<(String, i32)> = None;
    for path in paths {
        let mut reader = ParquetRecordBatchReaderBuilder::try_new(File::open(path)?)?.build()?;
        while let Some(batch) = reader.next() {
            let batch = batch?;
            let chroms = batch
                .column(0)
                .as_any()
                .downcast_ref::<LargeStringArray>()
                .ok_or("chrom type")?;
            let pos = batch
                .column(3)
                .as_any()
                .downcast_ref::<Int32Array>()
                .ok_or("pos type")?;
            for i in 0..batch.num_rows() {
                let key = (chroms.value(i).to_owned(), pos.value(i));
                if last.as_ref().is_some_and(|prev| prev >= &key) {
                    return Err(format!(
                        "non-increasing final profile coordinates at {}:{}",
                        key.0, key.1
                    )
                    .into());
                }
                last = Some(key);
            }
            stats.add_batch(&batch, args.reference.is_some())?;
            writer.write(&batch)?;
        }
    }
    writer.close()?;
    fs::rename(temp, output)?;
    Ok(stats)
}

fn main() -> AnyResult<()> {
    let start = Instant::now();
    let args = parse_args()?;
    if let Some(reference) = &args.reference {
        let mut fai = reference.as_os_str().to_os_string();
        fai.push(".fai");
        if !Path::new(&fai).exists() {
            faidx::build(reference).map_err(|e| format!("FASTA indexing failed: {e}"))?;
        }
    }
    fs::create_dir_all(&args.output_dir)?;
    let inputs = load_inputs(&args)?;
    let schema = profile_schema(args.reference.is_some());
    let stem = args
        .bam
        .file_stem()
        .ok_or("BAM has no filename")?
        .to_string_lossy();
    let out = args.output_dir.join(format!("{stem}_profile.parquet"));
    let tmp = args.output_dir.join(format!(".{stem}_rust_profile_work"));
    fs::create_dir_all(&tmp)?;
    eprintln!(
        "Profiling {} chunks with {} Rust workers",
        inputs.tasks.len(),
        args.threads
    );
    let workers_start = Instant::now();
    let paths = run_workers(&args, &inputs, &schema, &tmp)?;
    eprintln!("Workers: {:.2}s", workers_start.elapsed().as_secs_f64());
    let merge_start = Instant::now();
    let mut acc = merge_chunks(&paths, &args, &schema, &out)?;
    eprintln!("Merge: {:.2}s", merge_start.elapsed().as_secs_f64());
    let stats_start = Instant::now();
    acc.write_stats(&args, &inputs, &stem)?;
    eprintln!("Stats: {:.2}s", stats_start.elapsed().as_secs_f64());
    fs::remove_dir_all(&tmp)?;
    eprintln!(
        "Completed {} in {:.2}s",
        out.display(),
        start.elapsed().as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tasks_preserve_disjoint_lexicographic_intervals() {
        let regions = vec![
            Region {
                chrom: "chr1".into(),
                start: 0,
                end: 7,
            },
            Region {
                chrom: "chr10".into(),
                start: 2,
                end: 9,
            },
            Region {
                chrom: "chr2".into(),
                start: 1,
                end: 4,
            },
        ];
        let tasks = split_tasks(&regions, 5);
        let flattened: Vec<_> = tasks.iter().flat_map(|task| task.iter()).collect();
        assert!(tasks.len() <= 5);
        assert_eq!(flattened.iter().map(|r| r.end - r.start).sum::<i64>(), 17);
        for pair in flattened.windows(2) {
            assert!(
                pair[0].chrom < pair[1].chrom
                    || (pair[0].chrom == pair[1].chrom && pair[0].end <= pair[1].start)
            );
        }
    }

    #[test]
    fn reference_mask_codes_are_stable() {
        assert_eq!(
            [b'A', b'C', b'G', b'T', b'N', b'a'].map(ref_mask),
            [1, 2, 4, 8, 0, 1]
        );
    }
}
