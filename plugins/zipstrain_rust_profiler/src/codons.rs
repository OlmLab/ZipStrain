//! Read-phased codon counts from a coordinate-sorted BAM.
//!
//! This implements the extraction half of `zipstrain.codons.write_codon_counts`:
//! counts are per fragment, a mate overlap keeps the higher-quality triplet, equal-quality
//! disagreement discards the observation, and the mate cache expires as coordinates advance.
//! ZipStrain applies the retained-allele gate to the raw table written here.

use arrow_array::{ArrayRef, RecordBatch, StringArray, UInt32Array, UInt64Array};
use arrow_schema::{DataType, Field, Schema};
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{self, Read};
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::rc::Rc;
use std::sync::Arc;

type AnyResult<T> = Result<T, Box<dyn Error + Send + Sync>>;
type Codon = [u8; 3];

const BUCKET: i64 = 16384;
// Unmapped, secondary, QC-failed, duplicate, supplementary.
const EXCLUDED_FLAGS: u16 = 0xF04;
const WRITE_BATCH_ROWS: usize = 1 << 20;

struct Gene {
    id: u64,
    start: i64,
    end: i64,
    length: i64,
    strand: i64,
    phase: i64,
}

#[derive(Clone, Copy)]
struct Call {
    base: u8,
    quality: u8,
    query: i64,
    present: bool,
}

const ABSENT: Call = Call {
    base: 0,
    quality: 0,
    query: 0,
    present: false,
};

fn floor_div(a: i64, b: i64) -> i64 {
    let quotient = a / b;
    if a % b != 0 && ((a < 0) != (b < 0)) {
        quotient - 1
    } else {
        quotient
    }
}

/// Map an alignment onto reference coordinates. Returns the span's first position, or None
/// when no base aligns. `calls` is reused between reads.
fn aligned_calls(record: &bam::Record, start: i64, calls: &mut Vec<Call>) -> AnyResult<Option<(i64, i64)>> {
    calls.clear();
    let cigar = record.cigar();
    if cigar.is_empty() {
        return Err("Invalid CIGAR: *".into());
    }
    let sequence = record.seq();
    let qualities = record.qual();
    let length = sequence.len() as i64;
    let (mut query, mut reference) = (0i64, start);
    let (mut first, mut last) = (None, 0i64);
    for op in cigar.iter() {
        let n = op.len() as i64;
        if n < 1 {
            return Err("CIGAR operations must be positive".into());
        }
        match op {
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                if query + n > length {
                    return Err("CIGAR exceeds sequence length".into());
                }
                let offset = (reference - start) as usize;
                if calls.len() < offset {
                    calls.resize(offset, ABSENT);
                }
                for i in 0..n {
                    let index = (query + i) as usize;
                    calls.push(Call {
                        base: sequence[index],
                        quality: qualities[index],
                        query: query + i,
                        present: true,
                    });
                }
                first.get_or_insert(reference);
                last = reference + n - 1;
                reference += n;
                query += n;
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) => query += n,
            Cigar::Del(_) | Cigar::RefSkip(_) => reference += n,
            Cigar::HardClip(_) | Cigar::Pad(_) => (),
        }
    }
    if query != length {
        return Err("CIGAR query length does not match sequence".into());
    }
    Ok(first.map(|first| (first, last)))
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'A',
    }
}

/// Call `emit(index, codon, quality)` for each complete, in-frame triplet this read observes.
fn observe(
    calls: &[Call],
    start: i64,
    (lo, hi): (i64, i64),
    gene: &Gene,
    min_baseq: u8,
    mut emit: impl FnMut(u32, Codon, u8),
) {
    let at = |position: i64| -> Option<Call> {
        let offset = position - start;
        if offset < 0 {
            return None;
        }
        calls.get(offset as usize).copied().filter(|call| call.present)
    };
    let minus = gene.strand == -1;
    let origin = if minus { gene.end - gene.phase } else { gene.start + gene.phase };
    let first = if minus {
        floor_div(hi - origin, -3).max(0)
    } else {
        floor_div(lo - origin, 3).max(0)
    };
    for index in first..floor_div(gene.length - gene.phase, 3) {
        let position = |offset: i64| origin + gene.strand * (index * 3 + offset);
        if (!minus && position(0) > hi) || (minus && position(0) < lo) {
            break;
        }
        let (Some(a), Some(b), Some(c)) = (at(position(0)), at(position(1)), at(position(2))) else {
            continue;
        };
        if [a, b, c]
            .iter()
            .any(|call| !matches!(call.base, b'A' | b'C' | b'G' | b'T') || call.quality < min_baseq)
        {
            continue;
        }
        if b.query - a.query != gene.strand || c.query - b.query != gene.strand {
            continue;
        }
        let codon = if minus {
            [complement(a.base), complement(b.base), complement(c.base)]
        } else {
            [a.base, b.base, c.base]
        };
        emit(index as u32, codon, a.quality.min(b.quality).min(c.quality));
    }
}

fn load_genes(path: &PathBuf, names: &HashMap<Vec<u8>, usize>) -> AnyResult<(Vec<Gene>, Vec<HashMap<i64, Vec<u32>>>)> {
    let (mut genes, mut bins) = (Vec::new(), vec![HashMap::new(); names.len()]);
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        let [id, chrom, start, end, length, strand, phase] = fields[..] else {
            return Err(format!("Malformed gene row: {line}").into());
        };
        let Some(&tid) = names.get(chrom.as_bytes()) else {
            continue;
        };
        let gene = Gene {
            id: id.parse()?,
            start: start.parse()?,
            end: end.parse()?,
            length: length.parse()?,
            strand: strand.parse()?,
            phase: phase.parse()?,
        };
        for bucket in gene.start / BUCKET..=gene.end / BUCKET {
            let row: &mut HashMap<i64, Vec<u32>> = &mut bins[tid];
            row.entry(bucket).or_default().push(genes.len() as u32);
        }
        genes.push(gene);
    }
    Ok((genes, bins))
}

fn flush(writer: &mut ArrowWriter<File>, schema: &Arc<Schema>, counts: &mut HashMap<(u64, u32, Codon), i64>) -> AnyResult<()> {
    let (mut ids, mut indices, mut codons, mut values) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for ((id, index, codon), count) in counts.drain() {
        if count > 0 {
            ids.push(id);
            indices.push(index);
            codons.push(String::from_utf8_lossy(&codon).into_owned());
            values.push(count as u64);
        }
        if ids.len() == WRITE_BATCH_ROWS {
            write_rows(writer, schema, &mut ids, &mut indices, &mut codons, &mut values)?;
        }
    }
    write_rows(writer, schema, &mut ids, &mut indices, &mut codons, &mut values)
}

fn write_rows(
    writer: &mut ArrowWriter<File>,
    schema: &Arc<Schema>,
    ids: &mut Vec<u64>,
    indices: &mut Vec<u32>,
    codons: &mut Vec<String>,
    values: &mut Vec<u64>,
) -> AnyResult<()> {
    if ids.is_empty() {
        return Ok(());
    }
    let columns: Vec<ArrayRef> = vec![
        Arc::new(UInt64Array::from(std::mem::take(ids))),
        Arc::new(UInt32Array::from(std::mem::take(indices))),
        Arc::new(StringArray::from(std::mem::take(codons))),
        Arc::new(UInt64Array::from(std::mem::take(values))),
    ];
    writer.write(&RecordBatch::try_new(schema.clone(), columns)?)?;
    Ok(())
}

/// Count fragments' codons across the BAM and write the raw table. `allowed` applies the
/// profiler's read filters.
pub fn extract(
    bam_path: &PathBuf,
    genes_path: &PathBuf,
    output: &PathBuf,
    threads: usize,
    min_baseq: u8,
    allowed: impl Fn(&bam::Record) -> bool,
) -> AnyResult<()> {
    let mut reader = bam::Reader::from_path(bam_path)?;
    if threads > 1 {
        reader.set_threads(threads)?;
    }
    let names: HashMap<Vec<u8>, usize> = reader
        .header()
        .target_names()
        .iter()
        .enumerate()
        .map(|(tid, name)| (name.to_vec(), tid))
        .collect();
    let (genes, bins) = load_genes(genes_path, &names)?;
    let schema = Arc::new(Schema::new(vec![
        Field::new("gene_id", DataType::UInt64, false),
        Field::new("codon_index", DataType::UInt32, false),
        Field::new("codon", DataType::Utf8, false),
        Field::new("count", DataType::UInt64, false),
    ]));
    let properties = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build();
    let mut writer = ArrowWriter::try_new(File::create(output)?, schema.clone(), Some(properties))?;

    type Key = (Rc<[u8]>, u64, u32);
    let mut counts: HashMap<(u64, u32, Codon), i64> = HashMap::new();
    let mut active: HashMap<Key, (Option<Codon>, u8)> = HashMap::new();
    let mut expiry: BinaryHeap<Reverse<(i64, u64)>> = BinaryHeap::new();
    let mut expiring: HashMap<u64, Key> = HashMap::new();
    let mut sequence_number = 0u64;
    let (mut chrom, mut previous) = (None::<i32>, -1i64);
    let mut finished = HashSet::new();
    let (mut record, mut calls, mut candidates) = (bam::Record::new(), Vec::new(), Vec::new());
    while let Some(result) = reader.read(&mut record) {
        result?;
        if record.flags() & EXCLUDED_FLAGS != 0 || !allowed(&record) {
            continue;
        }
        let tid = record.tid();
        let start = record.pos() + 1;
        if chrom != Some(tid) {
            if finished.contains(&tid) {
                return Err("Codon profiling requires coordinate-sorted BAM".into());
            }
            if let Some(done) = chrom {
                finished.insert(done);
            }
            flush(&mut writer, &schema, &mut counts)?;
            active.clear();
            expiry.clear();
            expiring.clear();
            chrom = Some(tid);
            previous = -1;
        }
        if start < previous {
            return Err("Codon profiling requires coordinate-sorted BAM".into());
        }
        previous = start;
        while let Some(&Reverse((end, number))) = expiry.peek() {
            if end >= start {
                break;
            }
            expiry.pop();
            if let Some(key) = expiring.remove(&number) {
                active.remove(&key);
            }
        }
        if record.seq_len() == 0 || record.qual().first() == Some(&0xFF) {
            continue;
        }
        if start < 1 {
            return Err("SAM reference positions must be one-based".into());
        }
        let Some((lo, hi)) = aligned_calls(&record, start, &mut calls)? else {
            continue;
        };
        let scaffold = &bins[tid as usize];
        candidates.clear();
        for bucket in start / BUCKET..=hi / BUCKET {
            if let Some(found) = scaffold.get(&bucket) {
                candidates.extend_from_slice(found);
            }
        }
        if candidates.is_empty() {
            continue;
        }
        candidates.sort_unstable();
        candidates.dedup();
        let mut fragment = record.qname().to_vec();
        fragment.push(0);
        if let Ok(Aux::String(group)) = record.aux(b"RG") {
            fragment.extend_from_slice(group.as_bytes());
        }
        let fragment: Rc<[u8]> = fragment.into();
        for &gene_index in &candidates {
            let gene = &genes[gene_index as usize];
            observe(&calls, start, (lo, hi), gene, min_baseq, |index, codon, quality| {
                let key = (fragment.clone(), gene.id, index);
                match active.get(&key).copied() {
                    Some((old_codon, old_quality)) => {
                        if quality < old_quality || (quality == old_quality && old_codon == Some(codon)) {
                            return;
                        }
                        if let Some(old) = old_codon {
                            *counts.entry((gene.id, index, old)).or_insert(0) -= 1;
                        }
                        if quality == old_quality {
                            active.insert(key, (None, quality));
                            return;
                        }
                    }
                    None => {
                        expiry.push(Reverse((hi, sequence_number)));
                        expiring.insert(sequence_number, key.clone());
                        sequence_number += 1;
                    }
                }
                *counts.entry((gene.id, index, codon)).or_insert(0) += 1;
                active.insert(key, (Some(codon), quality));
            });
        }
    }
    flush(&mut writer, &schema, &mut counts)?;
    writer.close()?;
    Ok(())
}
