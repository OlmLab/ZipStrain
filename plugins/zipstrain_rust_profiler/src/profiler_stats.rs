use super::{allowed, writer, AnyResult, Args, Inputs};
use arrow_array::{
    ArrayRef, Float64Array, Int32Array, Int64Array, LargeStringArray, RecordBatch, UInt8Array,
};
use arrow_schema::{DataType, Field, Schema};
use rust_htslib::bam::{self, Read};
use std::collections::{BTreeMap, HashMap};
use std::sync::Arc;

#[derive(Default)]
struct Counts {
    sites: i64,
    bases: i64,
    depth_histogram: BTreeMap<i32, i64>,
    deep_sites: i64,
    hetero_sites: i64,
    hetero_denominator: i64,
    ref_sites: i64,
    ref_shared: i64,
    ref_consensus_shared: i64,
    sns: i64,
    snv: i64,
}

impl Counts {
    fn add(&mut self, bases: [i32; 4], reference: Option<u8>, genome: bool) {
        let depth: i32 = bases.iter().sum();
        self.sites += 1;
        self.bases += i64::from(depth);
        if genome {
            *self.depth_histogram.entry(depth).or_default() += 1;
            if depth >= 5 {
                self.deep_sites += 1;
            }
            if depth > 5 {
                self.hetero_denominator += 1;
                if (*bases.iter().max().unwrap_or(&0) as f64 / depth as f64) <= 0.8 {
                    self.hetero_sites += 1;
                }
            }
        }
        if let Some(mask) = reference {
            if mask > 0 && depth >= 5 {
                self.ref_sites += 1;
                let idx = match mask {
                    1 => 0,
                    2 => 1,
                    4 => 2,
                    8 => 3,
                    _ => return,
                };
                if bases[idx] > 0 {
                    self.ref_shared += 1;
                }
                if bases[idx] == *bases.iter().max().unwrap_or(&0) {
                    self.ref_consensus_shared += 1;
                }
                if genome {
                    let n = bases.iter().filter(|v| **v > 0).count();
                    if n == 1 && bases[idx] == 0 {
                        self.sns += 1;
                    }
                    if n >= 2 {
                        self.snv += 1;
                    }
                }
            }
        }
    }

    fn ref_ani(&self) -> Option<f64> {
        (self.ref_sites > 0).then(|| self.ref_shared as f64 / self.ref_sites as f64 * 100.0)
    }

    fn median(&self) -> Option<f64> {
        if self.sites == 0 {
            return None;
        }
        let value = |target: i64| -> f64 {
            let mut n = 0;
            for (depth, count) in &self.depth_histogram {
                n += count;
                if n > target {
                    return *depth as f64;
                }
            }
            0.0
        };
        Some((value((self.sites - 1) / 2) + value(self.sites / 2)) / 2.0)
    }

    fn std(&self) -> Option<f64> {
        if self.sites < 2 {
            return None;
        }
        let mean = self.bases as f64 / self.sites as f64;
        let sum: f64 = self
            .depth_histogram
            .iter()
            .map(|(v, n)| (*v as f64 - mean).powi(2) * *n as f64)
            .sum();
        Some((sum / (self.sites - 1) as f64).sqrt())
    }
}

#[derive(Default)]
pub(super) struct Accumulator {
    genomes: HashMap<String, Counts>,
    genes: HashMap<(String, String), Counts>,
}

impl Accumulator {
    pub(super) fn add_batch(&mut self, batch: &RecordBatch, with_ref: bool) -> AnyResult<()> {
        let strings = |n: usize| -> AnyResult<&LargeStringArray> {
            batch
                .column(n)
                .as_any()
                .downcast_ref::<LargeStringArray>()
                .ok_or_else(|| "wrong string column".into())
        };
        let ints = |n: usize| -> AnyResult<&Int32Array> {
            batch
                .column(n)
                .as_any()
                .downcast_ref::<Int32Array>()
                .ok_or_else(|| "wrong integer column".into())
        };
        let genomes = strings(1)?;
        let genes = strings(2)?;
        let a = ints(4)?;
        let c = ints(5)?;
        let g = ints(6)?;
        let t = ints(7)?;
        let reference = if with_ref {
            Some(
                batch
                    .column(8)
                    .as_any()
                    .downcast_ref::<UInt8Array>()
                    .ok_or("wrong reference column")?,
            )
        } else {
            None
        };
        for i in 0..batch.num_rows() {
            let bases = [a.value(i), c.value(i), g.value(i), t.value(i)];
            let mask = reference.map(|r| r.value(i));
            self.genomes
                .entry(genomes.value(i).to_owned())
                .or_default()
                .add(bases, mask, true);
            if genes.value(i) != "NA" {
                self.genes
                    .entry((genomes.value(i).to_owned(), genes.value(i).to_owned()))
                    .or_default()
                    .add(bases, mask, false);
            }
        }
        Ok(())
    }

    pub(super) fn write_stats(
        &mut self,
        args: &Args,
        inputs: &Inputs,
        stem: &str,
    ) -> AnyResult<()> {
        self.write_gene_stats(args, inputs, stem)?;
        self.write_genome_stats(args, inputs, stem)?;
        Ok(())
    }

    fn write_gene_stats(&self, args: &Args, inputs: &Inputs, stem: &str) -> AnyResult<()> {
        let with_ref = args.reference.is_some();
        let mut fields = vec![
            Field::new("genome", DataType::LargeUtf8, true),
            Field::new("gene", DataType::LargeUtf8, true),
            Field::new("length", DataType::Int64, true),
            Field::new("breadth", DataType::Float64, true),
            Field::new("coverage", DataType::Float64, true),
        ];
        if with_ref {
            fields.push(Field::new("ref_ani", DataType::Float64, true));
        }
        let schema = Arc::new(Schema::new(fields));
        let mut seen = BTreeMap::new();
        for (chrom, genes) in &inputs.genes {
            let genome = inputs.genome.get(chrom).map(String::as_str).unwrap_or("NA");
            for gene in genes {
                seen.insert(
                    (genome.to_owned(), gene.name.clone()),
                    i64::from(gene.end - gene.start + 1),
                );
            }
        }
        let mut genome = Vec::new();
        let mut gene = Vec::new();
        let mut length = Vec::new();
        let mut breadth = Vec::new();
        let mut coverage = Vec::new();
        let mut ref_ani = Vec::new();
        for ((g, n), len) in seen {
            let counts = self.genes.get(&(g.clone(), n.clone()));
            genome.push(g);
            gene.push(n);
            length.push(len);
            breadth.push(counts.map_or(0.0, |c| c.sites as f64 / len as f64));
            coverage.push(counts.map_or(0.0, |c| c.bases as f64 / len as f64));
            ref_ani.push(counts.and_then(Counts::ref_ani));
        }
        let mut arrays: Vec<ArrayRef> = vec![
            Arc::new(LargeStringArray::from(genome)),
            Arc::new(LargeStringArray::from(gene)),
            Arc::new(Int64Array::from(length)),
            Arc::new(Float64Array::from(breadth)),
            Arc::new(Float64Array::from(coverage)),
        ];
        if with_ref {
            arrays.push(Arc::new(Float64Array::from(ref_ani)));
        }
        let output = args.output_dir.join(format!("{stem}_gene_stats.parquet"));
        let mut w = writer(&output, schema.clone(), vec![])?;
        w.write(&RecordBatch::try_new(schema, arrays)?)?;
        w.close()?;
        Ok(())
    }

    fn write_genome_stats(&self, args: &Args, inputs: &Inputs, stem: &str) -> AnyResult<()> {
        let with_ref = args.reference.is_some();
        let mut fields = vec![
            Field::new("genome", DataType::LargeUtf8, true),
            Field::new("coverage", DataType::Float64, true),
            Field::new("coverage_median", DataType::Float64, true),
            Field::new("coverage_std", DataType::Float64, true),
            Field::new("breadth", DataType::Float64, true),
            Field::new("genome_length", DataType::Int64, true),
            Field::new("gap_mean", DataType::Float64, true),
            Field::new("gap_std", DataType::Float64, true),
            Field::new("5x_cov_sites", DataType::Int64, true),
            Field::new("heterogeneity", DataType::Float64, true),
            Field::new("ber", DataType::Float64, true),
            Field::new("fug", DataType::Float64, true),
            Field::new("reads_mapped", DataType::Int64, true),
        ];
        if with_ref {
            fields.extend([
                Field::new("ref_ani", DataType::Float64, true),
                Field::new("conANI_reference", DataType::Float64, true),
                Field::new("SNS_count", DataType::Int64, true),
                Field::new("SNV_count", DataType::Int64, true),
            ]);
        }
        let schema = Arc::new(Schema::new(fields));
        let gaps = scan_read_gaps(args, inputs)?;
        let mut keys: Vec<_> = self.genomes.keys().cloned().collect();
        keys.sort();
        let mut names = Vec::new();
        let mut mean_cov = Vec::new();
        let mut med = Vec::new();
        let mut std = Vec::new();
        let mut breadth = Vec::new();
        let mut lengths = Vec::new();
        let mut gap_mean = Vec::new();
        let mut gap_std = Vec::new();
        let mut deep = Vec::new();
        let mut het = Vec::new();
        let mut ber = Vec::new();
        let mut fug = Vec::new();
        let mut reads = Vec::new();
        let mut ref_ani = Vec::new();
        let mut conani = Vec::new();
        let mut sns = Vec::new();
        let mut snv = Vec::new();
        for name in keys {
            let c = &self.genomes[&name];
            let len = *inputs
                .genome_lengths
                .get(&name)
                .ok_or("genome length missing")?;
            let cov = c.bases as f64 / len as f64;
            let br = c.sites as f64 / len as f64;
            let gap = gaps.get(&name);
            names.push(name);
            mean_cov.push(cov);
            med.push(c.median());
            std.push(c.std());
            breadth.push(br);
            lengths.push(len);
            gap_mean.push(gap.and_then(|v| v.mean));
            gap_std.push(gap.and_then(|v| v.std));
            deep.push(c.deep_sites);
            het.push(
                (c.hetero_denominator > 0)
                    .then(|| c.hetero_sites as f64 / c.hetero_denominator as f64),
            );
            ber.push(br / (1.0 - (-0.883 * cov).exp()));
            fug.push(gap.and_then(|v| v.fug));
            reads.push(gap.map_or(0, |v| v.rn));
            ref_ani.push(c.ref_ani());
            conani.push(
                (c.ref_sites > 0)
                    .then(|| c.ref_consensus_shared as f64 / c.ref_sites as f64 * 100.0),
            );
            sns.push(c.sns);
            snv.push(c.snv);
        }
        let mut arrays: Vec<ArrayRef> = vec![
            Arc::new(LargeStringArray::from(names)),
            Arc::new(Float64Array::from(mean_cov)),
            Arc::new(Float64Array::from(med)),
            Arc::new(Float64Array::from(std)),
            Arc::new(Float64Array::from(breadth)),
            Arc::new(Int64Array::from(lengths)),
            Arc::new(Float64Array::from(gap_mean)),
            Arc::new(Float64Array::from(gap_std)),
            Arc::new(Int64Array::from(deep)),
            Arc::new(Float64Array::from(het)),
            Arc::new(Float64Array::from(ber)),
            Arc::new(Float64Array::from(fug)),
            Arc::new(Int64Array::from(reads)),
        ];
        if with_ref {
            arrays.extend([
                Arc::new(Float64Array::from(ref_ani)) as ArrayRef,
                Arc::new(Float64Array::from(conani)) as ArrayRef,
                Arc::new(Int64Array::from(sns)) as ArrayRef,
                Arc::new(Int64Array::from(snv)) as ArrayRef,
            ]);
        }
        let output = args.output_dir.join(format!("{stem}_genome_stats.parquet"));
        let mut w = writer(&output, schema.clone(), vec![])?;
        w.write(&RecordBatch::try_new(schema, arrays)?)?;
        w.close()?;
        Ok(())
    }
}

struct GapStats {
    rn: i64,
    mean: Option<f64>,
    std: Option<f64>,
    fug: Option<f64>,
}
#[derive(Default)]
struct GapScan {
    rn: i64,
    n_gaps: i64,
    sum: f64,
    sum_sq: f64,
    excess: f64,
}

fn scan_read_gaps(args: &Args, inputs: &Inputs) -> AnyResult<HashMap<String, GapStats>> {
    // Two scans avoid retaining every read start or gap in memory.
    let mut scans: HashMap<String, GapScan> = HashMap::new();
    for pass in 0..2 {
        let mut reader = bam::Reader::from_path(&args.bam)?;
        let header = reader.header().to_owned();
        let mut previous: HashMap<String, i64> = HashMap::new();
        for rec in reader.records() {
            let rec = rec?;
            if rec.flags() & 0x84 != 0 || !allowed(&rec, args) {
                continue;
            }
            let chrom = std::str::from_utf8(header.tid2name(rec.tid() as u32))?;
            let Some(regions) = inputs.regions.get(chrom) else {
                continue;
            };
            let left = rec.pos();
            let right = rec.cigar().end_pos();
            let index = regions.partition_point(|r| r.end <= left);
            if index == regions.len() || regions[index].start >= right {
                continue;
            }
            let genome = inputs.genome.get(chrom).map(String::as_str).unwrap_or("NA");
            let gap = previous
                .insert(chrom.to_owned(), left + 1)
                .map(|p| left + 1 - p);
            let scan = scans.entry(genome.to_owned()).or_default();
            if pass == 0 {
                scan.rn += 1;
                if let Some(g) = gap {
                    scan.n_gaps += 1;
                    scan.sum += g as f64;
                    scan.sum_sq += (g as f64).powi(2);
                }
            } else if let Some(g) = gap {
                let len = *inputs.genome_lengths.get(genome).unwrap_or(&0) as f64;
                let delta = (len / scan.rn as f64).round();
                if (g as f64) > delta {
                    scan.excess += g as f64 - delta;
                }
            }
        }
    }
    let mut out = HashMap::new();
    for (genome, scan) in scans {
        let length = *inputs.genome_lengths.get(&genome).unwrap_or(&0) as f64;
        let delta = (length / scan.rn as f64).round();
        let mean = (scan.n_gaps > 0).then(|| scan.sum / scan.n_gaps as f64);
        let std = if scan.n_gaps > 1 {
            let variance =
                (scan.sum_sq - scan.sum.powi(2) / scan.n_gaps as f64) / (scan.n_gaps - 1) as f64;
            Some(variance.max(0.0).sqrt())
        } else {
            None
        };
        let fug = if scan.excess > 0.0 && scan.rn > 1 && delta > 0.0 {
            Some(1.0 - scan.excess / (scan.rn - 1) as f64 / delta)
        } else {
            None
        };
        out.insert(
            genome,
            GapStats {
                rn: if fug.is_some() { scan.rn } else { 0 },
                mean: if fug.is_some() { mean } else { None },
                std: if fug.is_some() { std } else { None },
                fug,
            },
        );
    }
    Ok(out)
}
