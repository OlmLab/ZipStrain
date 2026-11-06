"""Read-observed codon evidence and explicit uncorrected selection summaries.

This module never reconstructs haplotypes by multiplying marginal frequencies.
See ``docs/v2.md`` for filters, estimators, exclusions, and runnable examples.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import permutations, product
from pathlib import Path
import heapq
import json
import re
import subprocess

import polars as pl

from zipstrain import gene_info as gi


# NCBI tables 1 and 11 have identical internal amino-acid assignments.
_BASE_ORDER = "TCAG"
_AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
GENETIC_CODE = dict(zip(("".join(c) for c in product(_BASE_ORDER, repeat=3)), _AMINO_ACIDS))
COUNT_SCHEMA = {"gene_id": pl.UInt64, "codon_index": pl.UInt32,
                "codon": pl.String, "count": pl.UInt64}
METHOD = "sns_snv_reference_opportunities_v2"
ALLELE_INTEGRATION_MODES = ("consensus", "weighted")


@lru_cache(maxsize=64)
def synonymous_sites(codon: str) -> float:
    """Return NG-style synonymous opportunities (0..3) for a sense codon.

    Each of the nine possible single-base mutations has weight 1/3. Invalid
    codons and stops are rejected. Tables 1 and 11 share this assignment.
    """
    if codon not in GENETIC_CODE or GENETIC_CODE[codon] == "*":
        raise ValueError(f"Expected a sense codon, received {codon!r}")
    return sum(GENETIC_CODE[codon[:i]+b+codon[i+1:]] == GENETIC_CODE[codon]
               for i in range(3) for b in "ACGT" if b != codon[i])/3.0


@lru_cache(maxsize=4096)
def codon_changes(source: str, target: str) -> tuple[float, float]:
    """Count synonymous/nonsynonymous changes across shortest sense-only paths.

    Average every valid nucleotide-change order equally, using the full source
    and target codons, including shared non-reference bases. Reject invalid or
    stop endpoints and pairs with no sense-only path. Identical codons return
    zero changes. This is a counting component, not a corrected rate estimator.
    """
    synonymous_sites(source)
    synonymous_sites(target)
    differences = [i for i in range(3) if source[i] != target[i]]
    paths = []
    for order in permutations(differences):
        current, synonymous = source, 0
        for index in order:
            mutated = current[:index]+target[index]+current[index+1:]
            if GENETIC_CODE[mutated] == "*":
                break
            synonymous += GENETIC_CODE[current] == GENETIC_CODE[mutated]
            current = mutated
        else:
            paths.append(synonymous)
    if not paths:
        raise ValueError("No sense-only substitution path")
    s = sum(paths)/len(paths)
    return s, len(differences)-s


def aligned_read_bases(position: int, cigar: str, sequence: str, qualities: str) -> dict:
    """Map one SAM alignment to 1-based reference bases and query indices.

    Return ``{position: (base, phred_quality, query_index)}``. Clipping,
    insertions, deletions, skips, padding, '=' and 'X' are honored. Reject
    malformed CIGARs, mismatched lengths, and missing qualities. Query indices
    let callers exclude triplets interrupted by an insertion.
    """
    tokens = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if position < 1:
        raise ValueError("SAM reference positions must be one-based")
    if not tokens or "".join(n+op for n, op in tokens) != cigar:
        raise ValueError(f"Invalid CIGAR: {cigar}")
    if qualities == "*" or len(qualities) != len(sequence):
        raise ValueError("Missing or mismatched base qualities")
    query, reference, result = 0, position, {}
    for length_text, op in tokens:
        length = int(length_text)
        if length < 1:
            raise ValueError("CIGAR operations must be positive")
        if op in "M=X":
            if query+length > len(sequence):
                raise ValueError("CIGAR exceeds sequence length")
            for offset in range(length):
                result[reference+offset] = (sequence[query+offset].upper(), ord(qualities[query+offset])-33, query+offset)
            reference += length
            query += length
        elif op in "IS":
            query += length
        elif op in "DN":
            reference += length
    if query != len(sequence):
        raise ValueError("CIGAR query length does not match sequence")
    return result


def read_codon_observations(bases: dict, gene: dict, min_baseq: int = 13):
    """Yield ``(codon_index, coding_triplet, minimum_quality)`` for a read/gene.

    Require three consecutive query bases covering a complete reference codon.
    Exclude ambiguous bases, low-quality calls, insertions within a triplet,
    deletions, and skips. Minus-strand genes are reverse-complemented. Partial
    terminal codons and the explicit initial phase are excluded.
    """
    if not bases:
        return
    origin = gene["start"]+gene["phase"] if gene["strand"] == 1 else gene["end"]-gene["phase"]
    lo, hi = min(bases), max(bases)
    for index in range(max(0, (hi-origin)//-3) if gene["strand"] == -1 else max(0, (lo-origin)//3),
                       (gene["length"]-gene["phase"])//3):
        positions = [origin + gene["strand"]*(index*3+offset) for offset in range(3)]
        if gene["strand"] == 1 and positions[0] > hi:
            break
        if gene["strand"] == -1 and positions[0] < lo:
            break
        if not all(p in bases for p in positions):
            continue
        calls = [bases[p] for p in positions]
        if any(b not in "ACGT" or q < min_baseq for b, q, _ in calls):
            continue
        if any(calls[i+1][2]-calls[i][2] != gene["strand"] for i in (0, 1)):
            continue
        codon = "".join(b for b, _, _ in calls)
        if gene["strand"] == -1:
            codon = codon.translate(str.maketrans("ACGT", "TGCA"))
        yield index, codon, min(q for _, q, _ in calls)


def write_codon_counts(bam_file, profile_file, gene_info, output_file, *,
                       min_baseq=13, min_mapq=0, min_read_ani=None,
                       read_inclusion="all-mapped") -> Path:
    """Collect actual full-read triplets from a coordinate-sorted BAM via samtools.

    Counts are per fragment: overlapping mates with the same QNAME/RG contribute
    once per gene/codon. Choose the higher minimum-quality triplet; equal-quality
    disagreement discards that observation. Never bridge an unobserved mate gap.
    The active mate cache expires as reference coordinates advance. Observations
    are flushed per scaffold, then gated against alleles retained in the position
    profile. Base qualities here are original SAM qualities, not BAQ-adjusted;
    the retained-position gate also applies the profiler's error/frequency policy.
    Secondary, supplementary, duplicate, QC-failed and unmapped reads are excluded.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq
    from zipstrain.profile import _sam_alignment_passes_profile_filters
    from zipstrain.utils import sha256_file

    if min_baseq < 0 or min_mapq < 0:
        raise ValueError("Quality thresholds must be nonnegative")
    metadata = gi.validate_profile_assets([profile_file], gene_info)
    genes = gi.read_gene_info(gene_info)
    bins = {}
    for gene in genes.iter_rows(named=True):
        for bucket in range(gene["start"]//16384, gene["end"]//16384+1):
            bins.setdefault((gene["chrom"], bucket), []).append(gene)
    metadata.update({"zipstrain_codon_evidence": "read_phased_v1",
                     "zipstrain_source_profile_sha256": sha256_file(profile_file),
                     "zipstrain_codon_quality": "original_phred_then_retained_position_gate",
                     "zipstrain_codon_filters": json.dumps({"min_baseq": min_baseq, "min_mapq": min_mapq,
                         "min_read_ani": min_read_ani, "read_inclusion": read_inclusion}, sort_keys=True),
                     "zipstrain_selection_method": METHOD})
    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    raw = output.with_suffix(".raw.parquet")
    schema = pa.schema([("gene_id", pa.uint64()), ("codon_index", pa.uint32()),
                        ("codon", pa.string()), ("count", pa.uint64())])
    counts, active, expiry = Counter(), {}, []
    chrom, previous, finished = None, -1, set()
    process = subprocess.Popen(["samtools", "view", str(bam_file)], stdout=subprocess.PIPE,
                               stderr=subprocess.DEVNULL, text=True)
    try:
        with pq.ParquetWriter(raw, schema, compression="zstd") as writer:
            for line in process.stdout:
                fields = line.rstrip("\n").split("\t")
                flag, start, mapq = int(fields[1]), int(fields[3]), int(fields[4])
                if flag & 0xF04:
                    continue
                if not _sam_alignment_passes_profile_filters(flag=flag, mapq=mapq, cigar=fields[5],
                    optional_fields=fields[11:], min_mapq=min_mapq, min_read_ani=min_read_ani,
                    read_inclusion=read_inclusion, rname=fields[2], rnext=fields[6]):
                    continue
                if fields[2] != chrom:
                    if fields[2] in finished:
                        raise ValueError("Codon profiling requires coordinate-sorted BAM")
                    if chrom is not None:
                        finished.add(chrom)
                    if counts:
                        writer.write_table(pa.Table.from_pylist([dict(gene_id=g, codon_index=i, codon=c, count=n)
                            for (g, i, c), n in counts.items() if n > 0], schema=schema))
                    counts, active, expiry = Counter(), {}, []
                    chrom, previous = fields[2], -1
                if start < previous:
                    raise ValueError("Codon profiling requires coordinate-sorted BAM")
                previous = start
                while expiry and expiry[0][0] < start:
                    _end, key = heapq.heappop(expiry)
                    active.pop(key, None)
                if fields[10] == "*":
                    continue
                bases = aligned_read_bases(start, fields[5], fields[9], fields[10])
                if not bases:
                    continue
                candidates = {}
                for bucket in range(start//16384, max(bases)//16384+1):
                    for gene in bins.get((chrom, bucket), []):
                        candidates[gene["gene_id"]] = gene
                rg = next((f for f in fields[11:] if f.startswith("RG:Z:")), "")
                for gene in candidates.values():
                    for index, codon, quality in read_codon_observations(bases, gene, min_baseq):
                        key = (fields[0], rg, gene["gene_id"], index)
                        old = active.get(key)
                        if old is not None:
                            old_codon, old_quality = old
                            if quality < old_quality or (quality == old_quality and codon == old_codon):
                                continue
                            if old_codon is not None:
                                counts[gene["gene_id"], index, old_codon] -= 1
                            if quality == old_quality:
                                active[key] = (None, quality)
                                continue
                        else:
                            heapq.heappush(expiry, (max(bases), key))
                        counts[gene["gene_id"], index, codon] += 1
                        active[key] = (codon, quality)
            if counts:
                writer.write_table(pa.Table.from_pylist([dict(gene_id=g, codon_index=i, codon=c, count=n)
                    for (g, i, c), n in counts.items() if n > 0], schema=schema))
        if process.wait() != 0:
            raise RuntimeError("samtools view failed during codon profiling")
        # A triplet must consist entirely of alleles retained by the profiler.
        observed = pl.scan_parquet(raw).join(genes.lazy(), on="gene_id", how="left")
        valid = observed
        for offset in range(3):
            pos = pl.when(pl.col("strand") == 1).then(pl.col("start")+pl.col("phase")+3*pl.col("codon_index")+offset).otherwise(
                pl.col("end")-pl.col("phase")-3*pl.col("codon_index")-offset).cast(pl.Int64)
            base = pl.col("codon").str.slice(offset, 1)
            base = pl.when(pl.col("strand") == 1).then(base).otherwise(base.replace_strict(dict(zip("ACGT", "TGCA"))))
            valid = valid.with_columns(pos.alias("pos"), base.alias("_base")).join(
                pl.scan_parquet(profile_file).select("chrom", pl.col("pos").cast(pl.Int64), "A", "C", "G", "T"),
                on=["chrom", "pos"], how="inner").filter(
                    pl.any_horizontal([(pl.col("_base") == b) & (pl.col(b) > 0) for b in "ACGT"])).drop("A", "C", "G", "T", "_base", "pos")
        valid.select(list(COUNT_SCHEMA)).sort(["gene_id", "codon_index", "codon"]).sink_parquet(
            output, compression="zstd", metadata=metadata)
    finally:
        if process.poll() is None:
            process.kill()
        process.wait()
        process.stdout.close()
        raw.unlink(missing_ok=True)
    return output


def read_codon_counts(path, gene_info) -> pl.LazyFrame:
    """Authenticate and validate observed-codon evidence before selection analysis.

    Reject marginal/reconstructed evidence, invalid triplets, zero counts,
    duplicate haplotype rows, unknown genes, and out-of-range codon indices.
    The adjacent reference-codon file must carry the same annotation identity.
    """
    gi.validate_profile_assets([path], gene_info)
    if pl.read_parquet_metadata(path).get("zipstrain_codon_evidence") != "read_phased_v1":
        raise ValueError("Selection analysis requires read-phased codon evidence")
    reference_file = Path(gene_info).resolve().with_name("reference_codons.parquet")
    gi.validate_profile_assets([reference_file], gene_info)
    from zipstrain.utils import sha256_file
    if sha256_file(reference_file) != pl.read_parquet_metadata(gene_info).get("zipstrain_reference_codons_sha256"):
        raise ValueError("Reference codon contents do not match the preprofile catalog")
    frame = pl.scan_parquet(path).select(list(COUNT_SCHEMA)).cast(COUNT_SCHEMA)
    if frame.filter(pl.col("count").is_null() | (pl.col("count") == 0) | pl.col("codon").is_null() |
                    ~pl.col("codon").str.contains(r"^[ACGT]{3}$")).limit(1).collect().height:
        raise ValueError("Invalid observed codon or count")
    if frame.select(pl.struct("gene_id", "codon_index", "codon").is_duplicated().any()).collect().item():
        raise ValueError("Duplicate observed codon rows")
    reference = pl.scan_parquet(reference_file).select("gene_id", "codon_index")
    if frame.join(reference, on=["gene_id", "codon_index"], how="anti").limit(1).collect().height:
        raise ValueError("Unknown gene ID or codon index")
    return frame


@lru_cache(maxsize=12288)
def site_changes(source: str, target: str, offset: int) -> tuple[float, float]:
    """Average the effect of one differing site across full-codon sense-only paths.

    This preserves the SNS/SNV site definitions from ``gene-overlap`` while
    fixing its reference-background substitution bug. Other shared nonreference
    bases remain present; multiple differing sites use all shortest sense-only
    orders. The three per-site results sum to :func:`codon_changes`.
    """
    synonymous_sites(source)
    synonymous_sites(target)
    if offset not in (0, 1, 2):
        raise ValueError("Codon offset must be 0, 1, or 2")
    if source[offset] == target[offset]:
        return 0.0, 0.0
    differences = [i for i in range(3) if source[i] != target[i]]
    effects = []
    for order in permutations(differences):
        current, synonymous = source, None
        for index in order:
            mutated = current[:index]+target[index]+current[index+1:]
            if GENETIC_CODE[mutated] == "*":
                break
            if index == offset:
                synonymous = float(GENETIC_CODE[current] == GENETIC_CODE[mutated])
            current = mutated
        else:
            effects.append(synonymous)
    if not effects:
        raise ValueError("No sense-only substitution path")
    s = sum(effects)/len(effects)
    return s, 1.0-s


def codon_groups(path, gene_info, *, min_cov=5, min_freq=0.0):
    """Yield sorted callable ``(key, reference_codon, Counter)`` groups in batches.

    Authenticate observed counts and reference definitions. Stage the lazy
    join/sort to a temporary Parquet and iterate Arrow batches, retaining at most
    one 64-haplotype codon group in Python. ``min_freq`` is an optional additional
    *codon* frequency floor, default zero; position-level allele filtering has
    already occurred and is not applied again at joint frequency.
    """
    import tempfile
    import pyarrow.parquet as pq

    if min_cov < 1 or not 0 <= min_freq <= 1:
        raise ValueError("Invalid codon coverage or frequency threshold")
    counts = read_codon_counts(path, gene_info)
    references = pl.scan_parquet(Path(gene_info).resolve().with_name("reference_codons.parquet"))
    with tempfile.TemporaryDirectory(prefix="zipstrain_codon_groups_") as tmp:
        staged = Path(tmp)/"groups.parquet"
        counts.join(references, on=["gene_id", "codon_index"], how="inner").sort(
            ["gene_id", "codon_index", "codon"]).sink_parquet(staged, compression="zstd")
        key, reference, observed = None, None, Counter()
        for batch in pq.ParquetFile(staged).iter_batches(batch_size=65536):
            for row in batch.to_pylist():
                new_key = row["gene_id"], row["codon_index"]
                if key is not None and new_key != key:
                    depth = sum(observed.values())
                    retained = Counter({c:n for c,n in observed.items() if n >= depth*min_freq})
                    if sum(retained.values()) >= min_cov and GENETIC_CODE.get(reference, "*") != "*":
                        yield key, reference, retained
                    observed = Counter()
                key, reference = new_key, row["reference_codon"]
                observed[row["codon"]] += row["count"]
        if key is not None:
            depth = sum(observed.values())
            retained = Counter({c:n for c,n in observed.items() if n >= depth*min_freq})
            if sum(retained.values()) >= min_cov and GENETIC_CODE.get(reference, "*") != "*":
                yield key, reference, retained


def allele_events(left: Counter, right: Counter | None, reference_base: str,
                  integration="consensus") -> tuple[str, list]:
    """Return SNS/SNV category and weighted differing allele choices for one site.

    Preserve ``gene-overlap`` semantics exactly: reference SNSs have one retained
    nonreference allele; pair SNSs have one allele on each side and differ.
    Polymorphic sites are SNVs. Consensus chooses the strongest differing choice,
    splitting exact ties equally. Weighted mode uses base frequency, or the
    product of independent sample frequencies, without renormalizing away matches.
    """
    if integration not in ALLELE_INTEGRATION_MODES:
        raise ValueError("allele_integration must be consensus or weighted")
    if right is None:
        category = "sns" if len(left) == 1 else "snv"
        candidates = [(reference_base, base, n) for base, n in left.items() if base != reference_base]
        denominator = sum(left.values())
    else:
        category = "sns" if len(left) == len(right) == 1 else "snv"
        candidates = [(a,b,n*m) for a,n in left.items() for b,m in right.items() if a != b]
        denominator = sum(left.values())*sum(right.values())
    if not candidates:
        return category, []
    if category == "sns":
        return category, [(a,b,1.0) for a,b,_ in candidates]
    if integration == "consensus":
        maximum = max(n for _,_,n in candidates)
        candidates = [(a,b,n) for a,b,n in candidates if n == maximum]
        return category, [(a,b,1/len(candidates)) for a,b,_ in candidates]
    return category, [(a,b,n/denominator) for a,b,n in candidates]


def selection_frame(totals: dict) -> pl.LazyFrame:
    """Normalize SNS/SNV change totals by callable reference opportunities.

    Each gene's totals are ``[SNS_N, SNS_S, SNV_N, SNV_S, S_sites, N_sites]``.
    Return only ``gene_id`` and the six agreed metrics. No multiple-hit correction
    is applied. Zero opportunity denominators and zero dS/pS ratios yield null.
    """
    rows = []
    for gene_id, (dn, ds, pn, ps, s_sites, n_sites) in totals.items():
        d_n, d_s = dn/n_sites if n_sites else None, ds/s_sites if s_sites else None
        p_n, p_s = pn/n_sites if n_sites else None, ps/s_sites if s_sites else None
        rows.append((gene_id, d_n, d_s, d_n/d_s if d_s and d_n is not None else None,
                     p_n, p_s, p_n/p_s if p_s and p_n is not None else None))
    return pl.DataFrame(rows, schema={"gene_id": pl.UInt64, **{c: pl.Float64 for c in gi.SELECTION_COLUMNS}}, orient="row").lazy()


def reference_selection(codon_counts, gene_info, *, min_cov=5, min_freq=0.0,
                        allele_integration="consensus") -> pl.LazyFrame:
    """Compute reference-relative SNS dN/dS and SNV pN/pS using the branch method.

    Callability and opportunities are per complete reference codon. Classification
    is per nucleotide site, so one polymorphic codon can contribute both fixed
    SNS and SNV events at different offsets. Consensus selects the most-supported
    nonreference allele (exact ties averaged); weighted mode uses all alternatives
    at their observed frequency. Consequences use full reference/observed codon
    paths conditional on the chosen allele, with contexts weighted by support.
    Invalid/stop consequences are excluded from change counts. Reference stop or
    ambiguous codons are excluded from opportunities. No Jukes–Cantor correction.
    """
    if allele_integration not in ALLELE_INTEGRATION_MODES:
        raise ValueError("allele_integration must be consensus or weighted")
    totals = {}
    for key, reference, observed in codon_groups(codon_counts, gene_info, min_cov=min_cov, min_freq=min_freq):
        total = totals.setdefault(key[0], [0.0]*6)
        sites = synonymous_sites(reference)
        total[4] += sites
        total[5] += 3-sites
        for offset in range(3):
            marginal = Counter()
            for haplotype, count in observed.items():
                marginal[haplotype[offset]] += count
            category, events = allele_events(marginal, None, reference[offset], allele_integration)
            for source_base, target_base, weight in events:
                for haplotype, count in observed.items():
                    if haplotype[offset] != target_base:
                        continue
                    try:
                        s,n = site_changes(reference, haplotype, offset)
                    except ValueError:
                        continue
                    event_weight = weight*count/marginal[target_base]
                    total[0 if category == "sns" else 2] += n*event_weight
                    total[1 if category == "sns" else 3] += s*event_weight
    return selection_frame(totals)


def pairwise_selection(codon_counts_1, codon_counts_2, gene_info, *, min_cov=5,
                       min_freq=0.0, allele_integration="consensus") -> pl.LazyFrame:
    """Compute pairwise SNS dN/dS and SNV pN/pS with the branch's integration rules.

    Merge sorted callable codons from both samples. Use reference opportunities
    over their intersection. Consensus selects the most-supported differing
    allele pair at each polymorphic site, averaging exact ties; weighted mode
    uses all differing pairs at their independent observed frequencies. Phased
    contexts are integrated conditionally, with full-codon path averaging. Thus
    GGA/GGT stays synonymous even against reference AAA. Identical polymorphic
    populations can have positive pN/pS: this measures polymorphic allele-pair
    consequences, not only substitutions separating populations. Raw proportions,
    no Jukes–Cantor correction; six output metrics, undefined ratios null.
    """
    if allele_integration not in ALLELE_INTEGRATION_MODES:
        raise ValueError("allele_integration must be consensus or weighted")
    left = iter(codon_groups(codon_counts_1, gene_info, min_cov=min_cov, min_freq=min_freq))
    right = iter(codon_groups(codon_counts_2, gene_info, min_cov=min_cov, min_freq=min_freq))
    a, b, totals = next(left, None), next(right, None), {}
    while a is not None and b is not None:
        if a[0] < b[0]:
            a = next(left, None)
            continue
        if b[0] < a[0]:
            b = next(right, None)
            continue
        key, reference, observed_a = a
        if reference != b[1]:
            raise ValueError("Reference codons do not match")
        observed_b = b[2]
        total = totals.setdefault(key[0], [0.0]*6)
        sites = synonymous_sites(reference)
        total[4] += sites
        total[5] += 3-sites
        for offset in range(3):
            marginal_a, marginal_b = Counter(), Counter()
            for haplotype, count in observed_a.items():
                marginal_a[haplotype[offset]] += count
            for haplotype, count in observed_b.items():
                marginal_b[haplotype[offset]] += count
            category, events = allele_events(marginal_a, marginal_b, reference[offset], allele_integration)
            for source_base, target_base, weight in events:
                for source, count_a in observed_a.items():
                    if source[offset] != source_base:
                        continue
                    for target, count_b in observed_b.items():
                        if target[offset] != target_base:
                            continue
                        try:
                            s,n = site_changes(source, target, offset)
                        except ValueError:
                            continue
                        event_weight = weight*count_a*count_b/(marginal_a[source_base]*marginal_b[target_base])
                        total[0 if category == "sns" else 2] += n*event_weight
                        total[1 if category == "sns" else 3] += s*event_weight
        a, b = next(left, None), next(right, None)
    return selection_frame(totals)
