"""Plain-Python selection estimator, kept as a test oracle for the vectorized kernel.

This is the original codon-by-codon implementation of the SNS/SNV rules. It is slow and
written to be read against the method description, so a change to the kernel in
``zipstrain.codons`` that alters results fails against an independent implementation.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import tempfile

import polars as pl
import pyarrow.parquet as pq

from zipstrain.codons import (ALLELE_INTEGRATION_MODES, GENETIC_CODE, read_codon_counts,
                              selection_frame, site_changes, synonymous_sites)


def codon_groups(path, gene_info, *, min_cov=5, min_freq=0.0):
    """Yield sorted callable ``(key, reference_codon, Counter)`` groups in batches.

    Authenticate observed counts and reference definitions. Stage the lazy
    join/sort to a temporary Parquet and iterate Arrow batches, retaining at most
    one 64-haplotype codon group in Python. ``min_freq`` is an optional additional
    *codon* frequency floor, default zero; position-level allele filtering has
    already occurred and is not applied again at joint frequency.
    """
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
