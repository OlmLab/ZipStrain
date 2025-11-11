from zipstrain import profile
import polars as pl
import pytest


a_chr1=     [ 1, 0, 1, 0, 0, 0, 4, 4, 2, 1]
t_chr1=     [ 0, 2, 3, 1, 0, 0, 1, 0, 0, 1]
c_chr1=     [10, 1, 2, 0, 0, 0, 1, 0, 0, 0]
g_chr1=     [ 0, 2, 0, 2, 3, 1, 0, 2, 7, 0]
a_chr2=     [ 0, 0, 4, 3, 2, 4, 0, 1, 2, 0, 0, 2, 3, 1, 0, 0, 0, 0, 0, 0]
a_chr2_1mut=[ 0, 0, 4, 3, 2, 4, 8, 1, 2, 0, 0, 2, 3, 1, 0, 0, 0, 0, 0, 0] # A single mutation at index (0->1)
t_chr2=     [ 1, 0, 0, 1, 0, 0, 7, 0, 0, 1, 0, 2, 3, 1, 3, 2, 2, 1, 1, 2]
t_chr2_1mut=[ 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 2, 3, 1, 3, 2, 2, 1, 1, 2] # A single mutation at index (19->20)
c_chr2=     [ 0, 0, 1, 2, 1, 1, 0, 3, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0]
c_chr2_1mut=[ 0, 0, 1, 2, 1, 1, 0, 3, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0] # A single mutation at index (18->19)
g_chr2=     [ 2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1]
g_chr2_1mut=[ 2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1]
a_chr3=     [ 0, 0, 1, 0, 0, 0,15, 0, 0, 0, 0, 0, 0,11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
t_chr3=     [ 1, 0, 0, 1, 0,17, 1, 0,10,11, 0, 2, 3, 1, 3, 2,12,10, 1, 6, 1, 0, 0, 1, 0, 0, 8, 5, 0, 1]
c_chr3=     [ 0, 0,19,21,18, 1, 2,10, 0, 0,12, 1, 0, 0, 0,13, 0, 0,11, 0, 1, 0, 1, 2, 1,20, 0, 1, 0, 0]
g_chr3=     [20,20, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 2, 0, 2, 3, 1, 1, 2, 7, 0]

gene_locs = (["NA","NA","gene1","gene1","gene1","gene1","NA","NA","NA","NA"])+ \
            (["NA","NA","gene2","gene2","gene2","gene2","gene2","NA","NA","NA",
              "NA","gene3","gene3","gene3","gene3","gene3","NA","NA","NA","NA"])+ \
            (["NA","NA","NA","gene1","gene1","gene1","gene1","NA","NA","NA",
              "NA","gene2","gene2","gene2","gene2","gene2","NA","NA","NA","NA",
              "NA","NA","gene3","gene3","gene3","gene3","gene3","NA","NA","NA"])


@pytest.fixture(scope="module")
def profile_1()->pl.LazyFrame:
    return pl.DataFrame({
        "chrom": ["chr1"]*10+["chr2"]*20+["chr3"]*30,
        "pos":list(range(0,len(a_chr1)))+list(range(0,len(a_chr2)))+list(range(0,len(a_chr3))),
        "gene": gene_locs,
        "A": a_chr1 + a_chr2 + a_chr3,
        "T": t_chr1 + t_chr2 + t_chr3,
        "C": c_chr1 + c_chr2 + c_chr3,
        "G": g_chr1 + g_chr2 + g_chr3,
    }).lazy()
    

@pytest.fixture(scope="module")
def stb()->pl.LazyFrame:
    return pl.DataFrame({
        "scaffold":["chr1","chr2","chr3"],
        "genome":["genome1","genome1","genome2"],
    }).lazy()
    
@pytest.mark.parametrize("min_cov,freq_threshold", [
    (5, 0.8),
    (10, 0.5)
])
def test_get_strain_hetrogeneity(profile_1, stb, min_cov, freq_threshold):
    het_profile=profile.get_strain_hetrogeneity(profile_1, stb, min_cov=min_cov, freq_threshold=freq_threshold).collect().rows_by_key("genome",unique=True,named=True)
    genome1_freq=[
        a_chr1 + a_chr2 ,
        t_chr1 + t_chr2 ,
        c_chr1 + c_chr2 ,
        g_chr1 + g_chr2 ,
    ]
    genome2_freq=[
        a_chr3,
        t_chr3,
        c_chr3,
        g_chr3,
    ]
    assert het_profile["genome1"][f"total_sites_at_{min_cov}_coverage"]==sum(1 for pos in range(len(genome1_freq[0])) 
                                                        if sum(list(zip(*genome1_freq))[pos]) >= min_cov)
    assert het_profile["genome2"][f"total_sites_at_{min_cov}_coverage"]==sum(1 for pos in range(len(genome2_freq[0])) 
                                                        if sum(list(zip(*genome2_freq))[pos]) >= min_cov)
    assert het_profile["genome1"]["heterogeneous_sites"]== sum(1 for pos in range(len(genome1_freq[0])) 
                                                        if sum(list(zip(*genome1_freq))[pos]) >= min_cov and 
                                                           max(list(zip(*genome1_freq))[pos])/sum(list(zip(*genome1_freq))[pos]) < freq_threshold)
    assert het_profile["genome2"]["heterogeneous_sites"]== sum(1 for pos in range(len(genome2_freq[0])) 
                                                        if sum(list(zip(*genome2_freq))[pos]) >= min_cov and 
                                                           max(list(zip(*genome2_freq))[pos])/sum(list(zip(*genome2_freq))[pos]) < freq_threshold)
    
    assert het_profile["genome1"]["strain_heterogeneity"]== sum(1 for pos in range(len(genome1_freq[0])) 
                                                        if sum(list(zip(*genome1_freq))[pos]) >= min_cov and 
                                                           max(list(zip(*genome1_freq))[pos])/sum(list(zip(*genome1_freq))[pos]) < freq_threshold) / \
                                                        sum(1 for pos in range(len(genome1_freq[0])) 
                                                        if sum(list(zip(*genome1_freq))[pos]) >= min_cov)
    
    assert het_profile["genome2"]["strain_heterogeneity"]== sum(1 for pos in range(len(genome2_freq[0])) 
                                                        if sum(list(zip(*genome2_freq))[pos]) >= min_cov and 
                                                           max(list(zip(*genome2_freq))[pos])/sum(list(zip(*genome2_freq))[pos]) < freq_threshold) / \
                                                        sum(1 for pos in range(len(genome2_freq[0])) 
                                                        if sum(list(zip(*genome2_freq))[pos]) >= min_cov)
