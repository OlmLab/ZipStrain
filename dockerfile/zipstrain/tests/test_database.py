from zipstrain import database,compare
import polars as pl
import pytest
import pathlib 
import random
import copy
import json
from pydantic import BaseModel, ValidationError
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
def scaffold_1()->str:
    return "chr1\nchr2\nchr3"

@pytest.fixture(scope="module")
def profile_2()->pl.LazyFrame:
    """Exactly the same as profile_1"""
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
def scaffold_2():
    return "chr1\nchr3"

@pytest.fixture(scope="module")
def profile_3()->pl.LazyFrame:
    """A profile with a single mutation in chr2 compared to profile_1"""
    return pl.DataFrame({
        "chrom": ["chr1"]*10+["chr2"]*20+["chr3"]*30,
        "pos":list(range(0,len(a_chr1)))+list(range(0,len(a_chr2)))+list(range(0,len(a_chr3))),
        "gene": gene_locs,
        "A": a_chr1 + a_chr2_1mut + a_chr3,
        "T": t_chr1 + t_chr2_1mut + t_chr3,
        "C": c_chr1 + c_chr2_1mut + c_chr3,
        "G": g_chr1 + g_chr2_1mut + g_chr3,
    }).lazy()


@pytest.fixture(scope="module")
def scaffold_3()->str:
    return "chr1\nchr2\nchr3"

@pytest.fixture(scope="module")
def profile_4_random()->pl.LazyFrame:
    """A profile with a single mutation in chr2 compared to profile_1"""
    return pl.DataFrame({
        "chrom": ["Scaffold_1"]*100+["Scaffold_2"]*200+["Scaffold_3"]*300,
        "pos":list(range(0,100))+list(range(0,200))+list(range(0,300)),
        "gene": ["NA"]*600,
        "A": random.choices(range(0,20),k=600),
        "T": random.choices(range(0,30),k=600),
        "C": random.choices(range(0,20),k=600),
        "G": random.choices(range(0,20),k=600),
    }).lazy()

@pytest.fixture(scope="module")
def scaffold_4()->str:
    return "Scaffold_1\nScaffold_2\nScaffold_3"

@pytest.fixture(scope="module")
def profile_5_random()->pl.LazyFrame:
    """A profile with a single mutation in chr2 compared to profile_1"""
    return pl.DataFrame({
        "chrom": ["Scaffold_1"]*90+["Scaffold_2"]*220+["Scaffold_3"]*400+["Scaffold_4"]*150,
        "pos":list(range(0,90))+list(range(0,220))+list(range(0,400))+list(range(0,150)),
        "gene": ["NA"]*860,
        "A": random.choices(range(0,20),k=860),
        "T": random.choices(range(0,30),k=860),
        "C": random.choices(range(0,20),k=860),
        "G": random.choices(range(0,20),k=860),
    }).lazy()

@pytest.fixture(scope="module")
def scaffold_5()->str:
    return "Scaffold_1\nScaffold_2\nScaffold_3\nScaffold_4"

@pytest.fixture(scope="module")
def stb()->pl.LazyFrame:
    return pl.DataFrame({
        "scaffold":["chr1","chr2","chr3"],
        "genome":["genome1","genome1","genome2"],
    }).lazy()

@pytest.fixture(scope="module")
def stb_4_5()->pl.LazyFrame:
    return pl.DataFrame({
        "scaffold":["Scaffold_1","Scaffold_2","Scaffold_3","Scaffold_4"],
        "genome":["genome1","genome1","genome2","genome2"],
    }).lazy()

@pytest.fixture(scope="module")
def null_model()->pl.LazyFrame:
    return pl.DataFrame({
        "cov":list(range(100)),
        "max_error_count":[int(i*0.1) for i in range(100)],
    }).lazy()

@pytest.fixture(scope="module")
def comps_lf_1_2(profile_1,profile_2,stb,null_model)->pl.LazyFrame:
    comp_1_2=compare.compare_genomes(
        mpile_contig_1=profile_1,
        mpile_contig_2=profile_2,
        scaffold_to_genome=stb,
        null_model=null_model,
        min_cov=2 ,
        min_gene_compare_len=10
    )
    comp_1_2=comp_1_2.with_columns(pl.lit("profile_1").alias("sample_1"), pl.lit("profile_2").alias("sample_2")).fill_null(0)
    return comp_1_2

@pytest.fixture(scope="module")
def comps_lf_31_32(profile_1,profile_2,profile_3,stb,null_model)->pl.LazyFrame:
    comp_1_3=compare.compare_genomes(
        mpile_contig_1=profile_1,
        mpile_contig_2=profile_3,
        scaffold_to_genome=stb,
        null_model=null_model,
        min_cov=2 ,
        min_gene_compare_len=10
    )
    comp_1_3=comp_1_3.with_columns(pl.lit("profile_1").alias("sample_1"), pl.lit("profile_3").alias("sample_2")).fill_null(0)
    
    comp_2_3=compare.compare_genomes(
        mpile_contig_1=profile_2,
        mpile_contig_2=profile_3,
        scaffold_to_genome=stb,
        null_model=null_model,
        min_cov=2 ,
        min_gene_compare_len=10
    )
    comp_2_3=comp_2_3.with_columns(pl.lit("profile_2").alias("sample_1"), pl.lit("profile_3").alias("sample_2")).fill_null(0)
    
    return pl.concat([comp_1_3,comp_2_3])

@pytest.fixture(scope="module")
def comps_lf_4_5(profile_4_random,profile_5_random,stb_4_5,null_model)->pl.LazyFrame:
    comp_4_5=compare.compare_genomes(
        mpile_contig_1=profile_4_random,
        mpile_contig_2=profile_5_random,
        scaffold_to_genome=stb_4_5,
        null_model=null_model,
        min_cov=5 ,
        min_gene_compare_len=50
    )
    comp_4_5=comp_4_5.with_columns(pl.lit("profile_4").alias("sample_1"), pl.lit("profile_5").alias("sample_2")).fill_null(0)
    return comp_4_5
    
@pytest.fixture(scope="module")
def profile_1_2_database(tmp_path_factory, profile_1, scaffold_1, profile_2, scaffold_2):
    db = database.ProfileDatabase()
    tmp_path=tmp_path_factory.mktemp("profiles")
    profile_1_dir = tmp_path / "profile_1.parquet"
    profile_1.sink_parquet(profile_1_dir)
    scaffold_1_dir = tmp_path / "profile_1.parquet.scaffold"
    scaffold_1_dir.write_text(scaffold_1)
    profile_2_dir = tmp_path / "profile_2.parquet"
    profile_2.sink_parquet(profile_2_dir)
    scaffold_2_dir = tmp_path / "profile_2.parquet.scaffold"
    scaffold_2_dir.write_text(scaffold_2)
    db.add_profile({
        "profile_name": "profile_1",
        "profile_location": str(profile_1_dir),
        "scaffold_location": str(scaffold_1_dir),
        "reference_db_id": "ref_1",
        "gene_db_id": "gene_ref_1"
    })
    db.add_profile({
        "profile_name": "profile_2",
        "profile_location": str(profile_2_dir),
        "scaffold_location": str(scaffold_2_dir),
        "reference_db_id": "ref_2",
        "gene_db_id": "gene_ref_2"
    })
    return db

@pytest.fixture(scope="module")
def profile_1_2_3_database(tmp_path_factory, profile_1, scaffold_1, profile_2, scaffold_2, profile_3, scaffold_3):
    db = database.ProfileDatabase()
    tmp_path=tmp_path_factory.mktemp("profiles")
    profile_1_dir = tmp_path / "profile_1.parquet"
    profile_1.sink_parquet(profile_1_dir)
    scaffold_1_dir = tmp_path / "profile_1.parquet.scaffold"
    scaffold_1_dir.write_text(scaffold_1)
    profile_2_dir = tmp_path / "profile_2.parquet"
    profile_2.sink_parquet(profile_2_dir)
    scaffold_2_dir = tmp_path / "profile_2.parquet.scaffold"
    scaffold_2_dir.write_text(scaffold_2)
    profile_3_dir = tmp_path / "profile_3.parquet"
    profile_3.sink_parquet(profile_3_dir)
    scaffold_3_dir = tmp_path / "profile_3.parquet.scaffold"
    scaffold_3_dir.write_text(scaffold_3)
    db.add_profile({
        "profile_name": "profile_1",
        "profile_location": str(profile_1_dir),
        "scaffold_location": str(scaffold_1_dir),
        "reference_db_id": "ref_1",
        "gene_db_id": "gene_ref_1"
    })
    db.add_profile({
        "profile_name": "profile_2",
        "profile_location": str(profile_2_dir),
        "scaffold_location": str(scaffold_2_dir),
        "reference_db_id": "ref_2",
        "gene_db_id": "gene_ref_2"
    })
    db.add_profile({
        "profile_name": "profile_3",
        "profile_location": str(profile_3_dir),
        "scaffold_location": str(scaffold_3_dir),
        "reference_db_id": "ref_3",
        "gene_db_id": "gene_ref_3"
    })
    return db

@pytest.fixture(scope="module")
def profile_4_5_database(tmp_path_factory, profile_4_random, scaffold_4, profile_5_random, scaffold_5):
    db = database.ProfileDatabase()
    tmp_path=tmp_path_factory.mktemp("profiles")
    profile_4_dir = tmp_path / "profile_4.parquet"
    profile_4_random.sink_parquet(profile_4_dir)
    scaffold_4_dir = tmp_path / "profile_4.parquet.scaffold"
    scaffold_4_dir.write_text(scaffold_4)
    profile_5_dir = tmp_path / "profile_5.parquet"
    profile_5_random.sink_parquet(profile_5_dir)
    scaffold_5_dir = tmp_path / "profile_5.parquet.scaffold"
    scaffold_5_dir.write_text(scaffold_4)
    db.add_profile({
        "profile_name": "profile_4",
        "profile_location": str(profile_4_dir),
        "scaffold_location": str(scaffold_4_dir),
        "reference_db_id": "ref_4",
        "gene_db_id": "gene_ref_4"
    })
    db.add_profile({
        "profile_name": "profile_5",
        "profile_location": str(profile_5_dir),
        "scaffold_location": str(scaffold_5_dir),
        "reference_db_id": "ref_5",
        "gene_db_id": "gene_ref_5"
    })
    return db




@pytest.fixture(scope="module")
def simple_genome_compare_config(stb,null_model,tmp_path_factory):
    tmp_path=tmp_path_factory.mktemp("genome_compare")
    stb_dir = tmp_path / "stb.tsv"
    stb.sink_csv(stb_dir, separator="\t")
    null_model_dir = tmp_path / "null_model.parquet"
    null_model.sink_parquet(null_model_dir)
    return database.GenomeComparisonConfig(
            reference_id="ref_1",
            gene_db_id="gene_ref_1",
            min_cov=5, 
            scope="all",
            min_gene_compare_len=100,
            stb_file_loc=str(stb_dir),
            null_model_loc=str(null_model_dir),
        )
    

def test_profile_create_new_empty_database():
    db = database.ProfileDatabase()
    assert db.db.collect().height == 0
    assert db.db.collect().columns == ["profile_name", "profile_location", "scaffold_location", "reference_db_id", "gene_db_id"]
    assert isinstance(db.db, pl.LazyFrame)
    db._validate_db()
    
def test_profile_add_profile_1_to_empty_database(tmp_path,profile_1,scaffold_1):
    db = database.ProfileDatabase()
    profile_1_dir = tmp_path / "profile_1.parquet"
    profile_1.sink_parquet(profile_1_dir)
    scaffold_1_dir = tmp_path / "profile_1.parquet.scaffold"
    scaffold_1_dir.write_text(scaffold_1)
    db.add_profile({
        "profile_name": "profile_1",
        "profile_location": str(profile_1_dir),
        "scaffold_location": str(scaffold_1_dir),
        "reference_db_id": "ref_1",
        "gene_db_id": "gene_ref_1"
    })
    df = db.db.collect()
    assert df.height == 1
    assert set(df["profile_name"].to_list()) == {"profile_1"}
    assert set(df["reference_db_id"].to_list()) == {"ref_1"}
    assert set(df["gene_db_id"].to_list()) == {"gene_ref_1"}
    assert all(isinstance(x, str) for x in df["profile_location"].to_list())
    assert all(isinstance(x, str) for x in df["scaffold_location"].to_list())

def test_profile_add_faulty_profile_1_to_empty_database(tmp_path,profile_1,scaffold_1):
    db = database.ProfileDatabase()
    profile_1_dir = tmp_path / "profile_1.parquet"
    profile_1.sink_parquet(profile_1_dir)
    scaffold_1_dir = tmp_path / "profile_1.parquet.scaffold"
    scaffold_1_dir.write_text(scaffold_1)
    
    bad_profile_scaffold_and_profile = {
        "profile_name": "faulty_profile",
        "profile_location": "/does/not/exist.parquet",
        "scaffold_location": "/does/not/exist.scaffold",
        "reference_db_id": "ref_bad",
        "gene_db_id": "gene_ref_bad"
    }
    
    bad_profile_missing_attr = {
        "profile_location":str(profile_1_dir),
        "scaffold_location": str(scaffold_1_dir),
        "reference_db_id": "ref_bad",
        "gene_db_id": "gene_ref_bad"
    }     
    bad_profile_extra_attr = {
        "extra_attr":"some stuf",
        "profile_name": "faulty_profile",
        "profile_location":str(profile_1_dir),
        "scaffold_location": str(scaffold_1_dir),
        "reference_db_id": "ref_bad",
        "gene_db_id": "gene_ref_bad"
    }     
    bad_profile_wrong_attr = {
        "profile_nam": "faulty_profile",
        "profile_location":str(profile_1_dir),
        "scaffold_location": str(scaffold_1_dir),
        "reference_db_id": "ref_bad",
        "gene_db_id": "gene_ref_bad"
    }
    bad_profile_empty_ref_db_id = {
        "profile_name": "faulty_profile",
        "profile_location":str(profile_1_dir),
        "scaffold_location": str(scaffold_1_dir),
        "reference_db_id": "",
        "gene_db_id": "gene_ref_bad"
    }     
    bad_profile_empty_gene_db_id = {
        "profile_name": "faulty_profile",
        "profile_location":str(profile_1_dir),
        "scaffold_location": str(scaffold_1_dir),
        "reference_db_id": "ref_bad",
        "gene_db_id": ""
    }

    #Fails for both scaffold and profile
    with pytest.raises(ValueError, match="The profile data provided is not valid: 2 validation errors"):
        db.add_profile(bad_profile_scaffold_and_profile)
    
    with pytest.raises(ValueError, match="The profile data provided is not valid: 2 validation error"):
        db.add_profile(bad_profile_wrong_attr)
    
    with pytest.raises(ValueError, match="The profile data provided is not valid: 1 validation error"):
        db.add_profile(bad_profile_missing_attr)
    
    with pytest.raises(ValueError, match="The profile data provided is not valid: 1 validation error"):
        db.add_profile(bad_profile_extra_attr)
    
    with pytest.raises(ValueError, match="The profile data provided is not valid: 1 validation error"):
        db.add_profile(bad_profile_empty_ref_db_id)
    
    with pytest.raises(ValueError, match="The profile data provided is not valid: 1 validation error"):
        db.add_profile(bad_profile_empty_gene_db_id)
    


def test_profile_read_from_csv(profile_1_2_database,tmp_path):
    profile_1_2_database.db.sink_csv(tmp_path/"profile_1_2.csv")
    new_db=database.ProfileDatabase.from_csv(tmp_path/"profile_1_2.csv")
    assert new_db.db.collect().height == 2
    assert new_db.db.collect().columns == ["profile_name", "profile_location", "scaffold_location", "reference_db_id", "gene_db_id"]
    assert isinstance(new_db.db, pl.LazyFrame)

    
def test_profile_test_write_database(profile_1_2_database,tmp_path):
    assert profile_1_2_database.db_loc is None
    profile_1_2_database.save_as_new_database(tmp_path/"db_12.parquet")
    assert profile_1_2_database.db_loc==pathlib.Path(tmp_path/"db_12.parquet")
    db=database.ProfileDatabase(tmp_path/"db_12.parquet")
    assert db.db.collect().height == 2
    assert db.db.collect().columns == ["profile_name", "profile_location", "scaffold_location", "reference_db_id", "gene_db_id"]
    assert isinstance(db.db, pl.LazyFrame)

def test_profile_update_database(profile_1_2_database,tmp_path):
    new_db=database.ProfileDatabase()
    new_db.save_as_new_database(tmp_path/"new_db.parquet")
    new_db_loaded=database.ProfileDatabase(db_loc=(tmp_path/"new_db.parquet"))
    assert new_db_loaded.db.collect().height == 0
    new_db_loaded.add_database(profile_1_2_database)
    new_db_loaded.update_database()
    updated_db_loaded=database.ProfileDatabase(db_loc=(tmp_path/"new_db.parquet"))
    assert updated_db_loaded.db.collect().height == 2

def test_profile_write_csv(profile_1_2_database,tmp_path):
    profile_1_2_database.to_csv(tmp_path/"tmp.csv")
    loaded_from_csv=database.ProfileDatabase.from_csv(tmp_path/"tmp.csv")
    assert loaded_from_csv.db.collect().sort("profile_name").equals(profile_1_2_database.db.collect(engine="streaming").sort("profile_name"))


def test_profile_add_profile_1_2_to_empty_database(profile_1_2_database):
    db = database.ProfileDatabase()
    db.add_database(profile_1_2_database)
    df = db.db.collect()
    assert df.height == 2
    assert set(df["profile_name"].to_list()) == {"profile_1", "profile_2"}
    assert set(df["reference_db_id"].to_list()) == {"ref_1", "ref_2"}
    assert all(isinstance(x, str) for x in df["profile_location"].to_list())
    assert all(isinstance(x, str) for x in df["scaffold_location"].to_list())
    db._validate_db()


def test_genome_compare_config_faulty()->None:
    """tests the response of GenomeComparisonConfig to wrong inputs"""
    with pytest.raises(ValidationError, match="3 validation errors for GenomeComparisonConfig"):
        database.GenomeComparisonConfig(
            reference_id="blah",
            min_cov="5das", # This should be an int
            scope="all",
            min_gene_compare_len=100,
        )
    ### But passing string is okay if they can be turned to their numeric counterparts
    database.GenomeComparisonConfig(
            reference_id="blah",
            min_cov="5", # This should be an int
            scope="all",
            min_gene_compare_len=100,
            stb_file_loc="Somefile",
            null_model_loc="Somefile"
        )
    ### An extra attribute that is not expected
    with pytest.raises(ValueError, match="1 validation error for GenomeComparisonConfig"):
        database.GenomeComparisonConfig(
            minimum_similarity="Wrongatte",
            reference_id="blah",
            min_cov="5", # This should be an int
            scope="all",
            min_gene_compare_len=100,
            stb_file_loc="Somefile",
            null_model_loc="Somefile"
        )

    ### A missing attr that must be provided
    with pytest.raises(ValueError, match="1 validation error for GenomeComparisonConfig"):
        database.GenomeComparisonConfig(
            min_cov="5", # This should be an int
            scope="all",
            min_gene_compare_len=100,
            stb_file_loc="Somefile",
            null_model_loc="Somefile"
        )
        
def test_genome_compare_config_compatibility(simple_genome_compare_config):
    """This test examines the compatibility of two config objects to see if the two Genome comparison objects can be merged"""

    config_2 = copy.copy(simple_genome_compare_config)
    config_2.min_cov=10 #Change min_cov to make it incompatible
    config_1 = simple_genome_compare_config
    config_3=copy.copy(simple_genome_compare_config)
    config_3.reference_id="ref_2" #Change reference_id to make it incompatible
    config_4 = copy.copy(simple_genome_compare_config)
    config_4.scope="genome1,genome2" #Change scope to something that is still compatible

    config_5 = copy.copy(simple_genome_compare_config)
    config_5.scope="genome1" #Change scope to something that is still compatible but with smaller scope from config_4
    config_6 = copy.copy(simple_genome_compare_config)
    config_6.scope="genome3" #Change scope to something that is incompatible with config_4
    
     #Tests
    assert not config_1.is_compatible(config_2) #Should fail because min_cov and min_gene_compare_len are different
    assert not config_1.is_compatible(config_3) #Should fail because reference_id is different
    assert config_1.is_compatible(config_4) #Should pass because scope is compatible (all is compatible with anything)
    assert config_4.is_compatible(config_5) #Should pass because scope is compatible (genome1 is in genome1,genome2)
    assert not config_4.is_compatible(config_6) #Should fail because scope is different (genome3 is not in genome1,genome2)
    
def test_genome_compare_config_io(tmp_path,simple_genome_compare_config):
    """This test examines the input/output functionality of the GenomeComparisonConfig class"""
    # Test serialization
    config_dict = simple_genome_compare_config.to_dict()
    assert config_dict == {
        "reference_id": "ref_1",
        "gene_db_id": "gene_ref_1",
        "min_cov": 5,
        "scope": "all",
        "min_gene_compare_len": 100,
        "null_model_p_value":simple_genome_compare_config.null_model_p_value,
        "stb_file_loc": simple_genome_compare_config.stb_file_loc,
        "null_model_loc": simple_genome_compare_config.null_model_loc,
    }
    # Test to_json
    json_path = tmp_path / "config.json"
    simple_genome_compare_config.to_json(json_path)
    assert json_path.exists()
    # Test from_json
    loaded_config = database.GenomeComparisonConfig.from_json(json_path)
    assert loaded_config == simple_genome_compare_config

def test_genome_compare_config_get_maximal_scope(simple_genome_compare_config):
    """This test examines the get_maximal_scope method of the GenomeComparisonConfig class"""
    other_config = copy.copy(simple_genome_compare_config)
    other_config.scope="genome1,genome2"
    maximal_scope = simple_genome_compare_config.get_maximal_scope_config(other_config)
    assert maximal_scope.scope == "genome1,genome2"
    other_config_2 = copy.copy(simple_genome_compare_config)
    other_config_2.scope="genome2,genome3"
    maximal_scope = other_config.get_maximal_scope_config(other_config_2)
    assert maximal_scope.scope == "genome2"
    #Order should not matter in scope
    other_config_2.scope="genome2,genome1"
    maximal_scope = other_config.get_maximal_scope_config(other_config_2)
    assert maximal_scope.scope == "genome1,genome2"

def test_genome_comparison_database_create_new_empty_database(profile_1_2_database,simple_genome_compare_config):
    db = database.GenomeComparisonDatabase(
        profile_db=profile_1_2_database,
        config=simple_genome_compare_config,
        comp_db_loc=None
    )
    assert db.comp_db.collect().height == 0
    assert db.comp_db.collect().columns == database.GenomeComparisonDatabase.COLUMN_NAMES
    assert isinstance(db.comp_db, pl.LazyFrame)
    assert db.comp_db_loc is None
    db._validate_db()

def test_genome_comparison_database_create_faulty_new_db(profile_1_2_database,simple_genome_compare_config,tmp_path):
    """This test examines the response of GenomeComparisonDatabase to the scenario where there are samples in the comparison database that are not in the profile database"""
    faulty_comp_db = pl.DataFrame({
        "genome":["genome1","genome1"],
        "total_positions":[100,200],
        "share_allele_pos":[90,180],
        "genome_pop_ani":[0.9,0.9],
        "max_consecutive_length":[10,20],
        "shared_genes_count":[5,10],
        "identical_gene_count":[3,6],
        "perc_id_genes":[60,60],
        "sample_1":["profile_1","profile_3"], # profile_3.parquet does not exist in profile database
        "sample_2":["profile_2","profile_2"],
    }).lazy() 
    faulty_comp_db_dir = tmp_path / "faulty_comp_db.parquet"
    faulty_comp_db.sink_parquet(faulty_comp_db_dir)

    db = database.GenomeComparisonDatabase(
        profile_db=profile_1_2_database,
        config=simple_genome_compare_config,
        comp_db_loc=faulty_comp_db_dir
    )
    
    with pytest.raises(ValueError, match="The following profile names are in the comparison database but not in the profile database: {'profile_3'}"):
        db._validate_db()

def test_genome_comparison_database_test_get_allprofiles_names(profile_1_2_database,simple_genome_compare_config,tmp_path):
    faulty_comp_db = pl.DataFrame({
        "genome":["genome1","genome1"],
        "total_positions":[100,200],
        "share_allele_pos":[90,180],
        "genome_pop_ani":[0.9,0.9],
        "max_consecutive_length":[10,20],
        "shared_genes_count":[5,10],
        "identical_gene_count":[3,6],
        "perc_id_genes":[60,60],
        "sample_1":["profile_1","profile_3"], # profile_3.parquet does not exist in profile database
        "sample_2":["profile_2","profile_2"],
    }).lazy() 
    faulty_comp_db_dir = tmp_path / "faulty_comp_db.parquet"
    faulty_comp_db.sink_parquet(faulty_comp_db_dir)

    db = database.GenomeComparisonDatabase(
        profile_db=profile_1_2_database,
        config=simple_genome_compare_config,
        comp_db_loc=faulty_comp_db_dir
    )
    assert set(db.get_all_profile_names()) == {"profile_1", "profile_2","profile_3"}


def test_genome_comparison_database_get_remaining_pairs(profile_1_2_3_database,profile_1_2_database,simple_genome_compare_config,tmp_path):
    comp_db = pl.DataFrame({
        "genome":["genome1"],
        "total_positions":[100],
        "share_allele_pos":[90],
        "genome_pop_ani":[0.9],
        "max_consecutive_length":[10],
        "shared_genes_count":[5],
        "perc_id_genes":[3],
        "sample_1":["profile_1"],
        "sample_2":["profile_2"],
    }).lazy() 
    comp_db_dir = tmp_path / "comp_db.parquet"
    comp_db.sink_parquet(comp_db_dir)

    db = database.GenomeComparisonDatabase(
        profile_db=profile_1_2_3_database,
        config=simple_genome_compare_config,
        comp_db_loc=comp_db_dir
    )
    
    empty_db=pl.DataFrame({
        "genome":[],
        "total_positions":[],
        "share_allele_pos":[],
        "genome_pop_ani":[],
        "max_consecutive_length":[],
        "shared_genes_count":[],
        "identical_gene_count":[],
        "perc_id_genes":[],
        "sample_1":[],
        "sample_2":[],
    },schema={
        "genome": pl.Utf8,
        "total_positions": pl.Int64,
        "share_allele_pos": pl.Int64,
        "genome_pop_ani": pl.Float64,
        "max_consecutive_length": pl.Int64,
        "shared_genes_count": pl.Int64,
        "identical_gene_count": pl.Int64,
        "perc_id_genes": pl.Float64,
        "sample_1": pl.Utf8,
        "sample_2": pl.Utf8,
    }).lazy() 
    empty_db_dir = tmp_path / "empty_db.parquet"
    empty_db.sink_parquet(empty_db_dir)
    db_empty = database.GenomeComparisonDatabase(
        profile_db=profile_1_2_3_database,
        config=simple_genome_compare_config,
        comp_db_loc=empty_db_dir
    )
    assert db_empty.get_remaining_pairs().collect().equals(pl.DataFrame({"profile_1":["profile_1","profile_1","profile_2"],"profile_2":["profile_2","profile_3","profile_3"]},schema={"profile_1":pl.Utf8,"profile_2":pl.Utf8}))
    remaining_pairs = db.get_remaining_pairs()
    assert remaining_pairs.collect().equals(pl.DataFrame({"profile_1":["profile_1","profile_2"],"profile_2":["profile_3","profile_3"]},schema={"profile_1":pl.Utf8,"profile_2":pl.Utf8}))
    assert not db.is_complete()
    db.profile_db=profile_1_2_database
    remaining_pairs = db.get_remaining_pairs()
    assert remaining_pairs.collect().is_empty()
    assert db.is_complete()


def test_genome_comparison_database_add_comp_database_valid(comps_lf_1_2,comps_lf_31_32,profile_1_2_3_database,simple_genome_compare_config,tmp_path):
    comp_db_1_2_dir = tmp_path / "comp_db_1_2.parquet"
    comps_lf_1_2.sink_parquet(comp_db_1_2_dir)
    comp_db_31_32_dir = tmp_path / "comp_db_31_32.parquet"
    comps_lf_31_32.sink_parquet(comp_db_31_32_dir)
    db_12 = database.GenomeComparisonDatabase(
        profile_db=profile_1_2_3_database,
        config=simple_genome_compare_config,
        comp_db_loc=comp_db_1_2_dir
    )
    db_3_12 = database.GenomeComparisonDatabase(
        profile_db=profile_1_2_3_database,
        config=simple_genome_compare_config,
        comp_db_loc=comp_db_31_32_dir
    )
    db_12._validate_db()
    assert db_12.comp_db.collect().height == 2
    assert not db_12.is_complete()
    assert db_3_12.comp_db.collect().height == 4
    assert not db_3_12.is_complete()
    db_12.add_comp_database(db_3_12)
    assert db_12.comp_db.collect().height == 6
    assert db_12.is_complete()
    
def test_genome_comparison_database_add_comp_database_invalid(profile_1_2_3_database,profile_4_5_database,comps_lf_1_2,comps_lf_31_32,comps_lf_4_5,tmp_path):
    config_1_2_3=database.GenomeComparisonConfig(
        reference_id="ref_1",
        gene_db_id="gene_ref_1",
        min_cov=5, 
        scope="all",
        min_gene_compare_len=100,
        stb_file_loc="Somefile",
        null_model_loc="Somefile"
    )
    config_4_5=database.GenomeComparisonConfig(
        reference_id="ref_2",
        gene_db_id="gene_ref_2",
        min_cov=5,
        scope="all",
        min_gene_compare_len=100,
        stb_file_loc="Somefile",
        null_model_loc="Somefile"
    )
    comp_db_1_2_dir = tmp_path / "comp_db_1_2.parquet"
    comps_lf_1_2.sink_parquet(comp_db_1_2_dir)
    comp_db_31_32_dir = tmp_path / "comp_db_31_32.parquet"
    comps_lf_31_32.sink_parquet(comp_db_31_32_dir)
    comp_db_4_5_dir = tmp_path / "comp_db_4_5.parquet"
    comps_lf_4_5.sink_parquet(comp_db_4_5_dir)
    db_12 = database.GenomeComparisonDatabase(
        profile_db=profile_1_2_3_database,
        config=config_1_2_3,
        comp_db_loc=comp_db_1_2_dir
    )
    db_3_12 = database.GenomeComparisonDatabase(
        profile_db=profile_1_2_3_database,
        config=config_1_2_3,
        comp_db_loc=comp_db_31_32_dir
    )
    db_4_5 = database.GenomeComparisonDatabase(
        profile_db=profile_4_5_database,
        config=config_4_5,
        comp_db_loc=comp_db_4_5_dir
    )
    db_12._validate_db()
    assert db_12.comp_db.collect().height == 2
    assert not db_12.is_complete()
    assert db_3_12.comp_db.collect().height == 4
    assert not db_3_12.is_complete()
    db_3_12.add_comp_database(db_12)
    assert db_3_12.comp_db.collect().height == 6
    assert db_3_12.is_complete()
    db_4_5._validate_db()
    with pytest.raises(ValueError, match="The comparison database provided is not compatible with the current comparison database."):
        db_3_12.add_comp_database(db_4_5)

def test_genome_comparison_database_save_comps_database(profile_1_2_3_database,comps_lf_1_2,comps_lf_31_32,tmp_path,simple_genome_compare_config):
    profile_1_2_3_database_dir = tmp_path / "profile_1_2_3.parquet"
    profile_1_2_3_database.save_as_new_database(profile_1_2_3_database_dir)
    comp_db_1_2_dir = tmp_path / "comp_db_1_2.parquet"
    comps_lf_1_2.sink_parquet(comp_db_1_2_dir)
    comp_db_31_32_dir = tmp_path / "comp_db_31_32.parquet"
    comps_lf_31_32.sink_parquet(comp_db_31_32_dir)
    db_12 = database.GenomeComparisonDatabase(
        profile_db=database.ProfileDatabase(db_loc=profile_1_2_3_database_dir),
        config=simple_genome_compare_config,
        comp_db_loc=comp_db_1_2_dir
    )
    db_3_12 = database.GenomeComparisonDatabase(
        profile_db=database.ProfileDatabase(db_loc=profile_1_2_3_database_dir),
        config=simple_genome_compare_config,
        comp_db_loc=comp_db_31_32_dir
    )
    with pytest.raises(ValueError, match="The output path must be different"):
        db_12.save_new_compare_database(comp_db_1_2_dir)

    db_12_dir = tmp_path / "db_12_new_db.parquet"
    
    db_12.save_new_compare_database(db_12_dir)
    
    loaded_db = database.GenomeComparisonDatabase(
        profile_db=database.ProfileDatabase(db_loc=profile_1_2_3_database_dir),
        config=simple_genome_compare_config,
        comp_db_loc=db_12_dir
    )
    assert loaded_db.comp_db.collect().equals(db_12.comp_db.collect())
    db_12.add_comp_database(db_3_12)
    assert db_12.comp_db.collect().height == 6
    db_12.update_compare_database()
    loaded_db_updated = database.GenomeComparisonDatabase(
        profile_db=database.ProfileDatabase(db_loc=profile_1_2_3_database_dir),
        config=simple_genome_compare_config,
        comp_db_loc=comp_db_1_2_dir
    )
    assert loaded_db_updated.comp_db.collect().height == 6
    assert loaded_db_updated.is_complete()
    
def test_genome_comparison_database_load_and_dump(profile_1_2_3_database,comps_lf_1_2,comps_lf_31_32,tmp_path,simple_genome_compare_config):
    profile_1_2_3_database_dir = tmp_path / "profile_1_2_3.parquet"
    profile_1_2_3_database.save_as_new_database(profile_1_2_3_database_dir)
    comp_db_1_2_dir = tmp_path / "comp_db_1_2.parquet"
    comps_lf_1_2.sink_parquet(comp_db_1_2_dir)
    comp_db_31_32_dir = tmp_path / "comp_db_31_32.parquet"
    comps_lf_31_32.sink_parquet(comp_db_31_32_dir)
    db_12 = database.GenomeComparisonDatabase(
        profile_db=database.ProfileDatabase(db_loc=profile_1_2_3_database_dir),
        config=simple_genome_compare_config,
        comp_db_loc=comp_db_1_2_dir
    )
    db_3_12 = database.GenomeComparisonDatabase(
        profile_db=database.ProfileDatabase(db_loc=profile_1_2_3_database_dir),
        config=simple_genome_compare_config,
        comp_db_loc=comp_db_31_32_dir
    )
    db_12.add_comp_database(db_3_12)
    db_12.update_compare_database()
    assert db_12.comp_db.collect().height == 6
    dump_path = tmp_path / "db_dump"
    db_12.dump_obj(dump_path)
    with open(dump_path) as f:
        serialized=json.load(f)
    assert serialized["comp_db_loc"] == str(comp_db_1_2_dir)
    assert serialized["config"] == simple_genome_compare_config.to_dict()
    assert serialized["profile_db_loc"] == str(profile_1_2_3_database_dir)
    db_obj=database.GenomeComparisonDatabase.load_obj(dump_path)
    assert db_obj.comp_db.collect().height == 6

def test_genome_comparison_database_to_complete_input_table(profile_1_2_3_database,simple_genome_compare_config):
    db = database.GenomeComparisonDatabase(
        profile_db=profile_1_2_3_database,
        config=simple_genome_compare_config,
        comp_db_loc=None
    )
    assert db.to_complete_input_table().collect().height == 3
    assert set(db.to_complete_input_table().collect_schema().keys()) == {"sample_name_1","profile_location_1","scaffold_location_1","sample_name_2","profile_location_2","scaffold_location_2"}
    assert db.to_complete_input_table().collect().filter(pl.col("sample_name_1")=="profile_1").height == 2
    assert db.to_complete_input_table().collect().filter(pl.col("sample_name_1")=="profile_2").height == 1
    assert db.to_complete_input_table().collect().filter(pl.col("sample_name_1")=="profile_3").height == 0
    assert db.to_complete_input_table().collect().filter(pl.col("sample_name_2")=="profile_1").height == 0
    assert db.to_complete_input_table().collect().filter(pl.col("sample_name_2")=="profile_2").height == 1
    assert db.to_complete_input_table().collect().filter(pl.col("sample_name_2")=="profile_3").height == 2



