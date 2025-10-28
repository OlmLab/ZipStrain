"""zipstrain.database
========================
This module provides classes and functions to manage profile and comparison databases for efficient data handling.
The ProfileDatabase class manages profiles, while the GenomeComparisonDatabase class handles comparisons between profiles.
See the documentation of each class for more details.

"""
from __future__ import annotations
import polars as pl
import pathlib 
import os
import tempfile
import json
import copy
from pydantic import BaseModel, Field, field_validator,ConfigDict
class ProfileItem(BaseModel):
    """
    This class describes all necessary attributes of a profile and makes sure they comply with the necessary formating.
    """
    model_config = ConfigDict(extra="forbid")
    profile_name: str = Field(description="An arbitrary name given to the profile (Usually sample name or name of the parquet file)")
    profile_location: str = Field(description="The location of the profile")
    scaffold_location: str = Field(description="The location of the scaffold")
    reference_db_id: str = Field(description="The ID of the reference database. This could be the name or any other identifier for the database that the reads are mapped to.")
    gene_db_id:str= Field(default="",description="The ID of the gene database in fasta format. This could be the name or any other identifier for the database that the reads are mapped to.")
    
    @field_validator("profile_location","scaffold_location")
    def check_file_exists(cls, v):
        if not os.path.exists(v):
            raise ValueError(f"The file {v} does not exist.")
        return v

    @field_validator("reference_db_id","gene_db_id")
    def check_reference_db_id(cls, v):
        if not v:
            raise ValueError("The reference_db_id and gene_db_id cannot be empty.")
        return v
    
    

class ProfileDatabase:
    """
    The profile database simply holds profile information. Does not need to be specific to a comparison database.
    The data behind a profile is stored in a parquet file. It is basically a table with the following columns:
    - profile_name: An arbitrary name given to the profile (Usually sample name or name of the parquet file)
    - profile_location: The location of the profile
    - scaffold_location: The location of the scaffold
    - reference_db_id: The ID of the reference database. This could be the name or any other identifier for the database that the reads are mapped to.
    - gene_db_id: The ID of the gene database in fasta format. This could be the name or any other identifier for the database that the reads are mapped to.
    """
    def __init__(self,
                 db_loc: str|None = None,
                 ):
        if db_loc is not None:
            self.db_loc = pathlib.Path(db_loc)
            self._db = pl.scan_parquet(self.db_loc)
        else:
            self._db=pl.LazyFrame({
                "profile_name": [],
                "profile_location": [],
                "scaffold_location": [],
                "reference_db_id": [],
                "gene_db_id": []
            }, schema={
                "profile_name": pl.Utf8,
                "profile_location": pl.Utf8,
                "scaffold_location": pl.Utf8,
                "reference_db_id": pl.Utf8,
                "gene_db_id": pl.Utf8
            })
            self.db_loc=None
            
    @property
    def db(self):
        return self._db

    def _validate_db(self,check_profile_exists: bool=True,check_scaffold_exists:bool=True)->None:
        """Simple method to see if the database has the minimum required structure."""

        ### Next check if the database has the required columns
        required_columns = ["profile_name","profile_location", "scaffold_location", "reference_db_id", "gene_db_id"]
        for col in required_columns:
            if col not in self.db.collect_schema().names():
                raise ValueError(f"Missing required column: {col}")

        if check_profile_exists:
            # Check if the profile exists in the database
            db_path_validated= self.db.select(pl.col("profile_location")).collect(engine="streaming").with_columns(
                (pl.col("profile_location").map_elements(lambda x: pathlib.Path(x).exists(),return_dtype=pl.Boolean)).alias("profile_exists")
            ).filter(~ pl.col("profile_exists"))
            if db_path_validated.height != 0:
                raise ValueError(f"There are {db_path_validated.height} profiles that do not exist: {db_path_validated['profile_location'].to_list()}")
            ### add log later
        if check_scaffold_exists:
            db_path_validated= self.db.select(pl.col("scaffold_location")).collect(engine="streaming").with_columns(
                (pl.col("scaffold_location").map_elements(lambda x: pathlib.Path(x).exists(),return_dtype=pl.Boolean)).alias("scaffold_exists")
            ).filter(~ pl.col("scaffold_exists"))
            if db_path_validated.height != 0:
                raise ValueError(f"There are {db_path_validated.height} scaffolds that do not exist: {db_path_validated['scaffold_location'].to_list()}")
            ### add log later
        
    def add_profile(self,
                    data: dict
                    ) -> None:
        """Add a profile to the database.
        Args:
            data (dict): The profile data to add.
        """
        try:
            profile_item = ProfileItem(**data)
            lf=pl.LazyFrame({
                "profile_name": [profile_item.profile_name],
                "profile_location": [profile_item.profile_location],
                "scaffold_location": [profile_item.scaffold_location],
                "reference_db_id": [profile_item.reference_db_id],
                "gene_db_id": [profile_item.gene_db_id]
            })
            self._db = pl.concat([self.db, lf]).unique()
            self._validate_db()
        except Exception as e:
            raise ValueError(f"The profile data provided is not valid: {e}")
        

    def add_database(self, profile_database: ProfileDatabase) -> None:
        """Merge the provided profile database into the current database.
        Args:
            profile_database (ProfileDatabase): The profile database to merge.
        """
        try:
            profile_database._validate_db()
        
        except Exception as e:
            raise ValueError(f"The profile database provided is not valid: {e}")

        self._db = pl.concat([self._db, profile_database.db]).unique()
    

    def save_as_new_database(self, output_path: str) -> None:
        """Save the database to a parquet file.
        Args:
            database (pl.LazyFrame): The database to save.
            output_path (str): The path to save the database to.
        """
        #The new database must be written to a new location
        if self.db_loc is not None and str(self.db_loc.absolute()) == str(pathlib.Path(output_path).absolute()):
            raise ValueError("The output path must be different from the current database location.")
        
        try:
            self.db.sink_parquet(output_path)
            self.db_loc=pathlib.Path(output_path)
        ### add log later
        except Exception as e:
            pass 
        ### add log later
    
    def update_database(self)->None:
        """Overwrites the database saved on the disk to the current database object
        """
        if self.db_loc is None:
            raise Exception("db_loc attribute is not determined yet!")
        try:
            self.db.collect(engine="streaming").write_parquet(self.db_loc)
        except Exception as e:
            raise Exception(f"Something went wrong when updating the database:{e}")
             
    
    @classmethod
    def from_csv(cls, csv_path: str) -> ProfileDatabase:
        """Create a ProfileDatabase instance from a CSV file with exactly same columns as the required columns for a profile database.
        Args:
            csv_path (str): The path to the CSV file.
        Returns:
            ProfileDatabase: The created ProfileDatabase instance.
        """
        lf=pl.scan_csv(csv_path).collect().lazy() # To avoid clash when using to_csv on same file
        prof_db=cls()
        prof_db._db=lf
        prof_db._validate_db()
        return prof_db
    
    def to_csv(self,output_dir:str)->None:
        """Writes the the current database object to a csv file"
        Args:
            output_dir (str): The path to save the CSV file.
        Returns:
            None
        """
        self.db.sink_csv(output_dir,engine="streaming")
        

    
    
     
class GenomeComparisonConfig(BaseModel):
    """
    This class defines object which have all necessary options to describe 
    Parameters used to compare profiles:
    - gene_db_id (str): The ID of the gene fasta database to use for the comparison. The file name is perfect.
    - reference_id (str): The ID of the reference fasta database to use for the comparison. The file name is perfect.
    - scope (str): The scope of the comparison- 'all' if all covered positions are desired. Otherwise, a bunch of genome names separated by commas.
    - min_cov (int): Minimum coverage a base on the reference fasta that must have in order to be compared.
    - null_model_p_value(float): P_value above which a base call is counted as sequencing error
    - min_gene_compare_len (int): Minimum length of a gene that needs to be covered at min_cov to be considered for gene similarity calculations
    - stb_file_loc (str): The location of the scaffold to bin file.
    - null_model_loc (str): The location of the null model file.
    """
    model_config = ConfigDict(extra="forbid")
    gene_db_id:str= Field(default="",description="An ID given to the gene fasta file used for profiling. IMPORTANT: Make sure that this is in agreement with gene database IDs in the Profile Database.")
    reference_id:str= Field(description="An ID given to the reference fasta file used for profiling. IMPORTANT: Make sure that this is in agreement with reference IDs in the Profile Database.")
    scope: str =Field(description="An ID given to the reference fasta file used for profiling. IMPORTANT: Make sure that this is in agreement with reference IDs in the Profile Database.")
    min_cov: int =Field(description="Minimum coverage a base on the reference fasta that must have in order to be compared.")
    min_gene_compare_len: int=Field(description="Minimum length of a gene that needs to be covered at min_cov to be considered for gene similarity calculations")
    null_model_p_value:float=Field(default=0.05,description="P_value above which a base call is counted as sequencing error")
    stb_file_loc:str=Field(description="The location of the scaffold to bin file.")
    null_model_loc:str=Field(description="The location of the null model file.")

    def is_compatible(self, other: GenomeComparisonConfig) -> bool:
        """
        Check if this comparison configuration is compatible with another. Two configurations are compatible if they have the same parameters, except for scope.
        Scope can be different as long as they are not disjoint. Also, all is compatible with any scope.
        Args:
            other (GenomeComparisonConfig): The other comparison configuration to check compatibility with.
        Returns:
            bool: True if the configurations are compatible, False otherwise.
        """
        attrs=self.__dict__
        for key in attrs:
            if key!="scope":
                if attrs[key] != getattr(other, key):
                    return False
        if other.scope != "all" and self.scope != "all":
            if (set(other.scope.split(",")).intersection(set(self.scope.split(","))) == set()):
                return False
        return True

    @classmethod
    def from_json(cls,json_file_dir:str)->GenomeComparisonConfig:
        """Create a GenomeComparisonConfig instance from a json file."""
        with open(json_file_dir, 'r') as f:
            config_dict = json.load(f)
        return cls(**config_dict)
    
    def to_json(self,json_file_dir:str)->None:
        """Writes the the current object to a json file"""
        with open(json_file_dir,"w") as f:
            json.dump(self.__dict__,f)
            
    def to_dict(self)->dict:
        """Returns the dictionary representation of the current object"""
        return copy.copy(self.__dict__)
    
    
    def get_maximal_scope_config(self, other: GenomeComparisonConfig) -> GenomeComparisonConfig:
        """
        Get a new GenomeComparisonConfig object with the maximal scope that is compatible with the two configurations.
        Args:
            other (GenomeComparisonConfig): The other comparison configuration to get the maximal scope with.
        Returns:
            GenomeComparisonConfig: The new comparison configuration with the maximal scope.
        """
        if not self.is_compatible(other):
            raise ValueError("The two comparison configurations are not compatible.")
        
        new_scope=None
        if other.scope == "all" and self.scope == "all":
            new_scope="all"
        
        elif other.scope == "all":
            new_scope=self.scope.split(",")
        
        elif self.scope == "all":
            new_scope=other.scope.split(",")
        
        else:
            new_scope=list(set(self.scope.split(",")).intersection(set(other.scope.split(","))))
        curr_config_dict=self.to_dict()
        curr_config_dict["scope"]=new_scope if new_scope=="all" else ",".join(sorted(new_scope))
        return GenomeComparisonConfig(**curr_config_dict)


class GenomeComparisonDatabase:
    """
    GenomeComparisonDatabase object holds a reference to a comparison parquet file. The methods in this class serve to provide
    functionality for working with the comparison data in an easy and efficient manner.
    The comparison parquet file the result of running compare, and optionally concatenating multiple compare parquet file from single comparisons.
    This parquet file must contain the following columns:
    - genome
    - total_positions
    - share_allele_pos
    - genome_pop_ani
    - max_consecutive_length
    - shared_genes_count
    - identical_gene_count
    - sample_1
    - sample_2
    A ComparisonDatabase object needs a ComparisonConfig object to specify the parameters used for the comparison.
    """
    COLUMN_NAMES = [
        "genome",
        "total_positions",
        "share_allele_pos",
        "genome_pop_ani",
        "max_consecutive_length",
        "shared_genes_count",
        "identical_gene_count",
        "perc_id_genes",
        "sample_1",
        "sample_2"
    ]

    def __init__(self,
                 profile_db: ProfileDatabase,
                 config: GenomeComparisonConfig,
                 comp_db_loc: str|None = None,
                 ):
        self.profile_db = profile_db
        self.config = config
        if comp_db_loc is not None:
            self.comp_db_loc = pathlib.Path(comp_db_loc)
            self._comp_db = pl.scan_parquet(self.comp_db_loc)
        else:
            self.comp_db_loc = None
            self._comp_db=pl.LazyFrame({
                "genome": [],
                "total_positions": [],
                "share_allele_pos": [],
                "genome_pop_ani": [],
                "max_consecutive_length": [],
                "shared_genes_count": [],
                "identical_gene_count": [],
                "perc_id_genes": [],
                "sample_1": [],
                "sample_2": []
            }, schema={
                "genome": pl.Utf8,
                "total_positions": pl.Int64,
                "share_allele_pos": pl.Int64,
                "genome_pop_ani": pl.Float64,
                "max_consecutive_length": pl.Int64,
                "shared_genes_count": pl.Int64,
                "identical_gene_count": pl.Int64,
                "perc_id_genes": pl.Float64,
                "sample_1": pl.Utf8,
                "sample_2": pl.Utf8
            })
            self.comp_db_loc=None

    @property
    def comp_db(self):
        return self._comp_db
    
    def _validate_db(self)->None:
        self.profile_db._validate_db()
        
        if set(self._comp_db.collect_schema()) != set(self.COLUMN_NAMES):
            raise ValueError(f"Your comparison database must provide these extra columns: { set(self.COLUMN_NAMES)-set(self._comp_db.collect_schema())}")
        #check if all profile names exist in the profile database
        profile_names_in_comp_db = set(self.get_all_profile_names())
        profile_names_in_profile_db = set(self.profile_db.db.select("profile_name").collect(engine="streaming").to_series().to_list())
        if not profile_names_in_comp_db.issubset(profile_names_in_profile_db):
            raise ValueError(f"The following profile names are in the comparison database but not in the profile database: {profile_names_in_comp_db - profile_names_in_profile_db}")
    
    def get_all_profile_names(self) -> set[str]:
        """
        Get all profile names that are in the comparison database.
        """
        return set(self.comp_db.select(pl.col("sample_1")).collect(engine="streaming").to_series().to_list()).union(
            set(self.comp_db.select(pl.col("sample_2")).collect(engine="streaming").to_series().to_list())
        )
    def get_remaining_pairs(self) -> pl.LazyFrame:
        """
        Get pairs of profiles that are in the profile database but not in the comparison database.
        """
        profiles = self.profile_db.db.select("profile_name")
        pairs=profiles.join(profiles,how="cross").rename({"profile_name":"profile_1","profile_name_right":"profile_2"}).filter(pl.col("profile_1")<pl.col("profile_2"))
        samplepairs = self.comp_db.group_by("sample_1", "sample_2").agg().with_columns(pl.min_horizontal(["sample_1", "sample_2"]).alias("profile_1"), pl.max_horizontal(["sample_1", "sample_2"]).alias("profile_2")).select(["profile_1", "profile_2"])

        remaining_pairs = pairs.join(samplepairs, on=["profile_1", "profile_2"], how="anti").sort(["profile_1","profile_2"])
        return remaining_pairs

    def is_complete(self) -> bool:
        """
        Check if the comparison database is complete, i.e., if all pairs of profiles in the profile database have been compared.
        """
        return self.get_remaining_pairs().collect(engine="streaming").is_empty()

    def add_comp_database(self, comp_database: GenomeComparisonDatabase) -> None:
        """Merge the provided comparison database into the current database.
        Args:
            comp_database (ComparisonDatabase): The comparison database to merge.
        """
        try:
            comp_database._validate_db()
        
        except Exception as e:
            raise ValueError(f"The comparison database provided is not valid: {e}")

        if not self.config.is_compatible(comp_database.config):
            raise ValueError("The comparison database provided is not compatible with the current comparison database.")
        
        self._comp_db = pl.concat([self._comp_db, comp_database.comp_db]).unique()
        self.config = self.config.get_maximal_scope_config(comp_database.config)
        
        
    def save_new_compare_database(self, output_path: str) -> None:
        """Save the database to a parquet file."""
        output_path = pathlib.Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # The new database must be written to a new location
        if self.comp_db_loc is not None and str(self.comp_db_loc.absolute()) == str(output_path.absolute()):
            raise ValueError("The output path must be different from the current database location.")

        self.comp_db.sink_parquet(output_path)
            
    
    def update_compare_database(self)->None:
        """Overwrites the comparison database saved on the disk to the current comparison database object
        """
        if self.comp_db_loc is None:
            raise Exception("comp_db_loc attribute is not determined yet!")
        try:
            tmp_path=pathlib.Path(tempfile.mktemp(suffix=".parquet",prefix="tmp_comp_db_",dir=str(self.comp_db_loc.parent)))
            self.comp_db.sink_parquet(tmp_path)
            os.replace(tmp_path,self.comp_db_loc)
            self._comp_db=pl.scan_parquet(self.comp_db_loc)
        except Exception as e:
            raise Exception(f"Something went wrong when updating the comparison database:{e}")
    
    def dump_obj(self, output_path: str) -> None:
        """Dump the current object to a json file.
        Args:
            output_path (str): The path to save the json file to.
        """
        obj_dict = {
            "profile_db_loc": str(self.profile_db.db_loc.absolute()) if self.profile_db.db_loc is not None else None,
            "config": self.config.to_dict(),
            "comp_db_loc": str(self.comp_db_loc.absolute()) if self.comp_db_loc is not None else None
        }
        with open(output_path, "w") as f:
            json.dump(obj_dict, f, indent=4)
    
    @classmethod
    def load_obj(cls, json_path: str) -> GenomeComparisonDatabase:
        """Load a GenomeComparisonDatabase object from a json file.
        Args:
            json_path (str): The path to the json file.
        Returns:
            GenomeComparisonDatabase: The loaded GenomeComparisonDatabase object.
        """
        with open(json_path, "r") as f:
            obj_dict = json.load(f)
        
        return cls(profile_db=ProfileDatabase(db_loc=obj_dict["profile_db_loc"]) , 
                   config=GenomeComparisonConfig(**obj_dict["config"]), 
                   comp_db_loc=obj_dict["comp_db_loc"])
    
    
    def to_complete_input_table(self)->pl.LazyFrame:
        """This method gives a table of all pairwise comparisons that is needed to make the comparison database complete. The table contains the following columns:
        - sample_name_1
        - sample_name_2
        - profile_location_1
        - scaffold_location_1
        - profile_location_2
        - scaffold_location_2
        Returns:
            pl.LazyFrame: The table of all pairwise comparisons needed to complete the comparison database.
        """
        lf=self.get_remaining_pairs().rename({"profile_1":"sample_name_1","profile_2":"sample_name_2"})
        return (lf.join(self.profile_db.db.select(["profile_name","profile_location","scaffold_location"]),left_on="sample_name_1",right_on="profile_name",how="left")
                .rename({"profile_location":"profile_location_1","scaffold_location":"scaffold_location_1"})
                .join(self.profile_db.db.select(["profile_name","profile_location","scaffold_location"]),left_on="sample_name_2",right_on="profile_name",how="left")
                .rename({"profile_location":"profile_location_2","scaffold_location":"scaffold_location_2"})
               )

