import polars as pl
import pathlib


def get_genome_accessions_from_sylph_output(sylph_output:pl.LazyFrame) -> pl.LazyFrame:
    pass

class LocalGenomeDB:
    LOCAL_GENOME_DB_NAME=".genome_db"
    def __init__(self, db_path:str):
        self.db_path = db_path
        if not pathlib.Path(self.db_path).exists():
            self._create_db()
        else:
            self._db=pl.read_parquet(self.db_path)
        self._is_synced=False
        self.to_be_fetched=None
    
    def _create_db(self):
        pl.DataFrame({
            "accession":[],
            "genome_name":[],
            "location":[],
            "exist":[],
        }).write_parquet(self.db_path)
        self._db=pl.read_parquet(self.db_path)
    
    def add_genome(self, accession:str, genome_name:str):
        lf = pl.read_parquet(self.db_path)
        new_entry = pl.DataFrame({
            "accession":[accession],
            "genome_name":[genome_name],
            "location":[None],
            "exist":[False],
        })
        lf = pl.concat([lf, new_entry])
        lf.write_parquet(self.db_path)
        self._db=pl.read_parquet(self.db_path)
        self._is_synced=False
    
    def sync(self):
        self._db=self._db.with_columns(
            pl.col("location").map_elements(lambda x: pathlib.Path(x).exists() if x is not None else False).alias("exist")
        )
        self.to_be_fetched = set(self._db.filter(~pl.col("exist"))["accession"].to_list())
        self._is_synced = True
    
    
    
    
    