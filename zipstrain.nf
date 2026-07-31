params.error_rate=0.001
params.max_total_reads=50000
params.p_threshold=0.05
params.mode = null
params.parallel_mode="batched"
params.min_cov=5
params.min_gene_compare_len=200
params.batch_size=10
params.compare_duckdb_memory_limit=""
params.compare_engine="polars"
params.compare_calculate="all"
// Gene boundaries are resolved from ranges at compare time, so gene comparison
// requires this table; it is no longer carried inside the profile.
params.gene_range_table=null
params.batch_compare_n_parallel=4
params.publish_mode="symlink"
params.compare_scope="all"
params.compare_ani_method="popani"
params.input_type="profile_table"
params.bowtie2_non_competitive_mapping=false
params.min_mapq=0
params.min_baseq=13
params.min_freq=0.01
params.min_read_ani=0.95
params.read_inclusion="paired"
params.sylph_db = null
params.sylph_db_link="http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r220-c200-dbv1.syldb"
params.genome_db_cache_dir="genome_cache"
params.prefetch_max_size="200g"
params.compression_level=6
def tableToDict(file, delimiter = ',') {
    /*
    * This function reads a CSV file and converts it into a dictionary (map) where the keys are the headers
    * and the values are lists of the corresponding column values.
    * @param file: The CSV file to read.
    * @param delimiter: The delimiter used in the CSV file (default is comma).
    * @return: A map where keys are headers and values are lists of column values.
    */



    def result = [:]
    def lines = file.readLines()
    
    // Get all headers
    def headers = lines[0].split(delimiter).collect { it.trim() }
    
    // Initialize empty arrays for each header
    headers.each { header ->
        result[header] = []
    }
    
    // Collect all values for each column
    lines[1..-1].each { line ->
        def values = line.split(delimiter).collect { it.trim() }
        headers.eachWithIndex { header, i ->
            result[header] << values[i]
        }
    }
    
    return result
}

def getProfilesTableColumn(input_table) {
    if (input_table.containsKey('profiles')) {
        return input_table['profiles']
    }
    throw new IllegalArgumentException("Profile-table input must include a 'profiles' column.")
}

def getProfileSampleNamesTableColumn(input_table) {
    if (input_table.containsKey('sample_name')) {
        return input_table['sample_name']
    }
    throw new IllegalArgumentException("Profile-table input must include a 'sample_name' column.")
}

def makeBatches(items, batchSize) {
    def size = batchSize.toString().toInteger()
    if (size < 1) {
        throw new IllegalArgumentException("--batch_size must be at least 1.")
    }

    def batches = []
    def current = []
    items.each { item ->
        if (current.size() == size) {
            batches << current
            current = []
        }
        current << item
    }
    if (current) {
        batches << current
    }
    return batches
}

process download_sylph_db{
    output:
    path "*.syldb", emit: sylph_db
    script:
    """
    wget ${params.sylph_db_link}
    """
}
process estimate_abundance_sylph{
    /*
    * This process estimates the abundance of bins using Sylph. It takes the reads and prepared database.
    * One task per sample.
    */
    publishDir "${params.output_dir}/sylph_abundance/", mode: params.publish_mode
    
    input:
    tuple val(sample_name), path(reads)
    path sylph_db

    output:
    path "${sample_name}_sylph_abundance.tsv", emit: abundance
    script:
    """
    num_reads=\$(ls ${reads} | wc -l)
    if [ \$num_reads -eq 2 ]; then
        r1=\$(ls ${reads} | head -n 1)
        r2=\$(ls ${reads} | tail -n 1)
        sylph profile ${sylph_db} -1 \$r1 -2 \$r2 -t ${task.cpus} > ${sample_name}_sylph_abundance.tsv
    else
        r1=\$(ls ${reads} | head -n 1)
        sylph profile ${sylph_db} -U \$r1 -t ${task.cpus} > ${sample_name}_sylph_abundance.tsv
    fi
    """
}

process merge_sylph_abundance_tables {
    /*
    * Merge per-sample Sylph abundance tables into one table.
    */
    publishDir "${params.output_dir}/sylph_abundance/", mode: params.publish_mode
    input:
    path abundance_tables
    output:
    path "sylph_abundance.tsv", emit: abundance
    script:
    """
    awk 'FNR==1 && NR!=1 {next} 1' ${abundance_tables} > sylph_abundance.tsv
    """
}

process build_db_from_Sylph{
    /*
    * This process builds a Bowtie2 database from the Sylph database. It takes the Sylph database as input and outputs the indexed files.
    */
    publishDir "${params.output_dir}/db_from_sylph/", mode: params.publish_mode
    containerOptions {
        if (workflow.containerEngine == 'docker') {
            return "--volume ${genome_cache_dir}:${genome_cache_dir}"
        } else if (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') {
            return "--bind ${genome_cache_dir}:${genome_cache_dir}"
        } else {
            return ''
        }
    }
    input:
    path sylph_abundance
    val genome_cache_dir
    output:
    path "reference_genomes.fna", emit:reference_genome
    path "reference_genomes.stb",emit:stb
    path "reference_genomes_gene.fasta",emit:reference_genome_genes
    script:
    """
    zipstrain utilities build-genome-db \
        --tool sylph \
        --abundance-table ${sylph_abundance} \
        --cache-dir ${genome_cache_dir}  \
        --output-dir .
    
    prodigal -i reference_genomes.fna -d reference_genomes_gene.fasta  -p meta

    """
}

process get_sequences_from_sra {
    /*
    * This process retrieves sequences from the SRA database using the fastq-dump tool.
    * It takes in a list of SRA IDs and outputs the corresponding FASTQ files.
    * @param sra_ids: SRA ID to retrieve.
    */
    
    publishDir "${params.output_dir}/sra_sequences", mode: params.publish_mode
    
    input:
    val sra_ids
    
    output:
    path "${sra_ids}/${sra_ids}*.fastq.gz",emit: fastq_files
    val sra_ids, emit: sra_ids
    tuple val(sra_ids), path("${sra_ids}/${sra_ids}*.fastq.gz"), emit: sample_reads
    
    script:
    """
    prefetch --max-size ${params.prefetch_max_size} ${sra_ids} 
    fasterq-dump --split-files --outdir ${sra_ids} ${sra_ids}
    pigz -p ${task.cpus} -${params.compression_level} ${sra_ids}/${sra_ids}*.fastq
    rm -rf ${sra_ids}/${sra_ids}.sra
    rm -f ${sra_ids}/${sra_ids}*.fastq
    """
}
process index_reference {
    /*
    * This process indexes the reference genome using BWA.
    * It takes in the reference genome and outputs the indexed files.
    */
    publishDir "${params.output_dir}", mode: params.publish_mode
    input:
    path reference_genome
    output:
    path "${reference_genome.name}.*.bt2*", emit: index_files
    path reference_genome, emit: reference_genome
    script:
    """
    bowtie2-build --threads ${task.cpus} ${reference_genome} ${reference_genome}

    """
}

process map_reads{
    /*
    * This process maps reads to a reference genome using BWA.
    * It takes in the reference genome and reads, and outputs the BAM file.
    */
    publishDir "${params.output_dir}", mode: params.publish_mode
    input:
    tuple val(sample_name), path(reads)
    path reference_genome
    path index_files
    output:
    path "${sample_name}.bam", emit: bamfile
    script:
    def competitiveness= (params.bowtie2_non_competitive_mapping) ? "-a" : ""
    if (reads.size() == 2) {
    """
    bowtie2 \\
            -x ${reference_genome} \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            ${competitiveness} \\
            --threads ${task.cpus} \\
            | samtools view -bS -F 4 - \\
            | samtools sort -@ ${task.cpus} -o ${sample_name}.bam -
    """}
    else {
        """
        bowtie2 \\
                -x ${reference_genome} \\
                -U ${reads[0]} \\
                ${competitiveness} \\
                --threads ${task.cpus} \\
                | samtools view -bS -F 4 - \\
                | samtools sort -@ ${task.cpus} -o ${sample_name}.bam -
        """
    }

}

process prepare_profile{

    publishDir "${params.output_dir}", mode: params.publish_mode
    input:
    path reference_genome
    path gene_fasta
    path stb_file
    output:
    path "genomes_bed_file.bed", emit: genome_bed
    path "gene_range_table.tsv", emit: gene_range_table
    path "null_model.parquet", emit: null_model
    path "profiling_contract.json", emit: profiling_contract
    script:
"""
zipstrain utilities prepare_profiling \\
        --reference-fasta ${reference_genome} \\
        --gene-fasta ${gene_fasta} \\
        --stb-file ${stb_file} \\
        --error-rate ${params.error_rate} \\
        --max-total-reads ${params.max_total_reads} \\
        --p-threshold ${params.p_threshold} \\
        --output-dir . 
"""

}


process profile_bam {
    /*
    * This process profiles BAM files into ZipStrain parquet profiles.
    * It takes in the BAM file and outputs profile/stat parquet files.
    */
    publishDir "${params.output_dir}", mode: params.publish_mode
    input:
    val sample_name
    path bamfile
    path reference_fasta
    path bed_file
    path stb_file
    path gene_range_table
    path null_model
    path profiling_contract
    output:
    path "${bamfile.baseName}_profile.parquet", emit: profile
    path "${bamfile.baseName}_genome_stats.parquet", emit: genome_stats
    path "${bamfile.baseName}_gene_stats.parquet", emit: gene_stats
    val sample_name, emit: sample_name
    script:
    """
    zipstrain utilities profile-single \\
                        --bam-file ${bamfile} \\
                        --reference-fasta ${reference_fasta} \\
                        --bed-file ${bed_file} \\
                        --gene-range-table ${gene_range_table} \\
                        --stb-file ${stb_file} \\
                        --null-model ${null_model} \\
                        --profiling-contract ${profiling_contract} \\
                        --max-concurrency ${task.cpus} \\
                        --min-mapq ${params.min_mapq} \\
                        --min-baseq ${params.min_baseq} \\
                        --min-freq ${params.min_freq} \\
                        --min-read-ani ${params.min_read_ani} \\
                        --read-inclusion ${params.read_inclusion} \\
                        --output-dir .
    """
}

process compare_fast_profiles_single {
    /*
    * This process compares fast profiles.
    * It takes in profile parquets and outputs the comparison results.
    */
    input:
    path profile_location_1
    path profile_location_2
    path stb
    path gene_range_table
    val pair_name
    output:
    path "${pair_name}_comparison.parquet", emit: comparison_results
    script:
    def add_scope = (params.compare_scope=="all") ? "" : "--scope ${params.compare_scope}"
    def add_duckdb_memory_limit = (params.compare_engine == "duckdb" && params.compare_duckdb_memory_limit) ? "--duckdb-memory-limit ${params.compare_duckdb_memory_limit}" : ""
    def add_calculate = params.compare_calculate ? "--calculate ${params.compare_calculate}" : "--calculate all"
    def add_gene_range = gene_range_table ? "--gene-range-table ${gene_range_table}" : ""
    def add_compare_engine = params.compare_engine ? "--engine ${params.compare_engine}" : "--engine polars"
    """
    zipstrain utilities single-compare  \
                        --profile-location-1 ${profile_location_1} \
                        --profile-location-2 ${profile_location_2} \
                        -s ${stb} \
                        -c ${params.min_cov} \
                        -l ${params.min_gene_compare_len} \
                        ${add_compare_engine} \
                        ${add_calculate} \
                        ${add_gene_range} \
                        --ani-method "${params.compare_ani_method}" \
                        ${add_duckdb_memory_limit} \
                        ${add_scope} \
                        -o ${pair_name}_comparison.parquet

    """
}
process compare_batched {
    publishDir "${params.output_dir}/batch_comparisons", mode: params.publish_mode
    afterScript """
    rm -rf comps pairs.txt
    rm -f ${profile_locations.collect{t->t.join(' ')}}
    rm -f ${stb}
    """
    input:
    path profile_locations
    val pairs
    path stb
    path gene_range_table


    output:
    path "Batch_*_comparisons.parquet", emit: comparison_results

    script:
    pairs_text = pairs.collect{p-> p.join('\t')}.join('\n')
    remove_profile_locations = profile_locations.join(' ')
    def add_scope = (params.compare_scope=="all") ? "" : "--scope ${params.compare_scope}"
    def add_duckdb_memory_limit = (params.compare_engine == "duckdb" && params.compare_duckdb_memory_limit) ? "--duckdb-memory-limit ${params.compare_duckdb_memory_limit}" : ""
    def add_calculate = params.compare_calculate ? "--calculate ${params.compare_calculate}" : "--calculate all"
    def add_gene_range = gene_range_table ? "--gene-range-table ${gene_range_table}" : ""
    def add_compare_engine = params.compare_engine ? "--engine ${params.compare_engine}" : "--engine polars"
    """
    echo -e "${pairs_text}" > pairs.txt
    cat pairs.txt | parallel --tmpdir . --colsep '\\t' -j ${params.batch_compare_n_parallel} 'zipstrain utilities single-compare \
                        --profile-location-1 {1} \
                        --profile-location-2 {2} \
                        -s ${stb} \
                        -c ${params.min_cov} \
                        -l ${params.min_gene_compare_len} \
                        ${add_compare_engine} \
                        ${add_calculate} \
                        ${add_gene_range} \
                        --ani-method "${params.compare_ani_method}" \
                        ${add_duckdb_memory_limit} \
                        ${add_scope} \
                        -o {1}_{2}_comparison.parquet'
    mkdir comps
    hash=\$(sha1sum pairs.txt | awk '{print \$1}')
    mv *_comparison.parquet comps/
    zipstrain utilities merge_parquet  --input-dir comps --output-file "Batch_\${hash}_comparisons.parquet"
    rm -rf comps
    rm -f pairs.txt
    rm -f ${remove_profile_locations}

    """


}

process merge_comparison_tables {
    /*
    * This process merges multiple comparison result files into a single file.
    * It takes in multiple comparison result files and outputs a merged file.
    */
    publishDir "${params.output_dir}", mode: params.publish_mode
    input:
    path comparison_files
    output:
    path "merged_comparisons.parquet", emit: merged_comparisons
    script:
    """
    zipstrain utilities merge_parquet --input-dir . --output-file merged_comparisons.parquet
    """
}
process fromSRAtoProfile{
    publishDir "${params.output_dir}/profiles", mode: params.publish_mode
    input:
    val sra_id
    path reference_genome
    path index_files
    path bed_file
    path gene_range_file
    path stb_file
    path null_model
    path profiling_contract
    output:
    path "${sra_id}_profile.parquet", emit: profiles
    path "${sra_id}_genome_stats.parquet", emit: genome_stats
    path "${sra_id}_gene_stats.parquet", emit: gene_stats
    val sra_id, emit: sample_name
    script:
    """
    prefetch --max-size ${params.prefetch_max_size} ${sra_id} 
    fasterq-dump --split-files --outdir ${sra_id} ${sra_id}
    rm -rf ${sra_id}/${sra_id}.sra
    num_seq_files=\$(ls ${sra_id}/*.fastq | wc -l)
    if [ \$num_seq_files -eq 2 ]; then
    bowtie2 -x ${reference_genome} -1 ${sra_id}/${sra_id}_1.fastq -2 ${sra_id}/${sra_id}_2.fastq --threads ${task.cpus} | samtools view -bS -F 4 - | samtools sort -@ ${task.cpus} -o ${sra_id}.bam -
    else
    bowtie2 -x ${reference_genome} -U ${sra_id}/${sra_id}*.fastq --threads ${task.cpus} | samtools view -bS -F 4 - | samtools sort -@ ${task.cpus} -o ${sra_id}.bam -
    fi
    zipstrain utilities profile-single \\
                        --bam-file ${sra_id}.bam \\
                        --reference-fasta ${reference_genome} \\
                        --bed-file ${bed_file} \\
                        --gene-range-table ${gene_range_file} \\
                        --stb-file ${stb_file} \\
                        --null-model ${null_model} \\
                        --profiling-contract ${profiling_contract} \\
                        --max-concurrency ${task.cpus} \\
                        --min-mapq ${params.min_mapq} \\
                        --min-baseq ${params.min_baseq} \\
                        --min-freq ${params.min_freq} \\
                        --min-read-ani ${params.min_read_ani} \\
                        --read-inclusion ${params.read_inclusion} \\
                        --output-dir .
    rm -rf ${sra_id}
    rm -f ${sra_id}.bam
    """
}

process prepare_profile_no_genes{

    publishDir "${params.output_dir}", mode: params.publish_mode
    input:
    path reference_genome
    path stb_file
    output:
    path "genomes_bed_file.bed", emit: genome_bed
    path "gene_range_table.tsv", emit: gene_range_table
    path "null_model.parquet", emit: null_model
    path "profiling_contract.json", emit: profiling_contract
    script:
"""
zipstrain utilities prepare_profiling \\
        --reference-fasta ${reference_genome} \\
        --stb-file ${stb_file} \\
        --error-rate ${params.error_rate} \\
        --max-total-reads ${params.max_total_reads} \\
        --p-threshold ${params.p_threshold} \\
        --output-dir .
"""

}

process fromSRAtoProfileBuildDb{
    publishDir "${params.output_dir}/profiles", mode: params.publish_mode
    input:
    val sra_id
    path sylph_db
    output:
    path "${sra_id}_profile.parquet", emit: profiles
    path "${sra_id}_genome_stats.parquet", emit: genome_stats
    path "${sra_id}_gene_stats.parquet", emit: gene_stats
    path "${sra_id}_sylph_abundance.tsv", emit: sylph_abundance
    val sra_id, emit: sample_name
    script:
    """
    prefetch --max-size ${params.prefetch_max_size} ${sra_id}
    fasterq-dump --split-files --outdir ${sra_id} ${sra_id}
    rm -rf ${sra_id}/${sra_id}.sra

    num_seq_files=\$(ls ${sra_id}/*.fastq | wc -l)
    if [ \$num_seq_files -eq 2 ]; then
        sylph profile ${sylph_db} -1 ${sra_id}/${sra_id}_1.fastq -2 ${sra_id}/${sra_id}_2.fastq -t ${task.cpus} > ${sra_id}_sylph_abundance.tsv
        bowtie_reads="-1 ${sra_id}/${sra_id}_1.fastq -2 ${sra_id}/${sra_id}_2.fastq"
    else
        sylph profile ${sylph_db} -U ${sra_id}/${sra_id}*.fastq -t ${task.cpus} > ${sra_id}_sylph_abundance.tsv
        bowtie_reads="-U ${sra_id}/${sra_id}*.fastq"
    fi

    zipstrain utilities build-genome-db \\
        --tool sylph \\
        --abundance-table ${sra_id}_sylph_abundance.tsv \\
        --cache-dir .genome_cache \\
        --output-dir .

    bowtie2-build --threads ${task.cpus} reference_genomes.fna reference_genomes.fna

    zipstrain utilities prepare_profiling \\
        --reference-fasta reference_genomes.fna \\
        --stb-file reference_genomes.stb \\
        --error-rate ${params.error_rate} \\
        --max-total-reads ${params.max_total_reads} \\
        --p-threshold ${params.p_threshold} \\
        --output-dir .

    bowtie2 -x reference_genomes.fna \$bowtie_reads --threads ${task.cpus} | samtools view -bS -F 4 - | samtools sort -@ ${task.cpus} -o ${sra_id}.bam -

    zipstrain utilities profile-single \\
        --bam-file ${sra_id}.bam \\
        --reference-fasta reference_genomes.fna \\
        --bed-file genomes_bed_file.bed \\
        --stb-file reference_genomes.stb \\
        --null-model null_model.parquet \\
        --profiling-contract profiling_contract.json \\
        --max-concurrency ${task.cpus} \\
        --min-mapq ${params.min_mapq} \\
        --min-baseq ${params.min_baseq} \\
        --min-freq ${params.min_freq} \\
        --min-read-ani ${params.min_read_ani} \\
        --read-inclusion ${params.read_inclusion} \\
        --output-dir .

    rm -rf ${sra_id}
    rm -f ${sra_id}.bam
    """
}
workflow
{
    if (!params.mode) {
        error "Set --mode to one of: map_reads, profile, from_sra_to_profile, compare"
    }

    if (params.mode == 'map_reads') {
        
        if (params.input_type=="sra")
        {
            table=tableToDict(file("${params.input_table}"))
            get_sequences_from_sra(Channel.fromList(table["Run"]))
            get_sequences_from_sra.out.sample_reads.set{sample_reads}

        }
        if (params.input_type=="local")
        {
            table=tableToDict(file("${params.input_table}"))
            if (table.containsKey("reads2")) {
                sample_reads = Channel.from(
                    [table["sample_name"], table["reads1"], table["reads2"]]
                        .transpose()
                        .collect { row -> tuple(row[0], [file(row[1]), file(row[2])]) }
                )
            }
            else {
                sample_reads = Channel.from(
                    [table["sample_name"], table["reads1"]]
                        .transpose()
                        .collect { row -> tuple(row[0], [file(row[1])]) }
                )
            }
        }
        if (!params.reference_genome)
        {
            if (params.sylph_db) {
                sylph_db = file(params.sylph_db)
            }
            else {
                download_sylph_db()
                download_sylph_db.out.sylph_db.set{ sylph_db }
            }
            estimate_abundance_sylph(sample_reads, sylph_db)
            merge_sylph_abundance_tables(estimate_abundance_sylph.out.abundance.collect())
            merge_sylph_abundance_tables.out.abundance.set{ abundance }
            build_db_from_Sylph(abundance,params.genome_db_cache_dir)
            build_db_from_Sylph.out.reference_genome.set{ reference_genome }
            index_reference(reference_genome)
            index_files = index_reference.out.index_files
        }
        else
        {
            reference_genome = file(params.reference_genome)
                if (params.index_files) {
            index_files = files(params.index_files)
                                        }
                else {
            index_reference(reference_genome)
            index_files = index_reference.out.index_files
                     }

        }

        map_reads(sample_reads,reference_genome,index_files)
    }
    if (params.mode == "from_sra_to_profile") {
        table=tableToDict(file("${params.input_table}"))
        sra_ids = Channel.fromList(table["Run"])
        if (!params.reference_genome) {
            if (params.sylph_db) {
                sylph_db = file(params.sylph_db)
            }
            else {
                download_sylph_db()
                download_sylph_db.out.sylph_db.set{ sylph_db }
            }
            fromSRAtoProfileBuildDb(sra_ids, sylph_db)
        }
        else {
            reference_genome = file(params.reference_genome)
            if (params.index_files) {
                index_files = file(params.index_files)
            }
            else {
                index_reference(reference_genome)
                index_files = index_reference.out.index_files
            }
            if (params.gene_file) {
                gene_file = file(params.gene_file)
                prepare_profile(reference_genome, gene_file, file(params.stb))
                genome_bed = prepare_profile.out.genome_bed
                gene_range_table = prepare_profile.out.gene_range_table
                null_model = prepare_profile.out.null_model
                profiling_contract = prepare_profile.out.profiling_contract
            }
            else {
                prepare_profile_no_genes(reference_genome, file(params.stb))
                genome_bed = prepare_profile_no_genes.out.genome_bed
                gene_range_table = prepare_profile_no_genes.out.gene_range_table
                null_model = prepare_profile_no_genes.out.null_model
                profiling_contract = prepare_profile_no_genes.out.profiling_contract
            }
            fromSRAtoProfile(sra_ids, reference_genome, index_files, genome_bed, gene_range_table, file(params.stb), null_model, profiling_contract)
        }
    }
    if (params.mode =='profile') {
        input_table = tableToDict(file(params.input_table))
        bamfiles = Channel.fromPath(input_table['bamfile'].collect{t->file(t)})
        sample_names = Channel.fromList(input_table['sample_name'])
        gene_file = params.gene_file ? file(params.gene_file) : null
        reference_genome = file(params.reference_genome)
        profile(bamfiles, sample_names, gene_file, reference_genome)
    }
    
    else if (params.mode == "compare")
    {
        input_table = tableToDict(file(params.input_table))
        def profile_pairs = []

        if (params.input_type=="profile_table"){
            profile_locations_list = getProfilesTableColumn(input_table).collect{t->file(t)}
            sample_names_list = getProfileSampleNamesTableColumn(input_table)
            profiles=[profile_locations_list,sample_names_list].transpose()
            (0..<profiles.size()).each { i ->
                (i+1..<profiles.size()).each { j ->
                    profile_pairs << (profiles[i] + profiles[j])
                }
            }
        }
        else if (params.input_type=="pair_table")
        {
            sample_1=input_table['sample_name_1']
            sample_2=input_table['sample_name_2']
            profile_1=input_table["profile_location_1"].collect{t->file(t)}
            profile_2=input_table["profile_location_2"].collect{t->file(t)}
            profile_pairs=([profile_1]+[sample_1]+[profile_2]+[sample_2]).transpose()
        }
        pair_channel = params.parallel_mode == "batched"
            ? Channel.fromList(makeBatches(profile_pairs, params.batch_size))
            : Channel.fromList(profile_pairs)
        stb = file(params.stb)
        // One comparison path: --compare_calculate picks the grain. Gene metrics
        // need ranges, since gene boundaries no longer live in the profile.
        if (params.compare_calculate?.toString()?.contains("gene") && !params.gene_range_table) {
            error "--compare_calculate '${params.compare_calculate}' asks for gene metrics, which require --gene_range_table (headerless TSV of gene, scaffold, start, end)."
        }
        // An empty collection is the standard optional-path pattern in Nextflow.
        // A real table is staged into every comparison task and tracked as an
        // input dependency; the empty collection keeps genome-only runs valid.
        gene_range_table = params.gene_range_table
            ? file(params.gene_range_table, checkIfExists: true)
            : []
        compare(pair_channel, stb, gene_range_table)
    }
}

workflow profile{
    take:
    bamfiles
    sample_names
    gene_file
    reference_genome
    main:
    if (gene_file) {
        prepare_profile(reference_genome, gene_file, file(params.stb))
        genome_bed = prepare_profile.out.genome_bed
        gene_range_table = prepare_profile.out.gene_range_table
        null_model = prepare_profile.out.null_model
        profiling_contract = prepare_profile.out.profiling_contract
    }
    else {
        prepare_profile_no_genes(reference_genome, file(params.stb))
        genome_bed = prepare_profile_no_genes.out.genome_bed
        gene_range_table = prepare_profile_no_genes.out.gene_range_table
        null_model = prepare_profile_no_genes.out.null_model
        profiling_contract = prepare_profile_no_genes.out.profiling_contract
    }
    profile_bam(sample_names, bamfiles, reference_genome, genome_bed, file(params.stb), gene_range_table, null_model, profiling_contract)
}


workflow compare
{
    take:
    profile_pairs
    stb
    gene_range_table
    main:
    if (params.parallel_mode=="single") {

    profile_pairs.multiMap{ v ->
        profile_1: v[0]
        profile_2: v[2]
        pair_name: v[1]+"_" + v[3]
    }.set{profile_pairs}
    compare_fast_profiles_single(profile_pairs.profile_1, profile_pairs.profile_2, stb, gene_range_table, profile_pairs.pair_name)
    merge_comparison_tables(compare_fast_profiles_single.out.comparison_results.collect() )
    }
    else if (params.parallel_mode=="batched") {
    profile_pairs.map{t->t.transpose()}.set{batch_t}
    batch_t.multiMap{ v ->
        unique_profiles: (v[0]+v[2]).unique().sort()
        pairs: [v[0].collect{t->t.name}, v[2].collect{t->t.name}].transpose()
    }.set{batch_pairs}

    compare_batched(batch_pairs.unique_profiles, batch_pairs.pairs, stb, gene_range_table)
    merge_comparison_tables(compare_batched.out.comparison_results.collect() )
    }
    else {
        error "Set --parallel_mode to either single or batched"
    }


}
