params.error_rate=0.001
params.max_total_reads=10000
params.p_threshold=0.05
params.mode = 'fast_compare' // 'fast_compare' or 'profile_genes' or
params.parallel_mode="multiple"
params.min_cov=5
params.min_gene_compare_len=200
params.batch_size=10
params.breadth_min_cov=1
params.bed_max_scaffold_length=500000
params.compare_duckdb_memory_limit=""
params.batch_compare_n_parallel=4
params.publish_mode="link"
params.compare_genome_scope="all"
params.compare_gene_scope="all:all"
params.compare_genome_calculate="all"
params.compare_ani_method="popani"
params.input_type="profile_table"
params.bowtie2_non_competitive_mapping=false
params.sylph_db = null
params.sylph_db_link="http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r220-c200-dbv1.syldb"
params.genome_db_cache_dir="genome_cache"
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
    publishDir "${params.output_dir}/sylph_abundance/", mode: 'copy'
    
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
    publishDir "${params.output_dir}/sylph_abundance/", mode: 'copy'
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
    publishDir "${params.output_dir}/db_from_sylph/", mode: 'link'
    containerOptions {
        switch( workflow.containerEngine ) {
            case 'docker':
                return "--volume ${genome_cache_dir}:${genome_cache_dir}"
            case 'apptainer':
                return "--bind ${genome_cache_dir}:${genome_cache_dir}"
            default:
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
    
    publishDir "${params.output_dir}/sra_sequences", mode: 'link'
    
    input:
    val sra_ids
    
    output:
    path "${sra_ids}/${sra_ids}*.fastq.gz",emit: fastq_files
    val sra_ids, emit: sra_ids
    tuple val(sra_ids), path("${sra_ids}/${sra_ids}*.fastq.gz"), emit: sample_reads
    
    script:
    """
    prefetch --max-size 200g ${sra_ids} 
    fasterq-dump --split-files --outdir ${sra_ids} ${sra_ids}
    gzip ${sra_ids}/${sra_ids}*.fastq
    rm -rf ${sra_ids}/${sra_ids}.sra
    """
}
process index_reference {
    /*
    * This process indexes the reference genome using BWA.
    * It takes in the reference genome and outputs the indexed files.
    */
    publishDir "${params.output_dir}", mode: 'link'
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
    publishDir "${params.output_dir}", mode: 'link'
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

    publishDir "${params.output_dir}", mode: 'link'
    input:
    path reference_genome
    path gene_fasta
    path stb_file
    output:
    path "genomes_bed_file.bed", emit: genome_bed
    path "gene_range_table.tsv", emit: gene_range_table
    path "genome_lengths.parquet", emit: genome_lengths
    script:
"""
zipstrain utilities prepare_profiling \\
        --reference-fasta ${reference_genome} \\
        --gene-fasta ${gene_fasta} \\
        --stb-file ${stb_file} \\
        --output-dir . 
"""

}


process profile_bam {
    /*
    * This process generates mpileup files from BAM files.
    * It takes in the BAM file and outputs the mpileup file.
    */
    publishDir "${params.output_dir}", mode: 'link'
    input:
    val sample_name
    path bamfile
    path bed_file
    path stb_file
    path gene_range_table
    path null_model
    output:
    path "${bamfile.baseName}_profile.parquet", emit: profile
    path "${bamfile.baseName}_genome_stats.parquet", emit: genome_stats
    path "${bamfile.baseName}_gene_stats.parquet", emit: gene_stats
    val sample_name, emit: sample_name
    script:
    """
    zipstrain utilities profile-single \\
                        --bam-file ${bamfile} \\
                        --bed-file ${bed_file} \\
                        --gene-range-table ${gene_range_table} \\
                        --stb-file ${stb_file} \\
                        --null-model ${null_model} \\
                        --max-concurrency ${task.cpus} \\
                        --output-dir .
    """
}

process compare_genome_fast_profiles_single {
    /*
    * This process compares fast profiles.
    * It takes in the mpileup files and outputs the comparison results.
    */
    input:
    path mpileup_file1
    path mpileup_file2
    path stb
    val pair_name
    output:
    path "${pair_name}_comparison.parquet", emit: comparison_results
    script:
    def add_genome_scope= (params.compare_genome_scope=="all") ? "" : "-g ${params.compare_genome_scope}"
    def add_duckdb_memory_limit = params.compare_duckdb_memory_limit ? "--duckdb-memory-limit ${params.compare_duckdb_memory_limit}" : ""
    def add_calculate = params.compare_genome_calculate ? "--calculate ${params.compare_genome_calculate}" : ""
    """
    zipstrain utilities single_compare_genome  \
                        --mpileup-contig-1 ${mpileup_file1} \
                        --mpileup-contig-2 ${mpileup_file2} \
                        -s ${stb} \
                        -c ${params.min_cov} \
                        -l ${params.min_gene_compare_len} \
                        --ani-method ${params.compare_ani_method} \
                        ${add_calculate} \
                        ${add_duckdb_memory_limit} \
                        ${add_genome_scope} \
                        -o ${pair_name}_comparison.parquet

    """
}
process compare_gene_fast_profiles_single {
    /*
    * This process compares fast profiles.
    * It takes in the mpileup files and outputs the comparison results.
    */
    input:
    path mpileup_file1
    path mpileup_file2
    path stb
    val pair_name
    output:
    path "${pair_name}_comparison.parquet", emit: comparison_results
    script:
    def add_duckdb_memory_limit = params.compare_duckdb_memory_limit ? "--duckdb-memory-limit ${params.compare_duckdb_memory_limit}" : ""
    """
    zipstrain utilities single_compare_gene  \
                        --mpileup-contig-1 ${mpileup_file1} \
                        --mpileup-contig-2 ${mpileup_file2} \
                        -s ${stb} \
                        -c ${params.min_cov} \
                        -l ${params.min_gene_compare_len} \
                        --ani-method ${params.compare_ani_method} \
                        ${add_duckdb_memory_limit} \
                        --scope ${params.compare_gene_scope} \
                        -o ${pair_name}_comparison.parquet

    """
}


process get_genome_breadth {
    /*
    * This process calculates the genome breadth.
    * It takes in the mpileup files and outputs the genome breadth.
    */
    publishDir "${params.output_dir}", mode: 'link'
    input:
    path profile
    path stb_file
    path bed_file
    output:
    path "${profile.baseName}_breadth.parquet", emit: genome_breadth
    script:
    """
    zipstrain utilities get_genome_lengths --stb-file ${stb_file} --bed-file ${bed_file} --output-file genome_length.parquet 
    zipstrain utilities genome_breadth_matrix --profile ${profile} \
                        --genome-length genome_length.parquet \
                        --stb ${stb_file} \
                        --min-cov ${params.breadth_min_cov} \
                        --output-file ${profile.baseName}_breadth.parquet
    """
}


process compare_genome_batched {
    publishDir "${params.output_dir}/batch_comparisons", mode: 'link'
    afterScript """
    rm -rf comps pairs.txt
    rm -f ${mpiles.collect{t->t.join(' ')}}
    rm -f ${stb}
    """
    input:
    path mpiles
    val pairs
    path stb


    output:
    path "Batch_*_comparisons.parquet", emit: comparison_results

    script:
    pairs_text = pairs.collect{p-> p.join('\t')}.join('\n')
    remove_mpiles = mpiles.join(' ')
    def add_genome_scope= (params.compare_genome_scope=="all") ? "" : "-g ${params.compare_genome_scope}"
    def add_duckdb_memory_limit = params.compare_duckdb_memory_limit ? "--duckdb-memory-limit ${params.compare_duckdb_memory_limit}" : ""
    def add_calculate = params.compare_genome_calculate ? "--calculate ${params.compare_genome_calculate}" : ""
    """
    echo -e "${pairs_text}" > pairs.txt
    cat pairs.txt | parallel --tmpdir . --colsep '\\t' -j ${params.batch_compare_n_parallel} 'zipstrain utilities single_compare_genome \
                        --mpileup-contig-1 {1} \
                        --mpileup-contig-2 {2} \
                        -s ${stb} \
                        -c ${params.min_cov} \
                        -l ${params.min_gene_compare_len} \
                        --ani-method ${params.compare_ani_method} \
                        ${add_calculate} \
                        ${add_duckdb_memory_limit} \
                        -o {1}_{2}_comparison.parquet' ${add_genome_scope}
    mkdir comps
    hash=\$(sha1sum pairs.txt | awk '{print \$1}')
    mv *_comparison.parquet comps/
    zipstrain utilities merge_parquet  --input-dir comps --output-file "Batch_\${hash}_comparisons.parquet"
    rm -rf comps
    rm -f pairs.txt
    rm -f ${remove_mpiles}

    """


}

process compare_gene_batched {
    publishDir "${params.output_dir}/batch_comparisons", mode: 'link'
    afterScript """
    rm -rf comps pairs.txt
    rm -f ${mpiles.collect{t->t.join(' ')}}
    rm -f ${stb}
    """
    input:
    path mpiles
    val pairs
    path stb
    output:
    path "Batch_*_comparisons.parquet", emit: comparison_results
    script:
    pairs_text = pairs.collect{p-> p.join('\t')}.join('\n')
    remove_mpiles = mpiles.join(' ')
    def add_gene_scope= (params.compare_gene_scope=="all") ? "" : "--scope ${params.compare_gene_scope}"
    def add_duckdb_memory_limit = params.compare_duckdb_memory_limit ? "--duckdb-memory-limit ${params.compare_duckdb_memory_limit}" : ""
    """
    echo -e "${pairs_text}" > pairs.txt
    cat pairs.txt | parallel --tmpdir . --colsep '\\t' -j ${params.batch_compare_n_parallel} 'zipstrain utilities single_compare_gene \
                        --mpileup-contig-1 {1} \
                        --mpileup-contig-2 {2} \
                        -s ${stb} \
                        -c ${params.min_cov} \
                        -l ${params.min_gene_compare_len} \
                        --ani-method ${params.compare_ani_method} \
                        ${add_duckdb_memory_limit} \
                        -o {1}_{2}_comparison.parquet' ${add_gene_scope}
    mkdir comps
    hash=\$(sha1sum pairs.txt | awk '{print \$1}')
    mv *_comparison.parquet comps/
    zipstrain utilities merge_parquet  --input-dir comps --output-file "Batch_\${hash}_comparisons.parquet"
    rm -rf comps
    rm -f pairs.txt
    rm -f ${remove_mpiles}

    """ 
}



process merge_comparison_tables {
    /*
    * This process merges multiple comparison result files into a single file.
    * It takes in multiple comparison result files and outputs a merged file.
    */
    publishDir "${params.output_dir}", mode: 'link'
    input:
    path comparison_files
    output:
    path "merged_comparisons.parquet", emit: merged_comparisons
    script:
    """
    zipstrain utilities merge_parquet --input-dir . --output-file merged_comparisons.parquet
    """
}
process merge_breadth_tables {
    /*
    * This process merges multiple genome breadth files into a single file.
    * It takes in multiple genome breadth files and outputs a merged file.
    */
    publishDir "${params.output_dir}"
    input:
    path breadth_files
    output:
    path "merged_breadth.parquet", emit: merged_breadth
    script:
    """
    mkdir -p breadth_tables
    mv ${breadth_files} breadth_tables/
    zipstrain utilities collect_breadth_tables --breadth-tables-dir breadth_tables --output-file merged_breadth.parquet
    """
}

process build_null_model {
publishDir "${params.output_dir}", mode: 'link'
    /*
    * This process builds a null model for the gene abundance data.
    * It takes in the gene file and scaffold table, and outputs the null model.
    */
    output:
    path "null_model.parquet", emit: model
    script:
    """
    zipstrain utilities build-null-model \
                      --error-rate ${params.error_rate} \
                      --max-total-reads ${params.max_total_reads} \
                      --p-threshold ${params.p_threshold} \
                      --output-file "null_model.parquet" 
    """
}

process fromSRAtoProfile{
    publishDir "${params.output_dir}/profiles"
    input:
    val sra_id
    path reference_genome
    path index_files
    path bed_file
    path gene_range_file
    path genome_length_file
    path stb_file
    path null_model
    output:
    path "${sra_id}_profile.parquet", emit: profiles
    path "${sra_id}_genome_stats.parquet", emit: genome_stats
    path "${sra_id}_gene_stats.parquet", emit: gene_stats
    val sra_id, emit: sample_name
    script:
    """
    prefetch --max-size 200g ${sra_id} 
    fasterq-dump --split-files --outdir ${sra_id} ${sra_id}
    gzip ${sra_id}/${sra_id}*.fastq
    rm -rf ${sra_id}/${sra_id}.sra
    num_seq_files=\$(ls ${sra_id}/*.fastq.gz | wc -l)
    if [ \$num_seq_files -eq 2 ]; then
    bowtie2 -x ${reference_genome} -1 ${sra_id}/${sra_id}_1.fastq.gz -2 ${sra_id}/${sra_id}_2.fastq.gz --threads ${task.cpus} | samtools view -bS -F 4 - | samtools sort -@ ${task.cpus} -o ${sra_id}.bam -
    else
    bowtie2 -x ${reference_genome} -U ${sra_id}/${sra_id}*.fastq.gz --threads ${task.cpus} | samtools view -bS -F 4 - | samtools sort -@ ${task.cpus} -o ${sra_id}.bam -
    fi
    zipstrain utilities profile-single \\
                        --bam-file ${sra_id}.bam \\
                        --bed-file ${bed_file} \\
                        --gene-range-table ${gene_range_file} \\
                        --stb-file ${stb_file} \\
                        --null-model ${null_model} \\
                        --max-concurrency ${task.cpus} \\
                        --output-dir .
    rm -rf ${sra_id}
    rm -f ${sra_id}.bam
    """
}
workflow
{

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
        gene_file = file(params.gene_file)
        reference_genome = file(params.reference_genome)
        if (params.index_files) {
            index_files = file(params.index_files)
        }
        else {
            index_reference(reference_genome)
            index_files = index_reference.out.index_files
        }
        prepare_profile(reference_genome, gene_file, file(params.stb))
        build_null_model()

        fromSRAtoProfile(sra_ids, reference_genome, index_files, prepare_profile.out.genome_bed, prepare_profile.out.gene_range_table, prepare_profile.out.genome_lengths, file(params.stb), build_null_model.out.model)
    }
    if (params.mode =='profile') {
        input_table = tableToDict(file(params.input_table))
        bamfiles = Channel.fromPath(input_table['bamfile'].collect{t->file(t)})
        sample_names = Channel.fromList(input_table['sample_name'])
        gene_file = file(params.gene_file)
        reference_genome = file(params.reference_genome)
        fast_profile(bamfiles, sample_names, gene_file, reference_genome)
    }
    
    else if (params.mode =="compare_genes")
    {
        input_table = tableToDict(file(params.input_table))
        
        if (params.input_type=="profile_table"){
            mpileup_files_list = input_table['mpileup_files'].collect{t->file(t)}
            sample_names_list = input_table['sample_names']
            profiles=[mpileup_files_list,sample_names_list].transpose()
            def profile_pairs = []
            for (int i = 0; i < profiles.size(); i++) {
                for (int j = i + 1; j < profiles.size(); j++) {
                    profile_pairs << (profiles[i] + profiles[j])
                    }
            }
            pair_channel=Channel.from(profile_pairs)
        }
        else if (params.input_type=="pair_table")
        {
            sample_1=input_table['sample_name_1']
            sample_2=input_table['sample_name_2']
            mpile_1=input_table["profile_location_1"].collect{t->file(t)}
            mpile_2=input_table["profile_location_2"].collect{t->file(t)}
            profile_pairs=([mpile_1]+[sample_1]+[mpile_2]+[sample_2]).transpose()
            pair_channel=Channel.from(profile_pairs)

        }
        stb = file(params.stb)
        compare_genes(pair_channel, stb)
    }
    
    
    else if (params.mode =='compare_genomes') {
        input_table = tableToDict(file(params.input_table))
        
        if (params.input_type=="profile_table"){
            mpileup_files_list = input_table['mpileup_files'].collect{t->file(t)}
            sample_names_list = input_table['sample_names']
            profiles=[mpileup_files_list,sample_names_list].transpose()
            def profile_pairs = []
            for (int i = 0; i < profiles.size(); i++) {
                for (int j = i + 1; j < profiles.size(); j++) {
                    profile_pairs << (profiles[i] + profiles[j])
                    }
            }
            pair_channel=Channel.from(profile_pairs)
        }
        else if (params.input_type=="pair_table")
        {
            sample_1=input_table['sample_name_1']
            sample_2=input_table['sample_name_2']
            mpile_1=input_table["profile_location_1"].collect{t->file(t)}
            mpile_2=input_table["profile_location_2"].collect{t->file(t)}
            profile_pairs=([mpile_1]+[sample_1]+[mpile_2]+[sample_2]).transpose()
            pair_channel=Channel.from(profile_pairs)

        }
        stb = file(params.stb)
        compare_genomes(pair_channel, stb)

}}
workflow profile_contigs
{
    input_table = tableToDict(file(params.input_table))
    contig_tables = file(params.contig_tables)
    sample_names = Channel.fromList(input_table['sample_name'])
    bamfiles = Channel.fromPath(input_table['bamfile'].collect{t->file(t)})
    get_mpileup_contigs(sample_names, contig_tables, bamfiles, file(params.gene_file))
}

workflow fast_profile{
    take:
    bamfiles
    sample_names
    gene_file
    reference_genome
    main:
    prepare_profile(reference_genome, gene_file, file(params.stb))
    build_null_model()
    profile_bam(sample_names, bamfiles, prepare_profile.out.genome_bed,file(params.stb), prepare_profile.out.gene_range_table, build_null_model.out.model)
}


workflow compare_genomes
{
    take:
    profile_pairs
    stb
    main:
    if (params.parallel_mode=="single") {

    profile_pairs.multiMap{ v ->
        mpile_1: v[0]
        mpile_2: v[2]
        pair_name: v[1]+"_" + v[3]
    }.set{profile_pairs}
    compare_genome_fast_profiles_single(profile_pairs.mpile_1, profile_pairs.mpile_2, stb, profile_pairs.pair_name)
    merge_comparison_tables(compare_genome_fast_profiles_single.out.comparison_results.collect() )
    }
    else if (params.parallel_mode=="batched") {
    batch = profile_pairs.collate(params.batch_size)
    batch.map{t->t.transpose()}.set{batch_t}
    batch_t.multiMap{ v ->
        unique_mpiles: (v[0]+v[2]).unique().sort()
        pairs: [v[0].collect{t->t.name}, v[2].collect{t->t.name}].transpose()
    }.set{batch_pairs}

    compare_genome_batched(batch_pairs.unique_mpiles, batch_pairs.pairs, stb)
    merge_comparison_tables(compare_genome_batched.out.comparison_results.collect() )
    }


}


workflow compare_genes{
    take:
    profile_pairs
    stb
    main:
    if (params.parallel_mode=="single") {

    profile_pairs.multiMap{ v ->
        mpile_1: v[0]
        mpile_2: v[2]
        pair_name: v[1]+"_" + v[3]
    }.set{profile_pairs}
    compare_gene_fast_profiles_single(profile_pairs.mpile_1, profile_pairs.mpile_2, stb, profile_pairs.pair_name)
    merge_comparison_tables(compare_gene_fast_profiles_single.out.comparison_results.collect() )
    }
    else if (params.parallel_mode=="batched") {
    batch = profile_pairs.collate(params.batch_size)
    batch.map{t->t.transpose()}.set{batch_t}
    batch_t.multiMap{ v ->
        unique_mpiles: (v[0]+v[2]).unique().sort()
        pairs: [v[0].collect{t->t.name}, v[2].collect{t->t.name}].transpose()
    }.set{batch_pairs}

    compare_gene_batched(batch_pairs.unique_mpiles, batch_pairs.pairs, stb)
    merge_comparison_tables(compare_gene_batched.out.comparison_results.collect() )
    }

}
