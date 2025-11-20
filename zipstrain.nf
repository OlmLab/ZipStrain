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
params.compare_memory_mode="light"
params.compare_chrom_batch_size=10000
params.batch_compare_n_parallel=4
params.publish_mode="link"
params.compare_genome_scope="all"
params.input_type="profile_table"


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
process get_sequences_from_sra {
    /*
    * This process retrieves sequences from the SRA database using the fastq-dump tool.
    * It takes in a list of SRA IDs and outputs the corresponding FASTQ files.
    * @param sra_ids: SRA ID to retrieve.
    */
    
    publishDir "${params.output_dir}/sra_sequences", mode: 'copy'
    
    input:
    val sra_ids
    
    output:
    path "${sra_ids}/${sra_ids}*.fastq.gz",emit: fastq_files
    val sra_ids, emit: sra_ids
    
    script:
    """
    prefetch ${sra_ids}
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
    bowtie2-build \\
            -f ${reference_genome} \\
            ${reference_genome}
    """
}

process map_reads{
    /*
    * This process maps reads to a reference genome using BWA.
    * It takes in the reference genome and reads, and outputs the BAM file.
    */
    publishDir "${params.output_dir}", mode: 'link'
    input:
    val sample_name
    path reference_genome
    path index_files
    path reads
    output:
    path "${sample_name}.bam", emit: bamfile
    script:
    if (reads.size() == 2) {
    """
    bowtie2 \\
            -x ${reference_genome} \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            --threads ${task.cpus} \\
            | samtools view -bS -F 4 - \\
            | samtools sort -@ ${task.cpus} -o ${sample_name}.bam -
    """}
    else {
        """
        bowtie2 \\
                -x ${reference_genome} \\
                -U ${reads[0]} \\
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
zipstrain profile prepare_profiling \\
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
    path gene_range_table
    output:
    path "${sample_name}.parquet", emit: profile
    path "${sample_name}.parquet.scaffolds", emit: covered_scaffolds
    val sample_name, emit: sample_name
    script:
    """
    zipstrain profile profile-single \\
                        --bam-file ${bamfile} \\
                        --bed-file ${bed_file} \\
                        --gene-range-table ${gene_range_table} \\
                        --num-workers ${task.cpus} \\
                        --output-dir .
    samtools idxstats ${bamfile} | awk '\$3 > 0 {print \$1}' >> ${sample_name}.parquet.scaffolds
    """
}

process compare_fast_profiles_single {
    /*
    * This process compares fast profiles.
    * It takes in the mpileup files and outputs the comparison results.
    */
    input:
    path mpileup_file1
    path mpileup_file2
    path scaffold_file1
    path scaffold_file2
    path null_model
    path stb
    val pair_name
    output:
    path "${pair_name}_comparison.parquet", emit: comparison_results
    script:
    """
    zipstrain compare single_compare_genome  \
                        --mpileup-contig-1 ${mpileup_file1} \
                        --mpileup-contig-2 ${mpileup_file2} \
                        --scaffolds-1 ${scaffold_file1} \
                        --scaffolds-2 ${scaffold_file2} \
                        --memory-mode ${params.compare_memory_mode} \
                        --chrom-batch-size ${params.compare_chrom_batch_size} \
                        -n ${null_model} \
                        -s ${stb} \
                        -c ${params.min_cov} \
                        -l ${params.min_gene_compare_len} \
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


process compare_fast_profiles_batched {
    publishDir "${params.output_dir}/batch_comparisons", mode: 'link'
    afterScript """
    rm -rf comps pairs.txt
    rm -f ${mpiles.collect{t->t.join(' ')}}
    rm -f ${scaffolds.collect{t->t.join(' ')}}
    rm -f ${null_model} ${stb}
    """
    input:
    path mpiles
    path scaffolds
    val pairs
    path null_model
    path stb


    output:
    path "Batch_*_comparisons.parquet", emit: comparison_results

    script:
    pairs_text = pairs.collect{p-> p.join('\t')}.join('\n')
    remove_mpiles = mpiles.join(' ')
    remove_scaffolds = scaffolds.join(' ')
    def add_genome_scope= (params.compare_genome_scope=="all") ? "" : "-g ${params.compare_genome_scope}"
    """
    echo -e "${pairs_text}" > pairs.txt
    cat pairs.txt | parallel --tmpdir . --colsep '\\t' -j ${params.batch_compare_n_parallel} 'zipstrain compare single_compare_genome \
                        --mpileup-contig-1 {1} \
                        --mpileup-contig-2 {2} \
                        --scaffolds-1 {1}.scaffolds \
                        --scaffolds-2 {2}.scaffolds \
                        --memory-mode ${params.compare_memory_mode} \
                        --chrom-batch-size ${params.compare_chrom_batch_size} \
                        -n ${null_model} \
                        -s ${stb} \
                        -c ${params.min_cov} \
                        -l ${params.min_gene_compare_len} \
                        -o {1}_{2}_comparison.parquet' ${add_genome_scope}
    mkdir comps
    hash=\$(sha1sum pairs.txt | awk '{print \$1}')
    mv *_comparison.parquet comps/
    zipstrain utilities merge_parquet  --input-dir comps --output-file "Batch_\${hash}_comparisons.parquet"
    rm -rf comps
    rm -f pairs.txt
    rm -f ${remove_mpiles}
    rm -f ${remove_scaffolds}

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
    output:
    path "${sra_id}.parquet", emit: profiles
    path "${sra_id}.parquet.scaffolds", emit: covered_scaffolds
    path "${sra_id}.parquet.breadth", emit: breadth
    val sra_id, emit: sample_name
    script:
    """
    prefetch ${sra_id}
    fasterq-dump --split-files --outdir ${sra_id} ${sra_id}
    gzip ${sra_id}/${sra_id}*.fastq
    rm -rf ${sra_id}/${sra_id}.sra
    num_seq_files=\$(ls ${sra_id}/*.fastq.gz | wc -l)
    if [ \$num_seq_files -eq 2 ]; then
    bowtie2 -x ${reference_genome} -1 ${sra_id}/${sra_id}_1.fastq.gz -2 ${sra_id}/${sra_id}_2.fastq.gz --threads ${task.cpus} | samtools view -bS -F 4 - | samtools sort -@ ${task.cpus} -o ${sra_id}.bam -
    else
    bowtie2 -x ${reference_genome} -U ${sra_id}/${sra_id}*.fastq.gz --threads ${task.cpus} | samtools view -bS -F 4 - | samtools sort -@ ${task.cpus} -o ${sra_id}.bam -
    fi
    zipstrain profile profile-single \\
                        --bam-file ${sra_id}.bam \\
                        --bed-file ${bed_file} \\
                        --gene-range-table ${gene_range_file} \\
                        --num-workers ${task.cpus} \\
                        --output-dir .
    samtools idxstats ${sra_id}.bam | awk '\$3 > 0 {print \$1}' > ${sra_id}.parquet.scaffolds
    zipstrain utilities genome_breadth_matrix --profile ${sra_id}.parquet \
                        --genome-length ${genome_length_file} \
                        --stb ${stb_file} \
                        --min-cov ${params.breadth_min_cov} \
                        --output-file ${profile.baseName}_breadth.parquet
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
            get_sequences_from_sra.out.fastq_files.map{t-> [t]}.set{reads}
            get_sequences_from_sra.out.sra_ids.set{sample_names}

        }
        if (params.input_type=="local")
        {
            table=tableToDict(file("${params.input_table}"))
            reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
            reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
            reads=reads_1.merge(reads_2)
            sample_names=Channel.fromList(table["sample_name"])
        }
        reference_genome = file(params.reference_genome)
        if (params.index_files) {
            index_files = files(params.index_files)
        }
        else {
            index_reference(reference_genome)
            index_files = index_reference.out.index_files
        }
        
        map_reads(sample_names,reference_genome,index_files,reads)
    }
    if (params.mode == "from_sra_to_profile") {
        table=tableToDict(file("${params.input_table}"))
        sra_ids = Channel.fromList(table["Run"])
        gene_file = file(params.gene_file)
        reference_genome = file(params.reference_genome)
        if (params.index_files) {
            index_files = files(params.index_files)
        }
        else {
            index_reference(reference_genome)
            index_files = index_reference.out.index_files
        }
        prepare_profile(reference_genome, gene_file, file(params.stb))

        fromSRAtoProfile(sra_ids, reference_genome, index_files, prepare_profile.out.genome_bed, prepare_profile.out.out.gene_range_table)
    }
    if (params.mode =='fast_profile') {
        input_table = tableToDict(file(params.input_table))
        bamfiles = Channel.fromPath(input_table['bamfile'].collect{t->file(t)})
        sample_names = Channel.fromList(input_table['sample_name'])
        gene_file = file(params.gene_file)
        reference_genome = file(params.reference_genome)
        fast_profile(bamfiles, sample_names, gene_file, reference_genome)
    }
    else if (params.mode =='fast_compare') {
        input_table = tableToDict(file(params.input_table))
        
        if (params.input_type=="profile_table"){
            mpileup_files_list = input_table['mpileup_files'].collect{t->file(t)}
            sample_names_list = input_table['sample_names']
            scaffolds_list = input_table['scaffolds'].collect{t->file(t)}
            profiles=[mpileup_files_list,sample_names_list,scaffolds_list].transpose()
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
            scaffold_1=input_table["scaffold_location_1"].collect{t->file(t)}
            scaffold_2=input_table["scaffold_location_2"].collect{t->file(t)}
            profile_pairs=([mpile_1]+[sample_1]+[scaffold_1]+[mpile_2]+[sample_2]+[scaffold_2]).transpose()
            pair_channel=Channel.from(profile_pairs)

        }
        stb = file(params.stb)
        fast_compare(pair_channel, stb)

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
    profile_bam(sample_names, bamfiles, prepare_profile.out.genome_bed, prepare_profile.out.gene_range_table)
    get_genome_breadth(profile_bam.out.profile,
                       file(params.stb),
                       make_bed_file.out.bed_file)
    merge_breadth_tables(get_genome_breadth.out.genome_breadth.collect())

}


workflow fast_compare
{
    take:
    profile_pairs
    stb
    main:

    if (!params.null_model) {
        build_null_model()
        null_model = build_null_model.out.model
    }
    else {
        null_model =file(params.null_model)
    }
    if (params.parallel_mode=="single") {

    profile_pairs.multiMap{ v ->
        mpile_1: v[0]
        mpile_2: v[3]
        scaffold_1: v[2]
        scaffold_2: v[5]
        pair_name: v[1]+"_" + v[4]
    }.set{profile_pairs}
    compare_fast_profiles_single(profile_pairs.mpile_1, profile_pairs.mpile_2, profile_pairs.scaffold_1, profile_pairs.scaffold_2, null_model, stb, profile_pairs.pair_name)
    merge_tables(compare_fast_profiles_single.out.comparison_results.collect() )
    }
    else if (params.parallel_mode=="batched") {
    batch = profile_pairs.collate(params.batch_size)
    batch.map{t->t.transpose()}.set{batch_t}
    batch_t.multiMap{ v ->
        unique_mpiles: (v[0]+v[3]).unique().sort()
        unique_scaffolds: (v[2]+v[5]).unique().sort()
        pairs: [v[0].collect{t->t.name}, v[3].collect{t->t.name}].transpose()
    }.set{batch_pairs}

    compare_fast_profiles_batched(batch_pairs.unique_mpiles, batch_pairs.unique_scaffolds, batch_pairs.pairs, null_model, stb)
    merge_tables(compare_fast_profiles_batched.out.comparison_results.collect() )
    }


}


