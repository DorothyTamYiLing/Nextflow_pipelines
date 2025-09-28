#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define parameters at the TOP of the script
params.raw = "./*{1,2}.fastq.gz"
params.outdir = "results"
params.reference = "./SRR5005407_contig_1.fasta"

process FASTQC_RAW {
    conda "/shared/team/conda/tyl205.bordetella-grou/nf-fastqc"
    container 'off'

    publishDir "${params.outdir}/01_rawfastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{html,zip}"

    script:
    """
    fastqc -t 10 ${reads[0]} ${reads[1]}
    """
}

process FASTP_TRIMMING {
    conda "/shared/team/conda/tyl205.bordetella-grou/fastp"
    container 'off'

    publishDir "${params.outdir}/02_trimmed", mode: 'copy'
    publishDir "${params.outdir}/03_fastp_reports", pattern: '*.{json,html}', mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1_trimmed.fastq.gz"), path("${sample_id}_2_trimmed.fastq.gz"), emit: trimmed_reads_gz
    path "*.{json,html}", emit: reports

    script:
    """
    fastp \
      -i ${reads[0]} -I ${reads[1]} \
      -o "${sample_id}_1_trimmed.fastq.gz" -O "${sample_id}_2_trimmed.fastq.gz" \
      --detect_adapter_for_pe --cut_front --cut_tail \
      --trim_front1 15 --trim_front2 15 \
      --trim_poly_g --trim_poly_x \
      --json "${sample_id}_fastp.json" \
      --html "${sample_id}_fastp_report.html"
    """
}

process GUNZIP {
    container 'off'
    
    publishDir "${params.outdir}/02_trimmed_unzipped", mode: 'copy'

    input:
    tuple val(sample_id), path(read1_gz), path(read2_gz)  

    output:
    tuple val(sample_id),path("${read1_gz.baseName}"), path("${read2_gz.baseName}"), emit: trimmed_fastq

    script:
    """
    gunzip -c "${read1_gz}" > "${read1_gz.baseName}"
    gunzip -c "${read2_gz}" > "${read2_gz.baseName}"
    """
}

process GEM_INDEXER {
    conda "/shared/team/conda/tyl205.bordetella-grou/gem3_mapper"
    container 'off'
    
    publishDir "${params.outdir}/00_index", mode: 'copy'

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}.gem", emit: gem_index
    path "${fasta_file.baseName}.info"

    script:
    """
    gem-indexer -h
    gem-indexer \
      -i "${fasta_file}"
    """
}

 process GEM3_MAPPING {
    conda "/shared/team/conda/tyl205.bordetella-grou/gem3_mapper"
    container 'off'
    
    publishDir "${params.outdir}/04_mapped", mode: 'copy'

//accept one channel only (i.e. the combined channel)
    input:
    tuple path(gem_index), val(sample_id), path(read1), path(read2)
    
    output:
tuple val(sample_id), path("${sample_id}.sam"), emit: sam_files

    script:
    """
    # Run gem-mapper with the already decompressed files
    /shared/team/Dorothy/map651PP1PP2samples_dir/gem3-mapper/bin/gem-mapper \
      -I "${gem_index}" \
      -1 "${read1}" \
      -2 "${read2}" \
      -o "${sample_id}.sam"
    """
}


process SAM2BAM {
    conda "/shared/team/conda/tyl205.bordetella-grou/bowtie2"
    container 'off'
    
    publishDir "${params.outdir}/05_bam", mode: 'copy'

    input:
    tuple val(sample_id), path(samfile)

    output:
    tuple val(sample_id), path("${sample_id}_GEM3_markdup_2652579to2654612.sorted.bam"), path("${sample_id}_GEM3_markdup_2652579to2654612.sorted.bam.bai"), emit: bam_files

    script:
    """
samtools view --threads 8 -bS "${samfile}" | \
samtools fixmate --threads 8 -m - - | \
samtools sort --threads 8 - | \
samtools markdup --threads 8 - "${sample_id}_GEM3_markdup.bam" && \
samtools index --threads 8 "${sample_id}_GEM3_markdup.bam"

#bam with metG range only
samtools view -b "${sample_id}_GEM3_markdup.bam" "contig_1:2652579-2654612" | \
samtools sort -o "${sample_id}_GEM3_markdup_2652579to2654612.sorted.bam" - && \
samtools index "${sample_id}_GEM3_markdup_2652579to2654612.sorted.bam"

rm "${sample_id}_GEM3_markdup.bam"
rm "${sample_id}_GEM3_markdup.bam.bai"
    """
}

process CONSENSUS {
    conda "/shared/team/conda/tyl205.bordetella-grou/bowtie2"
    container 'off'
    
    publishDir "${params.outdir}/06_consensus", mode: 'copy'

    input:
    tuple path(fasta_file), val(sample_id), path(bam), path(bam_bai)
    
    output:
    path "${sample_id}_metG_GEM3_40_forced_alt.bcf" 
    path "${sample_id}_metG_GEM3_40_forced_alt.bcf.csi" 
    path "${sample_id}_metG_consensus.fasta"

    script:
    """
    # Single pipeline with fewer intermediate files
bcftools mpileup -r contig_1:2652579-2654612 --annotate FORMAT/AD -f "$fasta_file" "$bam" | \
bcftools call -A --ploidy 1 -m | \
bcftools view -i 'FORMAT/AD[0:1] >= 40' | \
awk -F'\t' 'BEGIN{OFS="\\t"} {if(/^#/) {print} else {\$9="GT"; \$10="1/1"; print}}' | \ #This awk command modifies a VCF file to force all genotype (GT) fields to be "1/1" (homozygous alternate) for all variant positions. If the line starts with #, print the line unchanged
bcftools view -Ob -o "${sample_id}_metG_GEM3_40_forced_alt.bcf"

# Index and create consensus
tabix "${sample_id}_metG_GEM3_40_forced_alt.bcf"
samtools faidx "${fasta_file}" contig_1:2652579-2654612 | bcftools consensus "${sample_id}_metG_GEM3_40_forced_alt.bcf" > "${sample_id}_metG_consensus.fasta"
"""
}

workflow {
    // Create channels
    reads_ch = Channel.fromFilePairs(params.raw, checkIfExists: true)
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true)

    // Run processesm
    FASTQC_RAW(reads_ch)
    trimmed_gz_results = FASTP_TRIMMING(reads_ch)
 trimmed_gunzip_results = GUNZIP(trimmed_gz_results.trimmed_reads_gz)
 index_results = GEM_INDEXER(reference_ch)

     // COMBINE the index with each sample pair (combining two channels into one)
    mapping_input = index_results.gem_index.combine(trimmed_gunzip_results.trimmed_fastq)
 map_results = GEM3_MAPPING(mapping_input)
 bam_results= SAM2BAM(map_results.sam_files)

// COMBINE two channels
 reference_with_bams = reference_ch.combine(bam_results.bam_files)
consensus_results = CONSENSUS(reference_with_bams)
}

