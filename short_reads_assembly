#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.raw = "./*_{1,2}.fastq.gz"
params.outdir = "results"
params.reference = "./SRR5005407_contig_1.fasta"

process FASTQC_RAW {
    // Use working FastQC Docker container
    container 'quay.io/biocontainers/fastqc:0.11.9--0'
    
    // Specify platform for Apple Silicon
    containerOptions '--platform linux/amd64'    

    publishDir "${params.outdir}/01_rawfastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*.{html,zip}"
    
    script:
    """
    fastqc -t 10 ${reads}
    """
}


process FASTP_TRIMMING {
        // Use working FastQC Docker container
    container 'biocontainers/fastp:v0.20.1_cv1'
    
    // Specify platform for Apple Silicon
    containerOptions '--platform linux/amd64'  

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
container 'biocontainers/fastp:v0.20.1_cv1'  // Use a container with gunzip
    containerOptions '--platform linux/amd64'
    
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

process MINIMAP2_INDEX {
    container 'quay.io/biocontainers/minimap2:2.24--h7132678_1'
    containerOptions '--platform linux/amd64'

    publishDir "${params.outdir}/00_index", mode: 'copy'

    input:
    path reference

    output:
    path "${reference}.mmi", emit: minimap2_index

    script:
    """
    minimap2 -d ${reference}.mmi ${reference}
    """
}

process MINIMAP2 {
    container 'quay.io/biocontainers/minimap2:2.24--h7132678_1'
    containerOptions '--platform linux/amd64'

    publishDir "${params.outdir}/04_mapped", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2), path(index)

    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: sam_files

    script:
    """
    minimap2 -t 4 -a ${index} ${read1} ${read2} > ${sample_id}.sam
    """
}

process SAM2BAM {
container 'pegi3s/samtools_bcftools'
    containerOptions '--platform linux/amd64'
    
    publishDir "${params.outdir}/05_bam", mode: 'copy'

    input:
    tuple val(sample_id), path(samfile)

    output:
    tuple val(sample_id), path("${sample_id}_minimap_markdup_2652579to2654612.sorted.bam"), path("${sample_id}_minimap_markdup_2652579to2654612.sorted.bam.bai"), emit: bam_files

    script:
    """
samtools view --threads 8 -bS "${samfile}" | \
samtools fixmate --threads 8 -m - - | \
samtools sort --threads 8 - | \
samtools markdup --threads 8 - "${sample_id}_minimap_markdup.bam" && \
samtools index --threads 8 "${sample_id}_minimap_markdup.bam"

#bam with metG range only
samtools view -b "${sample_id}_minimap_markdup.bam" "contig_1:2652579-2654612" | \
samtools sort -o "${sample_id}_minimap_markdup_2652579to2654612.sorted.bam" - && \
samtools index "${sample_id}_minimap_markdup_2652579to2654612.sorted.bam"

rm "${sample_id}_minimap_markdup.bam"
rm "${sample_id}_minimap_markdup.bam.bai"
    """
}

process CONSENSUS {
container 'pegi3s/samtools_bcftools'
    containerOptions '--platform linux/amd64'
    
    publishDir "${params.outdir}/06_consensus", mode: 'copy'

    input:
    tuple path(fasta_file), val(sample_id), path(bam), path(bam_bai)
    
    output:
    path "${sample_id}_metG_minimap_40_forced_alt.bcf" 
    path "${sample_id}_metG_minimap_40_forced_alt.bcf.csi" 
    path "${sample_id}_metG_consensus.fasta"

    script:
    """
    # Single pipeline with fewer intermediate files
bcftools mpileup -r contig_1:2652579-2654612 --annotate FORMAT/AD -f "$fasta_file" "$bam" | \
bcftools call -A --ploidy 1 -m | \
bcftools view -i 'FORMAT/AD[0:1] >= 40' | \
awk -F'\t' 'BEGIN{OFS="\\t"} {if(/^#/) {print} else {\$9="GT"; \$10="1/1"; print}}' | \
bcftools view -Ob -o "${sample_id}_metG_minimap_40_forced_alt.bcf"

# Index and create consensus
bcftools index "${sample_id}_metG_minimap_40_forced_alt.bcf"
samtools faidx "${fasta_file}" contig_1:2652579-2654612 | bcftools consensus "${sample_id}_metG_minimap_40_forced_alt.bcf" > "${sample_id}_metG_consensus.fasta"
"""
}


workflow {
    // Create channels
    reads_ch = Channel.fromFilePairs(params.raw, checkIfExists: true)
reference_ch = Channel.fromPath(params.reference, checkIfExists: true)

    // Run processes
    FASTQC_RAW(reads_ch)
    trimmed_gz_results = FASTP_TRIMMING(reads_ch)
    trimmed_gunzip_results = GUNZIP(trimmed_gz_results.trimmed_reads_gz)

index_results = MINIMAP2_INDEX(reference_ch)

 mapping_input = trimmed_gunzip_results.trimmed_fastq.combine(index_results.minimap2_index)
 minimap_results = MINIMAP2(mapping_input)

 bam_results= SAM2BAM(minimap_results.sam_files)

// COMBINE two channels
 reference_with_bams = reference_ch.combine(bam_results.bam_files)
consensus_results = CONSENSUS(reference_with_bams)

}
