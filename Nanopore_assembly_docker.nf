#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = "results"

process SPLIT_FASTQ {
    tag "$sample_id"

    container 'quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0'
    containerOptions '--platform linux/amd64'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
 tuple val(sample_id), path("chunks/${sample_id}.part_*.fastq.gz"), emit: chunks

    script:
    """
seqkit split2 -s 100 --threads 4 -O chunks ${reads} --force
    """
}

process PORECHOP {
    tag "$sample_id:${chunk.baseName}"

    container 'quay.io/biocontainers/porechop:0.2.3_seqan2.1.1--py36h2d50403_3'
    containerOptions '--platform linux/amd64'

    cpus 4
    memory '8 GB'
    time '12h'

    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(chunk)

    output:
    tuple val(sample_id), path("${chunk.baseName}_rmadap.fastq")

    script:
    """
    porechop -i ${chunk} -o ${chunk.baseName}_rmadap.fastq --discard_middle -v 1
    """
}

process MERGE_TRIMMED {
    tag "$sample_id"

    container 'quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0'
    containerOptions '--platform linux/amd64'

    publishDir "${params.outdir}/merged", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_chunks) 
    //trimmed_chunks will be a list of FASTQ files belonging to the same sample.

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq")

    script:
    """
    cat ${trimmed_chunks} > ${sample_id}_trimmed.fastq
    """
}

process GZIP {
    container 'biocontainers/fastp:v0.20.1_cv1'  // Use a container with gunzip
    containerOptions '--platform linux/amd64'

    
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq")  

    output:
    tuple val(sample_id),path("${sample_id}_trimmed.fastq.gz")  

    script:
    """
    gzip -c "${sample_id}_trimmed.fastq" > "${sample_id}_trimmed.fastq.gz"
    """
}

process NANOPLOT {
    tag "$sample_id"

    container 'ttubb/nanoplot:latest'
    containerOptions '--platform linux/amd64'

    publishDir "${params.outdir}/nanoplot", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("nanoplot_${sample_id}"), emit: np_out

    script:
    """
    NanoPlot -t ${task.cpus} \
             --fastq ${reads} \
             -o nanoplot_${sample_id}
    """
}

process NANOFILT {
    tag "$sample_id"

    container 'mcfonsecalab/nanofilt'
    containerOptions '--platform linux/amd64'

    publishDir "${params.outdir}/nanofilt", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

output:
    tuple val(sample_id), path("${sample_id}_nanofilt.fastq"), emit: nanofilt_out


    script:
    """
    zcat ${reads} | NanoFilt --quality 7 --length 150 --headcrop 24 > ${sample_id}_nanofilt.fastq
    """
}

process FLYE {
    tag "$sample_id"

    container 'staphb/flye'
    containerOptions '--platform linux/amd64'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("flye_${sample_id}"), emit: flye_out

    script:
    """
flye --nano-corr ${reads} \
--genome-size 4.5m  --out-dir flye_${sample_id} \
--threads 8
    """
}


process MEDAKA {
    tag "$sample_id"

    container 'ontresearch/medaka'
    containerOptions '--platform linux/amd64'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id),path(reads),path(assembly)

    output:
    tuple val(sample_id), path("medaka_${sample_id}")

    script:
    """
medaka_consensus -i ${reads} \
-d ${assembly} \
-o medaka_${sample_id} -t 8 -m r103_sup_g507
    """
}


workflow {
//the sample_id in the channel tuple will be output ".fastq" or ".fastq.gz"
    reads_ch = Channel.fromPath("*.fastq.gz", checkIfExists: true)
                      .map { file ->
                          def id = file.name.replaceAll(/\.fastq(\.gz)?$/, "")
                          tuple(id, file)
                      }

reads_ch.view()


//After the SPLIT_FASTQ process (which splits fastq into chunks of X reads), create a tuple for each chunk that will be fed into porechop separately.

split_ch = SPLIT_FASTQ(reads_ch)
                    .flatMap { sample_id, files -> 
                        files.collect { file -> tuple(sample_id, file) } 
                    }

    split_ch.view()

    porechop_results=PORECHOP(split_ch)

    merged_results = porechop_results
                        .groupTuple()    // groups by sample_id, see notes below
                        .set { merged_in }

    trim_results=MERGE_TRIMMED(merged_in)

    gzip_Fastq=GZIP(trim_results)

    nanoplot_results = NANOPLOT(gzip_Fastq)

    nanofilt_results=NANOFILT(gzip_Fastq)

    flye_results=FLYE(nanofilt_results)

    flye_assembly = flye_results.flye_out.map { sample_id, dir ->
        tuple(sample_id, file("${dir}/assembly.fasta"))
    }

// COMBINE two channels
reads_with_assembly = nanofilt_results.join(flye_assembly)
    .map { sid, reads, assembly ->
        tuple(sid, reads, assembly)
    }

    reads_with_assembly.view()

 MEDAKA(reads_with_assembly)

}
