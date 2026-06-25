params.memory = "3g"
params.cpus = 10
params.outdir = "."

params.enable_conda = true

process VIRSTRAIN_CALL {

    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/04_serotype", mode: 'copy'

    conda "bioconda::virstrain"

    input:
        tuple val(sid), path(reads)
        each path(pathdb)

    output:
        tuple val(sid), path("${sid}.VirStrain_report.txt"), path(reads), emit: report_for_align
        path("${sid}.Mps_ps_depth.csv"), optional: true
        path("${sid}.Ops_ps_depth.csv"), optional: true
        path("${sid}.VirStrain_report.html"), optional: true
        //path("${sid}.VirStrain_report.txt"), optional: true
    script:
    if (params.mode == "PE")
    """
    virstrain \
        -i ${reads[0]} \
        -p ${reads[1]} \
        -d ${pathdb} \
        -o .
   ## rename VirStrain outputs if present from upstream staging
    [[ -f Mps_ps_depth.csv ]] && mv Mps_ps_depth.csv ${sid}.Mps_ps_depth.csv || true
    [[ -f Ops_ps_depth.csv ]] && mv Ops_ps_depth.csv ${sid}.Ops_ps_depth.csv || true
    [[ -f VirStrain_report.html ]] && mv VirStrain_report.html ${sid}.VirStrain_report.html || true
    [[ -f VirStrain_report.txt ]] && mv VirStrain_report.txt ${sid}.VirStrain_report.txt || true
    """
    else
    """
      virstrain \
        -i ${reads} \
        -d ${pathdb} \
        -o .
    [[ -f Mps_ps_depth.csv ]] && mv Mps_ps_depth.csv ${sid}.Mps_ps_depth.csv || true
    [[ -f Ops_ps_depth.csv ]] && mv Ops_ps_depth.csv ${sid}.Ops_ps_depth.csv || true
    [[ -f VirStrain_report.html ]] && mv VirStrain_report.html ${sid}.VirStrain_report.html || true
    [[ -f VirStrain_report.txt ]] && mv VirStrain_report.txt ${sid}.VirStrain_report.txt || true
    """
}

process VIRSTRAIN_ALIGN_BESTMATCH {

    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/04_serotype", mode: 'copy', saveAs: { filename ->
    (
        filename.endsWith('.fasta') ||
        filename.endsWith('.txt')   ||
        filename.endsWith('.csv')   ||
        filename.endsWith('.bam')   ||
        filename.endsWith('.bai')   ||
        filename.endsWith('.png')
    ) ? filename : null
}

    conda "bioconda::samtools bioconda::seqkit bioconda::mm2plus bioconda::bwa-mem2 bioconda::samplot bioconda::covtobed"

    input:
        tuple val(sid), path(vs_report), path(reads)
        each path(meta)

    output:
        tuple val(sid), path("${sid}*.fasta"), emit: serotype_contig_ordering
        tuple val(sid), path("${sid}*.fq.gz"), emit: viral_reads
        path("${sid}.serotype.txt"), emit: serotyper_res
        //path("${sid}.VirStrain_report.txt")
        path("*.bam"), optional: true
        path("*.bai"), optional: true
        path("*.png"), optional: true
        path("*.csv"), optional: true
        path("*.txt"), optional: true
        path("*.html"), optional: true

    script:
    if (params.mode == "PE")
    """
    ${projectDir}/bin/alignPE.sh \
        ${vs_report} \
        ${params.db}/*.aln \
        ${task.cpus} \
        ${reads[0].toRealPath()} \
        ${reads[1].toRealPath()} \
        ${sid} \
        ${meta}
    ## if alignPE.sh writes the serotype file already, keep it; otherwise derive a placeholder
    [[ -f ${sid}.serotype.txt ]] || awk 'NR==2 {print \$1}' ${vs_report} > ${sid}.serotype.txt
    """
    else
    """
    ${projectDir}/bin/alignSE.sh \
        ${vs_report} \
        ${params.db}/*.aln \
        ${task.cpus} \
        ${params.minimap_ext} \
        ${reads} \
        ${sid} \
        ${meta}
## if alignPE.sh writes the serotype file already, keep it; otherwise derive a placeholder
    [[ -f ${sid}.serotype.txt ]] || awk 'NR==2 {print \$1}' ${vs_report} > ${sid}.serotype.txt
    """
}