process CONCATREFS{

 input:
    path fasta_files

 output:
    path ("all.fasta"), emit: combinedrefs

script:
"""
cat ${params.referenceDir}/*.fasta > all.fasta
"""

}

process MAP_AND_SELECT_CONTIGS {
    cpus params.cpus
    tag "${sid}"
    memory '24 GB'
    publishDir "${params.outdir}/coinfection/01_contig_selection", mode: 'copy'

    input:
    tuple val(sid), path(reads)
    path ref_fasta

    output:
    tuple val(sid), path(reads), path("${sid}.allcontigs.bam"), path("${sid}.coverage.tsv"), path("${sid}.selected_contigs.tsv", optional: true)

    script:
    def minReads  = params.min_reads_coinf ?: 500
    def minCov    = params.coverageCoinfection ?: 30
    def minDepth  = params.meandepth ?: 5

    """
    set -euo pipefail

    minimap2 -ax ${params.minimap_ext} ${ref_fasta} ${reads} 2> ${sid}.minimap2.log \\
        | samtools sort -@ ${task.cpus} -o ${sid}.allcontigs.bam -

    samtools index ${sid}.allcontigs.bam

    samtools coverage -q 1 --min-depth 5 --ff UNMAP,QCFAIL ${sid}.allcontigs.bam > ${sid}.coverage.tsv

    awk -F'\\t' -v minReads=${minReads} -v minCov=${minCov} -v minDepth=${minDepth} '
        BEGIN { OFS="\\t" }
        NR==1 {
            print "sample_id","contig","numreads","coverage","meandepth"
            next
        }
        \$4 >= minReads && \$6 >= minCov && \$7 > minDepth {
            print "${sid}", \$1, \$4, \$6, \$7
        }
    ' ${sid}.coverage.tsv > ${sid}.selected_contigs.tmp.tsv

    # keep selected_contigs.tsv only if it has >2 data lines excluding header
    data_lines=\$(awk 'NR>1 && NF>0 {c++} END {print c+0}' ${sid}.selected_contigs.tmp.tsv)

    if [ "\$data_lines" -gt 2 ]; then
        mv ${sid}.selected_contigs.tmp.tsv ${sid}.selected_contigs.tsv
    else
        rm -f ${sid}.selected_contigs.tmp.tsv
    fi
    """
}

process REMAP_TO_SELECTED_CONTIG {
    cpus params.cpus
    tag "${sid}:${contig}"
    memory '24 GB'
    publishDir "${params.outdir}/coinfection/02_remap_per_contig", mode: 'copy'

    //conda "bioconda::seqkit bioconda::minimap2 bioconda::samtools"

    input:
    tuple val(sid), val(contig), path(reads)
    path ref_fasta

    output:
    tuple val(sid), val(contig), path("${sid}.${contig}.fa"), path("${sid}.${contig}.bam")

    script:
    """
    set -euo pipefail

    seqkit grep -nrp "${contig}" ${ref_fasta} > ${sid}.${contig}.fa
    samtools faidx ${sid}.${contig}.fa

    minimap2 -ax ${params.map_preset ?: 'lr:hq'} ${sid}.${contig}.fa ${reads} 2> ${sid}.${contig}.minimap2.log \
        | samtools sort -@ ${task.cpus} -o ${sid}.${contig}.bam -

    samtools index ${sid}.${contig}.bam
    """
}

process RUN_CLAIR3_PER_CONTIG {

    tag "${sid}:${contig}"
    cpus 20
    memory '48 GB'
    publishDir "${params.outdir}/coinfection/03_clair3", mode: 'copy'

    conda "bioconda::clair3=1.2.0 conda-forge::cudatoolkit conda-forge::cudnn conda-forge::python=3.10.0 conda-forge::numpy=1.26.4 conda-forge::tensorflow=2.14.* conda-forge::parallel"
    // Clair3 usually works better from container or dedicated env; adapt as needed

    input:
    tuple val(sid), val(contig), path(ref), path(bam)

    output:
    tuple val(sid), val(contig), path(ref), path(bam), path("${sid}.${contig}_clair3")

    script:
    """
    set -euo pipefail
    run_clair3.sh \\
        --bam_fn=${bam.toRealPath()} \\
        --ref_fn=${ref.toRealPath()} \\
        --threads=${task.cpus} \\
        --model_path=${params.clair3model} \\
        --output=${sid}.${contig}_clair3 \\
        ${params.clair3_ext ?: ''}

    """
}

process BUILD_CONTIG_CONSENSUS {

    tag "${sid}:${contig}"
    cpus 8
    memory '24 GB'
    publishDir "${params.outdir}/coinfection/04_consensus", mode: 'copy'

    //conda "bioconda::bcftools bioconda::samtools bioconda::tabix"

    input:
    tuple val(sid), val(contig), path(ref), path(bam), path(clair3dir)

    output:
    tuple val(sid), val(contig), path("${sid}.${contig}.fasta"), path("${sid}.${contig}.normalized.clair3.vcf.gz"), path("${sid}.${contig}.lowcov.bed"), path("${sid}.${contig}.clair3.log")

    script:
    def minDepth = params.meandepth ?: 5

    """
    set -euo pipefail
    covtobed -x ${params.coverageCoinfection} ${bam}  > ${sid}.${contig}.lowcov.bed

    bcftools sort ${clair3dir.toRealPath()}/merge_output.vcf.gz \
      | bcftools view -f 'PASS,.' \
      | bcftools norm --multiallelics -any --check-ref e --fasta-ref ${ref.toRealPath()} --old-rec-tag OLD_CLUMPED --atomize - \
      | bcftools norm --rm-dup exact --output-type z -o ${sid}.${contig}.normalized.clair3.vcf.gz

    tabix -p vcf ${sid}.${contig}.normalized.clair3.vcf.gz
    bcftools index ${sid}.${contig}.normalized.clair3.vcf.gz

    bcftools consensus \
        -f ${ref.toRealPath()} \
        --include 'FILTER="PASS"' \
        -p "${sid}_${contig} " \
        -m ${sid}.${contig}.lowcov.bed \
        --output ${sid}.${contig}.fasta \
        ${sid}.${contig}.normalized.clair3.vcf.gz

    cp ${clair3dir.toRealPath()}/run_clair3.log ${sid}.${contig}.clair3.log
    """
}