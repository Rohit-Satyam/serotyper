import java.nio.file.Paths
params.memory = "3g"
params.cpus = 1
params.outdir = "."


process FASTP{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/03_adapterTrimming", mode: 'copy'

    input:
        tuple val(sid), path(reads)

    output:
        tuple val(sid), path("*.fastq.gz")
        path("${sid}.fastp.json"), emit: fastp_logs
        path("${sid}.fastp.html")
        path "*.log"
        path "*.txt"
        path "versions.yml", emit: versions

    script:
    fq_1_paired = sid + '_trim_R1.fastq.gz'
    fq_2_paired = sid + '_trim_R2.fastq.gz'

    def fastp_ext = params.fastp_ext ? params.fastp_ext : ""

    def fastp_input = ""
    def trimmed_file = ""

    if ("${params.mode}" == "PE") {
        fastp_input = "--in1 ${reads[0]} --in2 ${reads[1]} --out1 ${fq_1_paired} --out2 ${fq_2_paired}"
        trimmed_file = fq_1_paired
    }
    else if ("${params.mode}" == "SE") {
        fastp_input = "--in1 ${reads} -o ${sid}.trim.fastq.gz"
        trimmed_file = "${sid}.trim.fastq.gz"
    }

    """
    fastp \
        ${fastp_input} \
        --thread ${task.cpus} \
        --json ${sid}.fastp.json \
        --html ${sid}.fastp.html \
        ${fastp_ext} \
        2> ${sid}.log

    raw=\$(zcat ${params.mode == "PE" ? reads[0] : reads} | wc -l | awk '{print \$1/4}')
    trimmed=\$(zcat ${trimmed_file} | wc -l | awk '{print \$1/4}')
    echo ${sid},\$raw,\$trimmed > ${sid}_fastp.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      fastp: \$(fastp --version 2>&1 | head -n 1)
    END_VERSIONS
    """
}

process POSTTRIMFASTQC{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/02_adapterTrimming/postTrimFASTQC", mode: 'copy'

    input:
        tuple val(sid), path(reads)

    output:
        path "*", emit: postfastqc
        path "versions.yml", emit: versions
    script:
    def fastqc_ext = params.fastqc_ext ? params.fastqc_ext : ''
    """
    fastqc -t ${task.cpus} $fastqc_ext ${reads}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      fastqc: \$(fastqc --version 2>&1 | head -n 1)
    END_VERSIONS
    """
}
