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


        script:
    fq_1_paired = sid + '_trim_R1.fastq.gz'
    fq_2_paired = sid + '_trim_R2.fastq.gz'

    def fastp_ext = params.fastp_ext ? params.fastp_ext : ""

if ("${params.mode}" == "PE")
    """
	fastp \
	--in1 ${reads[0]} \
	--in2 ${reads[1]}\
  --thread ${task.cpus} \
	--out1 $fq_1_paired \
	--out2 $fq_2_paired \
	--json ${sid}.fastp.json \
	--html ${sid}.fastp.html \
    ${fastp_ext} \
    2> ${sid}.log

    raw=\$(zcat ${reads[0]}  | wc -l | awk '{print \$1/4}')
    trimmed=\$(zcat $fq_1_paired | wc -l | awk '{print \$1/4}')
    echo ${sid},\$raw,\$trimmed > ${sid}_fastp.txt
    """
else if ("${params.mode}" == "SE")

"""
	fastp \
	--in1 ${reads} \
  --thread ${task.cpus} \
	-o ${sid}.trim.fastq.gz \
	--json ${sid}.fastp.json \
	--html ${sid}.fastp.html \
    2> ${sid}.log

    raw=\$(zcat ${reads}  | wc -l | awk '{print \$1/4}')
    trimmed=\$(zcat ${sid}.trim.fastq.gz | wc -l | awk '{print \$1/4}')
    echo ${sid},\$raw,\$trimmed > ${sid}_fastp.txt

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
    script:
    def fastqc_ext = params.fastqc_ext ? params.fastqc_ext : ''
    """
    fastqc -t ${task.cpus} $fastqc_ext ${reads}
    """
}
