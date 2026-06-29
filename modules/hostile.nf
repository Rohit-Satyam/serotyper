params.memory = "3g"
params.cpus = 1
params.outdir = "."

process HOSTILE{
	cpus params.cpus
	memory params.memory
	publishDir "${params.outdir}/02_dehosting", mode: 'copy'

    input:
    tuple val(sid), path(reads)

    output:
		tuple val(sid), path("*clean*.fastq.gz")


        script:

    def hostile_ext = params.hostile_ext ? params.hostile_ext : ""

if ("${params.mode}" == "PE")
    """
hostile clean --fastq1 ${reads[0]} --fastq2 ${reads[1]} --threads ${task.cpus} \
    ${hostile_ext} \
    2> ${sid}.log
    """
else

"""
hostile clean --fastq1 ${reads} --threads ${task.cpus} \
    ${hostile_ext} \
    2> ${sid}.log
"""
}
