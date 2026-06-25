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
        path "${sid}.log"
        path "versions.yml", emit: versions

        script:

    def hostile_ext = params.hostile_ext ? params.hostile_ext : ""

if ("${params.mode}" == "PE")
    """
hostile clean --fastq1 ${reads[0]} --fastq2 ${reads[1]} --threads ${task.cpus} \
    ${hostile_ext} \
    2> ${sid}.log

     cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      hostile: \$(hostile --version 2>&1 | head -n 1 || true)
    END_VERSIONS
    """
else

"""
hostile clean --fastq1 ${reads} --threads ${task.cpus} \
    ${hostile_ext} \
    2> ${sid}.log
 cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      hostile: \$(hostile --version 2>&1 | head -n 1 || true)
    END_VERSIONS
"""
}
