params.memory = "3g"
params.cpus = 1
params.outdir = "."

params.enable_conda = true

process ALIGNTOREFERENCE{
  cpus params.cpus
  memory params.memory
  publishDir "${params.outdir}/06_referenceAssembly/", mode: 'copy'

  input:
      tuple val(sid), path(reads)
      each path(refdir)

  output:
  path("${sid}_alignment_results.tsv")
  tuple val("${sid}"), path("bams/${sid}_*bam"), path("bams/reference.fasta"), path("bams/${sid}_lowcoverage.bed"), path("aligntrim")


  script:

"""
referenceBasedSerotyping.sh ${sid} ${reads} ${refdir.toRealPath()} ${task.cpus} ${params.minimap_ext} ${params.cov}
"""
}
