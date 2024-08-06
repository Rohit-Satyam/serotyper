params.memory = "3g"
params.cpus = 1
params.outdir = "."

params.enable_conda = true
process SPADES{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/04_assembly", mode: 'copy'


    input:
        tuple val(sid), path(reads)
        each path(pathdb)

    output:
        path "*.csv"
        path("*.txt")
        path("*.html")

    script:
    def spades_ext = params.spades_ext ? params.spades_ext : ""
if ("${params.mode}" == "PE")
  """
  rnaviralspades.py -1 ${reads[0]} -2 ${reads[1]} --threads ${task.cpus}
  """
else
  """
  rnaviralspades.py -1 ${reads} --threads ${task.cpus}
  """
}
