params.memory = "3g"
params.cpus = 1
params.outdir = "."

process YACRD{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/03_chimericReadRemoval", mode: 'copy'
  //conda "bioconda::yacrd bioconda::minimap2==2.18"
    input:
        tuple val(sid), path(reads)

    output:
        tuple val(sid), path("${sid}.scrubb.fastq.gz"), emit: scrubb
        path("*.txt")
         path "versions.yml", emit: versions

    script:
def scrubb_ext = params.scrubb_ext ? params.scrubb_ext : ''
def yacrd_ext = params.yacrd_ext ? params.yacrd_ext : ''
"""
## This has to be changed for pacbio and ONT
mm2plus -t ${task.cpus} ${scrubb_ext} ${reads} ${reads} > mapping.paf
yacrd -t ${task.cpus} -i mapping.paf -o ${sid}.yacrdsummary.txt ${yacrd_ext} scrubb -i ${reads} \
-o ${sid}.scrubb.fastq.gz

seqkit seq --threads ${task.cpus} --min-len 100 ${sid}.scrubb.fastq.gz  > ${sid}.scrubb.fastq
rm ${sid}.scrubb.fastq.gz
gzip ${sid}.scrubb.fastq

cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      yacrd: \$(yacrd --version 2>&1 | head -n 1 || true)
    END_VERSIONS
"""
}
