params.memory = "3g"
params.cpus = 10
params.outdir = "."

params.enable_conda = true
process VIRSTRAIN{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/04_serotype", mode: 'copy'
    //conda (params.enable_conda ? "bioconda::virstrain" : null)
    conda "bioconda::virstrain bioconda::samtools bioconda::seqkit bioconda::minimap2 bioconda::bwa-mem2 bioconda::samplot bioconda::covtobed"
 
    input:
        tuple val(sid), path(reads)
        each path(pathdb)
	each path(meta)

    output:
        path("${sid}.serotype.txt"), emit: serotyper_res
        tuple val("${sid}"), path("${sid}*.fasta"), emit:serotype_contig_ordering
        tuple val("${sid}"), path("${sid}*.fq.gz"), emit:viral_reads
        path("*.csv")
        path("*.txt")
        path("*.html")
        path("*.png")
        path("*.bam")
        path("*.bai")

    script:
  if ("${params.mode}" == "PE")
  """
  virstrain -i ${reads[0]} -p ${reads[1]} -d ${pathdb} -o .

  ## Getting the best match FASTA file for coverage estimation

  ${projectDir}/bin/alignPE.sh VirStrain_report.txt ${pathdb.toRealPath()}/*.aln ${task.cpus} \
  ${reads[0].toRealPath()} ${reads[1].toRealPath()} ${sid} ${meta}


  mv Mps_ps_depth.csv ${sid}.Mps_ps_depth.csv
  mv Ops_ps_depth.csv ${sid}.Ops_ps_depth.csv
  mv VirStrain_report.html ${sid}.VirStrain_report.html
  mv VirStrain_report.txt ${sid}.VirStrain_report.txt
  """
  else if ("${params.mode}" == "SE")

  """
  virstrain -i ${reads} -d ${pathdb} -o .

  ## Getting the best match FASTA file for coverage estimation

  ${projectDir}/bin/alignSE.sh VirStrain_report.txt ${pathdb.toRealPath()}/*.aln ${task.cpus} \
  ${params.minimap_ext} ${reads}  ${sid} ${meta}


  mv Mps_ps_depth.csv ${sid}.Mps_ps_depth.csv
  mv Ops_ps_depth.csv ${sid}.Ops_ps_depth.csv
  mv VirStrain_report.html ${sid}.VirStrain_report.html
  mv VirStrain_report.txt ${sid}.VirStrain_report.txt
  """
}
