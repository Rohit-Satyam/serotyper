params.memory = "3g"
params.cpus = 1
params.outdir = "."

params.enable_conda = true


process CLAIR3{
  cpus params.cpus
  memory params.memory
  //publishDir "${params.outdir}/06_referenceAssembly/"
  conda "bioconda::clair3=1.2.0 conda-forge::cudatoolkit conda-forge::cudnn conda-forge::python=3.10.0 conda-forge::numpy=1.26.4 conda-forge::tensorflow=2.14.* conda-forge::parallel bioconda::bcftools bioconda::tabix"

  input:
      tuple val(sid), path(bam), path(ref), path(bed), path(aligntrimbed)

  output:
      tuple val(sid),  path("${ref}"), path("${bed}"), path("${sid}_clair3")

  script:

"""
## ensuring we have bam index
run_clair3.sh --bam_fn=${bam.toRealPath()} --ref_fn=${ref.toRealPath()} --threads=${task.cpus} \
--model_path=${params.clair3model} --output=${sid}_clair3 \
${params.clair3_ext}
"""
}

process MAKEREFBASEDASSEMBLY{
  cpus params.cpus
  memory params.memory
  publishDir "${params.outdir}/06_referenceAssembly/", mode: 'copy'


  input:
      tuple val(sid),  path(ref), path(bed), path(clair3dir)

  output:
  path("${sid}.normalized.clair3.vcf.gz")
  file("${sid}.clair3.log")
  path("${sid}.fasta")
  path "*.csi"
  path "*.tbi"

  script:

  """

  ## Breaking multialleles if any
  bcftools sort ${clair3dir.toRealPath()}/merge_output.vcf.gz | bcftools view -f 'PASS,.' |  bcftools norm --multiallelics \
  -any --check-ref e --fasta-ref ${ref.toRealPath()} --old-rec-tag OLD_CLUMPED --atomize - | bcftools norm --rm-dup exact --output-type z -o  ${sid}.normalized.clair3.vcf.gz

  tabix -p vcf ${sid}.normalized.clair3.vcf.gz
  bcftools index ${sid}.normalized.clair3.vcf.gz


  ## Building the Assembly, masking low coverage regions
  bcftools consensus -f ${ref.toRealPath()}  --include 'FILTER="PASS"'  -p "${sid} " -m ${bed.toRealPath()} --output ${sid}.fasta ${sid}.normalized.clair3.vcf.gz

  cp ${clair3dir.toRealPath()}/run_clair3.log ${sid}.clair3.log
  """

}
