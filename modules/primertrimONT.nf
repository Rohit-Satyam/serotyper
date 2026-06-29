params.memory = "3g"
params.cpus = 1
params.outdir = "."

params.enable_conda = true


process ALIGNTRIM{
  cpus params.cpus
  memory params.memory
  publishDir "${params.outdir}/06_referenceAssembly/primerTrimmed_bams", mode: 'copy'
  conda "bioconda::align_trim conda-forge::csvkit conda-forge::xlsx2csv bioconda::covtobed bioconda::samtools"

  input:
      tuple val(sid), path(bam), path(ref), path(bed), path(aligntrimbed)

  output:
  tuple val(sid), path("${sid}_primertrimmed.bam"), path("${ref}"), path("${sid}_primertrimmedlowcoverage.bed"), path("${aligntrimbed}"), emit: trimmed
  path("${sid}_report.tsv"), emit: report
        
    shell:

  '''
# Read sheet name (first line) from the file aligntrimbed

SHEET=$(cat !{aligntrimbed})

# Convert Excel sheet to primer BED (tab-delimited)
xlsx2csv -n ${SHEET} !{params.primerfile} | csvformat -T > primer.bed

# Primer trimming
align_trim --samfile !{bam.toRealPath()} --report !{sid}_report.tsv --amp-depth-report !{sid}_depth_report.tsv --output trimmed.bam primer.bed !{params.aligntrim_ext}

## Recompute the coverage bed file for the resulting bam file
samtools sort trimmed.bam -@ !{task.cpus} -o !{sid}_primertrimmed.bam
samtools index !{sid}_primertrimmed.bam
covtobed -x !{params.cov} !{sid}_primertrimmed.bam  > !{sid}_primertrimmedlowcoverage.bed

 '''
}