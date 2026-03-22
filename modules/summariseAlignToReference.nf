params.memory = "3g"
params.cpus = 1
params.outdir = "."


process SUMMARIZEALIGNTOREFERENCE{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/06_referenceAssembly", mode: 'copy'
    input:
    path(summaryFiles)
    output:
        path ("all_samples_summary_referencebased.tsv"), emit: summaryreport2
    shell:
'''
head -n 1 !{summaryFiles[0]}  > all_samples_summary_referencebased.tsv
tail -n +2 -q !{summaryFiles} >> all_samples_summary_referencebased.tsv
'''
}
