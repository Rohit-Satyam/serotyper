params.memory = "3g"
params.cpus = 1
params.outdir = "."


process COMBINESEROTYPESUMMARY{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/04_serotype", mode: 'copy'
    input:
    path(summaryFiles)
    each path(meta)
    output:
        path ("all_samples_summary.tsv")
    shell:
'''
metadataHeader=$(head -n 1 !{meta} )
echo 'SampleName,Reference_Genome_Length,Horizontal_coverage_(>=20X),Base_Sequenced,'$metadataHeader | tr ', ' "\t" > all_samples_summary.tsv
cat !{summaryFiles} >> all_samples_summary.tsv
'''
}
