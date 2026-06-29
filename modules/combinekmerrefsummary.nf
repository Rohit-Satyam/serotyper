params.memory = "3g"
params.cpus = 1
params.outdir = "."


process COMBINEKMERREFSUMMARIES{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}", mode: 'copy'
    input:
        path(report1)
        path(report2)
    output:
        path ("Serotyper_report.tsv")
    shell:
"""
combineKmerRefSummaries.sh ${report1.toRealPath()} ${report2.toRealPath()}
"""
}
