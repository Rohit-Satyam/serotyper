params.memory = "3g"
params.cpus = 1
params.outdir = "."

process BLASTN{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/${name}/", mode: 'copy'


    input:
        tuple val(sid), path(fasta)
        val(name)

    output:
    path("${sid}.blastN_results.tsv")

    shell:

  '''
blastn -query !{fasta} -num_threads !{task.cpus} -db !{params.blastdb} \
-out temp1 \
-outfmt '6 qseqid sseqid pident length evalue bitscore qcovs' \
-max_target_seqs 1

awk '{print $2}' temp1 | parallel -j 1 "grep {} !{params.mafftMeta}" > temp2

paste temp1 temp2 > !{sid}.blastN_results.tsv

'''
}
