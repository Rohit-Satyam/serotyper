params.memory = "3g"
params.cpus = 1
params.outdir = "."

process SUMMARIZEBLASTN{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/${name}", mode: 'copy'


    input:
        path(summaryFiles)
        val(name)
    output:
        path ("all_samples_summary.tsv")

    shell:

'''
metadataHeader=$(head -n 1 !{params.mafftMeta} )

echo "QueryID,MatchedReferenceID,perc_identity,QuerylengthMatch,Evalue,bitscore,queryCoverage",$metadataHeader | tr ', ' "\t" > temp

cat !{summaryFiles} >> temp

# Split the QueryID column, keeping only the third part of the split.
#awk 'BEGIN {FS=OFS="\t"} NR==1 {print; next} {split($1, arr, "_"); $1=arr[3]; print}'  temp | \

# Remove columns 3 to 7.
cut --complement -f3-7 temp | \

# Remove duplicate rows.
awk '!seen[$0]++' > all_samples_summary.tsv
'''
}
