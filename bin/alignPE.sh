#!/bin/bash
set -euo pipefail
## exit immediately on exception

## Get the Virstrain best hit name
## pipefail is disabled here because `head -n 2` intentionally closes its
## read end early, which otherwise sends SIGPIPE (exit 141) to `grep '>'`
## whenever the report lists more than two candidate clusters.
set +o pipefail
first=$(grep '>' $1 | head -n 2 | grep Cluster | awk '{print $1}' | sed 's/>//g')
set -o pipefail
seqkit grep -p $first $2 | sed '/^>/! s/-//g' | awk '/^>/ {print (NR==1?"":"\n")$0; next} {printf "%s", $0} END {print ""}' | fold -w 60 > ${6}.${first}.fasta ## Ensuring we don't mess up dash in the header


bwa-mem2 index ${6}.${first}.fasta
samtools faidx ${6}.${first}.fasta
endCoord=$(awk '{print $2}' ${6}.${first}.fasta.fai)
bwa-mem2 mem -t $3 ${6}.${first}.fasta  $4 $5 2> $first.bwa.output.log | samtools sort --threads $3 - |samtools view -F 4 --threads $3 -bS -o $6.bam

## -f 2 (properly paired) ensures both mates are always present so R1/R2
## read counts stay in sync; without it a mapped singleton whose mate is
## unmapped lands in only one of the two files and breaks downstream
## assemblers expecting matched pairs.
samtools fastq -1 ${6}.R1.fq -2 ${6}.R2.fq -f 2 -n $6.bam
gzip ${6}.R1.fq
gzip ${6}.R2.fq

samtools index --threads $3 $6.bam
samplot plot -n $6 -b $6.bam -o $6.${first}.png -c ${first} -s 1 -e ${endCoord} --coverage_only -H 3 --max_coverage 100
covtobed -m 0 $6.bam > $6.covtobed.txt

Ref_genome_length=${endCoord}
COVERED_LENGTH=$(awk '$4 >= 20 {covered += $3 - $2} END {print covered}' $6.covtobed.txt)
perc_baseCovered=$(echo "scale=2; ($COVERED_LENGTH * 100/ $Ref_genome_length) " | bc)

## Making summary file for each sample with following columns

serotype=$(awk -F'\t' -v k="$first" '$1==k' "$7") ## update unanchored serotype lookup to exact match via awk
echo "$6,$Ref_genome_length,$perc_baseCovered,$COVERED_LENGTH,$serotype" | tr ',' "\t" >> $6.serotype.txt
