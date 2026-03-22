#!/bin/bash

# Input variables
sample_name=${1}  # Take sample name
fastq_file=${2}  # Your FASTQ file for alignment
referenceDir=${3}
cpus=${4}
minimap_ext=${5}
output_file="${1}_alignment_results.tsv"  # Output TSV file
coverage_threshold=${6}

## We will keep the best hit reference here and the bam file
mkdir bams


# To check if the referenceDir ends with .fa, .fasta, or .fas or user provided a directory with mulifasta
if [ -d "${referenceDir}" ]; then
    echo "The path is a directory. Saving all the fasta file to tmp"
    ls -1 ${referenceDir}/*.{fasta,fa,fas} > tmp 2>/dev/null
elif echo "${referenceDir}" | awk '/\.(fa|fasta|fas)$/ {exit 0} {exit 1}'; then
    echo "The input is a single fasta file"
    ls -1 ${referenceDir} > tmp
else
    echo "Invalid Path"
fi


# Temporary file to store read counts
temp_file="read_counts_temp.txt"
> $temp_file  # Clear the temporary file before use

# Perform alignment for each reference FASTA in the referenceDir folder
cat tmp | while read p; do
    # Extract reference fasta basename name (e.g., denv1, denv2)
    ref_name=$(basename "$p" | cut -f 1 -d'.')

    # Align the FASTQ sample against the current reference using minimap2 and store the BAM file
    minimap2 -ax ${minimap_ext} "$p" "$fastq_file" 2> "${ref_name}_minimap.output.log" \
    | samtools sort -@ ${cpus} -o "${sample_name}_${ref_name}.bam"

    samtools index -@ ${cpus} "${sample_name}_${ref_name}.bam" ## Because clair3 requires it.

    # Count the number of mapped reads in the BAM file. Flag 260 filters unmapped reads and non-primary alignments
    mapped_reads=$(samtools view -c -F 260 "${sample_name}_${ref_name}.bam")

    # Save the result to the temp file (Reference name and count of mapped reads)
    echo -e "$ref_name\t$mapped_reads" >> $temp_file

    # log the count for this reference (for debugging purposes)
    echo "Reference: $ref_name, Mapped reads: $mapped_reads"
done

# Prepare the TSV file with header
echo -e "Sample\t$(cat tmp | xargs -n 1 basename | cut -f 1 -d'.' | tr '\n' '\t')Most_Mapped_Reference" > "$output_file"

# Write the sample name as the first column in the TSV file
echo -n -e "$sample_name\t" >> "$output_file"

# Append the counts for each reference in the temp file to the TSV file
awk '{print $2}' $temp_file | tr '\n' '\t' >> "$output_file"

seqkit stats $fastq_file -T --basename | cut -f4,6,7,8 > seqkitstats

# Use sort to find the reference with the highest mapped reads
most_mapped_ref=$(sort -k2 -nr $temp_file | head -n1 | cut -f1)

# Get best reference for variant calling

best=$(cat tmp | grep $most_mapped_ref)
bestalignedbam=$(ls -1 *.bam | grep $most_mapped_ref)
covtobed -x ${coverage_threshold} ${bestalignedbam}  > bams/${sample_name}_lowcoverage.bed
cp ${best} bams/reference.fasta
cp ${bestalignedbam}* bams/
samtools faidx bams/reference.fasta    ## Because clair3 requires it.


# Append the reference with the most mapped reads to the bams column of the TSV
echo -e "$most_mapped_ref" >> "$output_file"

paste "$output_file" seqkitstats > temp.tt && mv temp.tt "$output_file"

## Echoing the name of the best match in a separate file to be used for alignment trimming
echo -e "$most_mapped_ref" >> aligntrim

# Clean up temporary file
rm $temp_file tmp

echo "Alignment complete. Results saved to $output_file."
