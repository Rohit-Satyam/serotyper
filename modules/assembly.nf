params.memory = "3g"
params.cpus = 1
params.outdir = "."


process RNAVIRALSPADES{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/05_rnaviralSpadesAssembly/", mode: 'copy'


    input:
        tuple val(sid), path(reads)

    output:
        tuple val("${sid}"), path("${sid}.scaffolds.fasta")
        tuple val("${sid}"), path("${sid}.ragtag.scaffold.fasta")
        path("*.gfa")
        path("*.log")
        path("*.png")
        path("*.html")
        path("*.tsv")

    script:
    def spades_ext = params.spades_ext ? params.spades_ext : ""
    if ("${params.mode}" == "PE")
  """
  rnaviralspades.py -1 ${reads[0]} -2 ${reads[1]} --threads ${task.cpus} \
  -o . ${spades_ext}

  mv scaffolds.fasta ${sid}.scaffolds.fasta
  mv assembly_graph_with_scaffolds.gfa ${sid}.assembly_graph_with_scaffolds.gfa
  mv assembly_graph_after_simplification.gfa ${sid}.assembly_graph_after_simplification.gfa
  mv spades.log ${sid}.spades.log

  ragtag.py scaffold -t ${task.cpus} -o temp  ${reads[2]} ${sid}.scaffolds.fasta
  #sed -i "s/_RagTag/_RagTag_${sid}/g" ${sid}.ragtag.scaffold.fasta
  awk -v seq="${sid}" '/^>/ {print ">" seq "." ++i; next} {print}'  temp/ragtag.scaffold.fasta > ${sid}.ragtag.scaffold.fasta

  pgv-pmauve ${reads[2]} ${sid}.ragtag.scaffold.fasta \
   -o pgmauve --block_cmap viridis --track_align_type left  \
   --show_scale_xticks --curve

  mv pgmauve/result.png ${sid}.mauveresult.png
  mv pgmauve/align_coords.tsv ${sid}.mauvealign_coords.tsv
  mv pgmauve/result.html ${sid}.mauveresult.html
  mv pgmauve/pgv-cli.log ${sid}.mauve.log
  """
  else
  """
  rnaviralspades.py -s ${reads[0]} --nanopore --threads ${task.cpus}\
  -o . ${spades_ext}

  mv scaffolds.fasta ${sid}.scaffolds.fasta
  mv assembly_graph_with_scaffolds.gfa ${sid}.assembly_graph_with_scaffolds.gfa
  mv assembly_graph_after_simplification.gfa ${sid}.assembly_graph_after_simplification.gfa
  mv spades.log ${sid}.spades.log

  ragtag.py scaffold -t ${task.cpus} -o temp  ${reads[1]} ${sid}.scaffolds.fasta
  awk -v seq="${sid}" '/^>/ {print ">" seq "." ++i; next} {print}'  temp/ragtag.scaffold.fasta > ${sid}.ragtag.scaffold.fasta

  pgv-pmauve ${reads[1]} ${sid}.ragtag.scaffold.fasta \
   -o pgmauve --block_cmap viridis --track_align_type left  \
   --show_scale_xticks --curve

  mv pgmauve/result.png ${sid}.mauveresult.png
  mv pgmauve/align_coords.tsv ${sid}.mauvealign_coords.tsv
  mv pgmauve/result.html ${sid}.mauveresult.html
  mv pgmauve/pgv-cli.log ${sid}.mauve.log
  """
}

process UNICYCLER{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/05_UnicyclerAssembly", mode: 'copy'


    input:
        tuple val(sid), path(reads)

    output:
    tuple val("${sid}"), path("${sid}.scaffolds.fasta")
    tuple val("${sid}"), path("${sid}.ragtag.scaffold.fasta")
    path("*.gfa")
    path("*.log")
    path("*.png")
    path("*.html")
    path("*.tsv")

    script:
  if ("${params.mode}" == "PE")
"""
  unicycler -1 ${reads[0]} -2 ${reads[1]} -o . --linear_seqs 1 --keep 0
  mv assembly.fasta ${sid}.scaffolds.fasta
  mv assembly.gfa ${sid}.assembly_graph_with_scaffolds.gfa
  mv unicycler.log ${sid}.unicycler.log

  ragtag.py scaffold -t ${task.cpus} -o temp  ${reads[2]} ${sid}.scaffolds.fasta
  awk -v seq="${sid}" '/^>/ {print ">" seq "." ++i; next} {print}'  temp/ragtag.scaffold.fasta > ${sid}.ragtag.scaffold.fasta

  pgv-pmauve ${reads[2]} ${sid}.ragtag.scaffold.fasta \
   -o pgmauve --block_cmap viridis --track_align_type left  \
   --show_scale_xticks --curve

  mv pgmauve/result.png ${sid}.mauveresult.png
  mv pgmauve/align_coords.tsv ${sid}.mauvealign_coords.tsv
  mv pgmauve/result.html ${sid}.mauveresult.html
  mv pgmauve/pgv-cli.log ${sid}.mauve.log
"""
else
"""
unicycler -l ${reads[0]} -o . --linear_seqs 1 --keep 0
mv assembly.fasta ${sid}.scaffolds.fasta
mv assembly.gfa ${sid}.assembly_graph_with_scaffolds.gfa
mv unicycler.log ${sid}.unicycler.log

ragtag.py scaffold -t ${task.cpus} -o temp  ${reads[1]} ${sid}.scaffolds.fasta
awk -v seq="${sid}" '/^>/ {print ">" seq "." ++i; next} {print}'  temp/ragtag.scaffold.fasta > ${sid}.ragtag.scaffold.fasta
"""
}

process FLYE{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/05_FlyeAssembly", mode: 'copy'


    input:
        tuple val(sid), path(reads)

    output:
    tuple val("${sid}"), path("${sid}.scaffolds.fasta")
    tuple val("${sid}"), path("${sid}.racon_ragtag_scaffold.fasta")
    path("*.gfa")
    path("*.log")
    path("*.png")
    path("*.html")
    path("*.tsv")

    script:
    def flye_ext = params.flye_ext ? params.flye_ext : ""
"""
flye --threads ${task.cpus} --genome-size ${params.genomesize} --out-dir . --scaffold  \
${flye_ext} ${reads[1]}

## One time polishing with Racon
minimap2 assembly.fasta ${reads[1]} \
 > minimap.racon.paf

racon -t ${task.cpus} ${reads[1]} \
minimap.racon.paf assembly.fasta \
> ${sid}.racon.consensus.fasta

mv assembly.fasta ${sid}.scaffolds.fasta
mv assembly_graph.gfa ${sid}.assembly_graph_with_scaffolds.gfa
mv flye.log ${sid}.flye.log

ragtag.py scaffold -t ${task.cpus} -o temp  ${reads[0]} ${sid}.racon.consensus.fasta
awk -v seq="${sid}" '/^>/ {print ">" seq "." ++i; next} {print}'  temp/ragtag.scaffold.fasta > ${sid}.racon_ragtag_scaffold.fasta

pgv-pmauve ${reads[0]} ${sid}.racon_ragtag_scaffold.fasta \
 -o pgmauve --block_cmap viridis --track_align_type left  \
 --show_scale_xticks --curve

mv pgmauve/result.png ${sid}.mauveresult.png
mv pgmauve/align_coords.tsv ${sid}.mauvealign_coords.tsv
mv pgmauve/result.html ${sid}.mauveresult.html
mv pgmauve/pgv-cli.log ${sid}.mauve.log

"""

}

process CANU{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/05_CanuAssembly", mode: 'copy'


    input:
        tuple val(sid), path(reads)

    output:
    tuple val("${sid}"), path("${sid}.contigs.fasta")
    tuple val("${sid}"), path("${sid}.racon_ragtag_scaffold.fasta")
    path("*.err")
    path("*.log")
    path("*.png")
    path("*.html")
    path("*.tsv")

    script:
    def canu_ext = params.canu_ext ? params.canu_ext : ""
"""

## Get the geometric mean for minimum read length because the FMDV data has read
## length between 250-750 bp and canu filters out anything below 1000bps

len=\$(${projectDir}/bin/getGeoLength.sh ${reads[1]} ${task.cpus})

canu -p ${sid} minReadLength=\$len minOverlapLength=\$len genomeSize=${params.genomesize} ${canu_ext} ${reads[1]}

## One time polishing with Racon
minimap2 ${sid}.contigs.fasta ${reads[1]} \
 > minimap.racon.paf

racon -t ${task.cpus} ${reads[1]} \
minimap.racon.paf ${sid}.contigs.fasta \
> ${sid}.racon.consensus.fasta


mv ${sid}.report ${sid}.canu.log
mv ${sid}.seqStore.err ${sid}.err


ragtag.py scaffold -t ${task.cpus} -o temp  ${reads[0]} ${sid}.racon.consensus.fasta
awk -v seq="${sid}" '/^>/ {print ">" seq "." ++i; next} {print}'  temp/ragtag.scaffold.fasta > ${sid}.racon_ragtag_scaffold.fasta

pgv-pmauve ${reads[0]} ${sid}.racon_ragtag_scaffold.fasta \
 -o pgmauve --block_cmap viridis --track_align_type left  \
 --show_scale_xticks --curve

mv pgmauve/result.png ${sid}.mauveresult.png
mv pgmauve/align_coords.tsv ${sid}.mauvealign_coords.tsv
mv pgmauve/result.html ${sid}.mauveresult.html
mv pgmauve/pgv-cli.log ${sid}.mauve.log
"""

}
