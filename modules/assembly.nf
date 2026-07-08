
params.memory = "3g"
params.cpus = 1
params.outdir = "."
process RNAVIRALSPADES{
    publishDir "${params.outdir}/05_rnaviralSpadesAssembly/", mode: 'copy'


    input:
        tuple val(sid), path(reads)

    output:
        tuple val("${sid}"), path("${sid}.scaffolds.fasta")
        tuple val("${sid}"), path("${sid}.ragtag.scaffold.fasta")
        path("*.gfa"), optional: true
        path("*.log"), optional: true
        path("*.png"), optional: true
        path("*.html"), optional: true
        path("*.tsv"), optional: true
        path "versions.yml", emit: versions

    script:
   def spades_ext = params.spades_ext ? params.spades_ext : ""
   def assemble_cmd = params.mode == 'PE'
        ? """
          rnaviralspades.py -1 ${reads[0]} -2 ${reads[1]} --threads ${task.cpus} -o . ${spades_ext}
          """
        : """
          rnaviralspades.py -s ${reads[0]} --nanopore --threads ${task.cpus} -o . ${spades_ext}
          """
  
  """
    ${assemble_cmd}
    mv scaffolds.fasta ${sid}.scaffolds.fasta
    [ -f assembly_graph_with_scaffolds.gfa ] && mv assembly_graph_with_scaffolds.gfa ${sid}.assembly_graph_with_scaffolds.gfa
    [ -f assembly_graph_after_simplification.gfa ] && mv assembly_graph_after_simplification.gfa ${sid}.assembly_graph_after_simplification.gfa
    [ -f spades.log ] && mv spades.log ${sid}.spades.log
    
    ragtag.py scaffold -t ${task.cpus} -o temp  ${reads[2]} ${sid}.scaffolds.fasta
    awk -v seq="${sid}" '/^>/ {print ">" seq "." ++i; next} {print}'  temp/ragtag.scaffold.fasta > ${sid}.ragtag.scaffold.fasta
    
    ## progressiveMauve (invoked by pgv-pmauve) is prone to segfaulting on
    ## some systems; tolerate its failure so a valid assembly is still
    ## published even when the comparison plot cannot be generated.
    pgv-pmauve ${reads[2]} ${sid}.ragtag.scaffold.fasta \
   -o pgmauve --block_cmap viridis --track_align_type left  \
   --show_scale_xticks --curve || echo "WARNING: pgv-pmauve failed; skipping mauve visualization" >&2

    [ -f pgmauve/result.png ]        && mv pgmauve/result.png ${sid}.mauveresult.png
    [ -f pgmauve/align_coords.tsv ]  && mv pgmauve/align_coords.tsv ${sid}.mauvealign_coords.tsv
    [ -f pgmauve/result.html ]       && mv pgmauve/result.html ${sid}.mauveresult.html
    [ -f pgmauve/pgv-cli.log ]       && mv pgmauve/pgv-cli.log ${sid}.mauve.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      rnaviralspades: \$(rnaviralspades.py --version 2>&1 | sed 's/.*SPAdes genome assembler v//; s/ .*//')
      ragtag: \$(ragtag.py --version 2>&1 | head -n 1 | sed 's/^RagTag //')
      pgv_pmauve: \$(pgv-pmauve --version 2>&1 | head -n 1 | sed 's/^v//')
    END_VERSIONS
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
    path("*.gfa"), optional: true
    path("*.log"), optional: true
    path("*.png"), optional: true
    path("*.html"), optional: true
    path("*.tsv"), optional: true
    path "versions.yml", emit: versions

    script:
    def assemble_cmd = params.mode == 'PE'
        ? "unicycler -1 ${reads[0]} -2 ${reads[1]} -o . --linear_seqs 1 --keep 0"
        : "unicycler -l ${reads[0]} -o . --linear_seqs 1 --keep 0"

    def downstream = params.mode == 'PE'
        ? """
          ragtag.py scaffold -t ${task.cpus} -o temp  ${reads[2]} ${sid}.scaffolds.fasta

          awk -v seq="${sid}" '/^>/ {print ">" seq "." ++i; next} {print}'  temp/ragtag.scaffold.fasta > ${sid}.ragtag.scaffold.fasta

          ## progressiveMauve (invoked by pgv-pmauve) is prone to segfaulting
          ## on some systems; tolerate its failure so a valid assembly is
          ## still published even when the comparison plot cannot be generated.
          pgv-pmauve ${reads[2]} ${sid}.ragtag.scaffold.fasta \\
          -o pgmauve --block_cmap viridis --track_align_type left  \\
          --show_scale_xticks --curve || echo "WARNING: pgv-pmauve failed; skipping mauve visualization" >&2

        [ -f pgmauve/result.png ]       && mv pgmauve/result.png ${sid}.mauveresult.png
        [ -f pgmauve/align_coords.tsv ] && mv pgmauve/align_coords.tsv ${sid}.mauvealign_coords.tsv
        [ -f pgmauve/result.html ]      && mv pgmauve/result.html ${sid}.mauveresult.html
        [ -f pgmauve/pgv-cli.log ]      && mv pgmauve/pgv-cli.log ${sid}.mauve.log
          """
        : """
          ragtag.py scaffold -t ${task.cpus} -o temp  ${reads[1]} ${sid}.scaffolds.fasta
          awk -v seq="${sid}" '/^>/ {print ">" seq "." ++i; next} {print}'  temp/ragtag.scaffold.fasta > ${sid}.ragtag.scaffold.fasta
          """
"""
${assemble_cmd}

mv assembly.fasta ${sid}.scaffolds.fasta
[ -f assembly.gfa ]    && mv assembly.gfa ${sid}.assembly_graph_with_scaffolds.gfa
[ -f unicycler.log ]   && mv unicycler.log ${sid}.unicycler.log

${downstream}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      unicycler: \$(unicycler --version 2>&1 | sed 's/^Unicycler v//')
      ragtag: \$(ragtag.py --version 2>&1 | head -n 1 | sed 's/^RagTag //')
      pgv_pmauve: \$(pgv-pmauve --version 2>&1 | head -n 1 | sed 's/^v//')
    END_VERSIONS
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
    path("*.gfa"), optional: true
    path("*.log"), optional: true
    path("*.png"), optional: true
    path("*.html"), optional: true
    path("*.tsv"), optional: true
    path "versions.yml", emit: versions

    script:
    def flye_ext = params.flye_ext ? params.flye_ext : ""
"""
flye --threads ${task.cpus} --genome-size ${params.genomesize} --out-dir . --scaffold  \
${flye_ext} ${reads[1]}

## One time polishing with Racon
mm2plus assembly.fasta ${reads[1]} \
 > minimap.racon.paf

racon -t ${task.cpus} ${reads[1]} \
minimap.racon.paf assembly.fasta \
> ${sid}.racon.consensus.fasta

mv assembly.fasta ${sid}.scaffolds.fasta
mv assembly_graph.gfa ${sid}.assembly_graph_with_scaffolds.gfa
mv flye.log ${sid}.flye.log

ragtag.py scaffold -t ${task.cpus} -o temp  ${reads[0]} ${sid}.racon.consensus.fasta
awk -v seq="${sid}" '/^>/ {print ">" seq "." ++i; next} {print}'  temp/ragtag.scaffold.fasta > ${sid}.racon_ragtag_scaffold.fasta

## progressiveMauve (invoked by pgv-pmauve) is prone to segfaulting on some
## systems; tolerate its failure so a valid assembly is still published even
## when the comparison plot cannot be generated.
pgv-pmauve ${reads[0]} ${sid}.racon_ragtag_scaffold.fasta \
 -o pgmauve --block_cmap viridis --track_align_type left  \
 --show_scale_xticks --curve || echo "WARNING: pgv-pmauve failed; skipping mauve visualization" >&2

[ -f pgmauve/result.png ]       && mv pgmauve/result.png ${sid}.mauveresult.png
[ -f pgmauve/align_coords.tsv ] && mv pgmauve/align_coords.tsv ${sid}.mauvealign_coords.tsv
[ -f pgmauve/result.html ]      && mv pgmauve/result.html ${sid}.mauveresult.html
[ -f pgmauve/pgv-cli.log ]      && mv pgmauve/pgv-cli.log ${sid}.mauve.log

cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      flye: \$(flye --version 2>&1 | head -n 1 | sed 's/^.* //')
      mm2plus: \$(mm2plus --version 2>&1 | head -n 1 || echo "unknown")
      racon: \$(racon --version 2>&1 | head -n 1 | sed 's/^v//')
      ragtag: \$(ragtag.py --version 2>&1 | head -n 1 | sed 's/^RagTag //')
      pgv_pmauve: \$(pgv-pmauve --version 2>&1 | head -n 1 | sed 's/^v//')
    END_VERSIONS
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
    path("*.err"), optional: true
    path("*.log"), optional: true
    path("*.png"), optional: true
    path("*.html"), optional: true
    path("*.tsv"), optional: true
    path "versions.yml", emit: versions

    script:
    def canu_ext = params.canu_ext ? params.canu_ext : ""
"""

## Get the geometric mean for minimum read length because the FMDV data has read
## length between 250-750 bp and canu filters out anything below 1000bps

len=\$(${projectDir}/bin/getGeoLength.sh ${reads[1]} ${task.cpus})

canu -p ${sid} minReadLength=\$len minOverlapLength=\$len genomeSize=${params.genomesize} ${canu_ext} ${reads[1]}

## One time polishing with Racon
mm2plus ${sid}.contigs.fasta ${reads[1]} \
 > minimap.racon.paf

racon -t ${task.cpus} ${reads[1]} \
minimap.racon.paf ${sid}.contigs.fasta \
> ${sid}.racon.consensus.fasta


mv ${sid}.report ${sid}.canu.log
mv ${sid}.seqStore.err ${sid}.err


ragtag.py scaffold -t ${task.cpus} -o temp  ${reads[0]} ${sid}.racon.consensus.fasta
awk -v seq="${sid}" '/^>/ {print ">" seq "." ++i; next} {print}'  temp/ragtag.scaffold.fasta > ${sid}.racon_ragtag_scaffold.fasta

## progressiveMauve (invoked by pgv-pmauve) is prone to segfaulting on some
## systems; tolerate its failure so a valid assembly is still published even
## when the comparison plot cannot be generated.
pgv-pmauve ${reads[0]} ${sid}.racon_ragtag_scaffold.fasta \
 -o pgmauve --block_cmap viridis --track_align_type left  \
 --show_scale_xticks --curve || echo "WARNING: pgv-pmauve failed; skipping mauve visualization" >&2

[ -f pgmauve/result.png ]       && mv pgmauve/result.png ${sid}.mauveresult.png
[ -f pgmauve/align_coords.tsv ] && mv pgmauve/align_coords.tsv ${sid}.mauvealign_coords.tsv
[ -f pgmauve/result.html ]      && mv pgmauve/result.html ${sid}.mauveresult.html
[ -f pgmauve/pgv-cli.log ]      && mv pgmauve/pgv-cli.log ${sid}.mauve.log

 cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      canu: \$(canu --version 2>&1 | head -n 1 | sed 's/^Canu //')
      mm2plus: \$(mm2plus --version 2>&1 | head -n 1 || echo "unknown")
      racon: \$(racon --version 2>&1 | head -n 1 | sed 's/^v//')
      ragtag: \$(ragtag.py --version 2>&1 | head -n 1 | sed 's/^RagTag //')
      pgv_pmauve: \$(pgv-pmauve --version 2>&1 | head -n 1 | sed 's/^v//')
    END_VERSIONS
"""

}
