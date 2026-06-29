params.memory = "3g"
params.cpus = 1
params.outdir = "."
import java.nio.file.Paths

process FASTQC{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/01_rawFastQC", mode: "copy"

    input:
        tuple val(sid), path(reads)

    output:
        path "*", emit: fastqc
def fastqc_ext = params.fastqc_ext ? params.fastqc_ext : ''
    """
##    fastqc -t ${task.cpus} $fastqc_ext ${reads[0]} ${reads[1]}
      fastqc -t ${task.cpus} $fastqc_ext ${reads}
    """
}

process MULTIQC {
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/${name}", mode: "copy"
    input:
    val(name)
    path(filepaths)
    val(filename)


    output:
    path "*"

    script:
    """
    multiqc --force --config ${projectDir}/bin/multiqc_config.yaml --filename ${filename} ${filepaths}
    """
}

process PYCOQC{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/01_longReadQC", mode: "copy"

    input:
        tuple val(sid), path(summaryfile)

    output:
        path "*"
    script:
    """
    pycoQC -f ${summaryfile.toRealPath()} -o ${sid}.pycoqc.html --min_pass_len 100

    NanoPlot --summary  ${summaryfile.toRealPath()} --loglength -o ${sid}.nanoplot_log_transformed
    """
}
