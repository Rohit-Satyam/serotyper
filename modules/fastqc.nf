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
    path "versions.yml", emit: versions
    script:
    """
    multiqc --force --config ${projectDir}/bin/multiqc_config.yaml --filename ${filename} ${filepaths}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      multiqc: \$(multiqc --version 2>&1 | head -n 1)
    END_VERSIONS
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
        path "versions.yml", emit: versions
    script:
    """
    pycoQC -f ${summaryfile.toRealPath()} -o ${sid}.pycoqc.html --min_pass_len 100

    NanoPlot --summary  ${summaryfile.toRealPath()} --loglength -o ${sid}.nanoplot_log_transformed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      pycoqc: \$(pycoQC --version 2>&1 | head -n 1 || true)
      nanoplot: \$(NanoPlot --version 2>&1 | head -n 1 || true)
    END_VERSIONS
    """
}
