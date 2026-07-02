#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
/*
 * ============================================================================
 * SeroTyper main workflow
 * ============================================================================
 * Purpose:
 *   Identify viral serotypes from Illumina PE or ONT SE reads using a mix of
 *   k-mer / strain matching, optional de novo assembly, optional reference-
 *   based assembly, and optional coinfection analysis.
 *
 * Design principles:
 *   - Fail early on invalid configuration / missing required parameters.
 *   - Tolerate expected per-sample failures caused by low-depth clinical data.
 *   - Keep channel contracts explicit and workflow branching readable.
 * ============================================================================
 */


// ----------------------------------------------------------------------------
// Help / usage
// ----------------------------------------------------------------------------

if( params.help ) {

    // Define colors using ANSI escape codes
    def RED = '\033[0;31m'
    def GREEN = '\033[0;32m'
    def YELLOW = '\033[0;33m'
    def BLUE = '\033[0;34m'
    def MAGENTA = '\033[0;35m'
    def CYAN = '\033[0;36m'
    def NC = '\033[0m' // No Color

    log.info """${GREEN}
* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------${NC} ${GREEN}

##
##                                           .                                               .oooooo.  oooo    oooo       .o.       ooooo     ooo  .oooooo..o ooooooooooooo
##                                         .o8                                              d'     `b  `888   .8P'       .888.      `888'     `8' d8P'    `Y8 8'   888   `8
##  .oooo.o  .ooooo.  oooo d8b  .ooooo.  .o888oo oooo    ooo oo.ooooo.   .ooooo.  oooo d8b d' .d"bd  8  888  d8'        .8"888.      888       8  Y88bo.           888
## d88(  "8 d88' `88b `888""8P d88' `88b   888    `88.  .8'   888' `88b d88' `88b `888""8P 8  8. 8  .d  88888[         .8' `888.     888       8   `"Y8888o.       888
## `"Y88b.  888ooo888  888     888   888   888     `88..8'    888   888 888ooo888  888     Y.  YoP"b'   888`88b.      .88ooo8888.    888       8       `"Y88b      888
## o.  )88b 888    .o  888     888   888   888 .    `888'     888   888 888    .o  888      8.      .8  888  `88b.   .8'     `888.   `88.    .8'  oo     .d8P      888
## 8""888P' `Y8bod8P' d888b    `Y8bod8P'   "888"     .8'      888bod8P' `Y8bod8P' d888b      YooooooP  o888o  o888o o88o     o8888o    `YbodP'    8""88888P'      o888o
##                                               .o..P'       888
##                                               `Y8P'       o888o

								Identifying Viral Serotypes using short / long reads						${NC}
${GREEN}
* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------${NC}
${CYAN}Usage:${NC}
        ${YELLOW}nextflow run main.nf --input "${params.input}" --outdir ${params.outdir}${NC}
${CYAN}Input:${NC}
        ${RED}#### Mandatory Arguments ####${NC}
    ${MAGENTA}* --input:${NC} Path to FastQ files. Default [${params.input}]
    ${MAGENTA}* --mode:${NC} If data is Paired-end pass "PE" else "SE". Only Illumina PE data and ONT nanopore data is supported. Default [${params.mode}]
	${MAGENTA}* --genomesize:${NC} Approximate size of the genome. This argument is used by assemblers to auto-optimize computational resources. Default [${params.genomesize}]
	${MAGENTA}* --mafftMeta:${NC} Absolute path to metadata file for serotyping. It should have 2 mandatory columns: Col1: Sequence accessions in the MSA and Col2: Serotype. Default [${params.mafftMeta}]
    ${MAGENTA}* --referenceDir:${NC} Absolute path to directory containing fasta file of all serotypes. Mandatory if you want to run reference-based genome assembly. Default [${params.referenceDir}]
    ${MAGENTA}* --db:${NC} Absolute path to directory containing Virstrain database. Mandatory if you want to run assembly free serotyping quickly. Default [${params.db}]
    ${MAGENTA}* --clair3model:${NC} Absolute path to directory containing clair3 model for Variant calling step. Mandatory if you want to run reference-based genome assembly. This depends on basecaller and sequencing kit used.
    Download appropriate model from https://github.com/HKU-BAL/clair3#pre-trained-models. Default [${params.clair3model}]
	${MAGENTA}* --assembler:${NC} Assembler to use for viral genome assembly.
			Possible options for Illumina PR: "unicycler", "spades" or "all".
			Possible options for ONT: "canu", "flye" or "all". Default [${params.assembler}]

	${RED}#### Optional Arguments${NC}
	${MAGENTA}* --summaryFile:${NC} ONT sequencing summary file to run pycoQC. Default [${params.summaryFile}]
	${MAGENTA}* --blastdb:${NC} Path of BlastN indexes with the prefix. Provide this if skipBlast is set to FALSE. Default [${params.summaryFiblastdble}]
    ${MAGENTA}* --primerfile:${NC} Path of Primer set containing .xlsx file to performing primer trimming before genome assembly. Default [${params.primerfile}]
    ${MAGENTA}* --coverageCoinfection:${NC} Percentage of genome coverage to confidently detect co-infection samples. Default [${params.coverageCoinfection}]
    ${MAGENTA}* --meandepth:${NC} Depth of coverage to confidently detect co-infection samples. Default [${params.meandepth}]
    ${MAGENTA}* --min_reads_coinf:${NC} Minimum read threshold to confidently detect co-infection samples. Default [${params.min_reads_coinf}]

	${RED}#### Parameters to pass additional Arguments to the tools ####${NC}
	${MAGENTA}* --aligntrim_ext:${NC} Additional arguments to pass to align_trim. Default [${params.aligntrim_ext}]
    ${MAGENTA}* --fastp_ext:${NC} Additional arguments to pass to FASTP. Default [${params.fastp_ext}]
    ${MAGENTA}* --fastqc_ext:${NC} Additional arguments to pass to FASTQC. Default [${params.fastqc_ext}]
	${MAGENTA}* --hostile_ext:${NC} Additional arguments to pass to Hostile for dehosting.
			Pass parameters such as "--index /path/to/customReference.fa" for ONT or to a bowtie2 indexes for non-human host. Default [${params.hostile_ext}]
	${MAGENTA}* --spades_ext:${NC} Additional arguments to pass to rnaviralspades.py. Default [${params.spades_ext}]
	${MAGENTA}* --minimap_ext:${NC} Additional arguments to pass to minimap2. Default [${params.minimap_ext}]
	${MAGENTA}* --flye_ext:${NC} Additional arguments to pass to flye assembler. Default [${params.flye_ext}]
	${MAGENTA}* --chimeric_ext:${NC} Additional arguments to pass to minimap2 as per YACRD instruction here: https://github.com/natir/yacrd. Default [${params.chimeric_ext}]
	${MAGENTA}* --yacrd_ext=:${NC} Additional arguments to pass to YACRD as per YACRD instruction here: https://github.com/natir/yacrd. Default [${params.yacrd_ext}]
	${MAGENTA}* --canu_ext=:${NC} Additional arguments to pass to Canu assembler. Default [${params.canu_ext}]

	${RED}#### Parameters to Skip certain Steps ####${NC}
	${MAGENTA}* --skipTrim:${NC} Set this "true" to skip Trimming Step. Default [${params.skipTrim}]
	${MAGENTA}* --skipAlignment:${NC} Set this "true" to skip Alignment Step. Default [${params.skipAlignment}]
	${MAGENTA}* --skipDehost:${NC} Set this "true" to skip Dehosting Step. Default [${params.skipDehost}]
	${MAGENTA}* --skipDenovoAssembly:${NC} Set this "true" to skip Assembly Step. Default [${params.skipDenovoAssembly}]
	${MAGENTA}* --skipBlast:${NC} Set this "true" to skip BlastN based serotype determination. Default [${params.skipBlast}]
	${MAGENTA}* --skipScrubbing:${NC} Set this "true" to skip Chimeric read detection and removal from ONT data using YACRD. Default [${params.skipScrubbing}]
    ${MAGENTA}* --skipPrimertrim:${NC} Set this "true" to skip Primer trimming step. Default [${params.skipPrimertrim}]
    ${MAGENTA}* --skipCoinfection:${NC} Set this "true" to skip running co-infection module. Default [${params.skipCoinfection}]

	${RED}#### Parameters to increase speed or limit computational resources #### ${NC}
	${MAGENTA}* --jobs:${NC} No. of jobs/ samples to process parallely. Default [${params.jobs}]
	${MAGENTA}* --cpus:${NC} No of threads to be used for running tools in the pipeline. Default [${params.cpus}]

${CYAN}Output:${NC}
	${MAGENTA}* --outdir:${NC} Path/Name of the output directory. Default [${params.outdir}]

"""

exit 0
}

include {FASTQC; PYCOQC} from './modules/fastqc'
include {FASTP; POSTTRIMFASTQC} from './modules/fastp'
include {MULTIQC as PRETRIM; MULTIQC as POSTTRIM; MULTIQC as SUMMARISEALL} from './modules/fastqc'
include {HOSTILE} from './modules/hostile'
include {VIRSTRAIN_CALL; VIRSTRAIN_ALIGN_BESTMATCH} from './modules/virstrain'
include {RNAVIRALSPADES;UNICYCLER;FLYE;CANU} from './modules/assembly'
include {COMBINESEROTYPESUMMARY} from './modules/combinedSerotypeSummary'
include {BLASTN; BLASTN as BLASTNUNI; BLASTN as BLASTNSPADE;  BLASTN as BLASTNFLYE} from './modules/blastn'
include {SUMMARIZEBLASTN; SUMMARIZEBLASTN as SUMMARIZEBLASTUNI; SUMMARIZEBLASTN as SUMMARIZEBLASTSPADE; SUMMARIZEBLASTN as SUMMARIZEBLASTFLYE} from './modules/summarizeBlastN'
include {ALIGNTOREFERENCE} from "./modules/reference_based_assembly"
include {CLAIR3; MAKEREFBASEDASSEMBLY} from "./modules/clair3"
include {SUMMARIZEALIGNTOREFERENCE} from "./modules/summariseAlignToReference"
include {YACRD} from './modules/yacrd'
include {ALIGNTRIM} from './modules/primertrimONT'
include {COMBINEKMERREFSUMMARIES} from './modules/combinekmerrefsummary'
include {CONCATREFS; MAP_AND_SELECT_CONTIGS; REMAP_TO_SELECTED_CONTIG; RUN_CLAIR3_PER_CONTIG;BUILD_CONTIG_CONSENSUS} from './modules/coinfection'
include {SOFTWARE_VERSIONS_HTML} from './modules/software_versions'
params.help= false
params.input = false
params.outdir= false

// ----------------------------------------------------------------------------
// Parameter validation helpers
// ----------------------------------------------------------------------------
def validateParams() {

    if (params.help) {
        printHelp()
        exit 0
    }

    // ------------------------------------------------------------------------
    // Core required parameters
    // ------------------------------------------------------------------------
    if (!params.input) {
        error "Missing required parameter: --input"
    }

    if (!(params.mode in ['PE', 'SE'])) {
        error "Invalid value for --mode: '${params.mode}'. Allowed values: PE or SE"
    }

    if (!params.genomesize) {
        error "Missing required parameter: --genomesize"
    }

    if (!params.mafftMeta) {
        error "Missing required parameter: --mafftMeta"
    }

    if (!params.assembler) {
        error "Missing required parameter: --assembler"
    }

    // ------------------------------------------------------------------------
    // Validate assembler compatibility with sequencing mode
    // ------------------------------------------------------------------------
    def peAssemblers = ['unicycler', 'spades', 'all']
    def seAssemblers = ['canu', 'flye', 'all']

    if (params.mode == 'PE' && !(params.assembler in peAssemblers)) {
        error "Invalid assembler '${params.assembler}' for mode PE. Allowed: ${peAssemblers.join(', ')}"
    }

    if (params.mode == 'SE' && !(params.assembler in seAssemblers)) {
        error "Invalid assembler '${params.assembler}' for mode SE. Allowed: ${seAssemblers.join(', ')}"
    }

    // ------------------------------------------------------------------------
    // Conditional parameter requirements
    // ------------------------------------------------------------------------
    if (!params.skipBlast && !params.blastdb) {
        error "--blastdb must be provided when --skipBlast is false"
    }

    if (!params.db) {
        error "Missing required parameter: --db (VirStrain database path/prefix)"
    }

    if (!params.skipReferenceAssembly && !params.referenceDir) {
        error "--referenceDir must be provided when --skipReferenceAssembly is false"
    }

    if (!params.skipCoinfection && !params.referenceDir) {
        error "--referenceDir must be provided when --skipCoinfection is false"
    }

    if (params.mode == 'SE' && !params.skipPrimertrim && !params.primerfile) {
        log.warn "Primer trimming is enabled but --primerfile is empty. Primer trimming module may fail unless it handles this internally."
    }

    if (params.skipAlignment) {
        log.warn "Parameter --skipAlignment is currently defined for compatibility/documentation, but it is not used in main.nf orchestration."
    }
}



workflow{

    // ------------------------------------------------------------------------
    // Shared file/path channels used in multiple places
    // ------------------------------------------------------------------------
    mafft_meta_ch = Channel.fromPath(params.mafftMeta, checkIfExists: true)
    virstrain_db_ch = Channel.fromPath(params.db, checkIfExists: true)
    ch_versions = Channel.empty()
    // ------------------------------------------------------------------------
    // Step 0: Read input FASTQ files
    // ------------------------------------------------------------------------
    if (params.input != false){
      if (params.mode == "PE"){
        Channel.fromFilePairs(params.input, checkIfExists: true ).set { input_fastqs }
        } else if (params.mode == "SE") {
          Channel.fromPath(params.input, checkIfExists: true ).map { file -> tuple(file.simpleName, file)}.set { input_fastqs }
	}
    }

    // ------------------------------------------------------------------------
    // Step 1: Initial QC
    //
    // PE:
    //   - Run FASTQC on raw reads
    //   - Aggregate raw FASTQC with MultiQC
    //
    // SE:
    //   - If a sequencing summary file is provided, run pycoQC
    // ------------------------------------------------------------------------
    if (params.mode == "PE") {
        FASTQC(input_fastqs)
        pretrim_input = FASTQC.out.fastqc.collect()
        PRETRIM("01_rawFastQC", pretrim_input, 'pre-trimming')
        ch_versions = ch_versions.mix(PRETRIM.out.versions)
    }
    else if (params.mode == "SE" && params.summaryFile !="") {
        summaryfile_ch = Channel.fromPath(params.summaryFile, checkIfExists: true).map { file -> tuple(file.simpleName, file) }
        PYCOQC(summaryfile_ch)
        ch_versions = ch_versions.mix(PYCOQC.out.versions)
    }

    // ------------------------------------------------------------------------
    // Step 2: Dehosting
    //
    // Expected behavior:
    //   Some low-depth samples may fail downstream after dehosting. Those are
    //   tolerated at the process level by workflow configuration.
    // ------------------------------------------------------------------------

    if (params.skipDehost) {
        dehost_reads_ch = input_fastqs
    }
    else {
        HOSTILE(input_fastqs)
        dehost_reads_ch = HOSTILE.out[0]
                ch_versions = ch_versions.mix(HOSTILE.out.versions)
    }

    // ------------------------------------------------------------------------
    // Step 3: Prepare QC'd reads for downstream analysis
    //
    // PE:
    //   optional FASTP trimming
    //
    // SE:
    //   optional YACRD scrubbing for chimeric read removal
    // ------------------------------------------------------------------------
    qc_reads_ch = dehost_reads_ch

    if (params.mode == "PE" && !params.skipTrim) {
        FASTP(dehost_reads_ch)
        qc_reads_ch = FASTP.out[0]
         ch_versions = ch_versions.mix(FASTP.out.versions)
    }
    else if (params.mode == "SE" && !params.skipScrubbing) {
        yacrd_out = YACRD(dehost_reads_ch)
        qc_reads_ch = yacrd_out.scrubb
        ch_versions = ch_versions.mix(YACRD.out.versions)
    }

    // ------------------------------------------------------------------------
    // Step 4: VirStrain serotype calling and best-match alignment
    //
    // Output assumptions:
    //   aligned_out.viral_reads              -> tuple(sample_id, reads/file(s))
    //   aligned_out.serotype_contig_ordering -> tuple(sample_id, reference_fasta)
    //   aligned_out.serotyper_res            -> result files for per-sample summary
    // ------------------------------------------------------------------------
    virstrain_out = VIRSTRAIN_CALL(qc_reads_ch, virstrain_db_ch)
    aligned_out   = VIRSTRAIN_ALIGN_BESTMATCH(virstrain_out.report_for_align, mafft_meta_ch)

    assemble_reads_ch = aligned_out.viral_reads
    // ------------------------------------------------------------------------
    // Step 5: Combine serotype summaries across samples
    // ------------------------------------------------------------------------
    COMBINESEROTYPESUMMARY(
        aligned_out.serotyper_res.collect(),
        Channel.fromPath(params.mafftMeta, checkIfExists: true)
    )




    // ------------------------------------------------------------------------
    // Step 6: Build assembly input channel
    //
    // Channel contract after grouping:
    //   tuple(sample_id, [reads_and_reference_files...])
    //
    // This combines:
    //   - serotype-filtered reads
    //   - serotype-specific reference / contig ordering file
    // ------------------------------------------------------------------------
    assemble_reads_ch
        .mix(aligned_out.serotype_contig_ordering)
        .flatMap { sample_id, files ->
            files instanceof List ? files.collect { [sample_id, it] } : [[sample_id, files]]
        }
        .groupTuple(by: 0)
        .set { assembly_ch }

    // ------------------------------------------------------------------------
    // Step 7: De novo assembly
    // ------------------------------------------------------------------------
    if (!params.skipDenovoAssembly) {

        switch ("${params.mode}:${params.assembler}") {

            case "PE:unicycler":
                UNICYCLER(assembly_ch)
                ch_versions = ch_versions.mix(UNICYCLER.out.versions)
                break

            case "PE:spades":
                RNAVIRALSPADES(assembly_ch)
                 ch_versions = ch_versions.mix(RNAVIRALSPADES.out.versions)
                break

            case "PE:all":
                UNICYCLER(assembly_ch)
                ch_versions = ch_versions.mix(UNICYCLER.out.versions)
                RNAVIRALSPADES(assembly_ch)
                 ch_versions = ch_versions.mix(RNAVIRALSPADES.out.versions)
                break

            case "SE:canu":
                CANU(assembly_ch)
                ch_versions = ch_versions.mix(CANU.out.versions)
                break

            case "SE:flye":
                FLYE(assembly_ch)
                ch_versions = ch_versions.mix(FLYE.out.versions)
                break

            case "SE:all":
                CANU(assembly_ch)
                ch_versions = ch_versions.mix(CANU.out.versions)
                FLYE(assembly_ch)
                ch_versions = ch_versions.mix(FLYE.out.versions)
                break

            default:
                error "Unsupported mode/assembler combination: mode=${params.mode}, assembler=${params.assembler}"
        }
    }

    // ------------------------------------------------------------------------
    // Step 7.5: Coinfection detection
    //
    // Notes:
    //   - Expected that many non-coinfected samples may not progress through all
    //     coinfection-specific steps.
    //   - This is why process-level failure tolerance remains enabled in config.
    // ------------------------------------------------------------------------
    if (!params.skipCoinfection) {

        refs_ch     = Channel.fromPath("${params.referenceDir}/*.fasta", checkIfExists: true).collect()
        all_refs_ch = CONCATREFS(refs_ch)

        mapped_and_selected = MAP_AND_SELECT_CONTIGS(qc_reads_ch, all_refs_ch.combinedrefs)

        /*
         * selected_contigs_ch emits:
         *   tuple(sample_id, selected_contig_id, reads)
         *
         * Samples with no selected TSV or empty selections are intentionally
         * dropped from this branch.
         */
        selected_contigs_ch = mapped_and_selected.mapandselectout.flatMap { sid, reads, bam, cov, selected_tsv ->
            if (!selected_tsv) {
                return []
            }

            selected_tsv.text.readLines()
                .findAll { it?.trim() }
                .drop(1)
                .collect { line ->
                    def fields = line.split('\t')
                    tuple(sid, fields[1], reads)
                }
        }

        remapped      = REMAP_TO_SELECTED_CONTIG(selected_contigs_ch, all_refs_ch.combinedrefs)
        clair3_out    = RUN_CLAIR3_PER_CONTIG(remapped.remap2selectedcontig)
        consensus_out = BUILD_CONTIG_CONSENSUS(clair3_out.runclair3percontig)
    }

    // ------------------------------------------------------------------------
    // Step 8: Reference-based assembly (SE / ONT only)
    // ------------------------------------------------------------------------
    if (params.mode == "SE" && !params.skipReferenceAssembly) {

        ALIGNTOREFERENCE(
            qc_reads_ch,
            Channel.fromPath(params.referenceDir, type: 'any', checkIfExists: true)
        )

        SUMMARIZEALIGNTOREFERENCE(ALIGNTOREFERENCE.out[0].collect())

        COMBINEKMERREFSUMMARIES(
            COMBINESEROTYPESUMMARY.out.summaryreport1,
            SUMMARIZEALIGNTOREFERENCE.out.summaryreport2
        )

        // Perform optional primer trimming before Clair3-based consensus creation
        if (params.skipPrimertrim) {
            CLAIR3(ALIGNTOREFERENCE.out[1])
            MAKEREFBASEDASSEMBLY(CLAIR3.out[0])
        }
        else {
            trim_res = ALIGNTRIM(ALIGNTOREFERENCE.out[1])
            CLAIR3(trim_res.trimmed)
            MAKEREFBASEDASSEMBLY(CLAIR3.out[0])

        }
        ch_versions = ch_versions.mix(CLAIR3.out.versions)
        ch_versions = ch_versions.mix(MAKEREFBASEDASSEMBLY.out.versions)
    }


    // ------------------------------------------------------------------------
    // Step 9: BLASTN-based confirmation / serotyping of assembly outputs
    // ------------------------------------------------------------------------
    if (!params.skipBlast && !params.skipDenovoAssembly) {

        switch ("${params.mode}:${params.assembler}") {

            case "PE:unicycler":
                BLASTN(UNICYCLER.out[1], "06_assemblyBLASTunicycler")
                SUMMARIZEBLASTN(BLASTN.out[0].collect(), "06_assemblyBLASTunicycler")
                ch_versions = ch_versions.mix(UNICYCLER.out.versions)
                break

            case "PE:spades":
                BLASTN(RNAVIRALSPADES.out[1], "06_assemblyBLASTspades")
                SUMMARIZEBLASTN(BLASTN.out[0].collect(), "06_assemblyBLASTspades")
                ch_versions = ch_versions.mix(RNAVIRALSPADES.out.versions)
                break

            case "PE:all":
                BLASTNUNI(UNICYCLER.out[1], "06_assemblyBLASTunicycler")
                SUMMARIZEBLASTUNI(BLASTNUNI.out[0].collect(), "06_assemblyBLASTunicycler")

                BLASTNSPADE(RNAVIRALSPADES.out[1], "06_assemblyBLASTspades")
                SUMMARIZEBLASTSPADE(BLASTNSPADE.out[0].collect(), "06_assemblyBLASTspades")
                ch_versions = ch_versions.mix(UNICYCLER.out.versions)
                ch_versions = ch_versions.mix(RNAVIRALSPADES.out.versions)
                ch_versions = ch_versions.mix(BLASTNUNI.out.versions)
                break

            case "SE:canu":
                BLASTN(CANU.out[1], "06_assemblyBLASTcanu")
                SUMMARIZEBLASTN(BLASTN.out[0].collect(), "06_assemblyBLASTcanu")
                ch_versions = ch_versions.mix(CANU.out.versions)
                break

            case "SE:flye":
                BLASTN(FLYE.out[1], "06_assemblyBLASTflye")
                SUMMARIZEBLASTN(BLASTN.out[0].collect(), "06_assemblyBLASTflye")
                ch_versions = ch_versions.mix(FLYE.out.versions)
                break

            case "SE:all":
                BLASTN(CANU.out[1], "06_assemblyBLASTcanu")
                SUMMARIZEBLASTN(BLASTN.out[0].collect(), "06_assemblyBLASTcanu")

                BLASTNFLYE(FLYE.out[1], "06_assemblyBLASTflye")
                SUMMARIZEBLASTFLYE(BLASTNFLYE.out[0].collect(), "06_assemblyBLASTflye")
                ch_versions = ch_versions.mix(CANU.out.versions)
                ch_versions = ch_versions.mix(FLYE.out.versions)
                ch_versions = ch_versions.mix(BLASTNFLYE.out.versions)
                break

            default:
                error "Unsupported BLAST mode/assembler combination: mode=${params.mode}, assembler=${params.assembler}"
        }
    }


SOFTWARE_VERSIONS_HTML(ch_versions.collect())
}
