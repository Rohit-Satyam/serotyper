#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


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
	${MAGENTA}* --assembler:${NC} Assembler to use for viral genome assembly.
			Possible options for Illumina PR: "unicycler", "spades" or "all".
			Possible options for ONT: "canu", "flye" or "all". Default [${params.assembler}]

	${RED}#### Optional Arguments${NC}
	${MAGENTA}* --summaryFile:${NC} ONT sequencing summary file to run pycoQC. Default [${params.summaryFile}]
	${MAGENTA}* --blastdb:${NC} Path of BlastN indexes with the prefix. Provide this if skipBlast is set to FALSE.

	${RED}#### Parameters to pass additional Arguments to the tools ####${NC}
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
	${MAGENTA}* --skipAssembly:${NC} Set this "true" to skip Assembly Step. Default [${params.skipAssembly}]
	${MAGENTA}* --skipBlast:${NC} Set this "true" to skip BlastN based serotype determination. Default [${params.skipBlast}]
	${MAGENTA}* --skipChimeraDetect:${NC} Set this "true" to skip Chimeric read detection and removal from ONT data using YACRD. Default [${params.skipChimeraDetect}]

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
include {VIRSTRAIN} from './modules/virstrain'
include {RNAVIRALSPADES;UNICYCLER;FLYE;CANU} from './modules/assembly'
include {COMBINESEROTYPESUMMARY} from './modules/combinedSerotypeSummary'
include {BLASTN; BLASTN as BLASTNUNI; BLASTN as BLASTNSPADE;  BLASTN as BLASTNFLYE} from './modules/blastn'
include {SUMMARIZEBLASTN; SUMMARIZEBLASTN as SUMMARIZEBLASTUNI; SUMMARIZEBLASTN as SUMMARIZEBLASTSPADE; SUMMARIZEBLASTN as SUMMARIZEBLASTFLYE} from './modules/summarizeBlastN'
include {YACRD} from './modules/yacrd'
params.help= false
params.input = false
params.outdir= false


workflow{

// Step 0: Reading the FASTQ files both single end and paired end
if (params.input != false){
	if (params.mode == "PE"){
		Channel.fromFilePairs(params.input, checkIfExists: true )
		.set { input_fastqs }
		} else if (params.mode == "SE") {
			Channel.fromPath(params.input, checkIfExists: true ).map { file -> tuple(file.simpleName, file)}
			.set { input_fastqs }
	}
}




//************ Illumina Shortread workflow *****************************

// Step 1 PE: Running FASTQC if Illumina
if (params.mode == "PE"){
  rawfqc_ch=FASTQC(input_fastqs)
  pretrim_input=FASTQC.out.fastqc.collect()
  PRETRIM("01_rawFastQC",pretrim_input,'pre-trimming')
} else if (params.mode == "SE" && params.summaryFile != ""){
  summaryfile_ch= Channel.fromPath(params.summaryFile, checkIfExists: true ).map { file -> tuple(file.simpleName, file)}
PYCOQC(summaryfile_ch)
}

// Step 2: Common dehosting Step
if (params.skipDehost == false){
  HOSTILE(input_fastqs)
  dehost_reads_ch = HOSTILE.out[0]
} else {
  dehost_reads_ch=input_fastqs
}


// Step 2 PE: Running Adapter Trimming and dehosting
if (params.skipTrim == false && params.mode == "PE"){
        FASTP(dehost_reads_ch)
        //assemble_reads_ch=FASTP.out[0]
        VIRSTRAIN(FASTP.out[0],Channel.fromPath(params.db, checkIfExists: true ))
        assemble_reads_ch=VIRSTRAIN.out.viral_reads
} else if (params.skipTrim == true && params.mode == "PE"){
        VIRSTRAIN(dehost_reads_ch,Channel.fromPath(params.db, checkIfExists: true ))
        assemble_reads_ch=VIRSTRAIN.out.viral_reads
}

// Step 2 SE: Decision to run yacrd to remove

if (params.mode == "SE" && params.skipChimeraDetect == false ){
  YACRD(dehost_reads_ch)
  VIRSTRAIN(YACRD.out[0],Channel.fromPath(params.db, checkIfExists: true ),Channel.fromPath(params.mafftMeta, checkIfExists: true ))
  assemble_reads_ch=VIRSTRAIN.out.viral_reads
} else if (params.mode == "SE" && params.skipChimeraDetect == true){
  VIRSTRAIN(dehost_reads_ch,Channel.fromPath(params.db, checkIfExists: true ),Channel.fromPath(params.mafftMeta, checkIfExists: true ))
  assemble_reads_ch=VIRSTRAIN.out.viral_reads
}

// Step 3: Combining the results from all the samples and
  COMBINESEROTYPESUMMARY(VIRSTRAIN.out.serotyper_res.collect(),Channel.fromPath(params.mafftMeta, checkIfExists: true ))

// Step 4: Now mixing FASTQ channels and fasta path identified as best strain of by virstrain

//mixing the channel output
  assemble_reads_ch.mix(VIRSTRAIN.out.serotype_contig_ordering)
    .flatMap { srr, files -> files instanceof List ? files.collect { [srr, it] } : [[srr, files]] } // Flatten the channel so that each file has its own tuple
    .groupTuple(by: 0) // Group the tuples by the SRR number
    .set{assembly_ch}   // Merge the files into a single list for each SRR number

//assembly_ch.view()
// Step 5: Assembly and scaffoling using ragtag because this way blastN searches will be less
if (params.skipAssembly == false && params.assembler=="unicycler"){
  UNICYCLER(assembly_ch)
} else if (params.skipAssembly == false && params.assembler=="spades"){
  RNAVIRALSPADES(assembly_ch)
} else if (params.skipAssembly == false && params.assembler=="flye") {
  FLYE(assembly_ch)
} else if (params.skipAssembly == false && params.assembler=="canu") {
  CANU(assembly_ch)
} else if (params.skipAssembly == false && params.assembler == "all" && params.mode == "PE"){
  UNICYCLER(assembly_ch)
  RNAVIRALSPADES(assembly_ch)
} else if (params.skipAssembly == false && params.assembler == "all" && params.mode == "SE"){
  CANU(assembly_ch)
  FLYE(assembly_ch)
}



// Step 6: BLASTN based serotyping and summarising
if (params.skipBlast == false && params.skipAssembly == false && params.assembler=="unicycler"){
  BLASTN(UNICYCLER.out[1],"06_assemblyBLASTunicycler")
  SUMMARIZEBLASTN(BLASTN.out[0].collect(),"06_assemblyBLASTunicycler")
} else if (params.skipBlast == false && params.skipAssembly == false && params.assembler=="spades"){
  BLASTN(RNAVIRALSPADES.out[1],"06_assemblyBLASTspades")
  SUMMARIZEBLASTN(BLASTN.out[0].collect(),"06_assemblyBLASTspades")
} else if (params.skipBlast == false && params.skipAssembly == false && params.assembler=="flye"){
  BLASTN(FLYE.out[1],"06_assemblyBLASTflye")
  SUMMARIZEBLASTN(BLASTN.out[0].collect(),"06_assemblyBLASTflye")
} else if (params.skipBlast == false && params.skipAssembly == false && params.assembler=="canu"){
  BLASTN(CANU.out[1],"06_assemblyBLASTcanu")
  SUMMARIZEBLASTN(BLASTN.out[0].collect(),"06_assemblyBLASTcanu")
} else if (params.skipBlast == false && params.skipAssembly == false && params.assembler=="all" && params.mode == "PE"){
  BLASTNUNI(UNICYCLER.out[1],"06_assemblyBLASTunicycler")
  SUMMARIZEBLASTUNI(BLASTNUNI.out[0].collect(),"06_assemblyBLASTunicycler")
  BLASTNSPADE(RNAVIRALSPADES.out[1],"06_assemblyBLASTspades")
  SUMMARIZEBLASTSPADE(BLASTNSPADE.out[0].collect(),"06_assemblyBLASTspades")
} else if (params.skipBlast == false && params.skipAssembly == false && params.assembler=="all" && params.mode == "SE"){
  BLASTN(CANU.out[1],"06_assemblyBLASTcanu")
  SUMMARIZEBLASTN(BLASTN.out[0].collect(),"06_assemblyBLASTcanu")
  BLASTNFLYE(FLYE.out[1],"06_assemblyBLASTflye")
  SUMMARIZEBLASTFLYE(BLASTNFLYE.out[0].collect(),"06_assemblyBLASTflye")
}





}
