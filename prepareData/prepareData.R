#!/usr/bin/env Rscript

library(getopt)
library(utility.scripts)

how_to <- function(){
  message(" ")
  message("This script filters and prepare raw files, so that they can be used")
  message("in the hrDetect pipeline.")
  message(" ")
  message("prepareData [OPTIONS]")
  message(" ")
  message("General options:")
  message("  --outdir (-o) DIR     Name of the output directory. If omitted a")
  message("                        name will be given automatically in the")
  message("                        current folder.")
  message("  --genomev (-g) GENOME Indicate the genome version to be used: hg19")
  message("                        or hg38")
  message("  --help (-h)           Show this explanation")
  message(" ")
  message("Input files options:")
  message("  --caveman (-C) FILE   File containing a tab separated table, with")
  message("                        sample names in the first column and caveman")
  message("                        files names in the second column.")
  message("  --pindel (-P) FILE    File containing a tab separated table, with")
  message("                        sample names in the first column and pindel")
  message("                        files names in the second column.")
  message("  --ascat (-A) FILE     File containing a tab separated table, with")
  message("                        sample names in the first column and ascat")
  message("                        files names in the second column.")
  message("  --canvas (-V) FILE    File containing a tab separated table, with")
  message("                        sample names in the first column and canvas")
  message("                        files names in the second column.")
  message("  --brass (-B) FILE     File containing a tab separated table, with")
  message("                        sample names in the first column and brass")
  message("                        files names in the second column.")
  message("  --strelka (-S) FILE   File containing a tab separated table, with")
  message("                        sample names in the first column and strelka")
  message("                        files names in the second column.")
  message("                        Files contain both SNVs and Indels.")
  message("  --strelkasnv (-N) FILE   File containing a tab separated table, with")
  message("                        sample names in the first column and strelka")
  message("                        files names in the second column.")
  message("                        Files contain only SNVs.")
  message("  --strelkaindels (-I) FILE   File containing a tab separated table, with")
  message("                        sample names in the first column and strelka")
  message("                        files names in the second column.")
  message("                        Files contain only Indels.")
  message("  --manta (-M) FILE     File containing a tab separated table, with")
  message("                        sample names in the first column and manta")
  message("                        files names in the second column.")
  message(" ")
  message("Input-specific options/info:")
  message(" ")
  message("CAVEMAN")
  message(" ")
  message("The files used as input should be from the Sanger Institute Caveman")
  message("script, use the *.annot.muts.vcf.gz files.")
  message("  --cavemanpass (-c)              Flag to filter PASS mutations.")
  message("  --cavemanclpm (-l) MAX_CLPM     maximum CLPM. Default is 0.")
  message("  --cavemanasmd (-a) MIN_ASMD     minimum ASMD. Default is 140.")
  message(" ")
  message("PINDEL")
  message(" ")
  message("The files used as input should be from the Sanger Institute Pindel")
  message("script, use the *.annot.vcf.gz files.")
  message("  --pindelpass (-p)              Flag to filter PASS mutations.")
  message("  --pindelqual (-q) QUAL         minimum QUAL. Default is 250.")
  message("  --pindelrep (-r) REP           maximum REP. Default is 9.")
  message(" ")
  message("ASCAT")
  message(" ")
  message("The files used as input should be from the Sanger Institute Ascat")
  message("script, use the *.copynumber.caveman.csv files.")
  message("This script converts an ascat comma separated values (CSV) file")
  message("into a tab separated values (TSV) file and adds headers.")
  message(" ")
  message("CANVAS")
  message(" ")
  message("The files used as input should be a vcf file containing Canvas")
  message("copy number calls. A PASS filter will be applied. Also only variants")
  message("of at least 10kb will be considered.")
  message(" ")
  message("BRASS")
  message(" ")
  message("The files used as input should be from the Sanger Institute Brass")
  message("script, use the *.annot.bedpe.gz files.")
  message("This script filters Brass BEDPE files. Filters applied:")
  message("  1. Remove if no assembly score is given")
  message("  2. Remove if the distance between breakpoints is <1000bp")
  message("  3. Remove if it is an inversion shorter than 5000bp and readpair")
  message("     count is less than 6")
  message(" ")
  message("STRELKA")
  message(" ")
  message("The files used as input should be from the Strelka algorithm. This")
  message("script will produce to separate files for SNV and Indels.")
  message("  --strelkapass (-s)         Flag to filter PASS mutations.")
  message("  --strelkasomevs (-e) EVS   minimum SomaticEVS. For high")
  message("                             specificity use 16, lower for more")
  message("                             sensitivity (default is disabled).")
  message(" ")
  message("MANTA")
  message(" ")
  message("The files used as input should be from the Manta algorithm. This")
  message("script will produce a BEDPE.")
  message("To run this script you should install the SomaticVariantAnnotation")
  message("R package using:")
  message("install_github(\"PapenfussLab/StructuralVariantAnnotation@cb584be9474318e7f91c2cd8fca99f1212cab021\")")
  message("  --mantapass (-m)         Flag to filter PASS mutations.")
  message("  --mantapr (-k) PR        minimum PR of ALT. Default value is PR=8,")
  message("                           increase it for more specificity, decrease it")
  message("                           for more sensitivity.")
  message(" ")
}

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'outdir',        'o', 1, "character",
  'genomev',        'g', 1, "character",
  'caveman',       'C', 1, "character",
  'pindel',        'P', 1, "character",
  'ascat',         'A', 1, "character",
  'canvas',        'V', 1, "character",
  'brass',         'B', 1, "character",
  'strelka',       'S', 1, "character",
  'strelkasnv',       'N', 1, "character",
  'strelkaindels',    'I', 1, "character",
  'manta',         'M', 1, "character",
  'cavemanpass',   'c', 0, "logical",
  'cavemanclpm',   'l', 1, "double",
  'cavemanasmd',   'a', 1, "double",
  'pindelpass',    'p', 0, "logical",
  'pindelqual',    'q', 1, "double",
  'pindelrep',     'r', 1, "double",
  'strelkapass',   's', 0, "logical",
  'strelkasomevs', 'e', 1, "double",
  'mantapass',     'm', 0, "logical",
  'mantapr',       'k', 1, "double",
  'help',          'h', 0, "logical"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  how_to()
  q(status=1,save = "no")
}

# set some reasonable defaults for the options that are needed,
# but were not specified.
cavemanFilterPASS = FALSE
pindelFilterPASS = FALSE
strelkaFilterPASS = FALSE
mantaFilterPASS = FALSE

if ( is.null(opt$outdir) ) {
  opt$outdir <- "prepareDataOutput"
}

opt$outdir <- paste0(opt$outdir,"/")

if ( is.null(opt$genomev ) ) {
  message("\nMissing GENOMEV. Quit.\n")
  how_to()
  q(status=1,save = "no")
}

if(!is.null(opt$ascat) & !is.null(opt$canvas)){
  message("\nPlease provide Ascat or Canvas, not both. Quit.\n")
  how_to()
  q(status=1,save = "no")
}

#Caveman inputs
if(!is.null(opt$caveman)){
  if ( !is.null(opt$cavemanpass ) ) {
    message("\nCaveman filter PASS: TRUE.\n")
    cavemanFilterPASS = TRUE
  }
  if ( is.null(opt$cavemanasmd  ) ) {
    opt$cavemanasmd = 140
    message("\nCaveman ASMD_MIN not specified. Set to ",opt$cavemanasmd,".\n")
  }
  if ( is.null(opt$cavemanclpm  ) ) {
    opt$cavemanclpm = 0
    message("\nCaveman CLPM_MAX not specified. Set to ",opt$cavemanclpm,".\n")
  }
}

#Pindel inputs
if(!is.null(opt$pindel)){
  if ( !is.null(opt$pindelpass ) ) {
    message("\nPindel filter PASS: TRUE.\n")
    pindelFilterPASS = TRUE
  }
  if ( is.null(opt$pindelqual  ) ) {
    opt$pindelqual = 250
    message("\nPindel QUAL not specified. Set to ",opt$pindelqual,".\n")
  }
  if ( is.null(opt$pindelrep  ) ) {
    opt$pindelrep = 9
    message("\nPindel REP not specified. Set to ",opt$pindelrep,".\n")
  }
}

#Strelka inputs
if(is.null(opt$strelkasomevs)){
  filterSomaticEVS <- FALSE
}else{
  filterSomaticEVS <- TRUE
}
if(!is.null(opt$strelka) | !is.null(opt$strelkasnv) | !is.null(opt$strelkaindels)){
  if ( !is.null(opt$strelkapass ) ) {
    message("\nStrelka filter PASS: TRUE.\n")
    strelkaFilterPASS = TRUE
  }
}
#Manta inputs

if(!is.null(opt$manta)){
  if ( !is.null(opt$mantapass ) ) {
    message("\nManta filter PASS: TRUE.\n")
    mantaFilterPASS = TRUE
  }
  if ( is.null(opt$mantapr      ) ) {
    opt$mantapr = 8
    message("\nManta Missing MINPR. Set to ",opt$mantapr,".\n")
  }
}

dir.create(opt$outdir,showWarnings = FALSE,recursive = TRUE)

#analysisTable for hrDetect pipeline
analysisTable <- data.frame()

#caveman
if(!is.null(opt$caveman)){
  if(file.exists(opt$caveman)){
    cavemantable <- read.table(opt$caveman,sep = "\t",stringsAsFactors = FALSE)
    for (i in 1:nrow(cavemantable)) {
      vcffile <- cavemantable[i,2]
      sample_name <- cavemantable[i,1]
      sample_dir <- paste0(opt$outdir,sample_name,"/")
      dir.create(sample_dir,showWarnings = FALSE,recursive = TRUE)
      vcfout <- paste0(sample_dir,sample_name,"_SNV_ASMD",opt$cavemanasmd,"_CLPM",opt$cavemanclpm,".vcf")
      cavemanVcfFilter(vcffile,vcfout,opt$genomev,sample_name,opt$cavemanclpm,opt$cavemanasmd,cavemanFilterPASS)
      analysisTable[sample_name,"SNV_vcf_files"] <- paste0(vcfout,".bgz")
    }
  }else{
    message("Caveman file does not exist: ",opt$caveman)
  }
}

#pindel
if(!is.null(opt$pindel)){
  if(file.exists(opt$pindel)){
    pindeltable <- read.table(opt$pindel,sep = "\t",stringsAsFactors = FALSE)
    for (i in 1:nrow(pindeltable)) {
      vcffile <- pindeltable[i,2]
      sample_name <- pindeltable[i,1]
      sample_dir <- paste0(opt$outdir,sample_name,"/")
      dir.create(sample_dir,showWarnings = FALSE,recursive = TRUE)
      vcfout <- paste0(sample_dir,sample_name,"_Indels_QUAL",opt$pindelqual,"_REP",opt$pindelrep,".vcf")
      pindelVcfFilter(vcffile,vcfout,opt$genomev,sample_name,opt$pindelqual,opt$pindelrep,pindelFilterPASS)
      analysisTable[sample_name,"Indels_vcf_files"] <- paste0(vcfout,".bgz")
    }
  }else{
    message("Pindel file does not exist: ",opt$pindel)
  }
}

#ascat
if(!is.null(opt$ascat)){
  if(file.exists(opt$ascat)){
    ascattable <- read.table(opt$ascat,sep = "\t",stringsAsFactors = FALSE)
    for (i in 1:nrow(ascattable)) {
      ascatin <- ascattable[i,2]
      sample_name <- ascattable[i,1]
      sample_dir <- paste0(opt$outdir,sample_name,"/")
      dir.create(sample_dir,showWarnings = FALSE,recursive = TRUE)
      tsvout <- paste0(sample_dir,sample_name,"_CNV.tsv")
      ascatToTsv(ascatin,tsvout,sample_name)
      analysisTable[sample_name,"CNV_tab_files"] <- tsvout
    }
  }else{
    message("Ascat file does not exist: ",opt$ascat)
  }
}

#canvas
if(!is.null(opt$canvas)){
  if(file.exists(opt$canvas)){
    canvastable <- read.table(opt$canvas,sep = "\t",stringsAsFactors = FALSE)
    for (i in 1:nrow(canvastable)) {
      canvasin <- canvastable[i,2]
      sample_name <- canvastable[i,1]
      sample_dir <- paste0(opt$outdir,sample_name,"/")
      dir.create(sample_dir,showWarnings = FALSE,recursive = TRUE)
      tsvout <- paste0(sample_dir,sample_name,"_CNV.tsv")
      canvasVcfToTsv(canvasVcf = canvasin,
                     tsvout = tsvout,
                     sample_name = sample_name,genomev = opt$genomev)
      analysisTable[sample_name,"CNV_tab_files"] <- tsvout
    }
  }else{
    message("Canvas file does not exist: ",opt$canvas)
  }
}

#brass
if(!is.null(opt$brass)){
  if(file.exists(opt$brass)){
    brasstable <- read.table(opt$brass,sep = "\t",stringsAsFactors = FALSE)
    for (i in 1:nrow(brasstable)) {
      bedpein <- brasstable[i,2]
      sample_name <- brasstable[i,1]
      sample_dir <- paste0(opt$outdir,sample_name,"/")
      dir.create(sample_dir,showWarnings = FALSE,recursive = TRUE)
      bedpeout <- paste0(sample_dir,sample_name,"_Rearr.bedpe")
      brassBedpeFilter(bedpein,bedpeout,sample_name)
      analysisTable[sample_name,"SV_bedpe_files"] <- bedpeout
    }
  }else{
    message("Brass file does not exist: ",opt$bedpe)
  }
}

#strelka
if(!is.null(opt$strelka)){
  if(file.exists(opt$strelka)){
    strelkatable <- read.table(opt$strelka,sep = "\t",stringsAsFactors = FALSE)
    for (i in 1:nrow(strelkatable)) {
      vcffile <- strelkatable[i,2]
      sample_name <- strelkatable[i,1]
      sample_dir <- paste0(opt$outdir,sample_name,"/")
      dir.create(sample_dir,showWarnings = FALSE,recursive = TRUE)
      if(filterSomaticEVS){
        vcfoutSNV <- paste0(sample_dir,sample_name,"_SNV_SomaticEVS",opt$strelkasomevs,".vcf")
        vcfoutIndels <- paste0(sample_dir,sample_name,"_Indels_SomaticEVS",opt$strelkasomevs,".vcf")
      }else{
        vcfoutSNV <- paste0(sample_dir,sample_name,"_SNV.vcf")
        vcfoutIndels <- paste0(sample_dir,sample_name,"_Indels.vcf")
      }
      strelkaVcfFilter(vcffile,vcfoutSNV,vcfoutIndels,sample_name,opt$genomev,opt$strelkasomevs,strelkaFilterPASS,filterSomaticEVS = filterSomaticEVS)
      analysisTable[sample_name,"SNV_vcf_files"] <- paste0(vcfoutSNV,".bgz")
      analysisTable[sample_name,"Indels_vcf_files"] <- paste0(vcfoutIndels,".bgz")
    }
  }else{
    message("Strelka file does not exist: ",opt$strelka)
  }
}

#strelka SNV
if(!is.null(opt$strelkasnv)){
  if(file.exists(opt$strelkasnv)){
    strelkasnvtable <- read.table(opt$strelkasnv,sep = "\t",stringsAsFactors = FALSE)
    for (i in 1:nrow(strelkasnvtable)) {
      vcffile <- strelkasnvtable[i,2]
      sample_name <- strelkasnvtable[i,1]
      sample_dir <- paste0(opt$outdir,sample_name,"/")
      dir.create(sample_dir,showWarnings = FALSE,recursive = TRUE)
      if(filterSomaticEVS){
        vcfoutSNV <- paste0(sample_dir,sample_name,"_SNV_SomaticEVS",opt$strelkasomevs,".vcf")
      }else{
        vcfoutSNV <- paste0(sample_dir,sample_name,"_SNV.vcf")
      }
      vcfoutIndels <- NULL
      strelkaVcfFilter(vcffile,vcfoutSNV,vcfoutIndels,sample_name,opt$genomev,opt$strelkasomevs,strelkaFilterPASS,filterSomaticEVS = filterSomaticEVS)
      analysisTable[sample_name,"SNV_vcf_files"] <- paste0(vcfoutSNV,".bgz")
    }
  }else{
    message("Strelka SNV file does not exist: ",opt$strelkasnv)
  }
}

#strelka Indels
if(!is.null(opt$strelkaindels)){
  if(file.exists(opt$strelkaindels)){
    strelkaindelstable <- read.table(opt$strelkaindels,sep = "\t",stringsAsFactors = FALSE)
    for (i in 1:nrow(strelkaindelstable)) {
      vcffile <- strelkaindelstable[i,2]
      sample_name <- strelkaindelstable[i,1]
      sample_dir <- paste0(opt$outdir,sample_name,"/")
      dir.create(sample_dir,showWarnings = FALSE,recursive = TRUE)
      vcfoutSNV <- NULL
      if(filterSomaticEVS){
        vcfoutIndels <- paste0(sample_dir,sample_name,"_Indels_SomaticEVS",opt$strelkasomevs,".vcf")
      }else{
        vcfoutIndels <- paste0(sample_dir,sample_name,"_Indels.vcf")
      }
      strelkaVcfFilter(vcffile,vcfoutSNV,vcfoutIndels,sample_name,opt$genomev,opt$strelkasomevs,strelkaFilterPASS,filterSomaticEVS = filterSomaticEVS)
      analysisTable[sample_name,"Indels_vcf_files"] <- paste0(vcfoutIndels,".bgz")
    }
  }else{
    message("Strelka Indels file does not exist: ",opt$strelkaindels)
  }
}

#manta
if(!is.null(opt$manta)){
  if(file.exists(opt$manta)){
    mantatable <- read.table(opt$manta,sep = "\t",stringsAsFactors = FALSE)
    for (i in 1:nrow(mantatable)) {
      vcffile <- mantatable[i,2]
      sample_name <- mantatable[i,1]
      sample_dir <- paste0(opt$outdir,sample_name,"/")
      dir.create(sample_dir,showWarnings = FALSE,recursive = TRUE)
      bedpeout <- paste0(sample_dir,sample_name,"_Rearr_PR",opt$mantapr,".bedpe")
      mantaVcfToBedpe(vcffile,bedpeout,sample_name,opt$genomev,PR_threshold = opt$mantapr,modeGS = "S",filterPASS = mantaFilterPASS)
      analysisTable[sample_name,"SV_bedpe_files"] <- bedpeout
    }
  }else{
    message("Manta file does not exist: ",opt$bedpe)
  }
}

#write analysisTable
analysisTable[,"sample"] <- rownames(analysisTable)
analysisTable <- cbind(data.frame(sample=rownames(analysisTable),stringsAsFactors = F),analysisTable)
write.table(analysisTable,file = paste0(opt$outdir,"analysisTable_hrDetect.tsv"),
            sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)


q(status = 0,save = "no")
