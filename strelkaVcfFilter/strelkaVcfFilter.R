#!/usr/bin/env Rscript

library(getopt)
library(utility.scripts)

how_to <- function(){
  message(" ")
  message("This script filters Strelka VCF files. Run this script as follows:")
  message(" ")
  message("strelkaVcfFilter [-i VCFIN] [-o VCFOUT] [-f] [-p MINSOMEVS] [-e GENOMEV]")
  message(" ")
  message("    -i VCFIN         Strelka VCF input file.")
  message("    -o VCFOUT        Name of the output VCF file. If omitted a name will be given automatically.")
  message("                     Indels and SNV will be output to different files.")
  message("    -f               Flag to filter PASS mutations.")
  message("    -s MINSOMEVS     minimum SomaticEVS. For high specificity use 16, lower for more sensitivity.")
  message("                     If not specificed, no SomaticEVS filter is applied (e.g. for germline)")
  message("    -e GENOMEV       Indicate the genome version to be used: hg19 or hg38")
  message("    -h               Show this explanation")
  message(" ")
}

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'input',         'i', 1, "character",
  'help',          'h', 0, "logical",
  'filterpass',    'f', 0, "logical",
  'output',        'o', 1, "character",
  'genomev',       'e', 1, "character",
  'minSomaticEVS', 's', 1, "double"
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
filterPASS = FALSE
filterSomaticEVS = FALSE
Phred_threshold <- NULL

if ( is.null(opt$input ) ) { 
  message("\nMissing VCFIN. Quit.\n")
  how_to()
  q(status=1,save = "no")
}
if ( is.null(opt$genomev ) ) { 
  message("\nMissing GENOMEV. Quit.\n")
  how_to()
  q(status=1,save = "no")
}
if ( !is.null(opt$filterpass ) ) { 
  message("\nFilter PASS: TRUE.\n")
  filterPASS = TRUE
}
if ( !is.null(opt$minSomaticEVS  ) ) { 
  filterSomaticEVS = TRUE
  Phred_threshold <- opt$minSomaticEVS
}
if ( is.null(opt$output) ) { 
  opt$output = paste0(basename(opt$input),"_SomaticEVS",opt$minSomaticEVS)    
}

vcffile <- opt$input
vcfoutSNV <- paste0(opt$output,"_SNV.vcf")
vcfoutIndels <- paste0(opt$output,"_Indels.vcf")
genomev <- opt$genomev
sample_name <- basename(vcffile)

strelkaVcfFilter(vcffile,vcfoutSNV,vcfoutIndels,sample_name,genomev,Phred_threshold,filterPASS,filterSomaticEVS)

q(status = 0,save = "no")
