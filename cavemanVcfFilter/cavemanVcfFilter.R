#!/usr/bin/env Rscript

library(getopt)
library(utility.scripts)

how_to <- function(){
  message(" ")
  message("This script filters Caveman VCF files. Run this script as follows:")
  message(" ")
  message("cavemanVcfFilter [-i VCFIN] [-o VCFOUT] [-f] [-c CLPM_MAX] [-a ASMD_MIN] [-e GENOMEV]")
  message(" ")
  message("    -i VCFIN         Caveman VCF input file.")
  message("    -o VCFOUT        Name of the output VCF file. If omitted a name will be given automatically.")
  message("    -f               Flag to filter PASS mutations.")
  message("    -c CLPM_MAX      maximum CLPM. Default is 0.")
  message("    -a ASMD_MIN      minimum ASMD. Default is 140.")
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
  'asmd.min',      'a', 1, "double",
  'clpm.max',      'c', 1, "double"
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
if ( is.null(opt$asmd.min  ) ) { 
  opt$asmd.min = 140
  message("\nASMD_MIN not specified. Set to ",opt$asmd.min,".\n")
}
if ( is.null(opt$clpm.max  ) ) { 
  opt$clpm.max = 0
  message("\nCLPM_MAX not specified. Set to ",opt$clpm.max,".\n")
}
if ( is.null(opt$output) ) { 
  opt$output = paste0(basename(opt$input),"_ASMD",opt$asmd.min,"_CLPM",opt$clpm.max,".vcf")    
}

vcffile <- opt$input
vcfout <- opt$output
genomev <- opt$genomev
sample_name <- basename(vcffile)
clpm_max <- opt$clpm.max
asmd_min <- opt$asmd.min

cavemanVcfFilter(vcffile,vcfout,genomev,sample_name,clpm_max,asmd_min)

q(status = 0,save = "no")
