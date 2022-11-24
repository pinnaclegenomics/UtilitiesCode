#!/usr/bin/env Rscript

library(getopt)
library(utility.scripts)

how_to <- function(){
  message(" ")
  message("This script filters Pindel VCF files. Run this script as follows:")
  message(" ")
  message("pindelVcfFilter [-i VCFIN] [-o VCFOUT] [-f] [-q QUAL] [-r REP] [-e GENOMEV]")
  message(" ")
  message("    -i VCFIN         Pindel VCF input file.")
  message("    -o VCFOUT        Name of the output VCF file. If omitted a name will be given automatically.")
  message("    -f               Flag to filter PASS mutations.")
  message("    -q QUAL          minimum QUAL. Default is 250.")
  message("    -r REP           maximum REP. Default is 9.")
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
  'qual.min',      'q', 1, "double",
  'rep.max',      'r', 1, "double"
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
if ( is.null(opt$qual.min  ) ) { 
  opt$qual.min = 250
  message("\nQUAL not specified. Set to ",opt$qual.min,".\n")
}
if ( is.null(opt$rep.max  ) ) { 
  opt$rep.max = 9
  message("\nREP not specified. Set to ",opt$rep.max,".\n")
}
if ( is.null(opt$output) ) { 
  opt$output = paste0(basename(opt$input),"_QUAL",opt$qual.min,"_REP",opt$rep.max,".vcf")    
}

vcffile <- opt$input
vcfout <- opt$output
genomev <- opt$genomev
sample_name <- basename(vcffile)
qualmin <- opt$qual.min
repmax <- opt$rep.max

pindelVcfFilter(vcffile,vcfout,genomev,sample_name,qualmin,repmax,filterPASS)

q(status = 0,save = "no")
