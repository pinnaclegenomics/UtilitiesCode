#!/usr/bin/env Rscript

library(getopt)
library(utility.scripts)

how_to <- function(){
  message(" ")
  message("This script converts an ascat comma separated values (CSV) file")
  message("into a tab separated values (TSV) file and adds headers.")
  message("Run this script as follows:")
  message(" ")
  message("ascatToTsv [-i ASCATIN] [-o TSVOUT]")
  message(" ")
  message("    -i ASCATIN       Ascat CSV file.")
  message("    -o TSVOUT        Name of the output TSV file. If omitted a name will be given automatically.")
  message("    -h               Show this explanation")
  message(" ")
}

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'input',         'i', 1, "character",
  'help',          'h', 0, "logical",
  'output',        'o', 1, "character"
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

if ( is.null(opt$input ) ) { 
  message("\nMissing ASCATIN. Quit.\n")
  how_to()
  q(status=1,save = "no")
}
if ( is.null(opt$output) ) { 
  opt$output = paste0(basename(opt$input),"_CNV.tsv")    
}

ascatin <- opt$input
tsvout <- opt$output
sample_name <- basename(ascatin)

ascatToTsv(ascatin,tsvout,sample_name)

q(status = 0,save = "no")
