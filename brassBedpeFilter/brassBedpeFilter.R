#!/usr/bin/env Rscript

library(getopt)
library(utility.scripts)

how_to <- function(){
  message(" ")
  message("This script filters Brass BEDPE files. Filters applied:")
  message("  1. Remove if no assembly score is given")
  message("  2. Remove if the distance between breakpoints is <1000bp")
  message("  3. Remove if it is an inversion shorter than 5000bp and readpair count is less than 6")
  message(" ")
  message("Run this script as follows:")
  message(" ")
  message("brassBedpeFilter [-i BEDPEIN] [-o BEDPEOUT]")
  message(" ")
  message("    -i BEDPEIN       Brass compressed Bedpe input file (.annot.bedpe.gz).")
  message("    -o BEDPEOUT      Name of the output bedpe file. If omitted a name will be given automatically.")
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
filterPASS = FALSE

if ( is.null(opt$input ) ) { 
  message("\nMissing BEDPEIN. Quit.\n")
  how_to()
  q(status=1,save = "no")
}
if ( is.null(opt$output) ) { 
  opt$output = paste0(basename(opt$input),"_filtered.bedpe")    
}

bedpein <- opt$input
bedpeout <- opt$output
sample_name <- basename(bedpein)

brassBedpeFilter(bedpein,bedpeout,sample_name)

q(status = 0,save = "no")
