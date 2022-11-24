#!/usr/bin/env Rscript

library(getopt)

how_to <- function(){
  message(" ")
  message("This script converts a Manta VCF into bedpe. Run this script as follows:")
  message(" ")
  message("mantaVcfToBedpe [-i VCFFILE] [-o BEDPEFILE] [-f] [-p MINPR] [-m MODEGS] [-e GENOMEV]")
  message(" ")
  message("    -i VCFFILE       Manta VCF input file")
  message("    -o BEDPEFILE     Name of the output BEDPE file. If omitted a name will be given automatically")
  message("    -f               Flag to filter PASS mutations.")
  message("    -p MINPR         minimum PR of ALT. For high specificity use 18, lower for more sensitivity.")
  message("                     If not specificed, default value is 18")
  message("    -m MODEGS        Use S for somatic (Manta was run with both normal and tumour)")
  message("                     and G for germline (Manta was run with only normal).")
  message("    -e GENOMEV       Indicate the genome version to be used: hg19 or hg38")
  message("    -h               Show this explanation")
  message(" ")
  message("To run this script you should install the SomaticVariantAnnotation R package using:")
  message("install_github(\"PapenfussLab/StructuralVariantAnnotation@cb584be9474318e7f91c2cd8fca99f1212cab021\")")
  message(" ")
}

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'input',      'i', 1, "character",
  'help',       'h', 0, "logical",
  'filterpass', 'f', 0, "logical",
  'output',     'o', 1, "character",
  'mode',       'm', 1, "character",
  'genomev',    'e', 1, "character",
  'minPR',      'p', 1, "double"
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
  message("\nMissing VCFFILE. Quit.\n")
  how_to()
  q(status=1,save = "no")
}
if ( is.null(opt$mode    ) ) { 
  message("\nMissing MODEGS. Please set to S (somatic) or G (germline). Quit.\n")
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
if ( is.null(opt$minPR      ) ) { 
  opt$minPR = 18
  message("\nMissing MINPR. Set to ",opt$minPR,".\n")
}
if ( is.null(opt$output) ) { 
  opt$output = paste0(basename(opt$input),"_PR",opt$minPR,".bedpe")    
}

vcffile <- opt$input
bedpefile <- opt$output
PR_threshold <- opt$minPR
modeGS <- opt$mode
genomev <- opt$genomev
sample_name <- basename(vcffile)

mantaVcfToBedpe(vcffile,bedpefile,sample_name,genomev,PR_threshold,modeGS,filterPASS)

q(status = 0,save = "no")
