# install.packages("devtools")
# install.packages("roxygen2")
# 
# devtools::create("utility.scripts")

# setwd("~/git/utility.scripts")

#install_github("PapenfussLab/StructuralVariantAnnotation@cb584be9474318e7f91c2cd8fca99f1212cab021")

usethis::use_package("VariantAnnotation")
usethis::use_package("SummarizedExperiment")
usethis::use_package("BiocGenerics")
usethis::use_package("GenomeInfoDb")
usethis::use_package("BSgenome")
# devtools::use_package("BSgenome.Hsapiens.UCSC.hg38")
# devtools::use_package("BSgenome.Hsapiens.1000genomes.hs37d5")
# devtools::use_package("BSgenome.Mmusculus.UCSC.mm10")
# devtools::use_package("foreach")
# devtools::use_package("doParallel")
# devtools::use_package("doMC")
# devtools::use_package("methods")
usethis::use_package("stats")
# devtools::use_package("gmp")
# devtools::use_package("plyr")
# devtools::use_package("RCircos")
# devtools::use_package("scales")
usethis::use_package("GenomicRanges")
usethis::use_package("IRanges")
usethis::use_package("getopt")
usethis::use_package("stringr")
devtools::load_all()



devtools::document()
devtools::install()

#test all
devtools::test()





