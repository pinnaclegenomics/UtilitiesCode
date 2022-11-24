
#' @export
cavemanVcfFilter <- function(vcffile,vcfout,genomev,sample_name,clpm_max,asmd_min,filterPASS){
  #need to remove caveman badly formatted header lines,
  #which break the writeVcf function in recent versions of VariantAnnotation
  tempvcf <- tempfile()
  system(paste0("zgrep -v vcfProcessLog ",vcffile," | gzip -c > ",tempvcf))

  message("Load Caveman file: ",sample_name)
  subs_VCF <- VariantAnnotation::readVcf(tempvcf,genome = genomev)
  e.snv <- VariantAnnotation::expand(subs_VCF)

  #PASS filter
  if(filterPASS){
    message("Caveman Filter PASS: ",sample_name)
    selected_snv <- VariantAnnotation::fixed(e.snv)[,"FILTER"]=="PASS"
    e.snv <- e.snv[selected_snv,]
  }

  message("Filter caveman SNV (CLPM/ASMD): ",sample_name)
  selected_snv <- VariantAnnotation::info(e.snv)$CLPM<=clpm_max & VariantAnnotation::info(e.snv)$ASMD>=asmd_min
  e.snv <- e.snv[selected_snv,]

  message("Caveman Number of mutations after filtering: ",nrow(e.snv))

  message("Saving Caveman to VCF file: ",vcfout)
  VariantAnnotation::writeVcf(e.snv,vcfout,index=TRUE)
}

#' @export
pindelVcfFilter <- function(vcffile,vcfout,genomev,sample_name,qualmin,repmax,filterPASS){

  message("Load Pindel file: ",sample_name)
  indels_VCF <- VariantAnnotation::readVcf(vcffile,genome = genomev)
  e.indels <- VariantAnnotation::expand(indels_VCF)

  #PASS filter
  if(filterPASS){
    message("Pindel Filter PASS: ",sample_name)
    selected_indels <- VariantAnnotation::fixed(e.indels)[,"FILTER"]=="PASS"
    e.indels <- e.indels[selected_indels,]
  }

  message("Filter Pindel Indels (QUAL/REP): ",sample_name)
  selected_indels <- VariantAnnotation::fixed(e.indels)[,"QUAL"]>=qualmin & VariantAnnotation::info(e.indels)$REP<=repmax
  e.indels <- e.indels[selected_indels,]

  message("Pindel Number of mutations after filtering: ",nrow(e.indels))

  message("Saving Indels to VCF file: ",vcfout)
  VariantAnnotation::writeVcf(e.indels,vcfout,index=TRUE)
}

#' @export
ascatToTsv <- function(ascatin,tsvout,sample_name){
  message("Load Ascat file: ",sample_name)
  ascat_table <- read.table(ascatin,sep = ",",stringsAsFactors = FALSE,header = FALSE)
  colnames(ascat_table) <- c("seg_no","Chromosome","chromStart","chromEnd",
                             "total.copy.number.inNormal","minor.copy.number.inNormal",
                             "total.copy.number.inTumour","minor.copy.number.inTumour")
  message("Saving CNV to tab file: ",tsvout)
  write.table(ascat_table,file = tsvout,sep = "\t",
              col.names = TRUE,row.names = FALSE,quote = FALSE)
}

#' @export
brassBedpeFilter <- function(bedpein,bedpeout,sample_name){

  library("R.utils")

  #brass
  message("Untar Brass BEDPE file: ",sample_name)
  bedpe_file <- bedpeout
  gunzip(filename = bedpein,destname = bedpe_file,remove = FALSE,overwrite = TRUE)
  #retrieve the header
  con=file(bedpe_file,open="r")
  bedpe_lines=readLines(con)
  close(con)
  header_nline <- grep("^# chr1",bedpe_lines)
  header_line <- bedpe_lines[header_nline]
  header_line <- substr(header_line,3,nchar(header_line))
  bedpe_header <- strsplit(header_line,split = "\t")[[1]]
  if(header_nline==length(bedpe_lines)){
    # there are no rearrangements here, just the header
    # prep empty table
    bedpe_table <- data.frame(matrix(NA,nrow = 0,ncol = length(bedpe_header),dimnames = list(c(),bedpe_header)),
                              check.names = F,stringsAsFactors = F)
  }else{
    #load table
    bedpe_table <- read.table(bedpe_file,sep = "\t",comment.char = "",header = FALSE,
                              check.names = FALSE,stringsAsFactors = FALSE,skip = grep("^# chr1",bedpe_lines))
    colnames(bedpe_table) <- bedpe_header
  }

  colnames(bedpe_table)[bedpe_header=="chr1"] <- "chrom1"
  colnames(bedpe_table)[bedpe_header=="chr2"] <- "chrom2"

  #filters
  message("Filter Rearr: ",sample_name)
  if(nrow(bedpe_table)>0){
    selected_rearr <- bedpe_table$assembly_score != "_" & (bedpe_table$bkdist == -1 | bedpe_table$bkdist >= 1000) & (bedpe_table$svclass != "inversion" | (bedpe_table$bkdist >= 5000 | bedpe_table$`readpair count` > 5))
    bedpe_table <- bedpe_table[selected_rearr,]
  }
  if(nrow(bedpe_table)>0){
    #remove germline
    #selected_rearr <- bedpe_table$sample == paste0(sample_name,"T")
    #selected_rearr <- ! sapply(strsplit(bedpe_table$sample,split = ","),function(x) normal_name %in% x)
    #if both normal and tumour files have the rearrangement (both listed and separated by a comma, then it is a germline mutation)
    selected_rearr <- ! sapply(strsplit(bedpe_table$sample,split = ","),function(x) length(x)>1)
    bedpe_table <- bedpe_table[selected_rearr,]
  }
  message("Saving brassfile: ",bedpeout)
  write.table(bedpe_table,file = bedpe_file,sep = "\t",
              quote = FALSE,col.names = TRUE,row.names = FALSE)
}

#' @export
strelkaVcfFilter <- function(vcffile,vcfoutSNV,vcfoutIndels,sample_name,genomev,Phred_threshold,filterPASS,filterSomaticEVS){
  message("Load Strelka file: ",sample_name)
  subs_VCF <- VariantAnnotation::readVcf(vcffile,genome = genomev)
  e.snv <- VariantAnnotation::expand(subs_VCF)

  #somaticEVS and PASS filter
  if(filterPASS){
    message("Strelka Filter PASS: ",sample_name)
    selected_snv <- VariantAnnotation::fixed(e.snv)[,"FILTER"]=="PASS"
    e.snv <- e.snv[selected_snv,]
  }

  if(filterSomaticEVS){
    message("Strelka Filter SomaticEVS >= ",Phred_threshold,": ",sample_name)
    selected_snv <- VariantAnnotation::info(e.snv)$SomaticEVS >= Phred_threshold
    e.snv <- e.snv[selected_snv,]
  }

  rd <- SummarizedExperiment::rowRanges(e.snv)
  wt <- as.character(rd$REF)
  mt <- as.character(rd$ALT)

  selectSNV <- nchar(wt) == 1 & nchar(mt) == 1
  selectIndels <- !selectSNV

  if(sum(selectSNV)>0 & !is.null(vcfoutSNV)){
    message("Saving Strelka SNV to VCF file: ",vcfoutSNV)
    VariantAnnotation::writeVcf(e.snv[selectSNV,],vcfoutSNV,index=TRUE)
  }
  if(sum(selectIndels)>0 & !is.null(vcfoutIndels)){
    message("Saving Strelka Indels to VCF file: ",vcfoutIndels)
    VariantAnnotation::writeVcf(e.snv[selectIndels,],vcfoutIndels,index=TRUE)
  }
}


#' @export
mantaVcfToBedpe <- function(vcffile,bedpefile,sample_name,genomev,PR_threshold,modeGS,filterPASS){

  library("VariantAnnotation")
  #install.packages("devtools")
  #library(devtools)
  #install_github("PapenfussLab/StructuralVariantAnnotation@cb584be9474318e7f91c2cd8fca99f1212cab021") #works with R >=3.5.1
  library(StructuralVariantAnnotation)
  library(stringr)

  #tested with vcf obtained from Illumina Manta package
  #note that for the Sanger classification of rearrangements to work,
  #we need to reverse + and - in the partner strand (strand2) in the implementation below

  vcfToBedpe <- function(sv_VCF,sample_name){
    #code copied from https://github.com/PapenfussLab/gridss/blob/master/example/somatic.R
    gr <- StructuralVariantAnnotation::breakpointRanges(sv_VCF)

    bedpe <- data.frame(
      chrom1=seqnames(gr),
      start1=start(gr) - 1,
      end1=end(gr),
      chrom2=seqnames(partner(gr)),
      start2=start(partner(gr)) - 1,
      end2=end(partner(gr)),
      sample=sample_name,
      name=names(gr),
      partner=names(partner(gr)),
      score=gr$QUAL,
      strand1=strand(gr),
      strand2=strand(partner(gr)),
      stringsAsFactors = FALSE
    )

    svtype <- c()
    svclass <- c()
    for (i in 1:nrow(bedpe)){
      if(bedpe[i,"chrom1"]!=bedpe[i,"chrom2"]){
        svclass <- c(svclass,"translocation")
        svtype <- c(svtype,"BND")
      }else if(bedpe[i,"strand1"]==bedpe[i,"strand2"]){
        svclass <- c(svclass,"inversion")
        svtype <- c(svtype,"INV")
      }else if(bedpe[i,"strand1"]=="+"){
        svclass <- c(svclass,"deletion")
        svtype <- c(svtype,"DEL")
      }else if(bedpe[i,"strand1"]=="-"){
        svclass <- c(svclass,"tandem-duplication")
        svtype <- c(svtype,"DUP")
      }
    }
    bedpe[,"TYPE"] <- svtype
    bedpe[,"svclass"] <- svclass

    #for some reason, the partner strand needs to be reversed to get the Sanger version of the strand correct
    tmp <- bedpe$strand2
    tmp[bedpe$strand2=="+"] <- "-"
    tmp[bedpe$strand2=="-"] <- "+"
    bedpe$strand2 <- tmp

    #need only one row for each pair
    res_bedpe <- NULL
    for (i in 1:nrow(bedpe)){
      addrow <- FALSE
      #first check that if the mates are on the same chromosome than 1 comes before 2
      #(optional, but you never know)
      if((bedpe$chrom1[i]==bedpe$chrom2[i] & bedpe$start1[i]<bedpe$start2[i]) | bedpe$chrom1[i]!=bedpe$chrom2[i]){
        if(is.null(res_bedpe)){
          addrow <- TRUE
        }else{
          #check if the mate was already added
          if(length(intersect(c(bedpe$name[i],bedpe$partner[i]),res_bedpe$name))==0) addrow <- TRUE
          #if (addrow==TRUE) message("already present ",bedpe$name[i])
        }
        if(addrow) res_bedpe <- rbind(res_bedpe,bedpe[i,,drop=FALSE])
      }
    }
    return(res_bedpe)
  }


  #SV
  message("Load Manta SV file: ",sample_name)
  e.sv <- VariantAnnotation::readVcf(vcffile,genome = genomev)
  # e.sv <- VariantAnnotation::expand(sv_VCF)

  # make sure these are Manta rows only
  # rr <- SummarizedExperiment::rowRanges(e.sv)
  selection <- grepl(rownames(e.sv),pattern = "^Manta")
  e.sv <- e.sv[selection,]

  if(filterPASS){
    message("Filter PASS: ",sample_name)
    # In case this is GEL data we need to keep the MGE10kb as well
    selected_sv <- VariantAnnotation::fixed(e.sv)[,"FILTER"]=="PASS" | VariantAnnotation::fixed(e.sv)[,"FILTER"]=="MGE10kb"
    e.sv <- e.sv[selected_sv,]
  }

  #Manta filters, and PASS filter
  # PR <- sapply(1:nrow(e.sv),function(j) VariantAnnotation::geno(e.sv)$PR[j,2][[1]][2])
  if(modeGS=="G"){
    PR <- sapply(VariantAnnotation::geno(e.sv)$PR[,1],function(j) j[2])
  }else if(modeGS=="S"){
    PR <- sapply(VariantAnnotation::geno(e.sv)$PR[,2],function(j) j[2])
  }else{
    message("Unknown mode: use either G (germline) or S (somatic). Quit.")
    how_to()
    q(status=1,save = "no")
  }

  message("Filter PR >= ",PR_threshold,": ",sample_name)
  selected_sv <- PR >= PR_threshold
  e.sv <- e.sv[selected_sv,]

  #useful plots
  # if (do_PLOTS){
  #   jpeg(filename = paste0("../results/",sample_name,"_Strelka_REARR_PR.jpg"),width = 600,height = 600,res = 120)
  #   par(mfrow = c(1,1))
  #   hist(PR,breaks = 100,
  #        main = "Histogram of PR for REARR",xlab = "PR")
  #   dev.off()
  # }
  #select Chromosomes only
  message("Filter Chromosome Only: ",sample_name)
  expected_seqnames <- c(1:22,"X","Y")
  select_chrom <- as.character(GenomeInfoDb::seqnames(e.sv)) %in% expected_seqnames
  if(sum(select_chrom)==0){
    expected_seqnames <- paste0("chr",expected_seqnames)
    select_chrom <- as.character(GenomeInfoDb::seqnames(e.sv)) %in% expected_seqnames
  }
  e.sv <- e.sv[select_chrom,]

  #save the filtered vcf
  #writeVcf(e.sv,paste0(outdirs[i],sample_name,"_PR",PR_threshold,"_SV.vcf.gz"),index=TRUE)

  #now convert the vcf to bedpe using https://github.com/ctsa/svtools as recommended by Illumina
  #https://github.com/Illumina/manta/tree/master/docs/userGuide#converting-manta-vcf-to-bedpe-format

  #using "PapenfussLab/StructuralVariantAnnotation"
  #sv_VCF <- readVcf(paste0(outdirs[i],sample_name,"_PR",PR_threshold,"_SV.vcf.bgz"),genome = "hg19")
  sv_bedpe <- vcfToBedpe(e.sv,sample_name)

  write.table(sv_bedpe,file = bedpefile,sep = "\t",
              quote = FALSE,row.names = FALSE,col.names = TRUE)
  message("BEDPE saved: ",bedpefile)

}

#' @export
canvasVcfToTsv <- function(canvasVcf,tsvout,sample_name,genomev){

  #CNV
  message("Load Canvas VCF file: ",sample_name)
  e.cnv <- VariantAnnotation::readVcf(canvasVcf,genome = genomev)
  # ID based selection using rownames needs to donne before expand
  # because if the ALT has multiple entries then the expansion will
  # have the effect of removing the rownames

  # make sure these are Canvas rows only
  selection <- grepl(rownames(e.cnv),pattern = "^Canvas")
  e.cnv <- e.cnv[selection,]

  # don't expand to avoid duplication of copy number variants
  # e.cnv <- VariantAnnotation::expand(sv_VCF)

  # filter PASS
  selected_cnv <- VariantAnnotation::fixed(e.cnv)[,"FILTER"]=="PASS"
  e.cnv <- e.cnv[selected_cnv,]

  #select Chromosomes only
  message("Filter Chromosome Only: ",sample_name)
  expected_seqnames <- c(1:22,"X","Y")
  select_chrom <- as.character(GenomeInfoDb::seqnames(e.cnv)) %in% expected_seqnames
  if(sum(select_chrom)==0){
    expected_seqnames <- paste0("chr",expected_seqnames)
    select_chrom <- as.character(GenomeInfoDb::seqnames(e.cnv)) %in% expected_seqnames
  }
  e.cnv <- e.cnv[select_chrom,]

  # CNV size of at least 10kb
  rr <- SummarizedExperiment::rowRanges(e.cnv)
  sv_size <- VariantAnnotation::info(e.cnv)$END - BiocGenerics::start(rr)
  selection <- sv_size >= 10000
  e.cnv <- e.cnv[selection,]

  # construct the CNV df
  chr <- sapply(as.character(GenomeInfoDb::seqnames(e.cnv)),function(x) ifelse(grepl(x,pattern = "^chr"),substr(x,4,5),x))
  rr <- SummarizedExperiment::rowRanges(e.cnv)
  sample_MCC <- VariantAnnotation::geno(e.cnv)$MCC[,2]
  sample_MCC[is.na(sample_MCC)] <- VariantAnnotation::geno(e.cnv)$CN[is.na(sample_MCC),2]
  sv_df <- data.frame(Chromosome = chr,
                      chromStart = BiocGenerics::start(rr),
                      chromEnd = VariantAnnotation::info(e.cnv)$END,
                      total.copy.number.inNormal = rep(2,nrow(e.cnv)),
                      minor.copy.number.inNormal = rep(1,nrow(e.cnv)),
                      total.copy.number.inTumour = VariantAnnotation::geno(e.cnv)$CN[,2],
                      minor.copy.number.inTumour = VariantAnnotation::geno(e.cnv)$CN[,2] - sample_MCC,
                      stringsAsFactors = F)

  message("Saving CNV to tab file: ",tsvout)
  write.table(sv_df,file = tsvout,sep = "\t",
              col.names = TRUE,row.names = FALSE,quote = FALSE)
}

