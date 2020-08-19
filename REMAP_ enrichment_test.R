# https://github.com/remap-cisreg/ReMapEnrich/blob/master/vignettes/basic_use.md 19/aug/2020 chat with walter

library(ReMapEnrich)
all_record <- dataBrowser(cell="K562")


                                                                      
cebpb_k562 <- loadPeaks(id = "MM1_HSA_K562_CEBPB", includeMotifOnly = TRUE)
cebpb_k562$id<-gsub("MM1_HSA_K562_|_peaks_with_motif_.+","",cebpb_k562$id,perl=TRUE)
colnames(cebpb_k562) <- c("chrom","chromStart","chromEnd", "name", "score")

################################################################################
cell_ID <- dplyr::pull(cell_TFBS,ID)
All_TF_cell_peaks <- TFregulomeR::loadPeaks(id = cell_ID, includeMotifOnly = includeMotifOnly) ;


catalog_ID <- lapply(cell_ID, function(ID){
  TFregulomeR::loadPeaks(id = ID, includeMotifOnly = includeMotifOnly) ;
}) 

 catalog_ID <- data.table::rbindlist(catalog_ID)

catalog_ID$id<-gsub("MM1_HSA_K562_|_peaks_with_motif_.+","",catalog_ID$id,perl=TRUE)
colnames(catalog_ID) <- c("chrom","chromStart","chromEnd", "name", "score")

######
# create GRanges
createGR_from_df <- function(df_object){
  bedData <- df_object
  # If the data frame has more than 6 columns then remove them.
  if(ncol(bedData) > 6)
    bedData <- bedData[,-c(7:ncol(bedData))]
  # If the data frame has less than 3 columns then throw an error.
  if(ncol(bedData)<3)
  stop("File has less than 3 columns")
  # If the strand is known in the data frame then replace it by the
  # granges equivalent.
  if ('strand' %in% colnames(bedData))
    bedData$strand <- gsub(pattern= "[^+-]+", replacement = '*',
                         x = bedData$strand)
  # Construct the grangesanges object depending on the number of columns.
  if (ncol(bedData) == 3) {
    granges <- with(bedData, GenomicRanges::GRanges(chrom, IRanges::IRanges(chromStart, chromEnd)))
  } else if (ncol(bedData)==4) {
    granges = with(bedData, GenomicRanges::GRanges(chrom, IRanges::IRanges(chromStart, chromEnd),
                                                 id = name))
  } else if (ncol(bedData)==5) {
    granges <- with(bedData, GenomicRanges::GRanges(chrom, IRanges::IRanges(chromStart, chromEnd),
                                                  id = name, score = score))
  } else if (ncol(bedData)==6) {
    granges <- with(bedData, GenomicRanges::GRanges(chrom, IRanges::IRanges(chromStart, chromEnd),
                                                  id = name, score = score, strand = strand))
  }
  return(granges)

}
########
cebpb_gr<-createGR_from_df(cebpb_k562)
catalog_ID_gr<-createGR_from_df(catalog_ID)

enrichment.df <- enrichment(cebpb_gr, catalog_ID_gr,shuffle=300)


