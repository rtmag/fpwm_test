# https://github.com/remap-cisreg/ReMapEnrich/blob/master/vignettes/basic_use.md 19/aug/2020 chat with walter

library(ReMapEnrich)
all_record <- dataBrowser(cell="K562")


                                                                      
cebpb_k562 <- loadPeaks(id = "MM1_HSA_K562_CEBPB", includeMotifOnly = TRUE)
                                                                      
cell_ID <- dplyr::pull(cell_TFBS,ID)
All_TF_cell_peaks <- TFregulomeR::loadPeaks(id = cell_ID, includeMotifOnly = includeMotifOnly) ;

####
catalog_ID <- lapply(cell_ID, function(ID){
  TFregulomeR::loadPeaks(id = ID, includeMotifOnly = includeMotifOnly) ;
}) 

 catalog_ID <- data.table::rbindlist(catalog_ID)

######


bed_to_granges <- function(file){
   df <- read.table(file,
                    header=F,
                    stringsAsFactors=F)
 
   if(length(df) > 6){
      df <- df[,-c(7:length(df))]
   }
 
   if(length(df)<3){
      stop("File has less than 3 columns")
   }
 
   header <- c('chr','start','end','id','score','strand')
   names(df) <- header[1:length(names(df))]
 
   if('strand' %in% colnames(df)){
      df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
   }
 
   library("GenomicRanges")
 
   if(length(df)==3){
      gr <- with(df, GRanges(chr, IRanges(start, end)))
   } else if (length(df)==4){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
   } else if (length(df)==5){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
   } else if (length(df)==6){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
   }
   return(gr)
}
