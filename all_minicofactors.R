library(TFregulomeR)
library(forkedTF)

all_record <- dataBrowser()

mm_table <- all_record[all_record$source=="MethMotif",]

mm_table$cell_tissue_name<- gsub(" ","",mm_table$cell_tissue_name)

for( i in 556:dim(mm_table)[1] ){
  miniCofactorReport( TF = mm_table$TF[i], cell = mm_table$cell_tissue_name[i],cobinding_threshold=.1 )
}

