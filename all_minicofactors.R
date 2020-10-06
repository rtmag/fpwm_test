library(TFregulomeR)

all_record <- dataBrowser()

mm_table <- all_record[all_record$source=="MethMotif",]

for( i in 555:dim(mm_table)[1] ){
  miniCofactorReport( TF = mm_table$TF[i], cell = mm_table$cell_tissue_name[i],cobinding_threshold=.1 )
}

