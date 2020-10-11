library(TFregulomeR)
library(forkedTF)

all_record <- dataBrowser()

mm_table <- all_record[all_record$source=="MethMotif",]

mm_table$cell_tissue_name<- gsub(" ","",mm_table$cell_tissue_name)

for( i in 556:dim(mm_table)[1] ){
  miniCofactorReport( TF = mm_table$TF[i], cell = mm_table$cell_tissue_name[i],cobinding_threshold=.1 )
}



###  TEST CASES
# HepG2 USF2 cofactor: TFE3

# HepG2 USF1 cofactor: TFE3

# HCT116 ELF1 cofactor: FOSL1

# A549 SREBF1 cofactor: JUNB
