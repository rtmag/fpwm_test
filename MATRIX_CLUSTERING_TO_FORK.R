library('TFregulomeR')
cebpa <- loadPeaks(id = "MM1_HSA_HepG2_CEBPA", includeMotifOnly = TRUE)
cebpa[,2] <- cebpa[,2]-50
cebpa[,3] <- cebpa[,3]+50
write.table(cebpa[,1:3],"HepG2_CEBPA.bed",sep="\t",quote=F,row.names=F,col.names=F)
