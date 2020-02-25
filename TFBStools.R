library(TFBSTools)
library(TFregulomeR)

CEBPB_pfm <- toTFBSTools(id = "MM1_HSA_K562_CEBPB")
ATF4_pfm <- toTFBSTools(id = "MM1_HSA_K562_ATF4")
CEBPD_pfm <- toTFBSTools(id = "MM1_HSA_K562_CEBPD")

# CEBPB + ATF4
CEBPB_ATF4 <- TFregulomeR::intersectPeakMatrix(peak_id_x = "MM1_HSA_K562_CEBPB", 
					  motif_only_for_id_x = TRUE, 
					  peak_id_y = "MM1_HSA_K562_ATF4", 
					  motif_only_for_id_y = TRUE)

motifTmp <- CEBPB_ATF4[1,1][[1]]

CEBPB_ATF4_pfm <- TFBSTools::PFMatrix(ID = motifTmp@id, 
			   name = motifTmp@MethMotif_x@MMmotif@id, 
			   strand = "*",
                           bg = motifTmp@MethMotif_x@MMmotif@background, 
			   profileMatrix = t(round(motifTmp@MethMotif_x@MMmotif@motif_matrix * motifTmp@MethMotif_x@MMmotif@nPeaks))
			  )

# CEBPB + CEBPD
CEBPB_CEBPD <- TFregulomeR::intersectPeakMatrix(peak_id_x = "MM1_HSA_K562_CEBPB", 
					  motif_only_for_id_x = TRUE, 
					  peak_id_y = "MM1_HSA_K562_CEBPD", 
					  motif_only_for_id_y = TRUE)

motifTmp <- CEBPB_CEBPD[1,1][[1]]

CEBPB_CEBPD_pfm <- TFBSTools::PFMatrix(ID = motifTmp@id, 
			   name = motifTmp@MethMotif_x@MMmotif@id, 
			   strand = "*",
                           bg = motifTmp@MethMotif_x@MMmotif@background, 
			   profileMatrix = t(round(motifTmp@MethMotif_x@MMmotif@motif_matrix * motifTmp@MethMotif_x@MMmotif@nPeaks))
			  )
