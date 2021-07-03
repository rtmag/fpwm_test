devtools::load_all("/home/roberto/forkedTF_May2021/forkedTF")

miniCofactorReport(TF = "JUND",cell = "HepG2", filterBy="fraction",           pdfName = "test_fraction_2.pdf" , server="ca")
miniCofactorReport(TF = "JUND",cell = "HepG2", filterBy="mapped.peaks.ratio", pdfName = "test_peaks_ratio.pdf", server="ca")
miniCofactorReport(TF = "JUND",cell = "HepG2", filterBy="q.significance",     pdfName = "test_qsignificance.pdf", server="ca")

miniCofactorReport(TF = "JUND",cell = "HepG2", filterBy="fraction",           pdfName = "test_fraction_2.pdf" , server="ca")
miniCofactorReport(TF = "JUND",cell = "HepG2", filterBy="mapped.peaks.ratio", pdfName = "test_peaks_ratio2.pdf", server="ca")
miniCofactorReport(TF = "JUND",cell = "HepG2", filterBy="q.significance",     pdfName = "test_qsignificance2.pdf", server="ca")


 scp -i ~/keys/motif_key.pem ubuntu@206.12.92.176:/home/roberto/*pdf .
