gene_annot[match("MRA_0223", gene_annot$Locus.tag),]
MRA_0223 = substr(genome, 259145, 260218)
dir.create(file.path("/media/rachel/Part2/OneDrive - south.du.ac.in/project_jnu/RNAT/2019/h37ra","unique_RNAT"))
writeLines(paste(">MRA_0223", MRA_0223, sep = "\n"), "/media/rachel/Part2/OneDrive - south.du.ac.in/project_jnu/RNAT/2019/h37ra/unique_RNAT/MRA_0223.fasta")

gene_annot[match("MRA_2380", gene_annot$Locus.tag),]
MRA_2380 = substr(genome, 2649656, 2651503)
writeLines(paste(">MRA_2380", MRA_2380, sep = "\n"), "/media/rachel/Part2/OneDrive - south.du.ac.in/project_jnu/RNAT/2019/h37ra/unique_RNAT/MRA_2380.fasta")
