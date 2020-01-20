library(readxl)
setwd("/media/rachel/Part2/OneDrive - south.du.ac.in/project_jnu/RNAT/2019/")
dir.create("h37rv")
# Read genome as a character sequence
genome = scan("genomes and annotations/h37rv.txt", what = "character", sep = NULL)
str(genome)

# Get genome annotation
gene_annot = readxl::read_xls("genomes and annotations/gff anf fa/ProteinTable166_159857_h37rv.xls")
gene_annot = read.csv(file = "genomes and annotations/genome_annotations_h37rv.csv")

# Order the table according to sense/antisense strand
ordered_gene_annot = data.frame(gene_annot[order(gene_annot$Strand),])
sense = subset(ordered_gene_annot, Strand == "+")
antisense = subset(ordered_gene_annot, Strand == "-")

# Get start stop positions for extraction
u_start_pos = sense$Start
u_stop_pos = sense$Stop

# Start position for sense strand is taken as 250 bases upstream of the start position
# whereas stop position is taken as 4 bases downstream of the start position (ATG + 2)
u_cut_start = u_start_pos - 35
u_cut_stop = u_start_pos + 4

# Get details for the genes
upstr = data.frame(sense[,3:12])
upstr = upstr[,-7]
upstr = upstr[,-7]
upstr = upstr[,-7]

upstr_seq = vector(mode = "character", length = nrow(sense))
for(i in 1:nrow(sense)){upstr_seq[i] = substr(genome, u_cut_start[i], u_cut_stop[i])}
u_seq = data.frame(upstr_seq)
upst_fin = cbind(upstr, upstr_seq)
write.csv(upst_fin, file = "./h37rv/h37rv-details-sense_35.txt", col.names = TRUE)
h37rv_fasta = data.frame(upst_fin$Locus.tag, upst_fin$upstr_seq)
write.csv(h37rv_fasta, file = "./h37rv/h37rv_fasta_sense_35.fasta", col.names = T)

# Start position for antisense strand is taken as 4 bases upstream of the start position (ATG + 2)
# whereas stop position is taken as 250 bases downstream of the start position 

d_start_pos = antisense$Stop
d_cut_start = d_start_pos - 4
d_cut_stop = d_start_pos + 35
dstr = data.frame(antisense[,3:12])
dstr = dstr[,-7]
dstr = dstr[,-8]
dstr_seq = vector(mode = "character", length = nrow(antisense))
for(i in 1:nrow(antisense)){dstr_seq[i] = substr(genome, d_cut_start[i], d_cut_stop[i])}
d_seq = data.frame(dstr_seq)
dst_fin = cbind(dstr, dstr_seq)
write.csv(dst_fin, file = "./h37rv/h37rv-details-antisense_35.txt", col.names = TRUE)
h37rv_fasta_downstream = data.frame(dst_fin$Locus.tag, dst_fin$dstr_seq)
write.csv(h37rv_fasta_downstream, file = "./h37rv/h37rv_fasta_antisense_35.fasta", col.names = TRUE)

# CONVERT ANTISENSE STRAND SEQUENCES TO REVERSE COMPLEMENT

setwd("h37rv/")
sequences = read.csv("h37rv.fasta", sep = "\t", header = F, stringsAsFactors = F)
colnames(sequences) = "seq"
seq_conformation = sequences
seq_conformation$conformation_20 = "open"
seq_conformation$conformation_25 = "open"
seq_conformation$conformation_37 = "open"
seq_conformation$conformation_40 = "open"
seq_conformation$conformation = "NULL"

flag_count_20 = data.frame(sequences)
flag_count_20$conf = "NULL"
flag_count_20$conf_25 = "NULL"
flag_count_20$conf_37 = "NULL"
flag_count_20$conf_40 = "NULL"

i=3
j=1
flag = "open"

for (i in seq(3,nrow(sequences),2)) 
{
  if(grepl("AGGAG|GGAGG",sequences[i+1,1],perl = T))
    {
      writeLines(paste(sequences[i,1],sequences[i+1,1], sep = "\n"), "temp_seq.fasta")
      system(command = "mfold SEQ=temp_seq.fasta NA=RNA T=40")
      system(command = "perl ../Ct2B.pl temp_seq.ct > dotb_test.txt")
      seq_dot_bra = read.csv(file = "dotb_test.txt", sep = "\t", row.names = NULL, header = F, stringsAsFactors = F)
      flag_mix = "NULL"
      for (k in 2:nrow(seq_dot_bra)) 
      {
        sd_start = unlist(gregexpr("GGAG", seq_dot_bra[1,1]))
        substr(seq_dot_bra[k,1], sd_start, sd_start+3)
        if(length(grep("[()]", substr(seq_dot_bra[k,1], sd_start, sd_start+3), perl = T)))
        {
         # write(sequences[i,], file = "RNAT_result_2019-Feb_40.txt", append = T)
         # writeLines(paste(sequences[i,1],sequences[i+1,1], sep = "\n"), con = paste("./rnat-sequences/fasta/fasta_40/",paste(unlist(strsplit(trimws(sequences[i,1]),">"))[[2]],".fasta", sep = ""), sep = ""))
         # write(x = paste(sequences[i,1],"closed", sep = "\t"), file = "./rnat-sequences/conformation/conformation_40.txt", append = T)
          flag = "closed"
        }
        flag_mix = append(flag_mix, flag)
        flag = "open"
      }
      flag_mix = flag_mix[-1]
      flag_count_20[i,5] = paste(unlist(flag_mix), collapse = "_")
      lm = length(levels(factor(flag_mix)))
      if (lm > 1) 
        {
          flag = "mix"
      }
      else
      flag = levels(factor(flag_mix))
      seq_conformation[i,5] = flag
  }
  flag = "open"
}

flag_count = flag_count_20[match(RNATs$seq, flag_count_20$seq),]
flag_count_prob = data.frame(flag_count$seq)
flag_count_prob$t20 = "NULL"
flag_count_prob$t25 = "NULL"
flag_count_prob$t37 = "NULL"
flag_count_prob$t40 = "NULL"
library(stringr)
for (i in seq(3,nrow(flag_count_prob),2)) {
  for (k in 2:5) {
    str_flag = strsplit(flag_count[i,k],split = "_")
    numer_l = str_count(str_flag,"open")
    denom_l = length(unlist(str_flag))
    flag_count_prob[i,k] = as.numeric(numer_l/denom_l)
  }
}
library(circlize)
hm_h37rv = flag_count_prob[seq(3,nrow(flag_count_prob),2),]
rownames(hm_h37rv) = hm_h37rv$flag_count.seq
hm_h37rv = hm_h37rv[,-1]
sapply(hm_h37rv, class)
cols.num <- c("t20", "t25", "t37", "t40")
hm_h37rv[cols.num] <- sapply(hm_h37rv[cols.num],as.numeric)
sapply(hm_h37rv, class)
write.csv(hm_h37rv, "/media/rachel/Part2/OneDrive - south.du.ac.in/project_jnu/RNAT/2019/heatmap/h37rv_hm_data_21_June_2019.csv")
f1 = colorRamp2(seq(min(hm_h37rv), max(hm_h37rv), length = 2), c("red", "black"), space = "LAB")
library(ComplexHeatmap)
Heatmap(hm_h37rv, cluster_columns = F, cluster_rows = T, col = f1, row_names_gp = gpar(fontsize = 4))



#CONDITION ONE
for (i in seq(3,nrow(seq_conformation),2)) 
{
  l = as.character(seq_conformation[i,2:5])
  f = length(levels(factor(l)))
  seq_conformation[i,6] = f
  seq_conformation[i+1,6] = f
 }

RNATs = subset(seq_conformation, conformation > 1)
write.csv(RNATs, file = "RNATs_h37rv.csv")

for (i in seq(3,nrow(sequences),2)) 
{
  if(grepl("AGGAG|GGAGG",sequences[i+1,1],perl = T))
  {
    write(paste(sequences[i,1],sequences[i+1,1], sep = "\t"), file = "/media/rachel/Part2/OneDrive - south.du.ac.in/project_jnu/RNAT/2019/h37rv/sequences_with_sd_region.txt", append = T)
  }
}

combined_hm = read.csv("./../heatmap/combined_hm_data_21_June_2019.csv", row.names = 1)
sapply(combined_hm, class)
f1 = colorRamp2(seq(min(combined_hm), max(combined_hm), length = 2), c("red", "black"), space = "LAB")
library(ComplexHeatmap)
Heatmap(combined_hm, cluster_columns = F, cluster_rows = T, col = f1, row_names_gp = gpar(fontsize = 4))
