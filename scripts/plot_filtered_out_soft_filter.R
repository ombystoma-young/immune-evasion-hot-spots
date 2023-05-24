# plot filtered out
setwd('work_dir/anti_defence/anti_defence_pipeline/')

library(gggenomes)
library(dplyr)

seq_1 <- read_seqs('ncbi_dataset/data/GCA_029376425.1/GCA_029376425.1_ASM2937642v1_genomic.fna')
seq_2 <- read_seqs('ncbi_dataset/data/GCF_003440975.1/GCF_003440975.1_ASM344097v1_genomic.fna')

g_1 <- read_feats('ncbi_dataset/data/GCA_029376425.1/genomic.gff')
g_2 <- read_feats('ncbi_dataset/data/GCF_003440975.1/genomic.gff')

seqs <- rbind(seq_1, seq_2)
genes <- merge(g_1, g_2, all.y = TRUE, all.x = TRUE)


p <- gggenomes(seqs=seqs, genes=genes) +
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=5, vjust = 2) + # label each sequence by this caption
  geom_gene(aes(fill = product), data=genes(.gene_types = c('CDS')), intron_shape=0, size=5) +  # add gene arrows
 # geom_gene_tag(aes(label=product), size=5, check_overlap=TRUE, vjust=-1, angle=25) +
  geom_gene_text(aes(label=locus_tag), data=genes(.gene_types = c('gene')), size=5, check_overlap=TRUE, vjust=3) +
  #scale_fill_brewer("Function", palette="Dark2", na.value="cornsilk3") +   # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=18))  # change font size and legend position
p 
