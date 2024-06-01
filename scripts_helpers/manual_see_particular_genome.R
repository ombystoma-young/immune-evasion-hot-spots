#GCA_003867335.1

library(gggenomes)
library(tidyverse)


setwd('work_dir/anti_defence/anti_defence_pipeline/')
genome <- 'GCA_020496145.1'
seq_1 <- read_seqs('ncbi_dataset/data/GCA_020496145.1/GCA_020496145.1_ASM2049614v1_genomic.fna')
g_1 <- read_feats(paste0('annotation/prokka/', genome, '/', genome, '.gff'))
p <- gggenomes(seqs=seq_1, genes=g_1) +
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=5, vjust = 2) + # label each sequence by this caption
  geom_gene(aes(fill = product), data=genes(.gene_types = c('CDS')), intron_shape=0, size=5) +  # add gene arrows
  # geom_gene_tag(aes(label=product), size=5, check_overlap=TRUE, vjust=-1, angle=25) +
  geom_gene_text(aes(label=locus_tag), data=genes(.gene_types = c('gene')), size=5, check_overlap=TRUE, vjust=3) +
  #scale_fill_brewer("Function", palette="Dark2", na.value="cornsilk3") +   # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=18))  # change font size and legend position
p 
