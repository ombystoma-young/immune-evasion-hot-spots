setwd('work_dir/anti_defence/anti_defence_pipeline/')

library(gggenomes)
library(dplyr)

seq_t3 <- read_seqs('ncbi_dataset/data/GCF_000841665.1/GCF_000841665.1_ViralProj14336_genomic.fna')
seq_t7 <- read_seqs('ncbi_dataset/data/GCF_000844825.1/GCF_000844825.1_ViralProj14460_genomic.fna')

# cat ncbi_dataset/data/GCF_000841665.1/genomic.gff ncbi_dataset/data/GCF_000844825.1/genomic.gff > temp_files/for_initial_pic.gff

vis <- c("CDS")

seqs <- rbind(seq_t3, seq_t7)
genes <- read_feats('temp_files/for_initial_pic.gff')

feat <- genes %>% filter(end < 1400) %>% 
          filter(type == 'CDS') %>% 
          filter(product != 'hypothetical protein')

antir <- genes %>% 
  filter(type == 'CDS') %>% 
  filter(locus_tag == 'T3p01' | locus_tag == 'T3p02' | locus_tag == 'T3p01' | locus_tag =='T7p01') 
        


p <- gggenomes(seqs=seqs, genes=genes, feats = feat) %>% 
  focus(.track_id = feats, .expand = 4620, .locus_id=seq_id, .overhang='drop') +
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(size=5) + # label each sequence by this caption
  geom_gene(aes(fill = product), data=genes(.gene_types = vis), intron_shape=0, size=5) +  # add gene arrows
  geom_gene_tag(aes(label=product), size=5, check_overlap=TRUE, vjust=-1, angle=25) +
  geom_gene_text(aes(label=locus_tag), data=genes(.gene_types = c('gene')), size=5, check_overlap=TRUE, vjust=3) +
  scale_fill_brewer("Function", palette="Dark2", na.value="cornsilk3") +   # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=18))  # change font size and legend position
p 
ggsave('pics/t7_and_t3.pdf', p, dpi=600, width = 21, height = 8)

# 
# feat <- genes %>% filter(end < 5950) %>% 
#   filter(type == 'CDS') %>% 
#   filter(product != 'hypothetical protein')
# 
# 
# p <- gggenomes(seqs=seqs, genes=genes, feats = feat) %>% 
#   focus(.track_id = feats, .expand = 0, .locus_id=seq_id, .overhang='drop') +
#   geom_seq() +  # draw contig/chromosome lines
#   geom_seq_label(size=5) + # label each sequence by this caption
#   geom_gene(data=genes(.gene_types = vis), intron_shape=0, size=5) +  # add gene arrows
#   geom_feat_tag(aes(label=product), size=5, check_overlap=TRUE, vjust=-1, angle=25) +
#   geom_gene_text(aes(label=locus_tag), data=genes(.gene_types = c('gene')), size=5, check_overlap=TRUE, vjust=3) +
#   scale_fill_brewer("Function", palette="Dark2", na.value="cornsilk3") +   # change fill, genes
#   theme(legend.position = 'bottom', axis.text.x = element_text(size=18))  # change font size and legend position
# p 
# 
# 
