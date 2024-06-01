library(gggenomes)
library(tidyverse)
setwd('PycharmProjects/anti_defence_pipeline/')

# g <- read_feats('data_autographiviridae_meta/annotation/concatenated_sele.gff')
g <- read_feats('data_autographiviridae_refseq/annotation/concatenated_with_clusters_phrogs.gff')

s <- read_seqs('data_autographiviridae_refseq/genomes/concatenated_genomes.fna')

sele <- data.frame(seq_id = c('NC_047751.1'))
# sele <- s[round(runif(50, min=1, max=nrow(s))), ]
g1 <- g %>% filter(seq_id %in% sele$seq_id)
s1 <- s %>% filter(seq_id %in% sele$seq_id)

early_coords <- read_table('data_autographiviridae_meta/upstreams/early_meta.bed', 
                           col_names = c('seq_id', 'start', 'end', 1:9))
early_coords <- early_coords %>% mutate(record_type = 'early_coords')


# g1 <- g1 %>% mutate(rnap = !is.na(pfam))
g1 <- g1 %>% mutate(rnap = str_detect(phro_gs, 'phrog_414'))
s1 <- s1 %>% filter(length < 50000)

test <- gggenomes(genes=g1, seqs=s1, feats = early_coords)  %>%
  flip(1) +
  # focus(.track_id = feats, .expand = 3000) +
  geom_seq() +
  geom_seq_label(aes(label=seq_id), size = 4) +
  geom_feat(aes(color=`record_type`), size=10, alpha=0.5) +
  geom_gene(aes(fill=phrog_category), intron_shape=0, size = 5) +
  # theme(legend.position = 'none') +
  theme(text = element_text(size = 12), axis.text = element_text(size = 15)) #+ 
  # scale_fill_manual(values = c('grey', 'tomato'))
test
ggsave('test_meta.pdf', test, limitsize = FALSE, height = 200, width = 20)
ggsave('test_meta_short.pdf', test, limitsize = FALSE, height = 30, width = 20)
