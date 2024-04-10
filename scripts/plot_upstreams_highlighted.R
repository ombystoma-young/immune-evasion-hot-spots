library(gggenomes)
library(tidyverse)
setwd('PycharmProjects/anti_defence_pipeline/')

g <- read_feats('data_autographiviridae_meta/annotation/concatenated_sele.gff')
t <- read.table('data_autographiviridae_meta/tdrs/best_tdrs.tsv', col.names = c('seq_id',
                                                      'start.tdr', 
                                                      'end.tdr',
                                                      'dist.tdr',
                                                      'start2.tdr',
                                                      'end2.tdr'))
s <- read_seqs('data_autographiviridae_meta/filtered_mags_flatten/concatenated_genomes.fna')
# sele <- read_table('metadata/filtered_upstreams_nuccore.id', col_names = c('seq_id'))
sele <- data.frame(seq_id = c('MGV-GENOME-0268257', 'MGV-GENOME-0297924', 'MGV-GENOME-0219685'))
g <- g %>% filter(seq_id %in% sele$seq_id)
s <- s %>% filter(seq_id %in% sele$seq_id)
t <- t %>% filter(seq_id %in% sele$seq_id)
t <- bind_rows(
    select(t, seq_id=seq_id, start=start.tdr, end=end.tdr, dist.tdr = dist.tdr),
    select(t, seq_id=seq_id, start=start2.tdr, end=end2.tdr, dist.tdr = dist.tdr))

early_coords <- read_table('data_autographiviridae_meta/upstreams/early_meta.bed', 
                           col_names = c('seq_id', 'start', 'end', 1:9))
early_coords <- early_coords %>% mutate(record_type = 'early_coords')
t <- t %>% mutate(record_type = 'tdrs')
# f <- early_coords
f <- bind_rows(
  select(early_coords, seq_id=seq_id, start=start, end=end, record=record_type),
  select(t, seq_id=seq_id, start=start, end=end, record=record_type))


g <- g %>% mutate(rnap = !is.na(pfam))
s <- s %>% filter(length < 50000)

test <- gggenomes(genes=g, seqs=s, feats = f)  %>%
  flip(1) +
  # focus(.track_id = feats, .expand = 3000) +
  geom_seq() +
  geom_seq_label(aes(label=seq_id), size = 4) +
  #geom_seq_label() + 
  geom_feat(aes(color=`record`), size=10, alpha=0.5) +
  geom_gene(aes(fill=rnap), intron_shape=0, size = 5) +
  # geom_gene_tag(aes(label=feat_id)) +
  theme(legend.position = 'none') +
  theme(text = element_text(size = 12), axis.text = element_text(size = 15)) + 
  scale_fill_manual(values = c('grey', 'tomato'))
test
ggsave('test_meta.pdf', test, limitsize = FALSE, height = 200, width = 20)
