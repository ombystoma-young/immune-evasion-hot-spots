library(gggenomes)
library(tidyverse)
setwd('PycharmProjects/anti_defence_pipeline/')

g <- read_feats('data/annotation/concatenated.gff')
t <- read.table('data/tdrs/best_tdrs.tsv', col.names = c('seq_id',
                                                      'start.tdr', 
                                                      'end.tdr',
                                                      'dist.tdr',
                                                      'start2.tdr',
                                                      'end2.tdr'))
s <- read_seqs('data/genomes/concatenated_genomes.fna')
sele <- read_table('metadata/filtered_upstreams_nuccore.id', col_names = c('seq_id'))
g <- g %>% filter(seq_id %in% sele$seq_id)
s <- s %>% filter(seq_id %in% sele$seq_id)
t <- t %>% filter(seq_id %in% sele$seq_id)
t <- bind_rows(
    select(t, seq_id=seq_id, start=start.tdr, end=end.tdr, dist.tdr = dist.tdr),
    select(t, seq_id=seq_id, start=start2.tdr, end=end2.tdr, dist.tdr = dist.tdr))

early_coords <- read_table('data/upstreams/early.bed', 
                           col_names = c('seq_id', 'start', 'end', 1:9))
early_coords <- early_coords %>% mutate(record_type = 'early_coords')
t <- t %>% mutate(record_type = 'tdrs')

f <- bind_rows(
  select(early_coords, seq_id=seq_id, start=start, end=end, record=record_type),
  select(t, seq_id=seq_id, start=start, end=end, record=record_type))


g <- g %>% mutate(rnap = ifelse(str_detect(phro_gs, 'phrog_414'), TRUE, FALSE))

test <- gggenomes(genes=g, seqs=s, feats = f) %>%
  focus(.track_id = feats, .expand = 2000) +
  geom_seq() +
  geom_seq_label(aes(label=seq_desc), size = 4) +
  #geom_seq_label() + 
  geom_feat(aes(color=`record`), size=10, alpha=0.5) +
  geom_gene(aes(fill=rnap), intron_shape=0, size = 5) +
  theme(legend.position = 'none') +
  theme(text = element_text(size = 12), axis.text = element_text(size = 15))

ggsave('test.pdf', test, limitsize = FALSE, height = 200, width = 20)
