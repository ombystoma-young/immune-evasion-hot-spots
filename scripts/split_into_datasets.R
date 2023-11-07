library(tidyverse)

# setwd('work_dir/anti_defence/anti_defence_pipeline/')
# tdrs_tsv <- 'data/tdrs/tdrs.tsv'
# intergenic_bed <- 'data/intergenics/all_intergenic_with_length.tsv'
# rnaps_gff <- 'data/annotation/concatenated_rnaps_only.gff'
# lengths_bed <- 'data/intergenics/chromosome_lengths.bed'
# 
# max.dist.tdr <- 10000
# max.dist.inter <- 5000
# min.dist.inter <- 1000
# max.gen <- 65000 # k
# min.len.inter <- 500  # k
# max.len.inter <- 3000# k


args <- commandArgs(trailingOnly=TRUE)

tdrs_tsv <- args[1]
intergenic_bed <- args[2]
rnaps_gff <- args[3]
lengths_bed <- args[4]
max.dist.tdr <- as.numeric(args[5])
max.dist.inter <- as.numeric(args[6])
min.dist.inter <- as.numeric(args[7])
max.gen <- as.numeric(args[8])  # k
min.len.inter <- as.numeric(args[9])  # k
max.len.inter <- as.numeric(args[10])  # k
output.txt <- args[11]
output.tdrs <- args[12]
output.igs <- args[13]


# read data
## TDRs
tdrs <- read.table(tdrs_tsv, 
                   col.names = c('seq_id', 'start', 'end', 'strand', 
                                 'start2', 'end2', 'map_match', 
                                 'map_length', 'de')) %>% 
  select(seq_id, start, end, start2, end2)

## rnaps
rnaps <- read.table(rnaps_gff, col.names = c('seq_id', 'source', 'type', 
                                             'start', 'end',
                                             'score', 'strand', 
                                             'phase', 'attribute')) %>% 
  select(seq_id, start, end, strand)

## intergenics
intergenics <- read.table(intergenic_bed, 
                          col.names = c('seq_id', 'start', 'end', 'length.ig'))
## chromosome lengths 
chr_lens <- read.table(lengths_bed, 
                       col.names = c('seq_id', 'start', 'seq_len')) %>% 
  select(-start)



# define outliers based on genome length
too_long_genome <- chr_lens %>% filter(seq_len > max.gen) %>% pull(seq_id) %>% unique()




# define distances

## tdrs
rnap_tdr_join <- rnaps %>% left_join(tdrs, by='seq_id',
                                     suffix = c('.rnap', '.tdr')) %>% 
                            left_join(chr_lens)

rnap_tdr_join <- rnap_tdr_join %>% 
  mutate(dist.tdr = case_when(
                    strand == '+' & 
                      start.rnap > end.tdr & end.rnap < start2 ~ 
                      start.rnap - end.tdr,
                    strand == '+' & 
                      end.rnap < start.tdr & end.rnap < start2 ~ 
                      seq_len - end2 + start.rnap, 
                    strand == '+' & 
                      start.rnap > end.tdr & start.rnap > end2 ~ 
                      start.rnap - end2,
                    strand == '-' &  
                      start.rnap > end.tdr & end.rnap < start2 ~ 
                      start2 - end.rnap,
                    strand == '-' & 
                      end.rnap < start.tdr & end.rnap < start2 ~ 
                      start.tdr - end.rnap,
                    strand == '-' & 
                      start.rnap > end.tdr & start.rnap > end2 ~ 
                      seq_len - end.rnap + start.tdr,
                    .default = NA
                            )
         )

rnap_tdr_join <- rnap_tdr_join %>%
  dplyr::group_by(seq_id) %>% 
  mutate(tdr.r = rank(dist.tdr, na.last='keep')) %>% 
  ungroup()

rnap_tdr_join <- rnap_tdr_join %>%   
  select(seq_id, start.tdr, end.tdr, start2,  end2,
                                 dist.tdr, tdr.r)
colnames(rnap_tdr_join) <- c('seq_id', 'start.tdr', 'end.tdr', 
                             'start2.tdr',  'end2.tdr',
                              'dist.tdr', 'r.tdr')


## intergenics 
rnap_igs_join <- rnaps %>% left_join(intergenics, by='seq_id',
                                             suffix = c('.rnap', '.ig')) %>% 
                           left_join(chr_lens)

rnap_igs_join <- rnap_igs_join %>% 
                                mutate(dist.ig = case_when(
                                  strand == '+' & 
                                    start.rnap > start.ig & start.rnap > end.ig ~ 
                                    start.rnap - end.ig,
                                  strand == '+' & 
                                    end.rnap <= start.ig & end.rnap < end.ig ~ 
                                    seq_len - end.ig + start.rnap,
                                  strand == '+' & 
                                    end.rnap < start.ig & start.rnap > end.ig ~ 
                                    start.rnap - end.ig,
                                  strand == '-' & 
                                    start.rnap > start.ig & start.rnap >= end.ig  ~ 
                                    seq_len - end.rnap + start.ig,
                                  strand == '-' & 
                                    end.rnap < start.ig & end.rnap < end.ig  ~ 
                                    start.ig - end.rnap,
                                  strand == '-' & 
                                    end.rnap < start.ig & start.rnap > end.ig  ~ 
                                    start.ig - end.rnap,
                                  .default=NA)) %>% 
                                select(seq_id, start.ig, end.ig, length.ig, dist.ig)


# define outliers based on length of intergenics
too_long_intergenic <- rnap_igs_join %>% 
                       filter(length.ig > max.len.inter) %>%
                       pull(seq_id) %>% 
                       unique()

# save best TDRs and intergenics
tdr_save <- rnap_tdr_join %>% 
  filter(r.tdr == 1) %>%
  select(seq_id, start.tdr, end.tdr, dist.tdr, start2.tdr, end2.tdr)

write.table(tdr_save, output.tdrs, sep='\t',col.names = F, row.names = F, quote = F)


igs_save <-  rnap_igs_join %>% 
  filter(length.ig > min.len.inter) %>%
  filter(dist.ig < max.dist.inter) %>%
  filter(dist.ig > min.dist.inter) %>%
  group_by(seq_id) %>%
  mutate(r.d.ig = rank(dist.ig, na.last='keep')) %>%
  mutate(len_max = max(length.ig, na.rm = TRUE)) %>% 
  mutate(r.d.ig_max = max(r.d.ig, na.rm = TRUE)) %>%
  mutate(r.ig = (length.ig * r.d.ig_max) /  ( len_max * r.d.ig )) %>% 
  slice_max(r.ig) %>% 
  select(seq_id, start.ig, end.ig, length.ig)

write.table(igs_save, output.igs, sep='\t', col.names = F, row.names = F, quote = F)


# filter datasets

## tdrs 
rnap_tdr_join_filter <- rnap_tdr_join %>% 
                               filter(r.tdr == 1) %>% 
                               filter(dist.tdr < max.dist.tdr) %>% 
                               filter(dist.tdr > (min.dist.inter + (end.tdr - start.tdr))) %>% 
                               filter(! seq_id %in% too_long_genome)
## intergenics
rnap_igs_join_filter <- rnap_igs_join %>% 
                           filter(length.ig > min.len.inter) %>%
                           filter(dist.ig < max.dist.inter) %>%
                           filter(dist.ig > min.dist.inter) %>%
                           filter(! seq_id %in% too_long_intergenic) %>% 
                           filter(! seq_id %in% too_long_genome) %>% 
                           group_by(seq_id) %>%
                           mutate(r.d.ig = rank(dist.ig, na.last='keep')) %>%
                           mutate(len_max = max(length.ig, na.rm = TRUE)) %>% 
                           mutate(r.d.ig_max = max(r.d.ig, na.rm = TRUE)) %>%
                           mutate(r.ig = (length.ig * r.d.ig_max) /  ( len_max * r.d.ig )) %>% 
                           slice_max(r.ig) %>% 
                           select(-r.d.ig_max) %>% 
                           select(-len_max) %>% 
                           select(-r.d.ig)


# join datasets

joined <- full_join(rnaps, chr_lens)
joined <- left_join(joined, rnap_tdr_join_filter)
joined <- left_join(joined, rnap_igs_join_filter)


# define datasets

# intergenic length within borders
joined <- joined %>% mutate(intergenic_length_okay = ifelse(!is.na(length.ig), TRUE, FALSE))

#  distance between TDR and RNAP within borders
joined <- joined %>% mutate(tdr_okay = ifelse(!is.na(dist.tdr), TRUE, FALSE))

joined <- joined %>%  mutate(dataset = case_when(
  !intergenic_length_okay ~ 'dataset_3',
  intergenic_length_okay & tdr_okay ~ 'dataset_1',
  intergenic_length_okay | tdr_okay  ~ 'dataset_2',
  .default = NA))        

good <- joined %>% filter(dataset != 'dataset_3') %>% pull(seq_id)
write.table(good, output.txt, col.names = F, row.names = F, quote = F)
