library(tidyverse)

# setwd('work_dir/anti_defence/anti_defence_pipeline/')
# tdrs_tsv <- 'data/tdrs/tdrs.tsv'
# intergenic_bed <- 'data/intergenics/all_intergenic_with_length.tsv'
# rnaps_gff <- 'data/annotation/concatenated_rnaps_only.gff'
# lengths_bed <- 'data/intergenics/chromosome_lengths.bed'
# 
# max.dist.tdr <- 4000
# max.dist.inter <- 5000
# min.dist.inter <- 1000
# max.gen <- 65000 # k
# min.len.inter <- 500  # k
# max.len.inter <- 3000# k


args <- commandArgs(trailingOnly=TRUE)


# intergenic_bed <- 'data_autographiviridae_meta/intergenics/all_intergenic_with_length.tsv'
# rnaps_gff <- 'data_autographiviridae_meta/annotation/concatenated_rnaps_only.gff'
# lengths_bed <- 'data_autographiviridae_meta/intergenics/chromosome_lengths.bed'
# max.dist.inter <- 5000
# min.dist.inter <- 1000
# max.gen <- 60000
# min.len.circ <- 35000
# min.len.inter <- 500
# max.len.inter <- 3000
# test_igs.tsv




intergenic_bed <- args[1]
rnaps_gff <- args[2]
lengths_bed <- args[3]
max.dist.inter <- as.numeric(args[4])
min.dist.inter <- as.numeric(args[5])
max.gen <- as.numeric(args[6])  # k
min.len.circ <- as.numeric(args[7])  # k
min.len.inter <- as.numeric(args[8])  # k
max.len.inter <- as.numeric(args[9])  # k
output.igs <- args[10]


# read data
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
## intergenics 
rnap_igs_join <- rnaps %>% left_join(intergenics, by='seq_id',
                                     suffix = c('.rnap', '.ig')) %>% 
                  left_join(chr_lens)

# add mask to remove RNAP downstream if length is less than threshold
rnap_igs_join <- rnap_igs_join %>% mutate(short_contig = (seq_len <= min.len.circ)) %>% 
    mutate(filter_target = case_when(
                    !short_contig ~ FALSE,
                    strand == '+' & end.ig >= end.rnap  ~ TRUE,
                    strand == '-' & end.ig <= start.rnap  ~ TRUE,
                            .default = FALSE)
                            )
rnap_igs_join <- rnap_igs_join %>% filter(!filter_target)

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

# save best intergenics
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
