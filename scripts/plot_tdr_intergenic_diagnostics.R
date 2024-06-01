if (!requireNamespace("gggenomes", quietly = TRUE))
  devtools::install_github("thackl/gggenomes")

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)
install.packages('latex2exp')
  
library(tidyverse)
library(gggenomes, verbose = T)
library(dplyr)
library(ggpubr)
library(latex2exp)

args <- commandArgs(trailingOnly=TRUE)

tdrs_tsv <- args[1]
intergenic_bed <- args[2]
rnaps_gff <- args[3]
lengths_bed <- args[4]
tdr_figure_name <- args[5]
intergenic_full_figure_name <- args[6]
intergenic_conditions_figure_name <- args[7]
setwd('PycharmProjects/anti_defence_pipeline/')
# tdrs_tsv <- 'data/tdrs/tdrs.tsv'
intergenic_bed <- 'data_autographiviridae/intergenics/all_intergenic_with_length.tsv'
rnaps_gff <- 'data_autographiviridae/annotation/concatenated_rnaps_only.gff'
lengths_bed <- 'data_autographiviridae/intergenics/chromosome_lengths.bed'

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
                          col.names = c('seq_id', 'start', 'end', 'length'))
## chromosome lengths 
chr_lens <- read.table(lengths_bed, 
                       col.names = c('seq_id', 'start', 'seq_len')) %>% 
            select(-start)

# define distances

## tdrs and rnaps
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
                                  ))

rnap_tdr_join <- rnap_tdr_join %>%
                                  dplyr::group_by(seq_id) %>% 
                                  mutate(tdr.r = rank(dist.tdr, na.last='keep')) %>% 
                                  ungroup()

rnap_tdr_join %>% ggplot(aes(y = dist.tdr, x=end.tdr-start.tdr, color=as.factor(tdr.r))) +
    geom_point(size=3, stroke=1, shape=21) + 
    theme_minimal(base_size=12) + 
    xlab('TDR size, bp') + 
    ylab('RNAP-TDR distance, bp') +
    labs(color='Rank of distances')
    
ggsave(tdr_figure_name, dpi=300)

# rnap_tdr_join <- rnap_tdr_join %>% filter(tdr.r == 1) %>% 
#     mutate(dataset = case_when(
#                        dist.tdr < 10000  & (start.tdr == 1 |  end2 == seq_len)  ~ '1',
#                        dist.tdr < 10000  ~ '2',
#                       .default = NA)) 
# dataset_1_2 <- rnap_tdr_join %>% 
#          filter(!is.na(dataset)) %>% select(seq_id, start.tdr, end.tdr, start2,  end2,
#                                             dist.tdr, dataset)
# 
# dataset_1_2_seqids <- dataset_1_2 %>% pull(seq_id) 


## intergenics and rnaps
rnap_intergenics_join <- rnaps %>% left_join(intergenics, by='seq_id',
                                     suffix = c('.rnap', '.ig')) %>% 
                                   left_join(chr_lens) %>% 
  filter(seq_len >= 35000)



rnap_intergenics_join <- rnap_intergenics_join %>% 
  mutate(dist.intergenic = case_when(
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
                                      .default=NA))


rnap_intergenics_join <- rnap_intergenics_join %>%
  dplyr::group_by(seq_id) %>% 
  mutate(ig.l.r = rank(length, na.last='keep')) %>% 
  mutate(ig.d.r = rank(dist.intergenic, na.last='keep')) %>%
  mutate(len_max_seq = max(length, na.rm = TRUE)) %>% 
  mutate(ig.d.r_max_seq = max(ig.d.r,na.rm = TRUE)) %>% 
  mutate(ig.r = (length * ig.d.r) /  ( len_max_seq  * ig.d.r_max_seq)) %>% 
  mutate(ig.r_inv = (length * ig.d.r_max_seq) /  ( len_max_seq  * ig.d.r)) %>% 
  ungroup()

rnap_intergenics_join <- rnap_intergenics_join %>% group_by(seq_id) %>%
  mutate(seq_max = max(ig.r, na.rm = TRUE)) %>%
  mutate(seq_max_inv = max(ig.r_inv, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(is_max_inv = ifelse(ig.r_inv == seq_max_inv, TRUE, FALSE)) %>% 
 mutate(is_max = ifelse(ig.r == seq_max, TRUE, FALSE))

a <- rnap_intergenics_join %>%
  ggplot(aes(x = dist.intergenic, y=length, color=ig.r)) +
  geom_point(size=2, stroke=2, shape=21, alpha=0.8) + 
  # geom_line(aes(group = seq_id)) +
  theme_minimal(base_size=12) + 
  ylab('Intergenic region size, bp') + 
  xlab('RNAP-intergenic region distance, bp') +
  # scale_y_log10() +
  labs(color=TeX(r'($\frac{L}{||L||}\cdot\frac{D}{||D||}$)')) +
  scale_color_continuous(type='viridis')

b <- rnap_intergenics_join %>%
  ggplot(aes(x = dist.intergenic, y=length, color=ig.r_inv)) +
  geom_point(size=2, stroke=2, shape=21, alpha=0.8) + 
  # geom_line(aes(group = seq_id)) +
  theme_minimal(base_size=12) + 
  ylab('Intergenic region size, bp') + 
  xlab('RNAP-intergenic region distance, bp') +
  # scale_y_log10() +
  labs(color=TeX(r'($\frac{L}{||L||}\cdot\frac{||D||}{D}$)')) +
  scale_color_continuous(type='viridis')


c <- rnap_intergenics_join %>%
  ggplot(aes(x = dist.intergenic, y=length, color=is_max)) +
  geom_point(size=3, stroke=2, shape=21, alpha=0.8) + 
  # geom_line(aes(group = seq_id)) +
  theme_minimal(base_size=12) + 
  ylab('Intergenic region size, bp') + 
  xlab('RNAP-intergenic region distance, bp') +
  # scale_y_log10() +
  labs(color=TeX(r"(Is $\max\left[\frac{L}{||L||}\cdot\frac{D}{||D||}\right]$)")) +
  scale_color_manual(values=c('grey', 'tomato'))


d <- rnap_intergenics_join %>%
  ggplot(aes(x = dist.intergenic, y=length, color=is_max_inv)) +
  geom_point(size=3, stroke=2, shape=21, alpha=0.8) + 
  # geom_line(aes(group = seq_id)) +
  theme_minimal(base_size=12) + 
  ylab('Intergenic region size, bp') + 
  xlab('RNAP-intergenic region distance, bp') +
  # scale_y_log10() +
  labs(color=TeX(r"(Is $\max\left[\frac{L}{||L||}\cdot\frac{||D||}{D}\right]$)")) +
  scale_color_manual(values=c('grey', 'tomato'))


ggarrange(a, b, c, d, nrow = 2, ncol = 2)

ggsave(intergenic_full_figure_name, dpi=300, width=10, height=10)


too_long_genome <- rnap_intergenics_join %>% filter(seq_len > 65000) %>% pull(seq_id) %>% unique()
too_long_intergenic <- rnap_intergenics_join %>% filter(length > 5000) %>% pull(seq_id) %>% unique()


rnap_intergenics_join <- rnap_intergenics_join %>%
  filter(!seq_id %in% too_long_genome) %>% 
  filter(!seq_id %in% too_long_intergenic) %>% 
  dplyr::group_by(seq_id) %>% 
  mutate(ig.l.r = rank(length, na.last='keep')) %>% 
  mutate(ig.d.r = rank(dist.intergenic, na.last='keep')) %>%
  mutate(len_max_seq = max(length, na.rm = TRUE)) %>% 
  mutate(ig.d.r_max_seq = max(ig.d.r,na.rm = TRUE)) %>% 
  mutate(ig.r = (length * ig.d.r) /  ( len_max_seq  * ig.d.r_max_seq)) %>% 
  mutate(ig.r_inv = (length * ig.d.r_max_seq) /  ( len_max_seq  * ig.d.r)) %>% 
  ungroup()

rnap_intergenics_join <- rnap_intergenics_join %>% 
  group_by(seq_id) %>%
  mutate(seq_max = max(ig.r, na.rm = TRUE)) %>%
  mutate(seq_max_inv = max(ig.r_inv, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(is_max_inv = ifelse(ig.r_inv == seq_max_inv, TRUE, FALSE)) %>% 
  mutate(is_max = ifelse(ig.r == seq_max, TRUE, FALSE))


a <- rnap_intergenics_join %>%
  ggplot(aes(x = dist.intergenic, y=length, color=ig.r)) +
  geom_point(size=2, stroke=2, shape=21, alpha=0.8) + 
  # geom_line(aes(group = seq_id)) +
  theme_minimal(base_size=12) + 
  ylab('Intergenic region size, bp') + 
  xlab('RNAP-intergenic region distance, bp') +
  # scale_y_log10() +
  labs(color=TeX(r'($\frac{L}{||L||}\cdot\frac{D}{||D||}$)')) +
  scale_color_continuous(type='viridis')

b <- rnap_intergenics_join %>%
  ggplot(aes(x = dist.intergenic, y=length, color=ig.r_inv)) +
  geom_point(size=2, stroke=2, shape=21, alpha=0.8) + 
  # geom_line(aes(group = seq_id)) +
  theme_minimal(base_size=12) + 
  ylab('Intergenic region size, bp') + 
  xlab('RNAP-intergenic region distance, bp') +
  # scale_y_log10() +
  labs(color=TeX(r'($\frac{L}{||L||}\cdot\frac{||D||}{D}$)')) +
  scale_color_continuous(type='viridis')

c <- rnap_intergenics_join %>%
  filter(!is.na(is_max)) %>% 
  ggplot(aes(x = dist.intergenic, y=length, fill=is_max)) +
  geom_point(size=3, stroke=1, color='black', shape=21, alpha=0.8) + 
  geom_vline(xintercept = 5000, color='black') +
  theme_classic(base_size=12) + 
  ylab('Intergenic region size, bp') + 
  xlab('RNAP-intergenic region distance, bp') +
  # scale_y_log10() +
  labs(fill=TeX(r"(Is $\max\left[\frac{L}{||L||}\cdot\frac{D}{||D||}\right]$)")) +
  scale_fill_manual(values=c('grey', 'tomato'))

c

d <- rnap_intergenics_join %>%
  filter(!is.na(is_max_inv)) %>% 
  ggplot(aes(x = dist.intergenic, y=length, fill=is_max_inv)) +
  geom_point(size=3, stroke=1, color='black', shape=21, alpha=0.8) + 
  geom_vline(xintercept = 5000, color='black') +
  # geom_line(aes(group = seq_id)) +
  theme_classic(base_size=12) + 
  ylab('Intergenic region size, bp') + 
  xlab('RNAP-intergenic region distance, bp') +
  # scale_y_log10() +
  labs(fill=TeX(r"(Is $\max\left[\frac{L}{||L||}\cdot\frac{||D||}{D}\right]$)")) +
  scale_fill_manual(values=c('#d9d9d9ff', '#8dd3c7ff'))
d
ggsave('data_autographiviridae/pics/intergenics.svg', width = 14, height = 6)

e <- rnap_intergenics_join %>%
  filter(is_max) %>%
  ggplot(aes(x = dist.intergenic, y=length, color=ig.r)) +
  geom_point(size=3, stroke=2, shape=21, alpha=0.8) + 
  # geom_line(aes(group = seq_id)) +
  theme_minimal(base_size=12) + 
  ylab('Intergenic region size, bp') + 
  xlab('RNAP-intergenic region distance, bp') +
  # scale_y_log10() +
  labs(color=TeX(r"($\frac{L}{||L||}\cdot\frac{D}{||D||}$)")) +
  scale_color_viridis_c()


f <- rnap_intergenics_join %>%
  filter(is_max_inv) %>%
  ggplot(aes(x = dist.intergenic, y=length, color=ig.r_inv)) +
  geom_point(size=3, stroke=2, shape=21, alpha=0.8) + 
  # geom_line(aes(group = seq_id)) +
  theme_minimal(base_size=12) + 
  ylab('Intergenic region size, bp') + 
  xlab('RNAP-intergenic region distance, bp') +
  # scale_y_log10() +
  labs(color=TeX(r"(\frac{L}{||L||}\cdot\frac{||D||}{D}$)")) +
  scale_color_viridis_c()


ggarrange(a, b, c, d, e, f, nrow = 3, ncol = 2)

ggsave(intergenic_conditions_figure_name, dpi=300, width=10, height=20)

# dataset3 <- rnap_intergenics_join %>% 
#             filter(dist.intergenic < 5000) %>% 
#             filter(dist.intergenic > 1000) %>% 
#             filter(length > 500) %>% 
#             group_by(seq_id) %>% 
#             mutate(len_max_seq = max(length, na.rm = TRUE)) %>% 
#             ungroup() %>% 
#             filter(length == len_max_seq)
# 
# dataset3 %>% pull(seq_id) %>% unique() %>%  length() ==
# dataset3 %>% pull(seq_id) %>% length()
# 
# dataset3 <- dataset3 %>%  
#           mutate(dataset=ifelse(!seq_id %in% dataset_1_2_seqids, '3', NA)) %>% select(seq_id, start.ig, end.ig, length, dist.intergenic,
#                                          dataset)
# 
# 
# datasets <- chr_lens %>% left_join(rnaps) %>% left_join(dataset_1_2, by='seq_id') %>% 
#   left_join(dataset3, by=c('seq_id')) %>% 
#   mutate(dataset = case_when(
#                           !is.na(dataset.x) ~ dataset.x,
#                           !is.na(dataset.y) ~ dataset.y,
#                           .default = NA
#                           ))
# 
# datasets %>% filter(!is.na(dataset)) %>% nrow()
# 




