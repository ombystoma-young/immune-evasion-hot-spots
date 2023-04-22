library(tidyverse)
library(gggenomes, verbose = T)
library(dplyr)
args <- commandArgs(TRUE)


in_df_path <- 'minimap2_out/all_phages.paf'
tdrs_tsv <- 'minimap2_out/TDRs_all.tsv'
figure_name <- 'pics/hist_length_TDRs.png'
mp_l_threshold <- 99
de_threshold <- 0.3

# print(mp_l_threshold)
# print(de_threshold)

tdrs_paf <- read_paf(in_df_path)
tdrs_paf <- tdrs_paf %>%
  dplyr::filter(seq_id == seq_id2)  %>%  
  dplyr::filter(start < start2) %>% 
  dplyr::filter(map_length > mp_l_threshold) %>% 
  dplyr::filter(de < de_threshold)

a <- tdrs_paf %>% mutate(state_100 = (100<=map_length)) %>% 
  mutate(state_300 = (300>=map_length)) %>% 
  mutate(state = state_100*state_300) %>% 
  ggplot() +
  geom_histogram(aes(y=..density.., x = map_length), fill='white', color='black') +
  ylab('Density') +
  xlab('TDR length') + 
  theme_minimal()

tdrs_paf <- tdrs_paf %>% dplyr::select(seq_id, start, end, strand, start2, end2, map_match, map_length, de) %>% 
  dplyr::filter(map_length > 70) %>% 
  dplyr::filter(map_length < 300) %>% 
  dplyr::filter(start2 - end > 0) %>%
  dplyr::filter( map_length - map_match == 0)



b <-  tdrs_paf %>%  ggplot(aes(map_length)) +
  geom_histogram(aes(y=..density..), color='black', fill='white') +
  ylab('Density') +
    xlab('TDR length') + 
  theme_minimal()

library(ggpubr)
c <- ggarrange(a, b, ncol=1, nrow=2)

write.table(tdrs_paf, tdrs_tsv, sep='\t', row.names = FALSE, quote = FALSE)
ggsave(filename=figure_name, plot=c, dpi=300, width=12, height=8)
