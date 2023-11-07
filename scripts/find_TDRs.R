if (!requireNamespace("gggenomes", quietly = TRUE))
    devtools::install_github("thackl/gggenomes")

library(tidyverse)
library(gggenomes, verbose = T)
library(dplyr)
library(ggpubr)

args <- commandArgs(trailingOnly=TRUE)


in_df_path <- args[1]
tdrs_tsv <- args[2]
figure_name <- args[3]
mp_l_threshold <- as.numeric(args[4])
de_threshold <- as.numeric(args[5])

print('The following arguments passed:')
print(in_df_path)
print(tdrs_tsv)
print(figure_name)
print(mp_l_threshold)
print(de_threshold)

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
  geom_histogram(aes(y=after_stat(density), x = map_length), fill='white', color='black') +
  ylab('Density') +
  xlab('TDR length') + 
  theme_minimal()

tdrs_paf <- tdrs_paf %>% dplyr::select(seq_id, start, end, strand, start2, end2, map_match, map_length, de) %>% 
  dplyr::filter(map_length > 70) %>% 
  dplyr::filter(map_length < 300) %>% 
  dplyr::filter(start2 - end > 0) %>%
  dplyr::filter( map_length - map_match == 0)



b <-  tdrs_paf %>%  ggplot(aes(map_length)) +
  geom_histogram(aes(y=after_stat(density)), color='black', fill='white') +
  ylab('Density') +
    xlab('TDR length') + 
  theme_minimal()

c <- ggarrange(a, b, ncol=1, nrow=2, labels=c('Full scale', 'Zoom on'))

write.table(tdrs_paf, tdrs_tsv, sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
ggsave(filename=figure_name, plot=c, dpi=300, width=12, height=8)
