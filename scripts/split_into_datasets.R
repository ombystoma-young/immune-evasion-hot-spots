library(tidyverse)

# setwd('work_dir/anti_defence/anti_defence_pipeline/')

upstreams_names <- c('chrm', 'upstream_start', 'upstream_end',
                     'dot', 'pos', 'strand',
                     'pol_start', 'pol_end',
                     'tdr_left_start', 'tdr_left_end',
                     'tdr_right_start', 'tdr_right_end')
upstreams_names_tr <- c('chrm', 'strand',
                        'pol_start', 'pol_end',
                        'tdr_left_start', 'tdr_left_end',
                        'tdr_right_start', 'tdr_right_end', 'u_length_')
intergenics_names <- c('chrm', 'i_start', 'i_end', 'i_length_')

found_upstreams <- read.table('upstream_search/upstream.bed',
                              fill=T, col.names = upstreams_names)
chr_length <- read.table(file ='promoters_search/chromosome_lengths.bed', 
                         col.names = c('chrm', 'chrm_start', 'chrm_end'))
intergenics <- read.table('promoters_search/all_intergenic_with_length.tsv', 
                          col.names = intergenics_names)
upstream_true_length <- read.table('upstream_search/tdr_pol_dist.tsv',
                                   fill=T, col.names = upstreams_names_tr) %>%  
                                            mutate(u_length_ = abs(u_length_))


found_upstreams <- found_upstreams %>% filter(is.na(tdr_left_start)) %>%
  select(-dot, -pos, -upstream_start, -upstream_end) %>% unique()

intergenics <- intergenics %>% mutate(i_length_ = abs(i_length_))

intergenics <- intergenics %>% 
  group_by(chrm) %>% slice_max(i_length_, n=1)


joined <- full_join(found_upstreams, upstream_true_length) 
joined <- left_join(joined, intergenics)
joined <- left_join(joined, chr_length) %>% select(-chrm_start)

# distance from RNAP to intergenic:
joined <- joined %>% 
  mutate(distance = case_when(strand == '+' & pol_start > i_end ~ pol_start - i_end,
                              strand == '+' & pol_start <= i_end ~ pol_start +
                                chrm_end - i_start,
                              strand == '-' & pol_start < i_start ~ i_start - pol_end,
                              strand == '-' & pol_start >= i_start ~ chrm_end - 
                                pol_start + i_start))
# plot distances
c <- joined %>%  
  ggplot(aes(x=distance)) +
  geom_histogram(aes(y=..density..), color= 'black', fill = 'white') + 
  theme_minimal() + 
  xlim(c(0,10000)) +
  ylab('Density') +
  xlab('Distance from the longest intergenic region to RNAP, bp') +
  theme(text = element_text(size=15))
c
ggsave('pics/intergenic_no_tdr_diagnostics_less_than10000_.png', width=10, height = 8, dpi=400)





# # define dataset

#intergenic more than 500:
joined <- joined %>% mutate(long_intergenic = ifelse(i_length_ > 500, TRUE, FALSE))

# has TDR and distance from TDR to RNAP <= 10000 
joined <- joined %>% mutate(has_tdr = !is.na(tdr_left_start)) %>% 
                     mutate(rna_tdr_less10 = case_when(u_length_ < 10001 & !is.na(tdr_left_start) ~ TRUE,
                                                       u_length_ > 10000 & !is.na(tdr_left_start) ~ FALSE,
                                                       .default = NA))

# TDR in intergenic:
joined <- joined %>% mutate(i_start  = ifelse(i_start > i_end, -i_start, i_start))
joined$left_in <- joined$tdr_left_end < joined$i_end & joined$tdr_left_start > joined$i_start
joined$right_in <- joined$tdr_right_end < joined$i_end & joined$tdr_right_start > joined$i_start
joined <- joined %>% mutate(TDR_in = left_in | right_in)

# TDR at the ends
joined <- joined %>% mutate(TDR_on_ends = case_when( tdr_left_start  == 1 & 
                                                       chrm_end - tdr_right_end < 5  ~ TRUE,
                                                      .default = FALSE))

# start splitting:

joined <- joined %>%  mutate(dataset = case_when(
  !long_intergenic ~ 'dataset_5',
   long_intergenic & has_tdr & rna_tdr_less10 & TDR_in & TDR_on_ends  ~ 'dataset_1',
   long_intergenic & has_tdr & rna_tdr_less10 & TDR_in & !TDR_on_ends  ~ 'dataset_2',
   long_intergenic & (!has_tdr | !rna_tdr_less10  |  !TDR_in)  & (distance < 5001) ~ 'dataset_3',
  long_intergenic & (!has_tdr | !rna_tdr_less10  |  !TDR_in)  & (distance > 5000) ~ 'dataset_4',
  .default = NA))                                                

sum(is.na(joined$dataset))

table(joined$dataset)

write.table(joined, 'define_datasets/joined.tsv',
            sep = '\t',
            col.names = T, row.names = F, quote=F)
