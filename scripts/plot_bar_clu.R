library(tidyverse)
library(ggpubr)
wd <- 'work_dir/anti_defence/anti_defence_pipeline/'
setwd(wd)
in_path <- 'data_autographiviridae_refseq/clans_info/res_table.tsv'

res_df <- read.table(in_path, sep='\t',
                     header=TRUE)

clans <- read.table('data_autographiviridae_refseq/clans_info/clans_wide_for_sankey.tsv',
                    sep='\t',
                    col.names = c('source', 'target'),
                    na.strings = "") 
res_df_narrow <- res_df %>% select(clu_netw,
                            num_reprs)
clans <- clans %>% left_join(res_df_narrow, by=c('source'='clu_netw'))

clans %>% group_by(target) %>% 
  summarize(s = sum(num_reprs)) %>%
  ggplot(aes(x=s)) +
  geom_histogram(
    # aes(y=after_stat(density)),
    color='black', fill='white', 
    # breaks = c(0, 5, 100, 500),
    closed='left') + 
  xlab('Clan size') +
  scale_x_log10() +
  ylab('# clans') +
  theme_classic2() +
  theme(text = element_text(size=16,
                            family="sans"),
        legend.position = 'bottom',
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent",
                                         colour = NA_character_),
        legend.box.background = element_rect(fill = "transparent",
                                             colour = NA_character_),
        legend.key = element_rect(fill = "transparent",
                                  colour = NA_character_))
ggsave('pics/hist_clan_size_log.png', height = 3, width = 4, dpi=300)


bar_df <- clans %>% group_by(target) %>% 
  summarize(s = sum(num_reprs)) %>% 
  left_join(res_df, by=c('target'='clan'))

in_path <- 'data_autographiviridae_refseq/clans_info/res_table_long.tsv'
res_df <- read.table(in_path, sep='\t',
                     header=TRUE, na.strings = '')
res_df <- res_df %>%  mutate(annot_no_hypo = case_when(
                                                  !is.na(prelim_info) ~ prelim_info,
                                                  annot == 'hypothetical protein' & is.na(pfam_dom) ~ NA,
                                                  !is.na(pfam_dom) ~ paste0('pfam: ', pfam_dom),
                                                  annot == 'unknown function' ~ NA,
                                                  .default = annot))

res_df_ <- res_df %>% 
          group_by(clan) %>%
          summarize(func = names(which.max(table(annot_no_hypo, 
                                                 exclude = if (length(na.omit(annot_no_hypo))>0) NA))),
                    Counts = n()) %>% 
          ungroup()

res_df_ <- res_df_ %>% 
           arrange(-Counts) %>% 
           filter(Counts > 30)
res_df_ <- res_df_ %>% mutate(func = ifelse(str_detect(func, 'kinase'), 'kinase', func))

res_df_ <- cbind(res_df_, ord = nrow(res_df_):1)
res_df_ <- res_df_ %>% mutate(has_func = !is.na(func))   
res_df_ <- res_df_ %>% mutate(func = if_else(is.na(func), 'hypo', func))

ggbarplot(res_df_, x = "ord", y = "Counts",
           fill = 'has_func',
           sorting = "descending",                        # Sort value in descending order
          palette='Dark2',
          font.label = list(size = 15,
                            vjust = -1),
          width = 0.6,     
          x.text.fill = TRUE,
          rotate = TRUE
) +
  scale_x_continuous(labels = res_df_$func,
                   breaks=res_df_$ord) +
  xlab('Clan') +
  ylab('Size') +
  theme(legend.position = 'bottom',
        axis.text.y = element_text(angle=0))

ggsave('pics/barplot_jan.pdf', width=8, height = 13)

