# plot tree of genome dataset, RNAP based
library(tidyverse)
library(patchwork)  # arrange multiple plots
library(ggtree)
library(treeio)
library(reshape2)
library(gggenomes)
library(ggnewscale)
library(stringr)
library(gridExtra)
library(stringi)
library(RColorBrewer)

col_pal <- brewer.pal(n = 12, name = "Set3")
setwd('work_dir/anti_defence/anti_defence_pipeline/')
s0 <- read_seqs("data_autographiviridae/genomes/concatenated_genomes.fna") %>% 
  select(seq_id, seq_desc)

hosts_tax <- read.table('metadata/lineage/autographiviridae_host_parsed_joined_lineage.tsv',
                        sep='\t', header=TRUE, na.strings = '')
hosts_tax <- hosts_tax %>% 
        mutate(order = ifelse(host_name == 'Enterobacteria', 'Enterobacterales', order)) %>% 
        mutate(family = ifelse(host_name == 'Enterobacteria', 'Enterobacteriaceae', family)) %>% 
        mutate(no.rank = ifelse(host_name == 'Enterobacteria', 'cellular organisms', no.rank))

hosts_tax <- hosts_tax %>% 
  mutate(order = ifelse(host_name == 'Pelagibacter', 'Pelagibacterales', order)) %>% 
  mutate(family = ifelse(host_name == 'Pelagibacter', 'Pelagibacteraceae', family)) %>% 
  mutate(no.rank = ifelse(host_name == 'Pelagibacter', 'cellular organisms', no.rank))

hosts_tax <- hosts_tax %>%
  mutate(no.rank = ifelse(host_name == 'Cyanophage', 'cellular organisms', no.rank))

hosts_tax <- hosts_tax %>% filter(!is.na(no.rank)) 



s3 <- s0 %>% left_join(hosts_tax, by=c('seq_id'='nuccore_id')) %>% select(seq_id, order)


virus_tax <- read.table('metadata/lineage/autographiviridae_phage_parsed_joined_lineage.tsv',
                        sep='\t', header=TRUE, na.strings = '')  


s2 <- s0 %>% left_join(virus_tax, by=c('seq_id'='nuccore_id')) %>% select(seq_id, subfamily) 


df <- read.table('data_autographiviridae/loci_similarity/has_target_adgs_htgs_bool.tsv',
                  sep='\t', header=TRUE) %>% 
      mutate(have_system = case_when(
        has_ocrs == 'True' & has_samases == 'True' ~ 'Both',
        has_ocrs == 'True' ~ 'Ocr',
        has_samases == 'True' ~ 'SAMase',
        .default = 'No'
        )) %>% 
      mutate(seq_id = str_replace_all(seq_id, '@', '_')) %>% 
      select(seq_id, have_system) %>% 
       mutate(have_system = factor(have_system, 
                                     levels = c('No', 'Ocr', 
                                                'SAMase', 'Both'), 
                                     ordered = TRUE))


df <- df %>% left_join(s3) %>% left_join(s2)# %>% unique()

rownames(df) <- df$seq_id 

df <- df %>% 
  select(-seq_id) %>% 
  mutate(order = factor(order))
df1 <- df %>% select(have_system)

df2 <- df %>% select(order)
in_other <- table(df2$order) / nrow(df2) * 100 < 1.5
in_other <- names(in_other[in_other])
  df2 <- df2 %>% mutate(order = case_when(order %in% in_other ~ 'Other',
                                 order == 'Enterobacteriaceae' ~ 'Enterobacterales',
                                 is.na(order) ~ 'Other',
                                 .default = order))
table(df2$order) / nrow(df2) * 100 < 1.5

df3 <- df %>% select(subfamily)
df3 <- df3 %>% mutate(subfamily = case_when(
                                            subfamily %in% c('Molineuxvirinae', 
                                                             'Colwellvirinae', 
                                                             'Studiervirinae') ~ subfamily,
                                            .default = 'Other'))

rnaps_tree <- read.newick("data_autographiviridae/known_adgs_analysis/trees/rnap_bootstrap_model_selection.iqtree.contree",
                          node.label='support')

rnaps_tree@data$support <- ifelse(rnaps_tree@data$support > 95, TRUE, NA)


color_all_df <- read.table('data_autographiviridae/loci_similarity/parsed_communities_info.tsv',
                 col.names = c('seq_id', 'group'), sep='\t')
rownames(color_all_df) <- color_all_df$seq_id 
t <- ggtree(rnaps_tree, 
            layout = 'circular',
            branch.length = 'none',
            linewidth=0.8,
            aes(color = group))   %<+% color_all_df +
  ggtree::geom_nodepoint(aes(subset=!is.na(support)), fill = 'white', alpha=1, shape=21) +
  geom_point2(aes(subset=(label %in% c('NC_001604.1', 'NC_003298.1'))), shape=23, size=2, fill='red') +
  scale_colour_hue(na.value='#000000') +
  theme(legend.position = 'none')

t


colors_ <- c('Both'=col_pal[1], 'SAMase'=col_pal[2], 'Ocr'=col_pal[3], 'No'='#ffffff')

colors_bac <- c('Burkholderiales'=col_pal[4], 'Vibrionales'= col_pal[5], 
                'Pseudomonadales'=col_pal[6], 
                'Enterobacterales'=col_pal[7], 'Other'=col_pal[9])

colors_phage <- c('Molineuxvirinae'=col_pal[8], 'Studiervirinae'=col_pal[10], 
                  'Colwellvirinae'=col_pal[11], 'Other'=col_pal[9])

  table(df2$order)[2]
p1 <- gheatmap(t, df1, offset=0, width=.1,
               colnames_position = 'top',
               custom_column_labels = c(''),  font.size = 1.5,
               colnames_angle=0, colnames_offset_y = 0.25) +
  scale_fill_manual(values=colors_, name="Known",
                    breaks = c('Both', 
                               'SAMase', 'Ocr', 'No')) 
p1
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, df2, offset=5.5, width=.1,
         custom_column_labels = c('Host'), font.size = 2.5,
         colnames_angle=0, colnames_offset_y = .25) +
  scale_fill_manual(values=colors_bac, name="Host",
                    breaks = c('Burkholderiales', 'Enterobacterales',
                               'Pseudomonadales', 
                               'Vibrionales',  'Other'),
                    na.value =col_pal[9]) 

p3 <- p2 + new_scale_fill()
p3 <- gheatmap(p3, df3, offset=11, width=.1, 
         custom_column_labels = c('Phage'), font.size = 2.5,
         colnames_angle=0, colnames_offset_y = .25) +
  scale_fill_manual(values=colors_phage, name="Virus Subfamily",
                    breaks = c('Studiervirinae', 'Molineuxvirinae', 
                               'Colwellvirinae',  'Other'), na.value =col_pal[9]) 

p3  
ggsave('data_autographiviridae/pics/rnap_tree_may.pdf', width=15, height = 9)


