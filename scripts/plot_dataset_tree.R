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

col_pal <- brewer.pal(n = 8, name = "Dark2")
setwd('work_dir/anti_defence/anti_defence_pipeline/')
s0 <- read_seqs("data_autographiviridae_refseq/genomes/concatenated_genomes.fna") %>% 
  select(seq_id, seq_desc)
s0 <- s0 %>% mutate(seq_desc = str_split(seq_desc, ",", simplify = T)[, 1]) %>% 
  mutate(seq_desc = str_split(seq_desc, " ", simplify = T)[, 1:6])

s0 <- s0 %>% mutate(first = seq_desc[,1],
              second = seq_desc[,2],
              third = seq_desc[,3]) %>% 
      mutate(host = case_when(
                  second == 'Reminis' ~ 'Ralstonia',
                  second == 'Titan-X' ~ 'Xanthomonas',
                  second == 'Phi' ~ NA,
                  first == 'MAG:' & 
                    (second == 'Prokaryotic' | second == 'Bacteriophage')  ~ NA,
                  first == 'MAG:' | first == 'Mutant' | 
                  first == 'UNVERIFIED:' | first == 'uncultured' ~ second,
                  first == 'Bacteriophage' ~ second,
                  first == 'Vibriophage' ~ 'Vibrio',
                  first == 'Caudovirales' | first == 'Caudoviricetes' ~ NA,
                  first == 'Uncultured' | first == 'Phage' ~ NA,
                  first == 'Cyanophage' ~ NA,
                 .default = first  
                  )) %>% select(seq_id, host)

col.nms <- c('taxid', 'emp', 'superkingdom',
             'A',  'B',  'C', 'D', 'E', 'F', 'G')
s1 <- read.table('metadata/host_lineage.tsv', sep=';', header=F, col.names = col.nms,
                 fill = TRUE) %>% filter(superkingdom == 'Bacteria')

s1 <- s1 %>% mutate(genus = case_when(
              G != '' ~ G,
              G == '' & `F` != '' ~ F,
              .default = E
                )) %>%
              mutate(order = case_when(
                G != '' ~ E,
                G == '' & `F` != '' ~ D,
                .default = C
              )) %>% 
            mutate(order = case_when(
              genus == 'Escherichia' ~ 'Escherichia',
              genus == 'Synechococcus' ~ 'Synechococcus',
              genus == 'Rhizobium' ~ 'Hyphomicrobiales',
              .default = order)) %>% 
            select(genus, order)


s3 <- left_join(s0, s1, by=c('host'='genus')) %>% 
  mutate(order = case_when(
    host %in% c('Enterobacteria', 'Lelliottia') ~ 'Enterobacterales',
    host %in% c('Burkholderia', 'Sphaerotilus',
                'Achromobacter', 'Curvibacter') ~ 'Burkholderiales',
    host %in% c('Acinetobacter') ~ 'Moraxellales',
    host %in% c('Pelagibacter') ~ 'Pelagibacterales',
    host %in% c('Roseobacter') ~ 'Rhodobacterales',
    host %in% c('Alteromonas', 'Shewanella',
                'Colwellia') ~ 'Alteromonadales',
    host %in% c('Stappia', 
                'Mesorhizobium',
                'Aquamicrobium',
                'Agrobacterium') ~ 'Hyphomicrobiales',
    host %in% c('Desulfovibrio') ~ 'Desulfovibrionales',
    host %in% c('Brevundimonas') ~ 'Caulobacterales',
    host %in% c('Tenacibaculum', 
                'Cellulophaga') ~ 'Flavobacteriales',
    host %in% c('Xylella') ~ 'Xanthomonadales',
    host %in% c('Marinomonas') ~ 'Oceanospirillales',
    host %in% c('Sphingomonas') ~ 'Sphingomonadales',
    .default = order
  ))

s3 <- s3 %>% select(seq_id, order)

s2 <- read.table('metadata/genome_lineage.tsv', sep=';', header=F,
                       fill = TRUE) %>% 
    mutate(subfamily = case_when(
                    str_detect(V8, 'sp.') ~ NA,
                    str_detect(V8, 'unclassified') ~ NA,
                    str_detect(V8, 'uncultured') ~ NA,
                          .default = V8
                        )) %>% select(V1, subfamily)
s4 <- read.table('metadata/genome_id2taxid.tsv') %>% 
      select(V2, V3)
s2 <- s2 %>% left_join(s4, by=c('V1'='V3')) %>% select(subfamily, V2) 
colnames(s2) <- c('subfamily', 'seqid')

# ocr 3, 107, 135, 297 (1-based)
ocr <- c(3814, 3136, 3132, 3030, 3020, 1911)
ocr <- paste0('clu=', ocr, ';', collapse = '|')
df2 <- ape::read.gff('data_autographiviridae_refseq/upstreams/early_with_clusters.gff') %>%
      filter(str_detect(attributes, ocr)) %>%
      mutate(have_system = 'Ocr') %>% 
      select(seqid, have_system) %>% unique()
# ard 256, 382, 488, 583 (1-based)
# df3 <- ape::read.gff('results/upstreams_with_clusters.gff') %>%
#   filter(str_detect(attributes, 'cluster_num=255;|cluster_num=381;|cluster_num=487;|cluster_num=582;')) %>% 
#   mutate(have_system = 'ArdA/ArdcN') %>% 
#   select(seqid, have_system) %>% unique()
# SAMase (1-based):
# samase_1_based <- c(5, 33, 55, 104, 196, 238, 289, 308, 372, 503, 626)

samase <- c(392, 516, 2538, 3167, 3357, 3486, 3509, 3704, 1677)
samase <- paste0('clu=', samase, ';', collapse = '|')
df4 <- ape::read.gff('data_autographiviridae_refseq/upstreams/early_with_clusters.gff') %>%
  filter(str_detect(attributes, samase)) %>% 
  mutate(have_system = 'SAMase') %>% 
  select(seqid, have_system) %>% unique()

df_none <- ape::read.gff('data_autographiviridae_refseq/upstreams/early_with_clusters.gff') %>%
  filter(! seqid %in% df2$seqid) %>% 
  filter(! seqid %in% df4$seqid) %>% 
  mutate(have_system = 'No') %>% 
  select(seqid, have_system) %>% unique()

df <- rbind(df2, df4, df_none) %>% 
  mutate(have_system = factor(have_system))
both <- df %>% select(seqid) %>% filter(duplicated(seqid))
df <- df %>% mutate(have_system = case_when(
                                  seqid %in% both$seqid ~ 'Both', 
                                  .default = have_system)) %>% 
        unique() %>% mutate(have_system = factor(have_system, 
                                                 levels = c('No', 'Ocr', 
                                                            'SAMase', 'Both'), 
                                                 ordered = TRUE))


df <- df %>% left_join(s3, by=c('seqid'='seq_id')) %>% left_join(s2) %>% unique()

rownames(df) <- df$seqid 

df <- df %>% select(-seqid) %>% 
  mutate(order = factor(order))
df1 <- df %>% select(have_system)

df2 <- df %>% select(order)
in_other <- table(df2$order) / nrow(df2) * 100 < 2.5 
in_other <- names(in_other[in_other == TRUE])
  df2 <- df2 %>% mutate(order = case_when(order %in% in_other ~ 'Other',
                                 order == 'Enterobacteriaceae' ~ 'Enterobacterales',
                                 is.na(order) ~ 'Other',
                                 .default = order))
table(df2$order) / nrow(df2) * 100 < 2.5

df3 <- df %>% select(subfamily)
in_other <- table(df3$subfamily) / nrow(df3) * 100 < 2.5 
in_other <- names(in_other[in_other == TRUE])
df3 <- df3 %>% mutate(subfamily = case_when(subfamily %in% in_other ~'Other', 
                                            is.na(subfamily) ~ 'Other',
                                            .default = subfamily))
table(df3$subfamily) / nrow(df3) * 100 < 2.5



# drops <- c('KJ183192.1', 'JQ780163.1', 'OP413828.1',
#            'MZ318361.1', 'MZ318362.1')
rnaps_tree <- read.tree("data_autographiviridae_refseq/known_proteins/trees/rnap_bootstrap_model_selection.iqtree.contree")
# rnaps_tree_with_drops <- drop.tip(rnaps_tree, drops)
a <- as.integer(rnaps_tree$node.label)
a[is.na(a)] <- 0
a <- ifelse(a >= 95, intToUtf8(9728), "")
rnaps_tree$node.label <- a

# rnaps_tree_with_drops <- groupOTU(rnaps_tree_with_drops, )
t <- ggtree(rnaps_tree, 
            layout = 'circular',
            branch.length = 'none') +
  #geom_hilight(node=478, fill="purple", alpha=0.3) +
  geom_point2(aes(subset=(label %in% c('NC_001604.1', 'NC_003298.1'))), shape=23, size=2, fill='red') +
  geom_nodelab(color='firebrick') +
  # geom_tippoint(aes(alpha = group), col = "red") +
  # scale_color_manual(values = c(0,1), aesthetics = "alpha") +
  theme(legend.position = 'none')

t
# write.table(df2, 'metadata/host_ids.tsv', col.names = F,
#             quote = F, sep='\t')
# colors_ <- c('Both'='#440154', 'SAMase'='#fde725',   'Ocr'='#35b779', 'No'='#ffffff')
# 
# colors_bac <- c('Escherichia'='#31aff5', 'Vibrionales'= '#faba39', 
#                 'Pseudomonadales'='#83ff52', 
#                 'Enterobacterales'='#440154', 'Other'='gray90')
# 
# colors_phage <- c('Molineuxvirinae'='#f1605d', 'Studiervirinae'='#721f81', 
#                   'Colwellvirinae'='#fcfdbf', 'Other'='gray90')

colors_ <- c('Both'=col_pal[6], 'SAMase'=col_pal[5], 'Ocr'=col_pal[7], 'No'='#ffffff')

colors_bac <- c('Escherichia'=col_pal[4], 'Vibrionales'= col_pal[2], 
                'Pseudomonadales'=col_pal[1], 
                'Enterobacterales'=col_pal[3], 'Other'=col_pal[8])

colors_phage <- c('Molineuxvirinae'=col_pal[6], 'Studiervirinae'=col_pal[5], 
                  'Colwellvirinae'=col_pal[7], 'Other'=col_pal[8])

  table(df1$have_system)
p1 <- gheatmap(t, df1, offset=0, width=.1,
               colnames_position = 'top',
               custom_column_labels = c(''),  font.size = 1.5,
               colnames_angle=0, colnames_offset_y = 0.25) +
  scale_fill_manual(values=colors_, name="Known",
                    breaks = c('Both', 
                               'SAMase', 'Ocr', 'No')) 

p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, df2, offset=5.5, width=.1,
         custom_column_labels = c('Host'), font.size = 2.5,
         colnames_angle=0, colnames_offset_y = .25) +
  scale_fill_manual(values=colors_bac, name="Host",
                    breaks = c('Enterobacterales', 'Escherichia',
                               'Pseudomonadales', 
                               'Vibrionales',  'Other')) 

p3 <- p2 + new_scale_fill()
p3 <- gheatmap(p3, df3, offset=11, width=.1, 
         custom_column_labels = c('Phage'), font.size = 2.5,
         colnames_angle=0, colnames_offset_y = .25) +
  scale_fill_manual(values=colors_phage, name="Virus Subfamily",
                    breaks = c('Studiervirinae', 'Molineuxvirinae', 
                               'Colwellvirinae',  'Other')) 

  
# p4 <- gridExtra::tableGrob(data.frame(table(df2$order)))
# 
# p5 <- gridExtra::tableGrob(data.frame(table(df3$subfamily)))

p3
ggsave('pics/rnap_tree_jan.svg', width=15, height = 9)



