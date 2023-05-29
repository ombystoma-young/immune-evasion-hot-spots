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

setwd('work_dir/anti_defence/anti_defence_pipeline/')
s0 <- read_seqs("blasted/phages_genomes_concat.fna") %>% select(seq_id,
                                                                seq_desc)
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

df2 <- ape::read.gff('results/upstreams_with_clusters.gff') %>%
      filter(str_detect(attributes, 'cluster_num=2;')) %>% 
      mutate(have_system = 'Ocr') %>% 
      select(seqid, have_system) %>% unique()
df4 <- ape::read.gff('results/upstreams_with_clusters.gff') %>%
  filter(str_detect(attributes, 'cluster_num=4;|cluster_num=33;')) %>% 
  mutate(have_system = 'SAMase') %>% 
  select(seqid, have_system) %>% unique()

df_none <- ape::read.gff('results/upstreams_with_clusters.gff') %>%
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



drops <- c('KJ183192.1', 'JQ780163.1', 'OP413828.1')
rnaps_tree <- read.tree("define_datasets/trees/polymerases_all.iqtree.treefile")
rnaps_tree_with_drops <- drop.tip(rnaps_tree, drops)
rnaps_tree_with_drops <- groupOTU(rnaps_tree_with_drops, c('NC_001604.1', 'NC_003298.1'))
t <- ggtree(rnaps_tree_with_drops, 
            layout = 'circular',
            branch.length = 'none') +
  geom_hilight(node=478, fill="purple", alpha=0.3) +
  geom_tippoint(aes(alpha = group), col = "red") +
  scale_color_manual(values = c(0,1), aesthetics = "alpha") +
  theme(legend.position = 'none')


# write.table(df2, 'metadata/host_ids.tsv', col.names = F,
#             quote = F, sep='\t')

colors_ <- c('Both'='#440154', 'SAMase'='#fde725', 'Ocr'='#35b779', 'No'='#ffffff')
p1 <- gheatmap(t, df1, offset=0.8, width=.1,
               colnames_position = 'top',
               custom_column_labels = c(''),  font.size = 1.5,
               colnames_angle=0, colnames_offset_y = 0.25) +
  scale_fill_manual(values=colors_, name="Known",
                    breaks = c('Both', 'Ocr', 
                               'SAMase',  'No')) 
  #geom_taxalink("NC_001604.1", "NC_003298.1", color="blue3") 
colors_bac <- c('Enterobacterales'='#440154', 'Pseudomonadales'='#1ae4b6', 
             'Vibrionales'='#faba39', 'Other'='gray90')
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, df2, offset=5.5, width=.1,
         custom_column_labels = c('Host'), font.size = 2.5,
         colnames_angle=0, colnames_offset_y = .25) +
  scale_fill_manual(values=colors_bac, name="Host Order",
                    breaks = c('Enterobacterales', 'Pseudomonadales', 
                               'Vibrionales',  'Other')) 

colors_phage <- c('Molineuxvirinae'='#f1605d', 'Studiervirinae'='#721f81', 
                'Colwellvirinae'='#fcfdbf', 'Other'='gray90')
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
ggsave('pics/tree_dataset_order_virus_subfamily.pdf', width=9, height = 9)


