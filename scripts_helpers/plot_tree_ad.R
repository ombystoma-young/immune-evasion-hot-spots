# plot tree of genome dataset, SAMase based
library(tidyverse)
library(patchwork)  # arrange multiple plots
library(ggtree)
library(treeio)
library(reshape2)
library(gggenomes)
library(ggnewscale)
library(stringr)
library(gridExtra)
library(RColorBrewer)

col_pal <- brewer.pal(n = 8, name = "Set3")

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
    genus == 'Escherichia' ~ 'Enterobacterales',
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
# write.table(s3, 'metadata/seq_id_to_host.tsv', sep='\t', quote=F, row.names = F)
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

# ocr 3, 100, 125, 260 (1-based)
ocr <-c(3814, 3136, 3132, 3030, 3020, 1911)
# new_clp <- c('GCA_020488535.1_00041',
#   'GCA_017903885.1_00036',
#   'GCF_000903135.1_00004',
#   'GCA_014071185.1_00045',
#   'GCF_002997865.1_00003')
# clu_107 <- read_feats('results/upstreams_with_clusters.gff') %>%
#   filter(cluster_num  == 99 | geom_id %in% new_clp)
# clu_135 <- read_feats('results/upstreams_with_clusters.gff') %>%
#   filter(cluster_num  == 124)

color_branches_ocr_df <- 
  gggenomes::read_feats('data_autographiviridae_refseq/upstreams/early_with_clusters.gff') %>% 
  filter(clu %in% ocr) %>% select(seq_id, clu)
color_branches_ocr_df <- color_branches_ocr_df %>% 
  mutate(clu = factor(clu))  

df2 <- read_feats('data_autographiviridae_refseq/upstreams/early_with_clusters.gff') %>%
  filter(clu %in% ocr) %>% 
  mutate(have_system = 'Ocr') %>% 
  select(seq_id, have_system) %>% unique()
colnames(df2) <- c('seqid', 'have_system')

samase <- c(392, 516, 2538, 3167, 3357, 3486, 3509, 3704, 1677)

color_branches_df <- 
  gggenomes::read_feats('data_autographiviridae_refseq/upstreams/early_with_clusters.gff') %>% 
  filter(clu %in% samase) %>% select(seq_id, clu)
color_branches_df <- color_branches_df %>% 
  mutate(clu = factor(clu))

# clu_308 <- read_feats('results/upstreams_with_clusters.gff') %>%
#   filter(clu == 269) %>%
#   select(seq_id) %>% unique()
# clu_289 <- read_feats('results/upstreams_with_clusters.gff') %>%
#   filter(clu == 251) %>%
#   select(seq_id) %>% unique()


df4 <- read_feats('data_autographiviridae_refseq/upstreams/early_with_clusters.gff') %>%
  filter(clu %in% samase) %>% 
  mutate(have_system = 'SAMase') %>% 
  select(seq_id, have_system) %>% unique()

colnames(df4) <- c('seqid', 'have_system')

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

df <- df %>% 
  mutate(order = factor(order))
df1 <- df %>% select(have_system) %>% as.matrix()
rownames(df1) <- df$seqid 

df2 <- df %>% select(order)


in_other <- table(df2$order) / nrow(df2) * 100 < 2.5 
in_other <- names(in_other[in_other == TRUE])
df2 <- df2 %>% mutate(order = case_when(order %in% in_other ~ 'Other',
                                        order == 'Enterobacteriaceae' ~ 'Enterobacterales',
                                        is.na(order) ~ 'Other',
                                        .default = order)) %>% as.matrix()
rownames(df2) <- df$seqid 
table(df2[, 1]) / nrow(df2) * 100 < 2.5

df3 <- df %>% select(subfamily)


in_other <- table(df3$subfamily) / nrow(df3) * 100 < 2.5 
in_other <- names(in_other[in_other == TRUE])
df3 <- df3 %>% mutate(subfamily = case_when(subfamily %in% in_other ~'Other', 
                                            is.na(subfamily) ~ 'Other',
                                            .default = subfamily)) %>% as.matrix()
rownames(df3) <- df$seqid 
table(df3[, 1]) / nrow(df3) * 100 < 2.5



# drops <- c('KJ183192.1', 'JQ780163.1', 'OP413828.1')
drops <- c('Svi3-7', 'Svi3-3', 'ORF1')
samase_tree <- read.tree("data_autographiviridae_refseq/known_proteins/trees/samase_bootstrap_model_selection.iqtree.contree")
# samase_tree <- drop.tip(samase_tree, drops)
samase_tree <- groupOTU(samase_tree, c('NC_003298.1'))
a <- as.integer(samase_tree$node.label)
a[is.na(a)] <- 0
a <- ifelse(a >= 95, intToUtf8(9728), "")
samase_tree$node.label <- a

t <- ggtree(samase_tree,
            branch.length = 'none',
            layout = 'circular',
            linewidth=1.5,
            aes(color = clu, label = clu)
            )  %<+% color_branches_df +
  geom_point2(aes(subset=(label == 'NC_003298.1')), shape=23, size=2.5, fill='red') +
  #geom_point2(aes(subset=(label %in% clu_289$seq_id)), shape=21, size=2.5, fill='#31aff5') +
  #geom_point2(aes(subset=(label %in% clu_308$seq_id)), shape=21, size=2.5, fill='#721f81') +
  # geom_hilight(node=478, fill="purple", alpha=0.3) +
  geom_nodelab(color='red')  + 
 # geom_tiplab(show.legend = FALSE, align = TRUE) +
  scale_color_brewer(na.value = "black", palette = 'Set3')

  # geom_tiplab(offset = 0.5, align=T, size=2) +
  # scale_color_manual(values = c(0,1), aesthetics = "alpha") 
  
t
# write.table(df2, 'metadata/host_ids.tsv', col.names = F,
#             quote = F, sep='\t')

# colors_ <- c('Both'='#440154', 'SAMase'='#fde725', 'Ocr'='#35b779', 'No'='#ffffff')
# 
# colors_bac <- c('Escherichia'='#31aff5', 'Vibrionales'= '#faba39', 
#                 'Pseudomonadales'='#83ff52', 
#                 'Enterobacterales'='#440154', 'Other'='gray90')
# 
# colors_phage <- c('Molineuxvirinae'='#f1605d', 'Studiervirinae'='#721f81', 
#                   'Colwellvirinae'='#fcfdbf', 'Other'='gray90')

colors_ <- c('Both'=col_pal[1], 'SAMase'=col_pal[2], 'Ocr'=col_pal[3], 'No'='#ffffff')

colors_bac <- c('Escherichia'=col_pal[4], 'Vibrionales'= col_pal[2], 
                'Pseudomonadales'=col_pal[1], 
                'Enterobacterales'=col_pal[3], 'Other'=col_pal[9])

colors_phage <- c('Molineuxvirinae'=col_pal[6], 'Studiervirinae'=col_pal[5], 
                  'Colwellvirinae'=col_pal[7], 'Other'=col_pal[8])


# p1 <- gheatmap(t, df1, offset=0.8, width=.1,
#                colnames_position = 'top',
#                custom_column_labels = c(''),  font.size = 1.5,
#                colnames_angle=0, colnames_offset_y = 0.25) +
#                # scale_fill_viridis_d(name="Known SAMase")
#   scale_fill_manual(values=colors_, name="Known",
#                     breaks = c('Both', 'Ocr',
#                                'SAMase',  'No'))

# p2 <- p1 + new_scale_fill()
p2 <- gheatmap(t, df2, offset=0, width=.15,
               custom_column_labels = c('Host'), font.size = 4,
               colnames_angle=0)  +
  scale_fill_manual(values=colors_bac, name="Host",
                    breaks = c('Enterobacterales', 'Escherichia',
                               'Pseudomonadales', 
                               'Vibrionales',  'Other')) 
p2
p3 <- p2 + new_scale_fill()
p3 <- gheatmap(p3, df3,  offset=6, width=.15, 
               custom_column_labels = c('Phage'), font.size = 4,
               colnames_angle=0) +
  scale_fill_manual(values=colors_phage, name="Virus Subfamily",
                    breaks = c('Studiervirinae', 'Molineuxvirinae', 
                               'Colwellvirinae',  'Other')) + ggtitle('SAMase')
p3
ggsave('pics/samase_tree_jan_may.svg', width=15, height = 9)


df2 <- cbind(df2, row.names(df2))
color_all_df <- color_branches_df %>% left_join(as.data.frame(df2), by=c('seq_id'='V2'))
t <- ggtree(samase_tree,
            branch.length = 'none',
            layout = 'circular',
            linewidth=1.2,
            aes(color = clu, label = clu)
)  %<+% color_all_df +
  geom_point2(aes(fill = order, subset=(label %in% color_all_df$seq_id)), color='black', shape=21, size=2.5) +
  #geom_point2(aes(subset=(label %in% clu_289$seq_id)), shape=21, size=2.5, fill='#31aff5') +
  #geom_point2(aes(subset=(label %in% clu_308$seq_id)), shape=21, size=2.5, fill='#721f81') +
  # geom_hilight(node=478, fill="purple", alpha=0.3) +
  geom_nodelab(color='firebrick')  + 
  # geom_tiplab(show.legend = FALSE, align = TRUE) +
  scale_color_brewer(na.value = "black", palette = 'Paired')+
  scale_fill_brewer(na.value = "black", palette = 'Paired')

t





# drops <- c('KJ183192.1', 'JQ780163.1', 'OP413828.1')
ocr_tree <- read.tree("data_autographiviridae_refseq/known_proteins/trees/ocr_bootstrap_model_selection.iqtree.contree")
# ocr_tree <- groupOTU(ocr_tree, c('NC_001604.1'))
a <- as.integer(ocr_tree$node.label)
a[is.na(a)] <- 0
a <- ifelse(a >= 95,  intToUtf8(9728), "")
ocr_tree$node.label <- a

t <- ggtree(ocr_tree,
            branch.length = 'none',
            layout='circular',
            linewidth=1.6,
            aes(color = clu, label = clu)
            )  %<+% color_branches_ocr_df +
  geom_point2(aes(subset=(label == 'NC_001604.1')), shape=23, size=2.5, fill='red') +
  # geom_point2(aes(subset=(label %in% clu_107$seq_id)), shape=21, size=2.5, fill='#35b779') +
  # geom_point2(aes(subset=(label %in% clu_135$seq_id)), shape=21, size=2.5, fill='#fde725') +
  # geom_hilight(node=478, fill="purple", alpha=0.3) +
  # geom_tippoint(aes(alpha = group), col = "red") +
   scale_color_manual(values = c(0,1), aesthetics = "alpha") +
 # geom_tiplab(offset = 0.22, align=T, size=2)+
 geom_nodelab(color='red') +
  theme(legend.position = 'none') +
  scale_color_brewer(na.value = "black", palette = 'Dark2')# + geom_tiplab2(size=2, hjust=8)
t

# write.table(df2, 'metadata/host_ids.tsv', col.names = F,
#             quote = F, sep='\t')

# p1 <- gheatmap(t, df1, offset=0.8, width=.1,
#                colnames_position = 'top',
#                custom_column_labels = c(''),  font.size = 1.5,
#                colnames_angle=0, colnames_offset_y = 0.25) +
#             scale_fill_viridis_d(name="Known Ocr")
#   # scale_fill_manual(values=colors_, name="Known",
  #                   breaks = c('Both', 'Ocr', 
  #                              'SAMase',  'No')) 
# p1
# p2 <- p1 + new_scale_fill()
p2 <- gheatmap(t, df2, offset=0, width=.15,
               custom_column_labels = c('Host'), font.size = 4,
               colnames_angle=0)  +
  scale_fill_manual(values=colors_bac, name="Host",
                    breaks = c('Enterobacterales', 'Escherichia',
                               'Pseudomonadales', 
                               'Vibrionales',  'Other')) 
p2
p4 <- p2 + new_scale_fill()
p4 <- gheatmap(p4, df3, offset=4.5, width=.15, 
               custom_column_labels = c('Phage'), font.size = 4,
               colnames_angle=0) +
  scale_fill_manual(values=colors_phage, name="Virus Subfamily",
                    breaks = c('Studiervirinae', 'Molineuxvirinae', 
                               'Colwellvirinae',  'Other')) + ggtitle('Ocr')

p4 
ggsave('pics/ocr_tree_may.svg', width=15, height = 9)

p3 + p4
