library(tidyverse)
library(gggenomes)
library(patchwork)  # arrange multiple plots
library(ggtree)
library(treeio)
library(viridis)


setwd('work_dir/anti_defence/anti_defence_pipeline/')
s0 <- read_seqs("blasted/phages_genomes_concat.fna")  # sequence
g0 <- read_feats("results/upstreams_with_clusters.gff")  # proteins 
f0 <- read.table("minimap2_out/TDRs_all.tsv", sep='\t', header=TRUE)  

f0 <- bind_rows(select(f0, seq_id=seq_id, start=start, end=end, de),
                      select(f0, seq_id=seq_id, start=start2, end=end2, de))
curated <- read.table('metadata/genomes_after_curation.tsv')
not_in_curated <- s0 %>% filter(! seq_id %in% curated$V1)

# drop_ <- c('OP413830.1',
#  'MW865291.1',
#  'NC_048034.1',
#  'NC_047707.1',
#  'NC_047706.1',
#  'LR797151.1',
#  'NC_023863.2',
#  'MK892742.1',
#  'LR796375.1',
#  'MT375528.1',
#  'MW202733.1',
#  'JQ780163.1',
#  'OX001577.1',
#  'OW991346.1',
#  'OP947227.1',
#  'OP947225.1',
#  'MN184885.1',
#  'NC_047935.1',
#  'LT961846.1',
#  'NC_047964.1',
#  'NC_029102.1',
#  'MZ428229.1')
# 
drops <- c('KJ183192.1', 'JQ780163.1', 'OP413828.1')
rnaps_tree <- read.tree("define_datasets/trees/polymerases_all.iqtree.treefile")
rnaps_tree_with_drops <- drop.tip(rnaps_tree, drops)
t <- ggtree(rnaps_tree_with_drops) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))
t

g0 <- g0 %>% mutate(cluster_num = as.integer(cluster_num)+1)
flips__ <- g0 %>% group_by(seq_id, strand) %>% 
  summarise(n()) %>% filter(seq_id != 'KJ183192.1') %>% 
  pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>% 
  ungroup() %>%
  filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>% 
  select(seq_id) %>% as.vector()

p <- gggenomes(seqs=s0, genes=g0, feats = f0)  %>%
  gggenomes::flip(flips__$seq_id) %>%  # flip sequences
  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='keep') + # focus on particular region +- 100 bp
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=10, nudge_y = -0.22) + # label each sequence by this caption
  geom_gene(aes(fill=`cluster_main_prod`), intron_shape=0, size=12) +  # add gene arrows
  #geom_gene(aes(fill=cluster_num), intron_shape=0, size=8) +  # add gene arrows
  geom_gene_text(aes(label=`cluster_num`), size=10) +  # add gene cluster text
  #geom_feat(color='red') +  # add TDRs
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=35))  # change font size and legend position

ggsave('pics/upstreams.pdf', t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5)) + 
                              theme(legend.position = 'none'), 
                            dpi=600,  height = 1.2 * 465, 
                          width = 50, limitsize = FALSE)





j_1 <- read.table('define_datasets/dataset_3_genomes_modified.txt')

seqs_ <-  s0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1) %>% 
  mutate(seq_desc = str_split(seq_desc, ",", simplify = T)[, 1])
genes_ <- g0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1) %>% mutate(cluster_num = as.integer(cluster_num)+1)
f0 <- f0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1)
flips__ <- genes_ %>% group_by(seq_id, strand) %>% 
  summarise(n()) %>%
  pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>% 
  ungroup() %>%
  filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>% 
  select(seq_id) %>% as.vector()

rnaps_tree <- read.tree("define_datasets/trees/polymerases_dataset_3.iqtree.treefile")
rnaps_tree_with_drops <- drop.tip(rnaps_tree, not_in_curated$seq_id) %>% drop.tip('JQ780163.1')
t <- ggtree(rnaps_tree_with_drops) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))
t

p <- gggenomes(seqs=seqs_, genes=genes_, feats = f0)  %>%
  gggenomes::flip(flips__$seq_id) %>%  # flip sequences
  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='keep') + # focus on particular region +- 100 bp
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=8, nudge_y = -0.19) + # label each sequence by this caption
  geom_gene(aes(fill=`cluster_main_prod`), intron_shape=0, size=8) +  # add gene arrows
  #geom_gene(aes(fill=cluster_num), intron_shape=0, size=8) +  # add gene arrows
  geom_gene_text(aes(label=`cluster_num`), size=8) +  # add gene cluster text
  #geom_feat(color='red') +  # add TDRs
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=18))  # change font size and legend position
t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5))
ggsave('pics/upstreams_dataset3.pdf', dpi=600,  height = 0.9 * nrow(seqs_), width = 40, limitsize = FALSE)



j_1 <- read.table('define_datasets/dataset_1_genomes_modified.txt')

seqs_ <-  s0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1) %>% 
  mutate(seq_desc = str_split(seq_desc, ",", simplify = T)[, 1])
genes_ <- g0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1) %>% mutate(cluster_num = as.integer(cluster_num)+1)
f0 <- f0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1)
flips__ <- genes_ %>% group_by(seq_id, strand) %>% 
  summarise(n()) %>%
  pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>% 
  ungroup() %>%
  filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>% 
  select(seq_id) %>% as.vector()

rnaps_tree <- read.tree("define_datasets/trees/polymerases_dataset_1.iqtree.treefile")
rnaps_tree_with_drops <- drop.tip(rnaps_tree, not_in_curated$seq_id) #%>% drop.tip('JQ780163.1')
t <- ggtree(rnaps_tree_with_drops) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))
t

p <- gggenomes(seqs=seqs_, genes=genes_, feats = f0)  %>%
  gggenomes::flip(flips__$seq_id) %>%  # flip sequences
  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='keep') + # focus on particular region +- 100 bp
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=8, nudge_y = -0.19) + # label each sequence by this caption
  geom_gene(aes(fill=`cluster_main_prod`), intron_shape=0, size=8) +  # add gene arrows
  #geom_gene(aes(fill=cluster_num), intron_shape=0, size=8) +  # add gene arrows
  geom_gene_text(aes(label=`cluster_num`), size=8) +  # add gene cluster text
  #geom_feat(color='red') +  # add TDRs
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=18))  # change font size and legend position
t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5))
ggsave('pics/upstreams_dataset1.pdf', dpi=600,  height = 0.9 * nrow(seqs_), width = 40, limitsize = FALSE)





j_1 <- read.table('define_datasets/dataset_2_genomes_modified.txt')

seqs_ <-  s0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1) %>% 
  mutate(seq_desc = str_split(seq_desc, ",", simplify = T)[, 1])
genes_ <- g0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1) %>% mutate(cluster_num = as.integer(cluster_num)+1)
f0 <- f0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1)
flips__ <- genes_ %>% group_by(seq_id, strand) %>% 
  summarise(n()) %>%
  pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>% 
  ungroup() %>%
  filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>% 
  select(seq_id) %>% as.vector()

rnaps_tree <- read.tree("define_datasets/trees/polymerases_dataset_2.iqtree.treefile")
rnaps_tree_with_drops <- drop.tip(rnaps_tree, not_in_curated$seq_id) #%>% drop.tip('JQ780163.1')
t <- ggtree(rnaps_tree_with_drops) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))
t

p <- gggenomes(seqs=seqs_, genes=genes_, feats = f0)  %>%
  #gggenomes::flip(flips__$seq_id) %>%  # flip sequences
  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='keep') + # focus on particular region +- 100 bp
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=8, nudge_y = -0.19) + # label each sequence by this caption
  geom_gene(aes(fill=`cluster_main_prod`), intron_shape=0, size=8) +  # add gene arrows
  #geom_gene(aes(fill=cluster_num), intron_shape=0, size=8) +  # add gene arrows
  geom_gene_text(aes(label=`cluster_num`), size=8) +  # add gene cluster text
  #geom_feat(color='red') +  # add TDRs
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=18))  # change font size and legend position
t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5))
ggsave('pics/upstreams_dataset2.pdf', dpi=600,  height = 1.2 * nrow(seqs_), width = 40, limitsize = FALSE)




j_1 <- read.table('define_datasets/dataset_4_genomes_modified.txt')

seqs_ <-  s0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1) %>% 
  mutate(seq_desc = str_split(seq_desc, ",", simplify = T)[, 1]) 
genes_ <- g0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1) %>% mutate(cluster_num = as.integer(cluster_num)+1)
f0 <- f0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1)
flips__ <- genes_ %>% filter(seq_id != 'KJ183192.1') %>%  group_by(seq_id, strand) %>% 
  summarise(n()) %>%
  pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>% 
  ungroup() %>%
  filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>% 
  select(seq_id) %>% as.vector() 

rnaps_tree <- read.tree("define_datasets/trees/polymerases_dataset_4.iqtree.treefile")
rnaps_tree_with_drops <- drop.tip(rnaps_tree, not_in_curated$seq_id) %>% drop.tip('KJ183192.1')
t <- ggtree(rnaps_tree_with_drops) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))
t

p <- gggenomes(seqs=seqs_, genes=genes_, feats = f0)  %>%
  gggenomes::flip(flips__$seq_id) %>%  # flip sequences
  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='keep') + # focus on particular region +- 100 bp
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=8, nudge_y = -0.19) + # label each sequence by this caption
  geom_gene(aes(fill=`cluster_main_prod`), intron_shape=0, size=8) +  # add gene arrows
  #geom_gene(aes(fill=cluster_num), intron_shape=0, size=8) +  # add gene arrows
  geom_gene_text(aes(label=`cluster_num`), size=8) +  # add gene cluster text
  #geom_feat(color='red') +  # add TDRs
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=18))  # change font size and legend position
t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5))
ggsave('pics/upstreams_dataset4.pdf', dpi=600,  height = 0.9 * nrow(seqs_), width = 40, limitsize = FALSE)

j_1 <- read.table('define_datasets/dataset_5_genomes_modified.txt')
seqs_ <-  s0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1) %>% 
  mutate(seq_desc = str_split(seq_desc, ",", simplify = T)[, 1])
genes_ <- g0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1) %>% mutate(cluster_num = as.integer(cluster_num)+1)
f0 <- f0 %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id %in% curated$V1)
flips__ <- genes_ %>% group_by(seq_id, strand) %>% 
  summarise(n()) %>%
  pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>% 
  ungroup() %>%
  filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>% 
  select(seq_id) %>% as.vector()

rnaps_tree <- read.tree("define_datasets/trees/polymerases_dataset_5.iqtree.treefile")
rnaps_tree_with_drops <- drop.tip(rnaps_tree, not_in_curated$seq_id) #%>% drop.tip('JQ780163.1')
t <- ggtree(rnaps_tree_with_drops) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))
t

p <- gggenomes(seqs=seqs_, genes=genes_, feats = f0)  %>%
  gggenomes::flip(flips__$seq_id) %>%  # flip sequences
  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='keep') + # focus on particular region +- 100 bp
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=8, nudge_y = -0.19) + # label each sequence by this caption
  geom_gene(aes(fill=`cluster_main_prod`), intron_shape=0, size=8) +  # add gene arrows
  #geom_gene(aes(fill=cluster_num), intron_shape=0, size=8) +  # add gene arrows
  geom_gene_text(aes(label=`cluster_num`), size=8) +  # add gene cluster text
  #geom_feat(color='red') +  # add TDRs
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=18))  # change font size and legend position
p
ggsave('pics/upstreams_dataset5.pdf', dpi=600,  height = 0.9 * nrow(seqs_), width = 40, limitsize = FALSE)














p <- gggenomes(seqs=s0, genes=g0, feats = f0) %>%
  pick_seqs(c('MN184885.1')) +  # change order of chromosomes (phages)
  # focus(.track_id = genes, .expand = 4620, .overhang='drop') +
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=5) + # label each sequence by this caption
  geom_gene(aes(fill=paste0(`cluster_num`, ': ', `cluster_main_prod`)), intron_shape=0, size=5) +  # add gene arrows
  geom_gene_tag(aes(label=`cluster_num`), nudge_y=0.1, size=6) +  # add gene cluster
  geom_feat(color='red') +  # add TDRs
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=18))  # change font size and legend position
p




# 
# suicide_function <- function(x){
#   s1 <- s0 %>% filter(seq_id == x)
#   g1 <- g0 %>% filter(seq_id == x)
#   f1 <- f0 %>% filter(seq_id == x)
#   p <- gggenomes(seqs=s1, genes=g1, feats = f1) %>%
#     # pick_seqs(x)  %>%   # change order of chromosomes (phages)
#     focus(.track_id = genes, .expand = 1000, .overhang='drop') +
#     geom_seq() +  # draw contig/chromosome lines
#     geom_seq_label(aes(label=seq_id), size=5) + # label each sequence by this caption
#     geom_gene(aes(fill=paste0(`cluster_num`, ': ', `cluster_main_prod`)), intron_shape=0, size=5) +  # add gene arrows
#     geom_gene_tag(aes(label=`cluster_num`), nudge_y=0.1, size=4) +  # add gene cluster
#     geom_feat(color='red') +  # add TDRs
#     scale_fill_discrete("Main cluster product") +  # change fill, genes
#     theme(axis.text.x = element_text(size=18))  # change font size and legend position
#     name_plot <- paste0('pics/upstreams/', s1[s1$seq_id == x, ]$seq_desc, '.png')
#     ggsave(name_plot, plot=p, dpi=150,
#            height = 5, width = 10)
# }
# 
# suicide_function('NC_047803.1')
# lapply(unique(g0$seq_id), suicide_function)


seqs <- as.vector(unique(g0 %>% filter(cluster_num == 68) %>% select(seq_id)))$seq_id  # choose sequences to plot
flips <- as.vector(unique(g0[g0$seq_id %in% seqs, ] %>%   # find sequences to flip (- strand -> + strand)
                            filter(strand == '-') %>% select(seq_id)))$seq_id
g0[g0$cluster_num == 14, ]$cluster_main_prod <-  'hypothetical protein'
g0[g0$cluster_num == 40, ]$cluster_main_prod <- 'protein kinase'
g0[g0$cluster_num == 5, ]$cluster_main_prod <- 'S-adenosyl-L-methionine hydrolase'
p <- gggenomes(seqs=s0, genes=g0, feats = f0) %>%
  gggenomes::pick_seqs(seqs) %>%   # choose chromosomes (phages)
  gggenomes::flip(flips) %>%  # flip sequences
  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='drop') + # focus on particular region +- 100 bp
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=6, nudge_y = -0.3) + # label each sequence by this caption
  geom_gene(aes(fill=`cluster_main_prod`), intron_shape=0, size=8) +  # add gene arrows
  geom_gene_text(aes(label=`cluster_num`), size=6) +  # add gene cluster text
  geom_feat(color='red') +  # add TDRs
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', legend.text = element_text(size=13), axis.text.x = element_text(size=15))  # change font size and legend position
p

ggsave('pics/both_ocr_samh.png', dpi=300,  height = 0.9714286 * length(seqs), width = 17)
    
