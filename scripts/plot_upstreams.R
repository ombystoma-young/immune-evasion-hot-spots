library(tidyverse)
library(gggenomes)
library(patchwork)  # arrange multiple plots
library(ggtree)
library(treeio)
library(viridis)


setwd('work_dir/anti_defence/anti_defence_pipeline/')
s0 <- read_seqs("data_autographiviridae_meta/filtered_mags_flatten/concatenated_genomes.fna")  # sequence
g0 <- read_feats("data_autographiviridae_meta/upstreams/early_with_clusters_phrogs.gff")  # proteins 
f0 <- read.table("data_autographiviridae_meta/tdrs/best_tdrs.tsv", sep='\t', header=FALSE)  

f0 <- bind_rows(select(f0, seq_id=V1, start=V2, end=V3),
                      select(f0, seq_id=V1, start=V6, end=V6))
curated <- read.table('metadata/filtered_upstreams_nuccore.id')

g0 <- g0 %>% mutate(seq_id = str_replace_all(seq_id, '@', '_'))
s0 <- s0 %>% mutate(seq_id = str_replace_all(seq_id, '@', '_'))
f0 <- f0 %>% mutate(seq_id = str_replace_all(seq_id, '@', '_'))
curated <- curated %>% mutate(V1 = str_replace_all(V1, '@', '_'))

not_in_curated <- s0 %>% filter(! seq_id %in% curated$V1)

rnaps_tree <- read.tree("data_autographiviridae_meta/known_adgs_analysis/trees/rnap_bootstrap_model_selection.iqtree.contree")
t <- ggtree(rnaps_tree) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))

flips__ <- g0 %>% group_by(seq_id, strand) %>% 
  summarise(n()) %>%
  pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>% 
  ungroup() %>%
  filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>% 
  select(seq_id) %>% as.vector()

g0 <- g0 %>% mutate(clan_col = ifelse(str_detect(clan, 'mono'), NA, clan))

p <- gggenomes(seqs=s0, genes=g0, feats = f0)  %>%
  gggenomes::flip(flips__$seq_id) %>%  # flip sequences
  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='keep') + # focus on particular region +- 100 bp
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_id), size=10, nudge_y = -0.22) + # label each sequence by this caption
  geom_gene(aes(fill=`clan_col`), intron_shape=0, size=12) +  # add gene arrows
  geom_gene_text(aes(label=`clu`), size=10) +  # add gene cluster text
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=35))  # change font size and legend position

ggsave('data_autographiviridae_meta/pics/upstreams_march.pdf', t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5)) + 
                              theme(legend.position = 'none'), 
                            dpi=600,  height = 1.2 * 465, 
                          width = 50, limitsize = FALSE)

plot_upstream <- function(cl_num){
 seqs <- as.vector(unique(g0 %>% filter(clu %in% cl_num) %>% select(seq_id)))$seq_id  # choose sequences to plot
flips <- as.vector(unique(g0[g0$seq_id %in% seqs, ] %>%   # find sequences to flip (- strand -> + strand)
                            filter(strand == '-') %>% select(seq_id)))$seq_id
drops <- as.vector(s0 %>% filter(!seq_id %in% seqs) %>% select(seq_id))$seq_id
rnaps_tree_with_drops <- drop.tip(rnaps_tree, drops)

t <- ggtree(rnaps_tree_with_drops) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))
if (length(flips) != 0){
p <- gggenomes(seqs=s0, genes=g0, feats = f0) %>%
  gggenomes::pick_seqs(seqs) %>%   # choose chromosomes (phages)
  gggenomes::flip(flips) %>%  # flip sequences
  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='drop') + # focus on particular region +- 100 bp
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=6, nudge_y = -0.3) + # label each sequence by this caption
  geom_gene(aes(fill=`clan`), intron_shape=0, size=8) +  # add gene arrows
  geom_gene_text(aes(label=`clu`), size=6) +  # add gene cluster text
  geom_feat(color='red') +  # add TDRs
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', legend.text = element_text(size=13), axis.text.x = element_text(size=15))  # change font size and legend position
} else {p <- gggenomes(seqs=s0, genes=g0, feats = f0) %>%
  gggenomes::pick_seqs(seqs) %>%   # choose chromosomes (phages)
  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='drop') + # focus on particular region +- 100 bp
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=6, nudge_y = -0.3) + # label each sequence by this caption
  geom_gene(aes(fill=`clan`), intron_shape=0, size=8) +  # add gene arrows
  geom_gene_text(aes(label=`clu`), size=6) +  # add gene cluster text
  geom_feat(color='red') +  # add TDRs
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', legend.text = element_text(size=13), axis.text.x = element_text(size=15))  # change font size and legend position
}
t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5)) 

path_ <- paste0('pics/upstreams/', '322_mono', '.svg')
ggsave(path_, dpi=200,  height = 1.29714286 * length(seqs), width = 25, limitsize = F)
}

  # lapply(c( 117, 189, 596, 625),  plot_upstream)
plot_upstream(c(2472))



with_ocr_samase <- g0 %>% filter(`cluster_num` %in% c(3, 100, 125, 5, 26, 54, 252, 270, 325, 440)) %>% select(seq_id) %>% unique()



drops <- c('KJ183192.1', 'JQ780163.1', 'OP413828.1')
rnaps_tree <- read.tree("define_datasets/trees/polymerases_all.iqtree.treefile")
rnaps_tree <- drop.tip(rnaps_tree, drops)
rnaps_tree <- drop.tip(rnaps_tree, with_ocr_samase$seq_id)
t <- ggtree(rnaps_tree) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))
t

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
  geom_gene_text(aes(label=`cluster_num`), size=10) +  # add gene cluster text
  #geom_feat(color='red') +  # add TDRs
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=35))  # change font size and legend position

ggsave('pics/upstreams03_without_ocr_samase.pdf', t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5)) + 
         theme(legend.position = 'none'), 
       dpi=600,  height = 1.2 * 280, 
       width = 55, limitsize = FALSE)




g0 <- g0 %>% filter(!(seq_id == 'MZ851152.1' & `end` > 35477)) 


p <- gggenomes(seqs=s0, genes=g0, feats = f0) %>%
  gggenomes::pick_seqs(c('MZ234024.1', 'MT862763.1', 'ON637250.1', 
                           'ON148527.1', 'ON602765.1', 'MW394389.1', #'OL744215.1',
                         'NC_047858.1', 'NC_028916.1', 
                         'MW671054.1', 'MT740748.1', 'MZ851152.1', 'ON604651.1')) %>%   # choose chromosomes (phages)
  gggenomes::flip('MT862763.1', 'ON602765.1', 'MW671054.1', 'MZ851152.1') %>%  # flip sequences
  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='drop') + # focus on particular region +- 100 bp
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=paste0(str_split_i(seq_desc, ',', 1), ' (', seq_id, ')')), size=9, nudge_y = -0.5) + # label each sequence by this caption
  geom_gene(aes(fill=`cluster_main_prod`), intron_shape=0, size=8) +  # add gene arrows
  geom_gene_text(aes(label=`cluster_num`), size=8, nudge_y=-0.3) +  # add gene cluster text
  geom_feat(color='red') +  # add TDRs
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', legend.text = element_text(size=13), axis.text.x = element_text(size=15))  # change font size and legend position
p
ggsave('pics/all_cand.svg', p, dpi=600, width = 20, height = 15)









# 
# plot_upstream <- function(clans){
#   seqs <- as.vector(unique(g0 %>% filter(clan %in% clans) %>% select(seq_id)))$seq_id  # choose sequences to plot
#   flips <- as.vector(unique(g0[g0$seq_id %in% seqs, ] %>%   # find sequences to flip (- strand -> + strand)
#                               filter(strand == '-') %>% select(seq_id)))$seq_id
#   drops <- as.vector(s0 %>% filter(!seq_id %in% seqs) %>% select(seq_id))$seq_id
#   rnaps_tree_with_drops <- drop.tip(rnaps_tree, drops)
#   
#   t <- ggtree(rnaps_tree_with_drops) + geom_tiplab(align=T, size=8) +
#     xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))
#   if (length(flips) != 0){
#     p <- gggenomes(seqs=s0, genes=g0, feats = f0) %>%
#       gggenomes::pick_seqs(seqs) %>%   # choose chromosomes (phages)
#       gggenomes::flip(flips) %>%  # flip sequences
#       gggenomes::focus(.track_id = genes, .expand = 100, .overhang='drop') + # focus on particular region +- 100 bp
#       geom_seq() +  # draw contig/chromosome lines
#       geom_seq_label(aes(label=seq_desc), size=6, nudge_y = -0.3) + # label each sequence by this caption
#       geom_gene(aes(fill=`clan`), intron_shape=0, size=8) +  # add gene arrows
#       geom_gene_text(aes(label=`clu`), size=6) +  # add gene cluster text
#       geom_feat(color='red') +  # add TDRs
#       scale_fill_discrete("Main cluster product") +  # change fill, genes
#       theme(legend.position = 'bottom', legend.text = element_text(size=13), axis.text.x = element_text(size=15))  # change font size and legend position
#   } else {p <- gggenomes(seqs=s0, genes=g0, feats = f0) %>%
#     gggenomes::pick_seqs(seqs) %>%   # choose chromosomes (phages)
#     gggenomes::focus(.track_id = genes, .expand = 100, .overhang='drop') + # focus on particular region +- 100 bp
#     geom_seq() +  # draw contig/chromosome lines
#     geom_seq_label(aes(label=seq_desc), size=6, nudge_y = -0.3) + # label each sequence by this caption
#     geom_gene(aes(fill=`clan`), intron_shape=0, size=8) +  # add gene arrows
#     geom_gene_text(aes(label=`clu`), size=6) +  # add gene cluster text
#     geom_feat(color='red') +  # add TDRs
#     scale_fill_discrete("Main cluster product") +  # change fill, genes
#     theme(legend.position = 'bottom', legend.text = element_text(size=13), axis.text.x = element_text(size=15))  # change font size and legend position
#   }
#   t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5)) 
#   
#   path_ <- paste0('pics/upstreams/', '481_mono', '.pdf')
#   ggsave(path_, dpi=200,  height = 1.49714286 * length(seqs), width = 25, limitsize = F)
# }
# 
# plot_upstream(c('481_mono'))
