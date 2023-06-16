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
rnaps_tree <- drop.tip(rnaps_tree, drops)
t <- ggtree(rnaps_tree) + geom_tiplab(align=T, size=8) +
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

plot_upstream <- function(cl_num){
 seqs <- as.vector(unique(g0 %>% filter(cluster_num == cl_num) %>% select(seq_id)))$seq_id  # choose sequences to plot
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
  geom_gene(aes(fill=`cluster_main_prod`), intron_shape=0, size=8) +  # add gene arrows
  geom_gene_text(aes(label=`cluster_num`), size=6) +  # add gene cluster text
  geom_feat(color='red') +  # add TDRs
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', legend.text = element_text(size=13), axis.text.x = element_text(size=15))  # change font size and legend position
} else {p <- gggenomes(seqs=s0, genes=g0, feats = f0) %>%
  gggenomes::pick_seqs(seqs) %>%   # choose chromosomes (phages)
  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='drop') + # focus on particular region +- 100 bp
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=6, nudge_y = -0.3) + # label each sequence by this caption
  geom_gene(aes(fill=`cluster_main_prod`), intron_shape=0, size=8) +  # add gene arrows
  geom_gene_text(aes(label=`cluster_num`), size=6) +  # add gene cluster text
  geom_feat(color='red') +  # add TDRs
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', legend.text = element_text(size=13), axis.text.x = element_text(size=15))  # change font size and legend position
}
t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5)) 

path_ <- paste0('pics/upstreams/', cl_num, '.pdf')
ggsave(path_, dpi=200,  height = 1.29714286 * length(seqs), width = 20, limitsize = F)
}

  lapply(c(38, 59, 91, 115, 117, 189, 596, 625),  plot_upstream)

