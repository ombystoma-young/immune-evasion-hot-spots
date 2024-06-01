library(tidyverse)
library(patchwork)  # arrange multiple plots
library(ggtree)
library(treeio)
library(gggenomes)
setwd('work_dir/anti_defence/anti_defence_pipeline')

s0 <- read_seqs("data_autographiviridae/genomes/co/phages_genomes_concat.fna")  # sequence
g0 <- read_feats("upstream_search/representative_genomes.gff")  # proteins 
f0 <- read.table("minimap2_out/TDRs_all.tsv", sep='\t', header=TRUE)  # domains, format: `chr(phage)	start	end	CDS_id	domain`

f0 <- bind_rows(select(f0, seq_id=seq_id, start=start, end=end, de),
                select(f0, seq_id=seq_id, start=start2, end=end2, de))


j_1 <- read.table('define_datasets/dataset_1_genomes_modified.txt')

seqs <- s0 %>% filter(seq_id %in% j_1$V1)
genes_ <- g0 %>% filter(seq_id %in% j_1$V1)
flips__ <- genes_ %>% group_by(seq_id, strand) %>% 
            summarise(n()) %>%
            pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>% 
            ungroup() %>%
            filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>% 
            select(seq_id) %>% as.vector()

rnaps_tree <- read.tree("define_datasets/trees/polymerases_dataset_1.iqtree.treefile")

t <- ggtree(rnaps_tree) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))
t



p <- gggenomes(seqs=seqs, genes=genes_, feats = f0)  %>%
  flip(flips__$seq_id) +
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=8, nudge_y = -0.19) + # label each sequence by this caption
  geom_gene(aes(fill=product), intron_shape=0, size=8) +  # add gene arrows
  #geom_gene(intron_shape=0, size=8) +  # add gene arrows
  # geom_gene_text(aes(label=product), size = 3) +  # add gene cluster
  geom_feat(color='red') +  # add TDRs
  #scale_fill_discrete("Product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=23),
        legend.text = element_text(size=23))  # change font size and legend position


d <- t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5))

ggsave('pics/have_tdrs_with_tree_my_ann_dataset1.pdf', d, dpi=600,  height = 0.9 * nrow(seqs), width = 40, limitsize = FALSE)



j_1 <- read.table('define_datasets/dataset_2_genomes_modified.txt')

seqs <- s0 %>% filter(seq_id %in% j_1$V1)
genes_ <- g0 %>% filter(seq_id %in% j_1$V1)
flips__ <- genes_ %>% group_by(seq_id, strand) %>% 
  summarise(n()) %>%
  pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>% 
  ungroup() %>%
  filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>% 
  select(seq_id) %>% as.vector()
rnaps_tree <- read.tree("define_datasets/trees/polymerases_dataset_2.iqtree.treefile")

t <- ggtree(rnaps_tree) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))
t



p <- gggenomes(seqs=seqs, genes=genes_, feats = f0) +#%>%
  #flip(flips__$seq_id) + 
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=8, nudge_y = -0.19) + # label each sequence by this caption
  geom_gene(aes(fill=product), intron_shape=0, size=8) +  # add gene arrows
  #geom_gene(intron_shape=0, size=8) +  # add gene arrows
  # geom_gene_text(aes(label=product), size = 3) +  # add gene cluster
  geom_feat(color='red') +  # add TDRs
  #scale_fill_discrete("Product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=23),
        legend.text = element_text(size=23))  # change font size and legend position

d <- t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5))

ggsave('pics/have_tdrs_with_tree_my_ann_dataset2.pdf', d, dpi=600,  height = 1.2 * nrow(seqs), width = 40, limitsize = FALSE)



j_1 <- read.table('define_datasets/dataset_3_genomes_modified.txt')

seqs <- s0 %>% filter(seq_id %in% j_1$V1)
genes_ <- g0 %>% filter(seq_id %in% j_1$V1)
flips__ <- genes_ %>% group_by(seq_id, strand) %>% 
  summarise(n()) %>%
  pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>% 
  ungroup() %>%
  filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>% 
  select(seq_id) %>% as.vector()
rnaps_tree <- read.tree("define_datasets/trees/polymerases_dataset_3.iqtree.treefile")
rnaps_tree_with_drops <- drop.tip(rnaps_tree, c('JQ780163.1'))
rnaps_tree_with_drops <- groupOTU(rnaps_tree_with_drops, "NC_001604.1")

cat('Номер таксона T7:', which(rnaps_tree_with_drops$tip.label == "NC_001604.1"))

t <- ggtree(rnaps_tree_with_drops, branch.length = 'none') + geom_tiplab(align=T, size=8) +
  geom_hilight(node = 617, fill = "gold") +
  geom_text2(aes(subset=!isTip, label=node), size = 8, col = "tomato") +
  xlim(0,90) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7)) 
t

p <- gggenomes(seqs=seqs, genes=genes_, feats = f0)   %>%
  flip(flips__$seq_id) + 
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=8, nudge_y = -0.19) + # label each sequence by this caption
  geom_gene(aes(fill=product), intron_shape=0, size=8) +  # add gene arrows
  #geom_gene(intron_shape=0, size=8) +  # add gene arrows
  # geom_gene_text(aes(label=product), size = 3) +  # add gene cluster
  geom_feat(color='red') +  # add TDRs
  #scale_fill_discrete("Product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=23),
        legend.text = element_text(size=23))  # change font size and legend position


d <- t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5))

ggsave('pics/have_tdrs_with_tree_my_ann_dataset3.pdf', d, dpi=600,  height = 0.9 * nrow(seqs), width = 40, limitsize = FALSE)



j_1 <- read.table('define_datasets/dataset_3_genomes_modified.txt')

seqs <- s0 %>% filter(seq_id %in% j_1$V1)
genes_ <- g0 %>% filter(seq_id %in% j_1$V1)
flips__ <- genes_ %>% group_by(seq_id, strand) %>% 
  summarise(n()) %>%
  pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>% 
  ungroup() %>%
  filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>% 
  select(seq_id) %>% as.vector()
rnaps_tree <- read.tree("define_datasets/trees/polymerases_dataset_3.iqtree.treefile")
rnaps_tree_with_drops <- drop.tip(rnaps_tree, c('JQ780163.1', 
                                                'NC_023863.2',
                                                'LR797151.1'))
#rnaps_tree_with_drops <- groupOTU(rnaps_tree_with_drops, "NC_001604.1")

cat('Номер таксона T7:', which(rnaps_tree_with_drops$tip.label == "NC_001604.1"))

t <- ggtree(rnaps_tree_with_drops, branch.length = 'none') + geom_tiplab(align=T, size=8) +
  geom_hilight(node = 613, fill = "gold") +
  geom_text2(aes(subset=!isTip, label=node), size = 8, col = "tomato") +
  xlim(0,90) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7)) 
t

p <- gggenomes(seqs=seqs, genes=genes_, feats = f0)   %>%
  flip(flips__$seq_id) + 
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=8, nudge_y = -0.19) + # label each sequence by this caption
  geom_gene(aes(fill=product), intron_shape=0, size=8) +  # add gene arrows
  #geom_gene(intron_shape=0, size=8) +  # add gene arrows
  # geom_gene_text(aes(label=product), size = 3) +  # add gene cluster
  geom_feat(color='red') +  # add TDRs
  #scale_fill_discrete("Product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=23),
        legend.text = element_text(size=23))  # change font size and legend position


d <- t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5))

ggsave('pics/have_tdrs_with_tree_my_ann_dataset3_exclude_large.pdf', d, dpi=600,  height = 0.9 * nrow(seqs), width = 40, limitsize = FALSE)





j_1 <- read.table('define_datasets/dataset_4_genomes_modified.txt')

seqs <- s0 %>% filter(seq_id %in% j_1$V1)
genes_ <- g0 %>% filter(seq_id %in% j_1$V1)
flips__ <- genes_ %>% group_by(seq_id, strand) %>% 
  summarise(n()) %>%
  pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>% 
  ungroup() %>%
  filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>% 
  select(seq_id) %>% as.vector()
rnaps_tree <- read.tree("define_datasets/trees/polymerases_dataset_4.iqtree.treefile")
#rnaps_tree_with_drops <- groupOTU(rnaps_tree_with_drops, "NC_001604.1")
rnaps_tree_with_drops <- drop.tip(rnaps_tree, c('KJ183192.1'))
#cat('Номер таксона T7:', which(rnaps_tree_with_drops$tip.label == "NC_001604.1"))

t <- ggtree(rnaps_tree_with_drops, branch.length = 'none') + geom_tiplab(align=T, size=10) +
  geom_hilight(node = 712, fill = "gold") +
  geom_text2(aes(subset=!isTip, label=node), size = 8, col = "tomato") +
  xlim(0,90) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7)) 
t

p <- gggenomes(seqs=seqs, genes=genes_, feats = f0)   %>%
  flip(flips__$seq_id) + 
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=10, nudge_y = -0.19) + # label each sequence by this caption
  geom_gene(aes(fill=product), intron_shape=0, size=10) +  # add gene arrows
  #geom_gene(intron_shape=0, size=8) +  # add gene arrows
  # geom_gene_text(aes(label=product), size = 3) +  # add gene cluster
  geom_feat(color='red') +  # add TDRs
  #scale_fill_discrete("Product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=23),
        legend.text = element_text(size=23))  # change font size and legend position


d <- t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5))

ggsave('pics/have_tdrs_with_tree_my_ann_dataset4.pdf', d, dpi=600,  height = 1 * nrow(seqs), width = 40, limitsize = FALSE)



j_1 <- read.table('define_datasets/dataset_4_genomes_modified.txt')

seqs <- s0 %>% filter(seq_id %in% j_1$V1)
genes_ <- g0 %>% filter(seq_id %in% j_1$V1)
flips__ <- genes_ %>% group_by(seq_id, strand) %>% 
  summarise(n()) %>%
  pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>% 
  ungroup() %>%
  filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>% 
  select(seq_id) %>% as.vector()
rnaps_tree <- read.tree("define_datasets/trees/polymerases_dataset_4.iqtree.treefile")
#rnaps_tree_with_drops <- groupOTU(rnaps_tree_with_drops, "NC_001604.1")
rnaps_tree_with_drops <- drop.tip(rnaps_tree, c('KJ183192.1',
                                                'OQ055246.1',
                                                'NC_018088.1',
                                                'MZ568826.1',
                                                'NC_020842.1'))
#cat('Номер таксона T7:', which(rnaps_tree_with_drops$tip.label == "NC_001604.1"))

t <- ggtree(rnaps_tree_with_drops, branch.length = 'none') + geom_tiplab(align=T, size=10) +
  geom_hilight(node = 705, fill = "gold") +
  geom_text2(aes(subset=!isTip, label=node), size = 8, col = "tomato") +
  xlim(0,90) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7)) 
t

p <- gggenomes(seqs=seqs, genes=genes_, feats = f0)   %>%
  flip(flips__$seq_id) + 
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=10, nudge_y = -0.19) + # label each sequence by this caption
  geom_gene(aes(fill=product), intron_shape=0, size=10) +  # add gene arrows
  #geom_gene(intron_shape=0, size=8) +  # add gene arrows
  # geom_gene_text(aes(label=product), size = 3) +  # add gene cluster
  geom_feat(color='red') +  # add TDRs
  #scale_fill_discrete("Product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=23),
        legend.text = element_text(size=23))  # change font size and legend position


d <- t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5))

ggsave('pics/have_tdrs_with_tree_my_ann_dataset4_exclude_large.pdf', d, dpi=600,  height = 1 * nrow(seqs), width = 40, limitsize = FALSE)


# OTHER

j_1 <- read.table('define_datasets/dataset_5_genomes_modified.txt')
j_1 <-  append(j_1$V1, c('KJ183192.1', 'JQ780163.1'))
seqs <- s0 %>% filter(seq_id %in% j_1)
genes_ <- g0 %>% filter(seq_id %in% j_1)
flips__ <- genes_ %>% group_by(seq_id, strand) %>% 
  summarise(n()) %>%
  pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>% 
  ungroup() %>%
  filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>% 
  select(seq_id) %>% as.vector()

p <- gggenomes(seqs=seqs, genes=genes_, feats = f0)  %>%
  flip(flips__$seq_id) + 
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=seq_desc), size=8, nudge_y = -0.19) + # label each sequence by this caption
  geom_gene(aes(fill=product), intron_shape=0, size=8) +  # add gene arrows
  #geom_gene(intron_shape=0, size=8) +  # add gene arrows
  # geom_gene_text(aes(label=product), size = 3) +  # add gene cluster
  geom_feat(color='red') +  # add TDRs
  #scale_fill_discrete("Product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=23),
        legend.text = element_text(size=23))  # change font size and legend position
p
ggsave('pics/have_tdrs_without_tree_my_ann_dataset5.pdf', p, dpi=600,  height = 0.9 * nrow(seqs), width = 40, limitsize = FALSE)





# Only trees

rnaps_tree <- read.tree("define_datasets/trees/polymerases_dataset_1.iqtree.treefile")

t <- ggtree(rnaps_tree, layout = 'circular') + geom_tiplab(align=T, size=2) +
  #geom_text2(aes(subset=!isTip, label=node), size = 2, col = "tomato") +
  xlim(0,4) #+ scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7)) 
t





