library(tidyverse)
library(gggenomes)

setwd('work_dir/anti_defence/anti_defence_pipeline/')


read.table('results/joined.tsv', header = T)

# BUILD PICTURE



seqs <- joined %>% filter(tdr_left_start < 1000) %>% 
  filter(TDR_in != 0) %>% 
  filter(u_length_ > 10000) %>% select(chrm)

s0 <- read_seqs("blasted/phages_genomes_concat.fna")  # sequence
g0 <- read_feats("upstream_search/representative_genomes.gff")  # proteins 
f0 <- read.table("minimap2_out/TDRs_all.tsv", sep='\t', header=TRUE)  # domains, format: `chr(phage)	start	end	CDS_id	domain`


f0 <- bind_rows(select(f0, seq_id=seq_id, start=start, end=end, de),
                select(f0, seq_id=seq_id, start=start2, end=end2, de))

sequences <- s0 %>% filter(seq_id %in% seqs$chrm)
flips__ <- as.vector(unique(g0[g0$seq_id %in% sequences$seq_id, ] %>% 
                              filter(strand == '-') %>% select(seq_id)))$seq_id

#g_stealed <- read_feats('metadata/not_all_genomic.gff')

# g_stealed <- g_stealed %>% filter(seq_id %in% seqs$chrm) %>% filter(`type` == 'CDS')
# sequences[sequences$seq_id %in% g_stealed$seq_id, ]$file_id <- 'not_all_genomic'
# seqs_from_our <- seqs %>% filter(!chrm %in% g_stealed$seq_id)
#genes_our <- g0 %>% filter(seq_id %in% seqs$chrm)
#genes_ <- merge(g_stealed, genes_our, all.x = T, all.y = T)

# rm(g0, g_stealed)



library(patchwork)  # arrange multiple plots
library(ggtree)
library(treeio)
rnaps_tree <- read.tree("upstream_search/tree_far_/polymerases_far.iqtree.treefile")

rnaps_tree_with_drops <- drop.tip(rnaps_tree, c( 'OP413829.1'
, 'MN794003.1'
, 'ON602754.1'
, 'MT197175.1'
, 'MZ634339.1'
, 'MN794000.1'
, 'MW250641.1'
, 'OP117450.1'
, 'ON080941.1'
, 'MZ099557.1'
, 'OP019135.1'
, 'NC_013638.1'
, 'MT104468.1'
, 'NC_047860.1'
, 'OK490494.1'
, 'MZ462995.1'
, 'ON602764.1'
, 'MN807240.1'
, 'MW023914.1'
, 'NC_047980.1'
, 'ON602742.1'
, 'ON995367.1'
, 'NC_048052.1'
, 'ON754979.1'
, 'OM937123.1'
, 'NC_048161.1'
, 'NC_048200.1'
, 'OW991346.1'
, 'MW671054.1'
, 'KY883649.1'
, 'KY883645.1'))






# write.table(genes_, 'genes_reload.tsv', sep='\t', row.names = FALSE, quote = FALSE)
# write.table(sequences, 'sequences_reload.tsv', sep='\t', row.names = FALSE, quote = FALSE)
# write.table(f0, 'feats_reload.tsv', sep='\t', row.names = FALSE, quote = FALSE)

t <- ggtree(rnaps_tree_with_drops) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))

p <- gggenomes(seqs=sequences, genes=g0, feats = f0) +#  %>%
 # flip(flips__) + 
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

ggsave('pics/have_tdrs_with_large_dist_tree_my_ann.pdf', d, dpi=600,  height = 0.9 * nrow(sequences), width = 40, limitsize = FALSE)


# sequences %>% filter(seq_id %in% rnaps_tree_with_drops[["tip.label"]]) %>% select(seq_id, seq_desc)
# ggtree(rnaps_tree_with_drops, layout='circular', branch.length = 'none') + geom_tiplab(align=T, size=5) +
#   geom_text2(aes(subset=!isTip, label=node), size = 3, col = "tomato")
# ggsave('pics/tree.pdf', dpi=600,  height = 15, width = 15)
# library(ape)
# 
# try <- extract.clade(rnaps_tree_with_drops, 182)
# # try <- drop.tip(rnaps_tree_with_drops, nodepath(rnaps_tree_with_drops, 107, 177))
# ggtree(try,branch.length = 'none', layout='circular') + geom_tiplab(align=T, size=5) +
#   geom_text2(aes(subset=!isTip, label=node), size = 3, col = "tomato")
# 
# 
# t <- ggtree(try) + geom_tiplab(align=T, size=8) +
#   xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))
# d <- t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5))
# ggsave('pics/have_tdrs_with_tree_182.pdf', d, dpi=600,  height = 0.5 * nrow(sequences), width = 50, limitsize = FALSE)
