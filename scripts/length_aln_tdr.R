# more statistics

library(tidyverse)
 
setwd('work_dir/anti_defence/anti_defence_pipeline/')


# PSI-BLAST
names_psi <- c('query_acc.ver', 'subject_acc.ver',
               'identity', 'alignment length',
               'mismatches', 'gap_opens', 'q._start',
               'q._end', 's._start', 's._end', 
               'evalue', 'bit_score', 'positives')
psi_blast_res <- read.table('psi_blast_more_genomes/psi_blast_round_three_attempt_3.csv',
                            sep=',', col.names = names_psi)
psi_blast_res %>% ggplot() +
  geom_histogram(aes(y=..density.., x = alignment.length, fill = ifelse(alignment.length>=750, 'tomato', 'grey95')), color='black') + 
  scale_fill_manual(labels = c("<750", ">=750"), values=c('grey95','tomato')) + 
  labs(fill = "Values") +
  xlab('Alignment length') +
  ylab('Density') +
  theme_minimal() +
  theme(text= element_text(size=13))


psi_blast_res %>% filter(alignment.length >= 750) %>% nrow()


# TDRs and intergenics
upstreams_names <- c('chrm', 'upstream_start', 'upstream_end',
                    'dot', 'pos', 'strand',
                    'pol_start', 'pol_end',
                    'tdr_left_start', 'tdr_left_end',
                    'tdr_right_start', 'tdr_right_end')
found_upstreams <- read.table('upstream_search/upstream.bed',
                              fill=T, col.names = upstreams_names)
found_upstreams_good <- found_upstreams %>% filter(!is.na(tdr_left_start))
length(unique(found_upstreams_good$chrm))

upstreams_names <- c('chrm', 'strand',
                     'pol_start', 'pol_end',
                     'tdr_left_start', 'tdr_left_end',
                     'tdr_right_start', 'tdr_right_end', 'u_length_')
upstream_true_length <- read.table('upstream_search/tdr_pol_dist.tsv',
                              fill=T, col.names = upstreams_names) %>%  mutate(u_length_ = abs(u_length_))

a <- upstream_true_length %>% 
  ggplot(aes(y=u_length_) ) +
  geom_histogram(aes(x=..density..), color= 'black', fill = 'white') + 
  theme_minimal() + 
  scale_y_log10() +
  scale_x_reverse() +
  xlab('Density') +
  ylab('Distance from RNAP to closest TDR, bp') +
  theme(text = element_text(size=15))
a


intergenics_names <- c('chrm', 'i_start', 'i_end', 'i_length_')
intergenics <- read.table('promoters_search/all_intergenic_with_length.tsv', 
                          col.names = intergenics_names)
intergenics <- intergenics %>% mutate(i_length_ = abs(i_length_))
longest_intergenics <- intergenics %>% 
  group_by(chrm) %>% slice_max(i_length_, n=1) 
b <- longest_intergenics %>%  
  ggplot(aes(x=i_length_)) +
  geom_histogram(aes(y=..density..), binwidth = 75, color= 'black', fill = 'white') + 
  theme_minimal() + 
  ylab('Density') +
  xlab('Length of longest intergenic region, bp') +
  theme(text = element_text(size=15))
b

longest_intergenics <- intergenics %>% filter(chrm %in% upstream_true_length$chrm) %>%
  group_by(chrm) %>% slice_max(abs(i_length_), n=1) 


c <- longest_intergenics %>%  
  ggplot(aes(x=i_length_)) +
  geom_histogram(aes(y=..density..), binwidth = 75, color= 'black', fill = 'white') + 
  theme_minimal() + 
  xlim(c(0,3336)) +
  ylab('Density') +
  xlab('Length of the longest intergenic region, bp') +
  theme(text = element_text(size=15))
c
library(ggpubr)
ggarrange(b, c, ncol = 1, nrow = 2, 
          common.legend = TRUE)  

ggsave('pics/intergenic_diagnostics.png', width=14, height = 10, dpi=400)



c <- c + scale_y_reverse() 
c

# filter(V1 %in% found_upstreams_good$V1) 



joined <- full_join(upstream_true_length, longest_intergenics)
joined <- joined %>% mutate(i_start  = ifelse(i_start > i_end, -i_start, i_start))
joined$left_in <- joined$tdr_left_end < joined$i_end & joined$tdr_left_start > joined$i_start
joined$right_in <- joined$tdr_right_end < joined$i_end & joined$tdr_right_start > joined$i_start
joined <- joined %>% mutate(TDR_in = factor(left_in + right_in))
  
d <- joined %>% ggplot(aes(x = i_length_, y = u_length_, color = TDR_in)) +
  geom_point(size = 3, aes(shape=ifelse(tdr_left_start<1000, TRUE, FALSE)), alpha = 0.75) + 
  scale_y_log10() + 
  scale_shape_manual('TDRs near ends  ', labels = c("No", "Yes"), values=c(15, 19)) + 
  scale_color_manual('TDRs in intergenic regions   ', labels = c("No", "1", '2'), values=c('tomato', 'forestgreen', 'blue')) + 
  theme_minimal() + 
  ylab('Distance from RNAP to closest TDR, bp') +
  xlab('Length of the longest intergenic region, bp') +
  theme(text = element_text(size=14))




a <- a + clean_theme()
c <- c + clean_theme()
ggarrange(a, d, NULL, c, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(1, 2), heights = c(2, 1),
          common.legend = TRUE)  

ggsave('pics/tdr_plus_intergenic_diagnostics.png', width=14, height = 10, dpi=400)
ggsave('pics/tdr_plus_intergenic_diagnostics.pdf', width=14, height = 10, dpi=400)
write.table(joined, 'results/joined.tsv', sep='\t', row.names = FALSE, quote = FALSE)


# BUILD PICTURE

library(gggenomes)

seqs <- joined %>% filter(tdr_left_start < 1000) %>% 
          filter(TDR_in != 0) %>% 
          filter(u_length_ < 10001) %>% select(chrm)

s0 <- read_seqs("blasted/phages_genomes_concat.fna")  # sequence
g0 <- read_feats("upstream_search/representative_genomes.gff")  # proteins 
f0 <- read.table("minimap2_out/TDRs_all.tsv", sep='\t', header=TRUE)  # domains, format: `chr(phage)	start	end	CDS_id	domain`


f0 <- bind_rows(select(f0, seq_id=seq_id, start=start, end=end, de),
                select(f0, seq_id=seq_id, start=start2, end=end2, de))

sequences <- s0 %>% filter(seq_id %in% seqs$chrm)
flips__ <- as.vector(unique(g0[g0$seq_id %in% sequences$seq_id, ] %>% 
                            filter(strand == '-') %>% select(seq_id)))$seq_id

g_stealed <- read_feats('metadata/not_all_genomic.gff')

g_stealed <- g_stealed %>% filter(seq_id %in% seqs$chrm) %>% filter(`type` == 'CDS')
sequences[sequences$seq_id %in% g_stealed$seq_id, ]$file_id <- 'not_all_genomic'
seqs_from_our <- seqs %>% filter(!chrm %in% g_stealed$seq_id)
genes_our <- g0 %>% filter(seq_id %in% seqs_from_our$chrm)
genes_ <- merge(g_stealed, genes_our, all.x = T, all.y = T)




library(patchwork)  # arrange multiple plots
library(ggtree)
library(treeio)
rnaps_tree <- read.tree("upstream_search/polymerases.mafft.treefile")

rnaps_tree_with_drops <- drop.tip(rnaps_tree, c( 'NC_025450.1',
 'OP413830.1',
 'OL744218.1',
 'OL744215.1',
 'NC_027330.1',
 'NC_048171.1',
 'MW015080.1',
 'NC_009531.1',
 'ON457556.1',
 'MT197176.1',
 'NC_019528.1',
 'NC_020483.1',
 'OX001577.1',
 'MF574151.1',
 'NC_048201.1',
 'MT711890.1',
 'NC_047819.1',
 'NC_028822.1',
 'NC_028863.1',
 'LT961846.1',
 'NC_047935.1',
 'MN101216.1',
 'OL799328.1',
 'ON602755.1',
 'OQ338185.1'
))

write.table(genes_, 'genes_reload.tsv', sep='\t', row.names = FALSE, quote = FALSE)
write.table(sequences, 'sequences_reload.tsv', sep='\t', row.names = FALSE, quote = FALSE)
write.table(f0, 'feats_reload.tsv', sep='\t', row.names = FALSE, quote = FALSE)

t <- ggtree(rnaps_tree_with_drops) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))

p <- gggenomes(seqs=sequences, genes=g0, feats = f0)  %>%
  flip(flips__) + 
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

ggsave('pics/have_tdrs_with_tree_my_ann.pdf', d, dpi=600,  height = 0.9 * nrow(sequences), width = 40, limitsize = FALSE)


sequences %>% filter(seq_id %in% rnaps_tree_with_drops[["tip.label"]]) %>% select(seq_id, seq_desc)
ggtree(rnaps_tree_with_drops, layout='circular', branch.length = 'none') + geom_tiplab(align=T, size=5) +
  geom_text2(aes(subset=!isTip, label=node), size = 3, col = "tomato")
ggsave('pics/tree.pdf', dpi=600,  height = 15, width = 15)
library(ape)

try <- extract.clade(rnaps_tree_with_drops, 182)
# try <- drop.tip(rnaps_tree_with_drops, nodepath(rnaps_tree_with_drops, 107, 177))
ggtree(try,branch.length = 'none', layout='circular') + geom_tiplab(align=T, size=5) +
  geom_text2(aes(subset=!isTip, label=node), size = 3, col = "tomato")


t <- ggtree(try) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))
d <- t + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5))
ggsave('pics/have_tdrs_with_tree_182.pdf', d, dpi=600,  height = 0.5 * nrow(sequences), width = 50, limitsize = FALSE)
