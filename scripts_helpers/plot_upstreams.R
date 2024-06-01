
library(tidyverse)
library(gggenomes)
library(patchwork)  # arrange multiple plots
library(ggtree)
library(treeio)
library(viridis)
setwd('work_dir/anti_defence/anti_defence_pipeline/')


community <- '93'
gff_path <- paste0("data_autographiviridae/loci_similarity/communities/target_with_clusters_phrogs_within_", community, "_community.gff")
curated_path <- paste0('data_autographiviridae/loci_similarity/communities/', community, '_community.ids')
out_path <- paste0('data_autographiviridae/pics/upstreams_similar_community_', community, '.pdf')


s0 <- read_seqs("data_autographiviridae/genomes/concatenated_genomes.fna")  # sequence
g0 <- read_feats(gff_path)  # proteins 

curated <- read.table(curated_path)
has_adg <- read.table('data_autographiviridae/loci_similarity/has_target_adgs_htgs_bool.tsv', header=TRUE,sep='\t',
                      na.strings='False') %>% 
  select(-clans)#%>% 
  # filter(has_adgs=='True') %>% pull(seq_id)
row.names(has_adg) <- has_adg$seq_id
has_adg <- has_adg %>% select(-seq_id, -has_adgs, -has_htgs)
colnames(has_adg) <- str_remove_all(colnames(has_adg), 'has_')
has_adg

not_in_curated <- s0 %>% filter(! seq_id %in% curated$V1)

rnaps_tree <- read.tree("data_autographiviridae/known_adgs_analysis/trees/rnap_fasttree.treefile")
rnaps_tree <- ape::drop.tip(rnaps_tree, not_in_curated$seq_id)

seqs <- g0 %>% pull(seq_id) %>% unique()
s0 <- s0 %>% filter(seq_id %in% seqs)

t <- ggtree(rnaps_tree) + #geom_tiplab(align=T, size=5) +
  # geom_tiplab()+
  # geom_point2(aes(subset=(label %in% has_adg)), shape=23, size=10, fill='red') +
  xlim(0,15) +
  theme_tree2()
t

g <- gheatmap(t, has_adg, offset=0.0, width=0.5, 
         colnames=FALSE, legend_title="Has",
         color='black') +
  scale_x_ggtree() +
  theme(axis.text = element_text(size=14, angle = 290, hjust = 0.15),
        legend.position = 'none')
g

flips <- g0[g0$seq_id %in% seqs, ] %>%   # find sequences to flip (- strand -> + strand)
                            filter(strand == '-') %>% pull(seq_id) %>% unique()

g0 <- g0 %>% mutate(clan_col = ifelse(str_detect(clan, 'mono'), NA, clan)) %>% 
            mutate(phrog_ann = str_replace(phrog_ann, 'hypothetical protein', ''))
s0 <- s0 %>% mutate(seq_desc = str_split_i(seq_desc, ',', 1))
p <- gggenomes(seqs=s0, genes=g0)  %>%
  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='keep') + # focus on particular region +- 100 bp
  geom_seq() +  # draw contig/chromosome lines
  geom_seq_label(aes(label=paste(seq_id, seq_desc, sep=', ')), size=10, nudge_y = -0.22) + # label each sequence by this caption
  geom_gene(aes(fill=`clan_col`), intron_shape=0, size=12) +  # add gene arrows
  geom_gene_text(aes(label=`clu`), size=7) +  # add gene cluster text
  scale_fill_discrete("Main cluster product") +  # change fill, genes
  theme(legend.position = 'bottom', axis.text.x = element_text(size=35))  # change font size and legend position
p
if (length(flips) != 0){
    p <- p %>% gggenomes::flip(flips)
    }   # flip sequences

g <-  g + scale_y_continuous(expand=c(0.000, 0.7, 0.000, 0.75)) 

ggsave(out_path, g + p %>% pick_by_tree(t) + plot_layout(widths = c(1,5)) + 
                              theme(legend.position = 'none'), 
                            dpi=600,  height = 2.2 * nrow(curated), 
                          width = 50, limitsize = FALSE)
