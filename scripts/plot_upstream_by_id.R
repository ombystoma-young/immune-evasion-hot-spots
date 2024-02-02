library(tidyverse)
library(gggenomes)
library(patchwork)  # arrange multiple plots
library(ggtree)
library(treeio)
library(viridis)


setwd('work_dir/anti_defence/anti_defence_pipeline/')
s0 <- read_seqs("data_autographiviridae_refseq/genomes/concatenated_genomes.fna")  # sequence
g0 <- read_feats("data_autographiviridae_refseq/upstreams/early_with_clusters_phrogs.gff")  # proteins 
f0 <- read.table("data_autographiviridae_refseq/tdrs/best_tdrs.tsv", sep='\t', header=FALSE)  

f0 <- bind_rows(select(f0, seq_id=V1, start=V2, end=V3),
                select(f0, seq_id=V1, start=V6, end=V6))
curated <- read.table('metadata/genomes_after_curation.txt')
not_in_curated <- s0 %>% filter(! seq_id %in% curated$V1)




plot_it <- function(interest){
              s <- s0 %>% filter(seq_id %in% interest)
              g <- g0 %>% filter(seq_id %in% interest)
              flips <- as.vector(unique(g[g$seq_id %in% s$seq_id, ] %>%   # find sequences to flip (- strand -> + strand)
                                          filter(strand == '-') %>% select(seq_id)))$seq_id
              g <- g %>% group_by(clan) %>% mutate(freq = n()) %>% ungroup()
              # 
              # flips__ <- g %>% group_by(seq_id, strand) %>%
              #   summarise(n()) %>%
              #   pivot_wider(names_from = strand, names_prefix = 'strand_', values_from = `n()`) %>%
              #   ungroup() %>%
              #   filter(is.na(`strand_+`) | `strand_-` > `strand_+`) %>%
              #   select(seq_id) %>% as.vector()
              if (length(flips) != 0){      
              p <- gggenomes(seqs=s, genes=g)  %>%
                gggenomes::flip(flips) %>%  # flip sequences
                gggenomes::focus(.track_id = genes, .expand = 100, .overhang='keep') + # focus on particular region +- 100 bp
                geom_seq() +  # draw contig/chromosome lines
                geom_seq_label(aes(label=seq_desc), size=5, nudge_y = -0.22) + # label each sequence by this caption
                # geom_gene(aes(fill=`clan`, color=`freq`), intron_shape=0, size=8) +  # add gene arrows
                geom_gene(aes(fill=`freq`), intron_shape=0, size=8) +  # add gene arrows
                #geom_gene(aes(fill=cluster_num), intron_shape=0, size=8) +  # add gene arrows
                geom_gene_text(aes(label=`clu`), size=5) +  # add gene cluster text
                #geom_feat(color='red') +  # add TDRs
                scale_fill_viridis_c("Clan") +  # change fill, genes
                theme(legend.position = 'bottom')  # change font size and legend position
              } else {
                p <- gggenomes(seqs=s, genes=g)  %>%
                  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='keep') + # focus on particular region +- 100 bp
                  geom_seq() +  # draw contig/chromosome lines
                  geom_seq_label(aes(label=seq_desc), size=5, nudge_y = -0.22) + # label each sequence by this caption
                  geom_gene(aes(fill=`clan`, color=`freq`), intron_shape=0, size=8) +  # add gene arrows
                  # geom_gene(aes(fill=`freq`), intron_shape=0, size=8) +  # add gene arrows
                  #geom_gene(aes(fill=cluster_num), intron_shape=0, size=8) +  # add gene arrows
                  geom_gene_text(aes(label=`clu`), size=5) +  # add gene cluster text
                  #geom_feat(color='red') +  # add TDRs
                  # scale_fill_discrete("Clan") +  # change fill, genes
                  theme(legend.position = 'bottom') 
              }
              return(p) 
}
# kinase 34
interest <- c('MZ234024.1',
              'MT862763.1',
              'ON637250.1')
plot_it(interest)
ggsave('pics/jan_kinase_loci.pdf', width=20, height=6)
# low pI in place of Ocr
interest <- c('NC_047858.1',
              'NC_028916.1')
plot_it(interest)
ggsave('pics/jan_low_pi_loci.pdf', width=20, height=6)
# ardA instead of ocr
interest <- c('MW671054.1', 'MT740748.1', 'MZ851152.1')
plot_it(interest)
ggsave('pics/jan_ardA_loci.pdf', width=20, height=6)



# community 8
interest <- c('MH179478.2',
'MT711887.2',
'NC_015208.1',
'NC_015264.1',
'NC_021062.1',
'NC_027292.1',
'NC_047965.1',
'NC_047997.1',
'NC_048200.1',
'NC_048201.1',
'OM471789.1',
'ON189046.1')
plot_it(interest)


interest <- c(
              # 'NC_047746.1', 
              'NC_047965.1',
              # 'NC_047997.1', 
              # 'MT711887.2',
              'NC_021062.1',
              # 'NC_027292.1'
              )
plot_it(interest)
ggsave('pics/jan_158.pdf', width=20, height=6)
