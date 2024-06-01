library(tidyverse)
library(gggenomes)
library(patchwork)  # arrange multiple plots
library(ggtree)
library(treeio)
library(viridis)


setwd('work_dir/anti_defence/anti_defence_pipeline/')
s0 <- read_seqs("data_autographiviridae/genomes/concatenated_genomes.fna")  # sequence
g0 <- read_feats("data_autographiviridae/target_loci/target_with_clusters.gff")  # proteins 

curated <- read.table('thesis_tables/plot_list.tsv', sep='\t',
                      header=TRUE) %>% pull(representative.with.known)

not_in_curated <- s0 %>% filter(! seq_id %in% curated)




plot_it <- function(interest){
              s <- s0 %>% filter(seq_id %in% interest)
              g <- g0 %>% filter(seq_id %in% interest)
              flips <- as.vector(unique(g[g$seq_id %in% s$seq_id, ] %>%   # find sequences to flip (- strand -> + strand)
                                          filter(strand == '-') %>% select(seq_id)))$seq_id
              g <- g %>% group_by(clan) %>% mutate(freq = n()) %>% ungroup()
                 if (length(flips) != 0){
              p <- gggenomes(seqs=s, genes=g)  %>%
                gggenomes::flip(flips) %>%  # flip sequences
                # pick_seqs(interest) %>% 
                gggenomes::focus(.track_id = genes, .expand = 100, .overhang='keep') +
                 # focus on particular region +- 100 bp
                geom_seq() +  # draw contig/chromosome lines
                geom_seq_label(aes(label=paste(seq_id, str_split_i(seq_desc, ',', 1), sep=', ')), size=5, nudge_y = -0.22) + # label each sequence by this caption
                geom_gene(aes(fill=`clan`), intron_shape=0, size=7) +  # add gene arrows
                # geom_gene(aes(fill=`freq`), intron_shape=0, size=8) +  # add gene arrows
                #geom_gene(aes(fill=cluster_num), intron_shape=0, size=8) +  # add gene arrows
                geom_gene_text(aes(label=`clu`), size=5) +  # add gene cluster text
                #geom_feat(color='red') +  # add TDRs
                # scale_colour_distiller(palette="Greys") +  # change color
                theme(legend.position = 'none')  # change font size and legend position
              } else {
                p <- gggenomes(seqs=s, genes=g)  %>%
                  gggenomes::focus(.track_id = genes, .expand = 100, .overhang='keep') + #%>% 
                  # pick_seqs(interest) + # focus on particular region +- 100 bp
                  geom_seq() +  # draw contig/chromosome lines
                  geom_seq_label(aes(label=seq_id), size=5, nudge_y = -0.22) + # label each sequence by this caption
                  geom_gene(aes(fill=`clan`), intron_shape=0, size=6) +  # add gene arrows
                  # geom_gene(aes(fill=`freq`), intron_shape=0, size=8) +  # add gene arrows
                  #geom_gene(aes(fill=cluster_num), intron_shape=0, size=8) +  # add gene arrows
                  geom_gene_text(aes(label=`clu`), size=4) +  # add gene cluster text
                  #geom_feat(color='red') +  # add TDRs
                  scale_colour_distiller(palette="Greys") +  # change color
                  # scale_fill_discrete("Clan") +  # change fill, genes
                  theme(legend.position = 'none') 
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


# in place of kinase noc
interest <- c('MN481365.1',
              'OP413828.1')
plot_it(interest)


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



# SAMase duplcation ?
interest <- c(
  'OK274245.1',
  'OK318991.1',
  'NC_023736.1',
  'OP313111.1',
  'MK290739.2'
)
plot_it(interest)
ggsave('data_autographiviridae/pics/samase_duplication.pdf', width=20, height=6)


interest <- c('NC_047751.1')



interest <- c(
  'AY264778.1',
  'NC_001604.1',
  'MZ375268.1'
)
plot_it(interest)
ggsave('data_autographiviridae/pics/kinases_shut_off_only.pdf', width=20, height=6)

interest <- c(
  'NC_031258.1',
'MZ501064.1',
'NC_047789.1',
'NC_048164.1'
)
plot_it(interest)
ggsave('data_autographiviridae/pics/kinases_shut_off_after_stop.pdf', width=20, height=6)


interest <- c(
  'AY264775.1',
  'MGV-GENOME-0262495',
  'NC_011045.1'
)
plot_it(interest)
ggsave('data_autographiviridae/pics/kinases_shut_off_after_stop_part_of_kinase.pdf', width=20, height=6)


interest <- c(
  'OM457002.1',
  'MGV-GENOME-0268676',
  'MGV-GENOME-0270012',
  'SAMN00792109_a1_ct16482@circular@Podoviridae__sp._cteXl482'
)
plot_it(interest)
ggsave('data_autographiviridae/pics/kinases_only_kinase.pdf', width=20, height=6)



interest <- c('SRS741334|NODE_2_length_43394_cov_11.226632',
              'OK499982.1',
              'NC_028863.1',
              'MGV-GENOME-0290541',
              'MN327636.1')
plot_it(interest)


plot_it(curated)
plot_it(c('Ma_2019_SRR341661_NODE_362_length_45429_cov_48.701701'))
ggsave('data_autographiviridae/pics/target_loci_representatives2.pdf', width=20, height=6/4)
