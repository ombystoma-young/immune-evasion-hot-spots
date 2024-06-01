library(tidyverse)
library(gggenomes)
library(patchwork)  # arrange multiple plots
library(ggtree)
library(treeio)
library(viridis)


setwd('work_dir/anti_defence/anti_defence_pipeline/')
s0 <- read_seqs("data_autographiviridae/genomes/concatenated_genomes.fna")  # sequence
g0 <- read_feats("data_autographiviridae/target_loci/target_with_clusters.gff")  # proteins 
f0 <- read.table("data_autographiviridae_refseq/tdrs/best_tdrs.tsv", sep='\t', header=FALSE)  

f0 <- bind_rows(select(f0, seq_id=V1, start=V2, end=V3),
                select(f0, seq_id=V1, start=V6, end=V6))
# curated <- read.table('data_autographiviridae/loci_similarity/community_representatives.tsv', sep='\t',
#                       header=TRUE) %>% filter(!str_detect(structure, ':')) %>% 
      # filter(str_detect(structure, 'community'))
curated <- read.table('thesis_tables/plot_list.tsv', sep='\t',
                      header=TRUE) %>% pull(representative.with.known)

not_in_curated <- s0 %>% filter(! seq_id %in% curated)




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

# PLOT REPRESENTATIVES


plot_it <- function(interest){
  s <- s0 %>% filter(seq_id %in% interest$max_k)
  s <- s %>% left_join(interest, by=c('seq_id' = 'max_k'))
  g <- g0 %>% filter(seq_id %in% interest$max_k)
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
      geom_seq_label(aes(label=structure), size=3, nudge_y = -0.22) + # label each sequence by this caption
      geom_gene(aes(fill=`clan`), intron_shape=0, size=5) +  # add gene arrows
      # geom_gene(aes(fill=`freq`), intron_shape=0, size=8) +  # add gene arrows
      #geom_gene(aes(fill=cluster_num), intron_shape=0, size=8) +  # add gene arrows
      geom_gene_text(aes(label=`clu`), size=3) +  # add gene cluster text
      #geom_feat(color='red') +  # add TDRs
      # scale_colour_distiller(palette="Greys") +  # change color
      theme(legend.position = 'none')  # change font size and legend position
  } else {
    p <- gggenomes(seqs=s, genes=g)  %>%
      gggenomes::focus(.track_id = genes, .expand = 100, .overhang='keep') + # focus on particular region +- 100 bp
      geom_seq() +  # draw contig/chromosome lines
      geom_seq_label(aes(label=seq_id), size=3, nudge_y = -0.22) + # label each sequence by this caption
      geom_gene(aes(fill=`clan`), intron_shape=0, size=5) +  # add gene arrows
      # geom_gene(aes(fill=`freq`), intron_shape=0, size=8) +  # add gene arrows
      #geom_gene(aes(fill=cluster_num), intron_shape=0, size=8) +  # add gene arrows
      geom_gene_text(aes(label=`clu`), size=5) +  # add gene cluster text
      #geom_feat(color='red') +  # add TDRs
      scale_colour_distiller(palette="Greys") +  # change color
      # scale_fill_discrete("Clan") +  # change fill, genes
      theme(legend.position = 'none') 
  }
  return(p) 
}

plot_it(curated)
ggsave('data_autographiviridae/pics/community_representatives.pdf', width=25, height=60, limitsize = FALSE)



# METAGENOMES 
s0 <- read_seqs("data_autographiviridae_meta/filtered_mags_flatten/concatenated_genomes.fna")  # sequence
g0 <- read_feats("data_autographiviridae_meta/upstreams/early_with_clusters_phrogs.gff")  # proteins 
# f0 <- read.table("data_autographiviridae_meta/tdrs/best_tdrs.tsv", sep='\t', header=FALSE)  

f0 <- bind_rows(select(f0, seq_id=V1, start=V2, end=V3),
                select(f0, seq_id=V1, start=V6, end=V6))
curated <- read.table('metadata/filtered_upstreams_nuccore.id')

g0 <- g0 %>% mutate(seq_id = str_replace_all(seq_id, '__', '_'))
s0 <- s0 %>% mutate(seq_id = str_replace_all(seq_id, '__', '_'))
f0 <- f0 %>% mutate(seq_id = str_replace_all(seq_id, '__', '_'))
curated <- curated %>% mutate(V1 = str_replace_all(V1, '__', '_'))

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


comm <- c(
  'Han_2018_ERR1398227_NODE_92_length_44625_cov_34.029616',
  'MGV-GENOME-0188590',
  'MGV-GENOME-0286699',
  'MGV-GENOME-0286724',
  'MGV-GENOME-0287854',
  'MGV-GENOME-0290208',
  'MGV-GENOME-0291317',
  'MGV-GENOME-0291506',
  'MGV-GENOME-0292351',
  'MGV-GENOME-0292584',
  'MGV-GENOME-0293792',
  'MGV-GENOME-0294219',
  'MGV-GENOME-0294808',
  'MGV-GENOME-0297924',
  'MGV-GENOME-0298320',
  'MGV-GENOME-0299782',
  'MGV-GENOME-0300114',
  'Ma_2019_SRR341661_NODE_362_length_45429_cov_48.701701',
  'Ma_2019_SRR341714_NODE_224_length_46619_cov_20.739262',
  'Rampelli_10023_NODE_30_length_45460_cov_4.603634',
  'SAMEA3951762_b1_ct5@circular@Caudovirales__sp._ctaq55',
  'SAMN00715211_b1_ct7@circular@Caudovirales__sp._ctjYH7',
  'SAMN00715264_a1_ct14598@circular@Podoviridae__sp._ctG0N598',
  'SAMN00990286_a1_ct8424@circular@Podoviridae__sp._ctAQR424',
  'SAMN03418267_a1_ct52266_vs1@Podoviridae__sp._ctjbH1',
  'SAMN05826883_a1_ct3263@circular@Podoviridae__sp._ctDQb263'
          )

rnaps_tree <- read.tree("data_autographiviridae_meta/known_adgs_analysis/trees/rnap_bootstrap_model_selection.iqtree.contree")
drops <- rnaps_tree$tip.label[!(rnaps_tree$tip.label %in% comm)]
rnaps_tree <- drop.tip(rnaps_tree, drops)
t <- ggtree(rnaps_tree) + geom_tiplab(align=T, size=8) +
  xlim(0,7) + scale_y_continuous(expand=c(0.01, 0.7, 0.01, 0.7))

ggsave( 'test_comm_6.pdf',  t + plot_it(comm) %>% pick_by_tree(t) + plot_layout(widths = c(1,5)) , height = 30, width=20)


ggsave( 'test_comm_3.pdf',  plot_it(comm), height = 23, width=17)


comm7 <- c('ERS743404|NODE_32_length_28217_cov_9.871387',
           'ERS795738|NODE_4_length_10144_cov_3.907058',
           'ERS795790|NODE_3_length_37118_cov_27.174421',
           'Han_2018_ERR1398246_NODE_859_length_15055_cov_4.942867',
           'MGV-GENOME-0189393',
           'MGV-GENOME-0207384',
           'MGV-GENOME-0262495',
           'MGV-GENOME-0263483',
           'MGV-GENOME-0264587',
           'MGV-GENOME-0264710',
           'MGV-GENOME-0264970',
           'MGV-GENOME-0265351',
           'MGV-GENOME-0265429',
           'MGV-GENOME-0268195',
           'MGV-GENOME-0268676',
           'MGV-GENOME-0269006',
           'MGV-GENOME-0269345',
           'MGV-GENOME-0269384',
           'MGV-GENOME-0269628',
           'MGV-GENOME-0270012',
           'MGV-GENOME-0270216',
           'MGV-GENOME-0270879',
           'MGV-GENOME-0270994',
           'MGV-GENOME-0271102',
           'MGV-GENOME-0271146',
           'MGV-GENOME-0271735',
           'MGV-GENOME-0271769',
           'MGV-GENOME-0271781',
           'MGV-GENOME-0272220',
           'MGV-GENOME-0272248',
           'MGV-GENOME-0272520',
           'MGV-GENOME-0272541',
           'MGV-GENOME-0272651',
           'MGV-GENOME-0272845',
           'MGV-GENOME-0272856',
           'MGV-GENOME-0273145',
           'MGV-GENOME-0273329',
           'MGV-GENOME-0274160',
           'MGV-GENOME-0274303',
           'MGV-GENOME-0274788',
           'MGV-GENOME-0275033',
           'MGV-GENOME-0275239',
           'MGV-GENOME-0276411',
           'MGV-GENOME-0276608',
           'Ma_2019_SRR341678_NODE_317_length_41243_cov_67.126032',
           'Rampelli_10017_NODE_250_length_40110_cov_15.836350',
           'Reyes_2015_F112_T1_NODE_5_length_48935_cov_36.137293',
           'SAMEA2737787_a1_ct23066@circular@Podoviridae__sp._ctBmy066',
           'SAMEA2737813_a1_ct27648@circular@Podoviridae__sp._ctydc648',
           'SAMEA2737831_b1_ct14_vs2@Caudovirales__sp._ctco714@linear',
           'SAMEA2737832_a1_ct25830@circular@Podoviridae__sp._ctA9Z830',
           'SAMEA3951678_b1_ct5@circular@Caudovirales__sp._ctiuh5',
           'SAMEA3951737_b1_ct4@circular@Caudovirales__sp._ctUc94',
           'SAMN00715263_a1_ct11339@circular@Podoviridae__sp._ctbS8339',
           'SAMN00791982_a1_ct9829@circular@Podoviridae__sp._cto1I829',
           'SAMN00792000_a1_ct10192@circular@Podoviridae__sp._ctL99192',
           'SAMN00792106_b1_ct73_vs02@Caudovirales__sp._ctsFK73@linear',
           'SAMN00792109_a1_ct16482@circular@Podoviridae__sp._cteXl482',
           'SAMN03418260_a1_ct40487@circular@Podoviridae__sp._ctrUk487',
           'SAMN10080877_a1_ct138@circular@Podoviridae__sp._ctxQr138',
           'SAMN10080901_a1_ct384@circular@Podoviridae__sp._ctOhC384',
           'SAMN10080901_a1_ct415@circular@Caudovirales__sp._cty8v415',
           'SAMN10080902_a1_ct31163_vs1@linear@Podoviridae__sp._ctzUX1',
           'SAMN10080912_a1_ct20335_vs1@linear@Podoviridae__sp._ctFfL1',
           'SRR7892430_NODE_27_length_11706_cov_57.387005',
           'SRR7892431_NODE_5_length_39514_cov_32.711067',
           'SRR7892456_NODE_33_length_7619_cov_29.800502',
           'Zuo_2017_SRR5677784_NODE_194_length_16760_cov_1674.172344')

ggsave( 'test_comm_7.pdf',  plot_it(comm7), height = 39, width=20)



comm18 <- c('MGV-GENOME-0150089',
            'MGV-GENOME-0150506',
            'MGV-GENOME-0199759',
            'MGV-GENOME-0224120',
            'MGV-GENOME-0243489',
            'MGV-GENOME-0253548',
            'MGV-GENOME-0257372',
            'MGV-GENOME-0267196',
            'MGV-GENOME-0271295',
            'MGV-GENOME-0275042',
            'MGV-GENOME-0276243',
            'MGV-GENOME-0277113',
            'MGV-GENOME-0282854',
            'MGV-GENOME-0283569',
            'MGV-GENOME-0287507',
            'MGV-GENOME-0295099',
            'MGV-GENOME-0296543',
            'MGV-GENOME-0296637',
            'MGV-GENOME-0299651',
            'MGV-GENOME-0322767',
            'SAMEA2580140_a1_ct46243_vs1@Podoviridae__sp._ctocX1',
            'SAMEA2737752_a1_ct13252_vs1@Podoviridae__sp._ctCJG1',
            'SAMN01915207_a1_ct21160@circular@Podoviridae__sp._ct3uM160',
            'SAMN05414963_a1_ct3461_vs1@Podoviridae__sp._ctuC31',
            'SAMN05826896_a1_ct11019_vs1@linear@Podoviridae__sp._ctmni1',
            'SAMN05827183_a1_ct4984@circular@Podoviridae__sp._ctWx7984',
            'SAMN10080928_a1_ct35515_vs1@linear@Podoviridae__sp._ctoh61',
            'SRS142883_a1_ct29027_vs1@Podoviridae__sp._ct1Q81',
            'SRS144000_a1_ct16781_vs1@Podoviridae__sp._ct5wZ1',
            'SRS2563908|NODE_15_length_37817_cov_14.857502',
            'SRS2563924|NODE_20_length_36675_cov_18.731486')

ggsave( 'test_comm_18.pdf',  plot_it(comm18), height = 25, width=15)


comm26 <- c('SAMN04262467_a1_ct4897_vs1@Podoviridae__sp._ctJcq1',
            'SRS011115_a1_ct154849_vs1@Podoviridae__sp._ct78r1',
            'SRS014404_a1_ct9375@circular@Podoviridae__sp._ctpUx375',
            'SRS014692_a1_ct5730_vs1@Podoviridae__sp._ctFTG1',
            'SRS014728_a1_ct18994_vs1@Podoviridae__sp._ct1y31',
            'SRS015644_a1_ct23904_vs1@Podoviridae__sp._ctG5D1',
            'SRS016569_a1_ct17846_vs1@Podoviridae__sp._cthux1',
            'SRS019027_a1_ct17128_vs1@Podoviridae__sp._ctOK81',
            'SRS019028_a1_ct91283_vs1@Podoviridae__sp._ctYqZ1',
            'SRS022464_a1_ct10786_vs1@Podoviridae__sp._ctv4K1',
            'SRS053603_a1_ct34958_vs1@Podoviridae__sp._ct0QC1',
            'SRS077508_a1_ct67868_vs1@Podoviridae__sp._ctZHi1',
            'SRS102808_a1_ct28624_vs1@Podoviridae__sp._ctVnB1',
            'SRS104830_a1_ct29152_vs1@Podoviridae__sp._ctq411',
            'SRS143678_a1_ct19922_vs1@Podoviridae__sp._ctPlt1')
ggsave( 'test_comm_26.pdf',  plot_it(comm26), height = 15, width=10)
