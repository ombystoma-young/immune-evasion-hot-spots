if (!requireNamespace("gggenomes", quietly = TRUE))
  devtools::install_github("thackl/gggenomes")

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)
install.packages('latex2exp')

library(tidyverse)
library(gggenomes, verbose = T)
library(dplyr)
library(ggpubr)
library(latex2exp)


args <- commandArgs(trailingOnly=TRUE)

# domains_tsv <- args[1]
# gff <- args[2]
# lengths_bed <- args[3]
# rnaps_gff <- args[4]
# selected <- args[5]
# tdr_figure_name <- args[6]
# intergenic_full_figure_name <- args[7]
# intergenic_conditions_figure_name <- args[8]
setwd('work_dir/anti_defence/anti_defence_pipeline/')
domains_tsv <- 'data_autographiviridae_refseq/domains/tables/apis_whole_with_descr.tsv'
gff <- 'data_autographiviridae_refseq/annotation/concatenated_with_clusters_phrogs.gff'
lengths_bed <- 'data_autographiviridae_refseq/intergenics/chromosome_lengths.bed'
rnaps_gff <- 'data_autographiviridae_refseq/annotation/concatenated_rnaps_only.gff'
selected <- 'metadata/filtered_upstreams_nuccore.id'

## rnaps
rnaps <- read.table(rnaps_gff, col.names = c('seq_id', 'source', 'type', 
                                             'start.rnap', 'end.rnap',
                                             'score', 'strand', 
                                             'phase', 'attribute')) %>% 
          select(seq_id, start.rnap, end.rnap, strand)

## chromosome lengths 
chr_lens <- read.table(lengths_bed, 
                       col.names = c('seq_id', 'start', 'seq_len')) %>% 
  select(-start)

# AntiDefense Genes
adgs <- read.table(domains_tsv, sep='\t', header=TRUE)

# Selection
selected <- read.table(selected, sep='\t', header=FALSE)

# read gff
gff_df <- read_feats(gff) %>% 
      filter(seq_id %in% selected$V1)

# filter rnaps only in selection
rnaps <- rnaps %>% filter(seq_id %in% selected$V1)


# join rnap df and gff df
rnap_join <- rnaps %>% left_join(gff_df, by='seq_id') %>% 
  left_join(chr_lens)

# calculate relative positions
# 1 way: directed 
rnap_join <- rnap_join %>%
  mutate(rel.dist = case_when(
    strand.x == '+' &
      start.rnap > start & start.rnap > end ~
      start.rnap - end,
    strand.x == '+' &
      end.rnap <= start & end.rnap < end ~
      seq_len - end + start.rnap,
    strand.x == '+' &
      end.rnap < start & start.rnap > end ~
      start.rnap - end,
    strand.x == '-' &
      start.rnap > start & start.rnap >= end  ~
      seq_len - end.rnap + start,
    strand.x == '-' &
      end.rnap < start & end.rnap < end  ~
      start - end.rnap,
    strand.x == '-' &
      end.rnap < start & start.rnap > end  ~
      start - end.rnap,
    .default=0))

rnap_join <- rnap_join %>%
  dplyr::group_by(seq_id) %>%
  mutate(rel.rank = rank(rel.dist, na.last='keep') - 1) %>%
  ungroup()
# 2 way: distinguish what is before and what is after

# 
# rnap_join <- rnap_join %>% 
#   mutate(rel.dist = case_when(
#     strand.x == '+' ~ start - start.rnap,
#     strand.x == '-' ~ start.rnap - start, 
#     .default=0))
# 
# rnap_join <- rnap_join %>%
#   dplyr::group_by(seq_id) %>% 
#   mutate(dist.rank = rank(rel.dist, na.last='keep')) %>%
#   ungroup()
# 
# rnap_ranks <- rnap_join %>% 
#   filter(rel.dist == 0) %>% 
#   select(seq_id, dist.rank) 
# colnames(rnap_ranks) <- c('seq_id', 'rnap.rank')
# 
# rnap_join <- 
#   rnap_join %>%
#   left_join(rnap_ranks) %>% 
#   mutate(rel.rank = dist.rank - rnap.rank)

# add information about adgs
rnap_join <- rnap_join %>% left_join(adgs, by=c('geom_id'='target_name')) %>% 
  mutate(is_adg = !is.na(query_name)) %>% 
  mutate(is_against_rm = str_detect(Defense.systems, '(RM)')) %>% 
  mutate(is_against_Thoeris = str_detect(Defense.systems, 'Thoeris')) %>% 
  mutate(is_against_ta = str_detect(Defense.systems, '(TA)')) %>% 
  mutate(is_against_recbcd = str_detect(Defense.systems, 'RecBCD'))

rnap_join[is.na(rnap_join$is_against_rm), ]$is_against_rm <- FALSE
rnap_join[is.na(rnap_join$is_against_Thoeris), ]$is_against_Thoeris <- FALSE
rnap_join[is.na(rnap_join$is_against_ta), ]$is_against_ta <- FALSE
rnap_join[is.na(rnap_join$is_against_recbcd), ]$is_against_recbcd <- FALSE


table(rnap_join$Defense.systems)

# group information about ADGs
group_adgs <- function(keyword, def_sys){
      df <- rnap_join %>% 
                        group_by(rel.rank) %>% 
                        summarise(num_ADGs = sum({{keyword}}),
                                  num_genes = n()) %>% 
                        mutate(`Is against` = paste0(def_sys))

      return(df)
}
grouped_adgs <- rbind(
                      #group_adgs(is_adg, 'all ADGs'), 
                      group_adgs(is_against_rm, 'RM or BREX'),
                      group_adgs(is_against_Thoeris, 'Thoeris'),
                      group_adgs(is_against_ta, 'TA'),
                      group_adgs(is_against_recbcd, 'recBCD'))
grouped_adgs %>% ggplot(aes(x=rel.rank, 
                            y=num_ADGs / num_genes,
                            color = `Is against`)) + 
                  geom_line(linewidth=1.5, alpha=1) + 
                  scale_color_brewer(palette='Dark2') +
                  theme_classic2() + 
                  scale_x_reverse() +
                  xlab('Relative position (ORFs)') +
                  ylab('Fraction of anti-defense') +
                  theme(text = element_text(size=18,
                                            family="sans"),
                        legend.position = 'bottom',
                        panel.background = element_rect(fill = "transparent",
                                                        colour = NA_character_), # necessary to avoid drawing panel outline
                        panel.grid.major = element_blank(), # get rid of major grid
                        panel.grid.minor = element_blank(), # get rid of minor grid
                        plot.background = element_rect(fill = "transparent",
                                                       colour = NA_character_), # necessary to avoid drawing plot outline
                        legend.background = element_rect(fill = "transparent",
                                                         colour = NA_character_),
                        legend.box.background = element_rect(fill = "transparent",
                                                             colour = NA_character_),
                        legend.key = element_rect(fill = "transparent",
                                                  colour = NA_character_))
ggsave('pics/freq_adgs.pdf', dpi=300, height = 6, width=8,  bg = "transparent")                  
  