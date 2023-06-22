# ADD DOMAINS TO TABLE OF CLUSTERS
library(tidyverse)
library(reshape2)
library(stringr)

read_CD <- function(path){
  cdd <- read_table(path)
  cdd <- cdd %>% mutate(feat_id = str_split_i(type, ">", 2)) %>% 
    mutate(feat_id = str_split_i(feat_id, "\\[", 1))  %>% 
    select(feat_id, Short, name)
  return(cdd)
}

setwd('work_dir/anti_defence/anti_defence_pipeline/')

phrog_index <- read.csv('../../Tools/databases/phrogs_v_4/PHROG_index.csv')

clusters <- read.table('results/upstream_proteins_clu_wide_seq_sorted_prokka.tsv', header= FALSE,
                       sep='\t') 
colnames(clusters) <- c('clu_parent', 'count', 'products', 'sequence')

pfam <- read.table('domain_tables/pfam.tsv')
phrog <- read.table('domain_tables/phrog.tsv', header = F)
cdd <- read_CD('domain_tables/cdd.txt')
tigrfam <- read_CD('domain_tables/tigrfam.txt')
conj_upstreams <- read.table('domain_tables/upstream_domains.tsv')

length(unique(pfam$V1))
length(unique(phrog$V2))

# add domain description
phrog <- phrog %>% left_join(phrog_index, by=c('V1'='X.phrog'))

# select features
phrog <- phrog %>% select(V1, V2, Annotation, Category)
pfam <- pfam %>% select(V1, V3, V4)
conj_upstreams <- conj_upstreams %>% select(V1, V3)

# rename columns
colnames(phrog) <- c('phrOG', 'feat_id', 'phrOG_annotation', 'phrOG_category')
colnames(pfam) <- c('feat_id', 'pfam_annotation', 'pfam')
colnames(conj_upstreams) <- c('feat_id', 'conj_pl_domains')
colnames(cdd) <- c('feat_id', 'cdd', 'cdd_name')
colnames(tigrfam) <- c('feat_id', 'tigr', 'tigrfam_name')


# squeeze columns
pfam <- pfam %>% group_by(feat_id) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))

conj_upstreams <- conj_upstreams %>% group_by(feat_id) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))

phrog <- phrog %>% group_by(feat_id) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))

cdd <- cdd %>% group_by(feat_id) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))

tigrfam <- tigrfam %>% group_by(feat_id) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))


clusters <- clusters %>% left_join(phrog, by=c('clu_parent'='feat_id')) %>% 
  left_join(pfam, by=c('clu_parent'='feat_id')) %>% 
  left_join(conj_upstreams, by=c('clu_parent'='feat_id')) %>% 
  left_join(cdd, by=c('clu_parent'='feat_id')) %>%                     
  left_join(tigrfam, by=c('clu_parent'='feat_id')) 

write.table(clusters, 'results/upstream_proteins_clu_wide_seq_sorted_prokka_domains.tsv', quote = F,
              sep='\t', col.names = T, row.names = T)
