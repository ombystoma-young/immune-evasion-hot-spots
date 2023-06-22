# ADD PFAM and PHROG results to gff
library(tidyverse)
library(reshape2)
library(gggenomes)
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

up_gff <- read_feats('results/upstreams_with_clusters.gff')

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
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ","))) %>% 
  mutate(phrOG_annotation = str_replace_all(phrOG_annotation, ';', ','))

cdd <- cdd %>% group_by(feat_id) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))

tigrfam <- tigrfam %>% group_by(feat_id) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))


up_gff <- up_gff %>% left_join(phrog) %>% 
                     left_join(pfam) %>% 
                     left_join(conj_upstreams) %>% 
                     left_join(cdd) %>%                     
                     left_join(tigrfam) 

write_gff3(up_gff, 'results/upstreams_clusters_domains.gff')
