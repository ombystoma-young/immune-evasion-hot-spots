library(tidyverse)
library(data.table)
library(gggenomes)
setwd('work_dir/anti_defence/anti_defence_pipeline')

datasets <- paste0('dataset_', 1:5)
s0 <- read_seqs("blasted/phages_genomes_concat.fna")  # sequence

remove_no_annotated_rnap <- read.table('metadata/genomes_not_found', col.names = 'genome') %>% filter(genome != 'OM716005.1')
soft_filter <- read.table('metadata/genomes_to_drop.tsv', skip = 1, sep=",")

s0 <- s0 %>% filter(!seq_id %in% remove_no_annotated_rnap$genome) %>% 
      filter(!seq_desc %in% soft_filter$V2)

read_datasets <- function(dataset){
  path <- file.path('define_datasets', paste0(dataset, '_genomes_modified.txt'), fsep = .Platform$file.sep)
  df <- read.table(path, sep='\t', col.names = 'genome')
  df$dataset <- dataset
  colnames(df) <- c('genome', 'dataset')
  return(df)
}

genomes_datasets <- lapply(datasets, read_datasets)
genomes_datasets <- rbindlist(genomes_datasets)
genomes_datasets <- unique(genomes_datasets, by='genome')
genomes_datasets <- genomes_datasets %>% add_row(genome = "OM716005.1", dataset = "dataset_3")
table(genomes_datasets$dataset)
nrow(genomes_datasets) +  s0 %>% filter(! seq_id %in% genomes_datasets$genome) %>%
 select(seq_id) %>% mutate(genome = seq_id) %>% 
 mutate(dataset = 'dataset_1') %>% select(genome, dataset) %>% nrow()

make_drop_dataset <- function(dataset){
  path <- file.path('metadata', paste0(dataset, '.drop'), fsep = .Platform$file.sep)
  df <- read.table(path, sep='\t')
  if (ncol(df) == 1) df$dataset <- dataset
  else df <- df[c(2,1)]
  colnames(df) <- c('genome', 'dataset')
  return(df)
}
genomes_to_drop <- lapply(datasets, make_drop_dataset)
genomes_to_drop <- rbindlist(genomes_to_drop) %>% filter(genome != 'NC_001604.1')

# number of genomes to drop in each dataset
table(genomes_to_drop$dataset)

genomes_after_curation <- genomes_datasets %>% filter(! genome %in% genomes_to_drop$genome)

#number of genomes after curation
table(genomes_after_curation$dataset)

write.table(genomes_after_curation, 
            file = file.path('metadata', 'genomes_after_curation.tsv', fsep = .Platform$file.sep),
            col.names = FALSE,
            row.names = FALSE,
            sep='\t',
            quote = FALSE
            )


# seq_description <- read.table('metadata/dataset_4_good.id', 
#                               sep='\t', header = F, 
#                               col.names = 'descr')

# 
# j_1 <- read.table('define_datasets/dataset_4_genomes_modified.txt')
# 
# a <- s0  %>% filter(seq_id %in% j_1$V1) %>% filter(!seq_desc %in% seq_description$descr)
# 
# write.table(data.frame('dataset_4', a$seq_id),
#           file = 'metadata/dataset_4.drop',
#           sep='\t', row.names = F, quote = F,
#           col.names = F)
# 
# j_1 <- read.table('define_datasets/dataset_5_genomes_modified.txt')
# a <- s0  %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id != 'OP413828.1') %>% 
#             filter(seq_id != 'NC_001604.1')
# write.table(data.frame('dataset_5', a$seq_id),
#             file = 'metadata/dataset_5.drop',
#             sep='\t', row.names = F, quote = F,
#             col.names = F)
# 
# j_1 <- read.table('metadata/drop_before_upstream_search.id')
# nrow(s0  %>% filter(seq_id %in% j_1$V2))
# 
# df <- read.table('define_datasets/joined.tsv', header=T) %>% filter(!chrm %in% j_1$V2)
# df_3 <- df %>% filter(dataset == 'dataset_3')
# 
# df_3 <- df_3 %>% filter(i_start > 0)
# df_3 
# 
# cl.nms <- c('chrm', 'upstream_start', 'upstream_end', 'dot',
#             'pos', 'strand', 'pol_start', 'pol_end',
#             'tdr_left_start', 'tdr_left_end', 
#             'tdr_right_start', 'tdr_right_end')
# 
# upstrms <- read.table('upstream_search/upstream.bed', 
#            header=FALSE,
#            fill = TRUE,
#            col.names = cl.nms) %>% filter(!chrm %in% j_1$V2)
# 
# upstrms <- left_join(upstrms, df_3) %>% filter(dataset == 'dataset_3')













