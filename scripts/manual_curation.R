library(tidyverse)

setwd('work_dir/anti_defence/anti_defence_pipeline')


seq_description <- read.table('metadata/dataset_4_good.id', 
                              sep='\t', header = F, 
                              col.names = 'descr')

s0 <- read_seqs("blasted/phages_genomes_concat.fna")  # sequence
j_1 <- read.table('define_datasets/dataset_4_genomes_modified.txt')

a <- s0  %>% filter(seq_id %in% j_1$V1) %>% filter(!seq_desc %in% seq_description$descr)

write.table(data.frame('dataset_4', a$seq_id),
          file = 'metadata/dataset_4.drop',
          sep='\t', row.names = F, quote = F,
          col.names = F)

j_1 <- read.table('define_datasets/dataset_5_genomes_modified.txt')
a <- s0  %>% filter(seq_id %in% j_1$V1) %>% filter(seq_id != 'OP413828.1') %>% 
            filter(seq_id != 'NC_001604.1')
write.table(data.frame('dataset_5', a$seq_id),
            file = 'metadata/dataset_5.drop',
            sep='\t', row.names = F, quote = F,
            col.names = F)

j_1 <- read.table('metadata/drop_before_upstream_search.id')
nrow(s0  %>% filter(seq_id %in% j_1$V2))

df <- read.table('define_datasets/joined.tsv', header=T) %>% filter(!chrm %in% j_1$V2)
df_3 <- df %>% filter(dataset == 'dataset_3')

df_3 <- df_3 %>% filter(i_start > 0)
df_3 

cl.nms <- c('chrm', 'upstream_start', 'upstream_end', 'dot',
            'pos', 'strand', 'pol_start', 'pol_end',
            'tdr_left_start', 'tdr_left_end', 
            'tdr_right_start', 'tdr_right_end')

upstrms <- read.table('upstream_search/upstream.bed', 
           header=FALSE,
           fill = TRUE,
           col.names = cl.nms) %>% filter(!chrm %in% j_1$V2)

upstrms <- left_join(upstrms, df_3) %>% filter(dataset == 'dataset_3')













