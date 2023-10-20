library(tidyverse)
library(gggenomes)
library(pheatmap)
setwd('work_dir/anti_defence/anti_defence_pipeline/')

df <- read_feats('results/upstreams_clusters_domains.gff')
df <- df %>% select(seq_id, cluster_num) %>% mutate(cluster_num = as.factor(as.numeric(cluster_num) + 1))
# norm_df <- df %>% group_by(cluster_num) %>% 
#         summarise(Count = n()) %>% arrange(desc(Count), desc)
# df_table <- df %>% pivot_wider(id_cols = 'seq_id', names_from = 'cluster_num', values_from = 'cluster_num')
df_small <- df %>% filter(as.numeric(cluster_num) < 200)
crosstab <- df_small %>%  group_by(seq_id, cluster_num) %>%  tally() %>%  spread(cluster_num, n) %>% replace(., is.na(.), 0) 
crosstab2 <- as.matrix(crosstab[, -1])
pheatmap(crosstab2, cluster_rows = TRUE, cluster_cols = FALSE,
         filename = "pics/test_heatmap_all2.pdf",
         width=30, height = 30)
# 
# 
# pca_res <- prcomp(crosstab2, center = TRUE, scale. = TRUE) 
# names(pca_res)
# summary(pca_res)
# plot(pca_res, type='l')
# biplot(pca_res)
# 
# library(embed)






merged = merge(df, df, by ='seq_id')
crosstab <- merged %>%  group_by(cluster_num.x, cluster_num.y) %>%  tally() %>%  spread(cluster_num.x, n) #%>% replace(., is.na(.), 0) 
crosstab <- as.matrix(crosstab[, -1])

# rownames(crosstab) <-  as.numeric(crosstab$cluster_num.y)
# colnames(crosstab) <-  as.numeric(colnames(crosstab))
# crosstab <- crosstab 
pheatmap(crosstab, cluster_rows = FALSE, cluster_cols = FALSE,  type = "lower",
         labels_col = FALSE)

df %>%  filter(seq_id %in% (df %>% filter(cluster_num=='0') %>% pull(seq_id))) %>%  filter(cluster_num == '31' ) %>% nrow()


