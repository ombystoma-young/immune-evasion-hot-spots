library(tidyverse)
library(ggpubr)
library(RColorBrewer)

setwd('work_dir/anti_defence/anti_defence_pipeline/')

df <- read.delim('data_autographiviridae_refseq/clans_info/res_table.tsv',
                 sep='\t')

df <- df %>%  group_by(new_clan) %>% 
        summarize(N = sum(num_reprs),
                  prelims = paste0(prelim_info, collapse = " "),
                  annotations = paste0(annotations, collapse = ' '))
summary_df <- df %>% filter(N > 7
                            ) %>% arrange(-N) 
summary_df <- cbind(summary_df, ord = 1:nrow(summary_df))
summary_df %>% 
        ggplot(aes(x=ord, y=N)) +
        geom_bar(stat='identity') +
    # scale_fill_viridis_d(na.value='#666666') +
  scale_x_continuous(labels = paste(summary_df$new_clan, summary_df$prelims, sep=':'), 
                     breaks=summary_df$ord) +
  theme(axis.text.x = element_text(angle = 315, vjust = 0.3, hjust=0.05),
        axis.title.x = element_blank(),
        legend.position = 'none')

# -----
df_na <- df %>% dplyr::select(clupar, counts, Product, In.place.of) %>% 
                filter(is.na(Product)) %>% 
                mutate(Product = 'Hypo') %>%
                select(Product, In.place.of, counts)
colnames(df_na) <- c( "Product", "In.place.of", "Counts")
df_no_na <- df %>% select(clupar, counts, Product, In.place.of) %>% 
                   filter(!is.na(Product))
c <- df_no_na %>% dplyr::group_by(Product, In.place.of) %>%
            summarise(Counts = sum(counts)) 
summary_df <- rbind(c, df_na) %>% 
  arrange(-Counts) %>% filter(Counts > 10)
summary_df <- cbind(summary_df, ord = 1:nrow(summary_df))



summary_df <- summary_df %>% 
            mutate(In.place.of = ifelse(In.place.of == 'Ocr',  '\u2666', In.place.of)) %>% 
            mutate(In.place.of = ifelse(In.place.of == '0.7-Kinase',  '\u2660', In.place.of)) %>% 
            mutate(In.place.of = ifelse(In.place.of == 'T7gp0.4',  '\u2663', In.place.of)) %>% 
            mutate(In.place.of = ifelse(In.place.of == 'T7gp0.5',  '\u2665', In.place.of))

summary_df %>% ggplot(aes(y=Counts, x=ord, fill=Product)) +
  geom_bar(stat='identity') + theme_minimal() +
  scale_fill_viridis_d(na.value='#666666') +
  scale_x_continuous(labels = summary_df$Product, 
                  breaks=summary_df$ord) +
  theme(axis.text.x = element_text(angle = 315, vjust = 0.3, hjust=0.05),
        axis.title.x = element_blank(),
        legend.position = 'none')


ggbarplot(summary_df, x = "ord", y = "Counts",
           fill = 'Product',
          # palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
           sorting = "descending",                        # Sort value in descending order
          ggtheme = theme_pubr(),                        # ggplot2 theme
          label = summary_df$In.place.of,
          font.label = list(size = 15, 
                            vjust = -1),
          dot.size = 5,                                 # Large dot size
          x.text.fill = TRUE
) +
  #theme_cleveland()
  scale_x_continuous(labels = summary_df$Product, 
                    breaks=summary_df$ord) +
  theme(axis.text.x = element_text(angle = 315, vjust = 0.3, hjust=0.05),
        axis.title.x = element_blank(),
        legend.position = 'none')
ggsave('pics/barplot_clusters.pdf', dpi=300, width=15, height = 8)
