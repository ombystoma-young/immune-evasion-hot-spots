library(tidyverse)
library(ggplot2)
library(ggpubr)

args<-commandArgs(TRUE)

in_df_path <- args[1]
figure_name <- args[2]
drop_file_name <- args[3]

# read table
df <- read.csv(in_df_path, sep='\t', na.strings = "None",
               col.names = c('Strain', 'AssemblyID', 'GC_content', 'Length'))

# quantile to drop outliers, length
qs <- quantile(df$Length, probs=c(0.005, 0.995), na.rm=T)

# specify things to filter
table_gc <- df %>% filter(is.na(GC_content))
table_length <- df %>% filter(Length < qs[1])
df_drop <- unique(rbind(table_gc,  table_length))

# tables of outliers
table.gc <- ggtexttable(table_gc, rows = NULL,
            theme = ttheme(colnames.style = colnames_style(color = "black", fill = "gray95"),
                           tbody.style = tbody_style(color = "black", fill = 'white')))

table.length <- ggtexttable(table_length, rows = NULL,
                        theme = ttheme(colnames.style = colnames_style(color = "black", fill = "gray95"),
                                       tbody.style = tbody_style(color = "black", fill = 'white')))

# hists of length and gc-content
hist_length <- df %>% ggplot(aes(Length)) +
        geom_histogram(aes(y=..count../sum(..count..)), color='black', fill='white') +
        geom_vline(xintercept=qs[1], color='tomato') +
        ylab('Density') +
        theme_minimal()

hist_gc <- df %>% ggplot(aes(GC_content)) +
  geom_histogram(aes(y=..count../sum(..count..)), color='black', fill='white') +
  ylab('Density') +
  theme_minimal()

# arrange plots
p <- ggarrange(hist_length, hist_gc, ncol = 2)
g <- ggarrange(table.length, table.gc, ncol=2)
ggarrange(p,g, ncol=1, heights = c(3,1))

# save
ggsave(figure_name, dpi=300, width=15, height=8)
write.csv(df_drop, drop_file_name, sep='\t')