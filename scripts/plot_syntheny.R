library(tidyverse)
library(gggenomes)

g0 <- read_feats("work_dir/anti_defence/anti_defence_pipeline/search_upstream/representative_genomes.gff")
s0 <- read_seqs("work_dir/anti_defence/anti_defence_pipeline/clusterization/phages_genomes_concat_clu.fna")

g1 <- g0 %>% filter(seq_id == "MN101216.1" | seq_id == 'NC_047969.1' | seq_id == 'ON995367.1' | seq_id == 'OX001577.1'| seq_id == 'NC_047980.1') %>% 
              filter(product == "T7 RNA polymerase")

s1 <- s0 %>% filter(seq_id == "MN101216.1" | seq_id == 'NC_047969.1' | 
                      seq_id == 'ON995367.1' | seq_id == 'OX001577.1' | seq_id == 'NC_047980.1')

# l0 <- read_paf("work_dir/crispr_cas_IE_context/genomes_completeness/loci_context.paf")
tirs_paf <- read_paf("work_dir/anti_defence/anti_defence_pipeline/minimap2_out/all_phages.paf") %>%
  filter(seq_id == seq_id2 & start < start2 & map_length > 99 & de < 0.1) %>%  
  filter(seq_id == "MN101216.1" | seq_id == 'NC_047969.1' | seq_id == 'ON995367.1' | seq_id == 'OX001577.1' | seq_id == 'NC_047980.1')


tirs_paf <- bind_rows(select(tirs_paf, seq_id=seq_id, start=start, end=end, de),
  select(tirs_paf, seq_id=seq_id2, start=start2, end=end2, de))

p3 <- gggenomes(genes=g1, seqs=s1, feats=tirs_paf) +
  geom_seq() +
  #geom_seq_label(aes(label=seq_desc), size = 4) +
  geom_seq_label() + 
  geom_feat(color='black', size=150) +
  geom_gene(aes(fill=product), intron_shape=0, size = 5) +
  theme(legend.position = 'none') +
  theme(text = element_text(size = 12), axis.text = element_text(size = 15))
p3
