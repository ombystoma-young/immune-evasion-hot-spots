library(ape)
library(tidyverse)
library(dendextend)
library(phytools)


setwd('work_dir/anti_defence/anti_defence_pipeline/')

tree_rnap_file <- 'data_autographiviridae_refseq/known_proteins/trees/rnap_bootstrap_model_selection.iqtree.contree'
tree_samase_file <- 'data_autographiviridae_refseq/known_proteins/trees/samase_bootstrap_model_selection.iqtree.contree'
tree_ocr_file <- 'data_autographiviridae_refseq/known_proteins/trees/ocr_bootstrap_model_selection.iqtree.contree'
tree_kinase_file <- 'data_autographiviridae_refseq/known_proteins/trees/kinase_bootstrap_model_selection.iqtree.contree'
tree_rnap <- read.tree(tree_rnap_file)
tree_samase <- read.tree(tree_samase_file)

drop_rnap <- tree_rnap$tip.label[!tree_rnap$tip.label %in% tree_samase$tip.label]
tree_rnap <- drop.tip(tree_rnap, drop_rnap)

association <- cbind(tree_rnap$tip.label, tree_rnap$tip.label)


a <- cophylo(tree_rnap, 
             tree_samase, 
             rotate=TRUE,
             tangle='both')
c <- cophylo(tree_rnap, 
             tree_samase, 
             rotate=FALSE,
             tangle='both')

tree_rnap <- read.tree(tree_rnap_file)
tree_ocr <- read.tree(tree_ocr_file)

drop_rnap <- tree_rnap$tip.label[!tree_rnap$tip.label %in% tree_ocr$tip.label]
tree_rnap <- drop.tip(tree_rnap, drop_rnap)

association <- cbind(tree_rnap$tip.label, tree_rnap$tip.label)

b <- cophylo(tree_rnap, 
             tree_ocr, 
             rotate=TRUE,
             tangle='both')
d <- cophylo(tree_rnap, 
             tree_ocr, 
             rotate=FALSE,
             tangle='both')


pdf(width=5, height=8)
plot(c, link.lwd=2, 
     link.lty="dashed", link.col='#E7298A',
     pts=FALSE, ftype="off")
plot(a, link.lwd=2, 
     link.lty="dashed", link.col='#E7298A',
     pts=FALSE, ftype="off")
plot(d, link.lwd=2, 
     link.lty="dashed", link.col='#7570B3',
     pts=FALSE, ftype="off")
plot(b, link.lwd=2, 
     link.lty="dashed", link.col='#7570B3',
     pts=FALSE, ftype="off")
a$trees[[1]]$edge.length<-NULL
a$trees[[2]]$edge.length<-NULL
b$trees[[1]]$edge.length<-NULL
b$trees[[2]]$edge.length<-NULL
plot(a, link.lwd=2, 
     link.lty="dashed", link.col='#E7298A',
     pts=FALSE, ftype="off")
plot(b, link.lwd=2, 
     link.lty="dashed", link.col='#7570B3',
     pts=FALSE, ftype="off")


a1pie <- ifelse(as.numeric(a$trees[[1]]$node.label) > 95, 1, 0)
a2pie <- ifelse(as.numeric(a$trees[[2]]$node.label) > 95, 1, 0)
b1pie <- ifelse(as.numeric(b$trees[[1]]$node.label) > 95, 1, 0)
b2pie <- ifelse(as.numeric(b$trees[[2]]$node.label) > 95, 1, 0)

plot(a, link.lwd=2, 
     link.lty="dashed", link.col='#E7298A',
     pts=FALSE, ftype="off")
nodelabels.cophylo(pie=a1pie,
                   cex=0.15, piecol=c('white', 'grey30'), 
                   frame='none')
nodelabels.cophylo(pie=a2pie,
                   cex=0.15, piecol=c('white', 'grey30'),
                   frame='none', which = 'right')

plot(b, link.lwd=2, 
     link.lty="dashed", link.col='#7570B3',
     pts=FALSE, ftype="off")
nodelabels.cophylo(pie=b1pie,
                   cex=0.15, piecol=c('white', 'grey30'), 
                   frame='none')
nodelabels.cophylo(pie=b2pie,
                   cex=0.15, piecol=c('white', 'grey30'),
                   frame='none', which = 'right')
dev.off()



tree_rnap <- read.tree(tree_rnap_file)
tree_kinase <- read.tree(tree_kinase_file)

drop_rnap <- tree_rnap$tip.label[!tree_rnap$tip.label %in% tree_kinase$tip.label]
tree_rnap <- drop.tip(tree_rnap, drop_rnap)
a <- cophylo(tree_rnap, 
             tree_kinase, 
             rotate=TRUE,
             tangle='both')
a1pie <- ifelse(as.numeric(a$trees[[1]]$node.label) > 95, 1, 0)
a2pie <- ifelse(as.numeric(a$trees[[2]]$node.label) > 95, 1, 0)
plot(a, link.lwd=2, 
     link.lty="dashed", link.col='#E7298A',
     pts=FALSE, ftype="off")
nodelabels.cophylo(pie=a1pie,
                   cex=0.15, piecol=c('white', 'grey30'), 
                   frame='none')
nodelabels.cophylo(pie=a2pie,
                   cex=0.15, piecol=c('white', 'grey30'),
                   frame='none', which = 'right')
