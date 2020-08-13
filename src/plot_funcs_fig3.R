#!/usr/bin/env Rscript

## Function to subset phylogenetic tree from larger tree
prune.tree <- function(larger_tree,metadata,outfile){
  
  # read in tree
  tree <- read.tree(larger_tree)
  
  # get tip names from metadata file
  n <- read.table(metadata, sep="\t", header=TRUE)
  n <- n %>% filter(paper=="yes")
  n$newick_label <- substr(n$strain,9,100)
  
  # get list of tree taxa to keep
  tips_to_keep <- n$newick_label
  
  # prune tree to keep only desired tips
  pruned.tree <- drop.tip(tree, setdiff(tree$tip.label, tips_to_keep))
  write.tree(pruned.tree,outfile)
  
}

#larger_tree <- "data/tree-na-20200609.nwk"
#md_tree <- "data/tree-extracted-20200609.nwk"
#all_metadata <- "data/tree-na-20200609-metadata.tsv"
#md_metadata <- "data/paper_metadata.txt"
#clinical_metadata <- "data/JHHS_sequenced-v2.txt"


## Function to make tanglegram-like tree figure
make.tanglegram.trees <- function(larger_tree,md_tree,all_metadata,md_metadata,clinical_metadata){
  
  # get full tree
  na_tree <- read.tree(larger_tree)
  
  # get subset tree
  md_tree <- read.tree(md_tree)
  
  # get metadata for full tree
  na_meta <- read.table(all_metadata,sep="\t",header=TRUE,fill=TRUE)
  na_meta <- na_meta %>% select(Strain,Clade)
  
  # get metadata for MD tree
  md_meta <- read.table(md_metadata, sep="\t", header=TRUE)
  md_meta <- md_meta %>% filter(paper=="yes")
  md_meta$newick_label <- substr(md_meta$strain,9,100) # extract tip label string
  md_meta <- md_meta %>% select(newick_label,nextstrain_clade,id)
  
  # get lists of strains corresponding to each clade
  clade_19a <- as.vector(na_meta[na_meta$Clade=="19A",]$Strain)
  clade_19b <- as.vector(na_meta[na_meta$Clade=="19B",]$Strain)
  clade_20a <- as.vector(na_meta[na_meta$Clade=="20A",]$Strain)
  clade_20b <- as.vector(na_meta[na_meta$Clade=="20B",]$Strain)
  clade_20c <- as.vector(na_meta[na_meta$Clade=="20C",]$Strain)
  
  # get node numbers for each clade
  node_19a <- ape::getMRCA(na_tree,tip=clade_19a)
  node_19b <- ape::getMRCA(na_tree,tip=clade_19b)
  node_20a <- ape::getMRCA(na_tree,tip=clade_20a)
  node_20b <- ape::getMRCA(na_tree,tip=clade_20b)
  node_20c <- ape::getMRCA(na_tree,tip=clade_20c)
  
  # get node numbers for all tips in paper
  md_meta$node <- NA
  tipnode <- seq_along(na_tree$tip.label)
  names(tipnode) <- na_tree$tip.label
  md_meta$node <- tipnode[md_meta$newick_label] ## convert the tip label to tip node number
  
  # group clades in tree
  grouped_na_tree <- groupClade(na_tree, .node=c(node_19a,node_19b,node_20a,node_20b,node_20c))
  
  # plot tree with branches colored by clade
  p_na <- ggtree(grouped_na_tree,aes(color=group)) +
    scale_color_manual(values=c("#7fbf7b","#1b7837","#e7d4e8","#af8dc3","#762a83")) +
    geom_point2(aes(subset=(node %in% md_meta$node)),shape=21,fill="white",size=1.5)
  
  ## move to plotting pruned tree of MD samples only
  
  # get clinical metadata
  s <- read.table(clinical_metadata, sep="\t", header=TRUE)
  s$newick_label <- substr(s$strain,9,100)
  s_paper <- s %>% filter(ID %in% md_meta$id)
  s_severity <- s_paper %>% select(ID,newick_label,ANYADMISSION,OXYGEN,ICU,VENT)
  
  # get clade name into this data frame
  md_clade <- md_meta %>% select(ID=id,nextstrain_clade)
  md_clade$ID <- as.numeric(as.matrix(md_clade$ID))
  s_severity <- left_join(s_severity,md_clade,by="ID")
  
  # turn severity into a score
  #s_severity <- s_severity %>% mutate(ANYADMISSION=ifelse(ANYADMISSION=="Yes",1,0)) %>%
  #  mutate(OXYGEN=ifelse(OXYGEN=="Yes",1,0)) %>% mutate(ICU=ifelse(ICU=="Yes",1,0)) %>%
  #  mutate(VENT=ifelse(VENT=="Yes",1,0)) %>% mutate(score=ANYADMISSION+OXYGEN+ICU+VENT)
  
  # turn severity into a size
  size_sm <- 2
  size_lg <- 3
  s_severity <- s_severity %>% mutate(sev=ifelse(VENT=="Yes",size_lg,
                                                 ifelse(ICU=="Yes",size_lg,
                                                        ifelse(OXYGEN=="Yes",size_sm,
                                                               ifelse(ANYADMISSION=="Yes",size_sm,NA)))))
  
  # get node numbers for this tree
  s_severity$node <- NA
  tipnode <- seq_along(md_tree$tip.label)
  names(tipnode) <- md_tree$tip.label
  s_severity$node <- tipnode[s_severity$newick_label] ## convert the tip label to tip node number
  
  # drop the ID for ggtree plotting
  s_severity <- s_severity %>% select(-ID)
  
  # group clades in tree
  grouped_md_tree <- groupClade(md_tree, .node=c(154,136,122,116))
  p_md <- ggtree(grouped_md_tree,aes(color=group)) %<+% s_severity +
    scale_color_manual(values=c("#e7d4e8","#762a83","#af8dc3","#1b7837","#7fbf7b")) +
    geom_point2(aes(subset=(node %in% s_severity[is.na(s_severity$sev),"node"])),
                shape=21,fill="white",size=1.25) +
    geom_tippoint(aes(size=sev)) +
    #geom_tippoint(aes(size=sev,fill=group),shape=21,color="black",stroke=0.8) +
    #scale_fill_manual(values=c("#e7d4e8","#762a83","#af8dc3","#1b7837","#7fbf7b")) +
    scale_size_identity() +
    scale_x_reverse()
  
  # put the plots together into a tanglegram-like object
  grid.arrange(p_na,p_md,ncol=2)
  
}
