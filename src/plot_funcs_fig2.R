#!/usr/bin/env Rscript

# load necessary packages
library(plyr)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(tidyr)
library(ggforce)
library(extrafont)
library(textshape)
library(gplots)
library(ggdendro)
library(ggnewscale)
library(scales)
library(grid)
library(dendextend)

## Function to make histograms stratified by different metadata variables
make_complete_genome_hists <- function(seq_metadata){
  
  ### stratify completeness histograms by different variables ###
  
  ### CT HISTOGRAM ###
  
  # clean up ct values
  # ignore unknowns, nans, and zero values
  # ignore any samples from failed runs
  ct_subset <- seq_metadata %>% filter(!is.na(diagnostic_ct_value)) %>%
    filter(diagnostic_ct_value!="0")

  # split up data into complete/not complete
  ct_df <- ct_subset
  ct_df_complete <- ct_df %>% filter(genome=="Yes")
  ct_df_fail <- ct_df %>% filter(genome!="Yes")
  
  # break data into bins
  hist_all <- hist(ct_df$diagnostic_ct_value, breaks = 7, plot = FALSE)
  hist_complete <- hist(ct_df_complete$diagnostic_ct_value, breaks = hist_all$breaks, plot = FALSE)
  hist_fail <- hist(ct_df_fail$diagnostic_ct_value, breaks = hist_all$breaks, plot = FALSE)
  
  # melt for plotting with ggplot
  full_df <- data.frame(bin=seq(1,length(hist_all$counts)),pass=hist_complete$counts,fail=hist_fail$counts)
  full_df <- reshape2::melt(full_df,id.vars="bin")
  full_df$variable <- factor(full_df$variable, levels = c("fail","pass"))
  
  # plot histogram
  ct_hist <- ggplot(full_df, aes(x=bin, y=value, fill=variable, color=variable)) +
    geom_bar(stat="identity", width=0.85) +
    scale_fill_manual(values=c("white", "black")) +
    scale_color_manual(values=c("black", "black")) +
    labs(x="Ct value", y="Number of samples") +
    #ggtitle("Successful genomes by Ct") +
    scale_x_reverse(breaks=seq(0.5,6.5,2), labels=seq(10,40,10), expand=c(0.05,0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          #plot.title = element_text(hjust = 0.5, size=12),
          legend.title = element_blank(),
          legend.position = "None") +
    geom_hline(yintercept=seq(0,40,10), color="white", size=0.2)
  
  ### AGE HISTOGRAM ###
  
  # stratify by age group
  ct_df <- seq_metadata %>% filter(!is.na(age))
  #ct_df$age <- as.character(as.matrix(ct_df$age))
  #ct_df <- mutate(ct_df, age = ifelse(grepl(">", age, fixed = TRUE),95,age))
  #ct_df$age <- as.numeric(as.matrix(ct_df$age))
  bins <- c("<10","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80-89","90+")
  ct_df$age_bins <- cut(ct_df$age,breaks = c(-1,9,19,29,39,49,59,69,79,89,99),labels = bins,right = TRUE)
  
  ct_df_complete <- ct_df %>% filter(genome=="Yes")
  ct_df_fail <- ct_df %>% filter(genome!="Yes")
  
  full_df <- data.frame(bin=bins,
                        pass=(ct_df_complete %>% group_by(age_bins) %>% dplyr::count() %>% 
                                ungroup() %>% complete(age_bins))$n,
                        fail=(ct_df_fail %>% group_by(age_bins) %>% dplyr::count() %>% 
                                ungroup() %>% complete(age_bins))$n)
  full_df <- reshape2::melt(full_df,id.vars="bin")
  full_df$variable <- factor(full_df$variable, levels = c("fail","pass"))
  full_df[is.na(full_df)] <- 0
  
  age_hist <- ggplot(full_df, aes(x=bin, y=value, fill=variable, color=variable)) +
    geom_bar(stat="identity", width=0.85) +
    scale_fill_manual(values=c("white", "black")) +
    scale_color_manual(values=c("black", "black")) +
    labs(x="Age", y="Number of samples") +
    scale_x_discrete(expand = c(0.08, 0), breaks = bins[seq(2,length(bins),2)],
                     labels = bins[seq(2,length(bins),2)]) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0,10,20,30,40), labels = c(0,10,20,30,40)) +
    theme_classic() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 30, vjust=0.8),
          legend.title = element_blank(),
          legend.position = "None") +
    geom_hline(yintercept=seq(0,40,10), color="white", size=0.2)
  
  ### SEX HISTOGRAM ###
  
  # stratify by sex
  ct_df <- seq_metadata %>% filter(sex!="unknown")

  ct_df_complete <- ct_df %>% filter(genome=="Yes")
  ct_df_fail <- ct_df %>% filter(genome!="Yes")
  
  full_df <- data.frame(bin=c("F","M"),
                        pass=(ct_df_complete %>% group_by(sex) %>% dplyr::count() %>% ungroup() %>% complete(sex))$n,
                        fail=(ct_df_fail %>% group_by(sex) %>% dplyr::count() %>% ungroup() %>% complete(sex))$n)
  full_df <- reshape2::melt(full_df,id.vars="bin")
  full_df$variable <- factor(full_df$variable, levels = c("fail","pass"))
  full_df[is.na(full_df)] <- 0
  
  sex_hist <- ggplot(full_df, aes(x=bin, y=value, fill=variable, color=variable)) +
    geom_bar(stat="identity", width=0.85) +
    scale_fill_manual(values=c("white", "black")) +
    scale_color_manual(values=c("black", "black")) +
    labs(x="Sex", y="Number of samples") +
    scale_x_discrete(expand = c(0.5, 0),labels=c("Female","Male")) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0,20,40,60,80,100), labels = c(0,20,40,60,80,100)) +
    theme_classic() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          #axis.text.x = element_text(angle = 30, vjust=0.8),
          legend.title = element_blank(),
          legend.position = "None") +
    geom_hline(yintercept=seq(0,80,20), color="white", size=0.2)
  
  ### DATE HISTOGRAM ###
  
  # stratify by date
  ct_df <- seq_metadata %>% filter(!is.na(collection_date))
  ct_df$collection_date <- as.numeric(ct_df$collection_date)
  
  ct_df_complete <- ct_df %>% filter(genome=="Yes")
  ct_df_fail <- ct_df %>% filter(genome!="Yes")
  
  hist_all <- hist(ct_df$collection_date, breaks = 6, plot = FALSE)
  hist_complete <- hist(ct_df_complete$collection_date, breaks = hist_all$breaks, plot = FALSE)
  hist_fail <- hist(ct_df_fail$collection_date, breaks = hist_all$breaks, plot = FALSE)
  
  full_df <- data.frame(bin=seq(1,length(hist_all$counts)),pass=hist_complete$counts,fail=hist_fail$counts)
  full_df <- melt(full_df,id.vars="bin")
  full_df$variable <- factor(full_df$variable, levels = c("fail","pass"))
  
  date_hist <- ggplot(full_df, aes(x=bin, y=value, fill=variable, color=variable)) +
    geom_bar(stat="identity", width=0.85) +
    scale_fill_manual(values=c("white", "black")) +
    scale_color_manual(values=c("black", "black")) +
    labs(x="Date", y="Number of samples") +
    scale_x_continuous(expand = c(0.01, 0), breaks = seq(1.5,6.5,1),
                       labels = format(as.Date(hist_all$breaks,origin="1970-01-01"),"%b-%d")[2:7]) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0,10,20,30,40,50), labels = c(0,10,20,30,40,50)) +
    theme_classic() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 30, vjust=0.7),
          legend.title = element_blank(),
          legend.position = "None") +
    geom_hline(yintercept=seq(0,50,10), color="white", size=0.2)
  
  plots <- list(ct_hist,age_hist,sex_hist,date_hist)
  return(plots)
  
}

## Function to make heatmap and dendogram
make.dendro.heatmap <- function(vars_file,outdir=NA,jhu_newick=NA,reorder_dendro=FALSE){
  
  # Load vars file
  allvars <- read.table(vars_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Add column for which bucket the position falls into
  allvars$posBucket = 60 * floor(allvars$pos / 60)
  
  # Add columns for ambiguity and masking
  allvars$isN = (allvars$alt != 'A' & allvars$alt != 'C' & allvars$alt != 'G' & allvars$alt != 'T')
  allvars$col = paste(allvars$isN, allvars$Masked)
  
  # Get the true variant calls grouped by position bucket
  grouped <- allvars %>% filter(pos >= 55 & pos <= 29836 & col == "FALSE No") %>% group_by(sample, posBucket) %>% tally()
  # Fill NA's with 0's
  grouped[is.na(grouped)] = 0
  
  # Get the masked and ambiguous positions grouped by position bucket
  maskeddf <- allvars %>% filter(pos >= 55 & pos <= 29836 & Masked == "Yes") %>% group_by(sample, posBucket) %>% tally()
  ambigdf <- allvars %>% filter(pos >= 55 & pos <= 29836 & col == "TRUE No") %>% group_by(sample, posBucket) %>% tally()
  # Fill NAs with 0's
  maskeddf[is.na(maskeddf)] = 0
  ambigdf[is.na(ambigdf)] = 0
  
  # Make a version with clade included to get updated leaf labels
  #allvars$cladeSample = paste("", allvars$clade, ": ", allvars$sample, sep = "")
  allvars$cladeSample = paste("", allvars$nextstrainclade, ": ", allvars$sample, sep = "")
  
  # Get all samples with counts of SNPs in each window
  allSitesWithClade <- allvars %>% filter(pos >= 55 & pos <= 29836) %>% group_by(cladeSample, posBucket) %>% tally(col == "FALSE No")
  allSitesWithClade[is.na(allSitesWithClade)] = 0
  
  # Convert variant counts to a matrix for clustering
  varmatrix = daply(.data=allSitesWithClade,
                    .variables=c("cladeSample","posBucket"),
                    .fun=function(x) sum(x$n))
  varmatrix[is.na(varmatrix)] = 0
  
  # Generate dendrogram with clades included
  varmatrix.cladedendro <- as.dendrogram(hclust(d = dist(x = varmatrix)))
  varmatrix.cladedendro <- sort(varmatrix.cladedendro, type = "labels", decreasing = FALSE)
  
  # reorder dendrogram is specified
  if (reorder_dendro==TRUE){
    
    # Load JHU tree in order to get node order
    jhu_tree <- read.tree(jhu_newick)
    # get order of tips in this tree
    is_tip_jhutree <- jhu_tree$edge[,2] <= length(jhu_tree$tip.label)
    tip_order_jhutree <- jhu_tree$tip.label[jhu_tree$edge[is_tip_jhutree,2]]
    
    # rename tips to match dendrogram
    fix.tree.names <- function(seq_name){
      
      seqext <- substr(seq_name,5,nchar(seq_name)-5)
      seqext <- str_replace(seqext,"..-HP","MDHP")
      return(seqext)
      
    }
    
    tip_order_jhutree <- sapply(tip_order_jhutree,fix.tree.names,USE.NAMES = FALSE)
    allvars_ordered <- left_join(data.frame(sample=tip_order_jhutree),allvars,by="sample")
    tip_order_jhutree <- rev(unique(allvars_ordered$cladeSample))
    
    # order dendrogram by jhu_tree tip order
    #varmatrix.cladedendro <- 
    #  dendextend::rotate(varmatrix.cladedendro, order(match(tip_order_jhutree,rownames(varmatrix))))
    
  }
  
  # Get the order of rows in the dendrogram
  varmatrix.order <- order.dendrogram(varmatrix.cladedendro)
  #varmatrix.order <- order(match(tip_order_jhutree,rownames(varmatrix)))
  newrows <- rownames(varmatrix)[varmatrix.order]
  #newrows
  
  # Generate a version of the dendrogram with no labels to incorporate it into the final plot later 
  dendro.cladeplot_nolabels <- ggdendrogram(data = varmatrix.cladedendro,rotate = TRUE) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"))
  #outfile <- paste(outdir, "cladedendro_nolabels.png", sep="")
  #ggsave(outfile, device = "png", width = 6, height = 18)
  
  # Merge this back with the variants
  df <- data.frame(cladeSample = unlist(newrows), index = rep(seq_along(newrows), lapply(newrows, length)))
  df_full <- left_join(df, allvars, by = "cladeSample", sort = F)
  find.idx <- function(val)df_full$index[match(val, df_full$sample)]
  
  # Reorder all of the data frames
  grouped$index = find.idx(grouped$sample)
  grouped <- grouped[order(grouped$index),]
  maskeddf$index = find.idx(maskeddf$sample)
  maskeddf <- maskeddf[order(maskeddf$index),]
  ambigdf$index = find.idx(ambigdf$sample)
  ambigdf <- ambigdf[order(ambigdf$index),]
  
  # Merge the SNPs with masked and ambiguous sites
  total <- merge(grouped, maskeddf, by = c("posBucket", "sample"), all = TRUE) %>% 
    merge(ambigdf, by = c("posBucket", "sample"), all = TRUE)
  total$index = find.idx(total$sample)
  total <- total[order(total$index),]
  total$sample.ordered <- factor(total$sample, levels=unique(total$sample), ordered=TRUE)
  total[is.na(total)] = 0
  
  #clademap <- distinct(allvars %>% select(clade, sample))
  clademap <- distinct(allvars %>% select(nextstrainclade, sample))
  total$id <- 1:nrow(total) 
  withclade <- left_join(total, clademap, by = "sample", sort = F)
  withclade <- withclade[order(withclade$id), ]
  withclade$sample.ordered <- factor(withclade$sample, levels=unique(withclade$sample), ordered=TRUE)
  
  withclade$type = "REF"
  withclade$type <- ifelse(withclade$n.x > 0, "SNP", withclade$type)
  withclade$type <- ifelse(withclade$n > 0, "Ambiguous Site", withclade$type)
  withclade$type <- ifelse(withclade$n.y > 0, "Masked Region", withclade$type)
  withclade$type <- factor(withclade$type, levels = c("SNP", "Ambiguous Site", "Masked Region", "REF", "B.1"))
  
  heatmap.plot <- ggplot(data = withclade, aes(x = posBucket, y = sample.ordered)) + theme_bw() +
    geom_tile(data = withclade, aes(fill = type)) +
    scale_x_continuous(limits = c(0,29903), expand = c(0, 0), breaks = c(2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000, 22500, 25000, 27500)) +
    scale_fill_manual(values = c("indianred2", "darkblue", "lightblue")) +
    #xlab("Genome Position") +
    #ylab("Sample") +
    theme(legend.key.size = unit(1,"line")) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          #legend.position = "top",
          legend.title = element_blank(),
          #legend.text = element_text(size = 4),
          legend.text = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"))
  
  #outfile <- paste(outdir, "clustersortedheatmap.png", sep="")
  #ggsave(outfile, device = "png", width = 18, height = 18)
  
  # Plot a single column with colors for clades
  colplot <- ggplot(data = withclade, aes(x = 1, y = sample.ordered)) + theme_bw() +
    geom_tile(data = withclade, aes(fill = nextstrainclade), color = "white") +
    labs(fill="Clade") +
    scale_fill_manual(values=c("#be1e2d","#fbb040","#8dc63f","#92278f","#00aeef")) +
    theme(axis.line=element_blank(),
          axis.title.x = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.text = element_text(size = 8),
          legend.title=element_text(size=8),
          #legend.key = element_rect(fill = NA),
          legend.spacing.x = unit(0.05, "cm"),
          plot.background=element_blank(),
          plot.margin = unit(c(0,0,0.27,0.2), "cm")) +
    guides(fill = guide_legend(ncol = 1,keywidth=0.75,keyheight=0.35,
                               default.unit="cm",override.aes = list(size = 4)))
  
  # Save legend for separate plotting later
  col_leg <- get_legend(colplot)
  
  # Update plot to not include legend
  colplot <- colplot + theme(legend.position = "none")
  #outfile <- paste(outdir, "cladecolumn.png", sep="")
  #ggsave(outfile, device = "png", width = 1, height = 18)
  
  plots <- list(heatmap.plot, colplot, dendro.cladeplot_nolabels,col_leg)
  return(plots)
  
}

## Function to plot total number of samples containing each variant
make_varplot <- function(alignment,ref_genome){
  
  # read in alignment of all sequences
  align <- read.alignment(alignment,"fasta")
  align <- as.matrix(align)
  
  # load reference genome
  ref <- read.fasta(ref_genome)
  ref <- as.vector(ref$MN908947.3)
  
  # function to get count of non-ref, non-N alleles per position
  count.vars <- function(pos,aln,ref){
    aln_col <- aln[,pos]
    aln_col <- as.vector(aln_col)
    varcount <- length(aln_col[(aln_col!="n") & (aln_col!=ref[pos])])
    return(varcount)
  }
  
  counts <- sapply(1:29903,count.vars,align,ref)
  counts[counts<1] <- NA
  
  df <- cbind(1:29903,counts)
  colnames(df) <- c("pos","count")
  df <- as.data.frame(df)
  
  # get key snps to set up annotations
  #key_snps <- c("C241T","C1059T\nThr265Ile","C3037T\nPhe924Phe","C14408T\nPro314Leu",
  #              "A23403G\nAsp614Gly","G25563T\nGln57His")
  key_snps <- c("-\n(C241T)","T265I\n(C1059T)","F924F\n(C3037T)","P314L\n(C14408T)",
                "D614G\n(A23403G)","Q57H\n(G25563T)")
  key_pos <- c(241,1059,3037,14408,23403,25563)
  key <- cbind(key_pos,key_snps)
  colnames(key) <- c("pos","text")
  key <- as.data.frame(key)
  key$pos <- as.numeric(as.matrix(key$pos))
  
  # merge data frames to get labels in the right place
  df <- left_join(df,key,by="pos")
  
  # set up genome annotations
  cds_starts <- c(266,13468,21563,25393,26245,26523,27202,27394,27894,28274,29558)
  cds_stops <- c(13483,21555,25384,26220,26472,27191,27387,27759,28259,29533,29674)
  cds_names <- c("orf1a","orf1b","S","orf3a","E","M","orf6","orf7a","orf8","N","orf10")
  
  num_genes <- length(cds_names)
  cds <- data.frame(xmin=cds_starts,xmax=cds_stops,ymin=rep(-10,num_genes),ymax=rep(-2,num_genes))
  #cds <- data.frame(xmin=cds_starts,xmax=cds_stops,
  #                  ymin=c(rep(c(-10,-18),num_genes/2),-10),
  #                  ymax=c(rep(c(-2,-10),num_genes/2),-2))
  
  # make the plot
  varplot <- ggplot() +
    geom_bar(data=df, aes(x=pos, y=count),stat="identity",width=50,fill="indianred2") +
    geom_text(data=df,aes(x=pos,y=count,label=text),vjust=1.5,hjust=-0.1,size=3) +
    scale_y_continuous(limits=c(-10,114),expand = c(0,0),breaks=c(19,38,57,76,95,114)) +
    #scale_y_continuous(expand = expand_scale(mult = c(0,0.1))) +
    scale_x_continuous(expand = c(0,0), breaks = seq(0,29903,5000)) +
    xlab("SARS-CoV-2 genome position") +
    ylab("Number of genomes") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(colour = "black"),
          #axis.title.x=element_blank(),
          axis.title.y=element_text(size = 10),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8)) +
    geom_hline(yintercept = 0,color="black",size=0.5) +
    geom_segment(aes(x = 0, y = 0, xend = 0, yend = 114),size=0.5) +
    geom_rect(data=cds[seq(1,11,2),],aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "lightblue", size = 0.2, color = "black", alpha = 0.8) +
    geom_rect(data=cds[seq(2,11,2),],aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "lightgreen", size = 0.2, color = "black", alpha = 0.5)
  
  for (i in c(1,3,2,6,10)){
    varplot <- varplot + annotate("text", label=cds_names[i], x = cds_starts[i] + 
                                    (cds_stops[i]-cds_starts[i])/2, y = -6, size = 3)
  }
  
  #for (i in c(2,6,10)){
  #  varplot <- varplot + annotate("text", label=cds_names[i], x = cds_starts[i] + 
  #                                  (cds_stops[i]-cds_starts[i])/2, y = -14, size = 3)
  #}
  
  return(varplot)
  
}