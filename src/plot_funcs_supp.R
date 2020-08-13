#!/usr/bin/env Rscript

# load necessary packages
library(dplyr)
library(ggplot2)
library(ggExtra)
library(tidyr)
library(ggforce)
library(extrafont)
library(ggpubr)
library(gridExtra)
#library(Cairo)

make_ont_ill_plots <- function(postfile_file, postfilt_summary, vars_file){

	# Read in true variant table
	allvars <- read.table(vars_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
	allvars$isN = (allvars$alt != 'A' & allvars$alt != 'C' & allvars$alt != 'G' & allvars$alt != 'T')

	af_data <- data.frame()

	# read in text file
	df <- read.table(postfilt_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
	df <- filter(df, illumina_support != ".")
	#df <- df %>% select(sample,pos)
	af_data <- rbind(af_data,df)
	af_data$sample <- gsub("-", "", af_data$sample)

	unique(af_data$sample)
	af_data <- subset(af_data, af_data$sample %in% allvars$sample)
	unique(af_data$sample)

	# loop through the postfilt summary files
	# to determine what samples are 'no' and should not be included

	sample_status <- data.frame()

	# read in text file
	df <- read.table(postfilt_summary, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

	# keep only the relevant columns
	df <- df %>% select(Sample,Status)
	df <- df %>% dplyr::rename(sample = Sample, overall_status = Status)

	# combine all sample data into one data frame
	sample_status <- rbind(sample_status,df)

	# clean up the data frame

	# combine allele frequency data with status data
	af_data <- merge(af_data,sample_status,by = "sample",all.x = TRUE)

	# exclude samples marked as no
	#af_data <- filter(af_data, overall_status != "No")

	# make sure the values in the allele frequency columns are numeric
	af_data$ont_AF = as.numeric(af_data$ont_AF)
	af_data$illumina_AF = as.numeric(af_data$illumina_AF)

	# add a column for low confidence flags
	maybe_flags = c("depth_flag","ntc_flag","sb_flag")
	af_data$maybe_flags = apply(af_data[maybe_flags],1,function(x) any(x!="."))

	# add a column for called or not called
	called_vars <- allvars %>% select(sample, pos, isN)
	called_vars$called <= "Called"
	called_vars$called <- ifelse(called_vars$isN, "Ambiguous", "Called")
	af_data <- merge(af_data,called_vars,by=c("sample","pos"),all.x = TRUE)
	af_data$called[is.na(af_data$called)] <- "Uncalled"

	# summary column for confidence and called
	af_data <- af_data %>% mutate(confidence=ifelse(called=="Uncalled","uncalled",
		                                            ifelse(maybe_flags==TRUE,"low","high")))
	af_data <- af_data %>% mutate(homopolymer=ifelse(pos==10712 | pos == 21575, "True", homopolymer))

	# add a column for the absolute value of the difference in allele frequencies
	af_data$freq_diff = af_data$ont_AF - af_data$illumina_AF

	af_data$samplepos = paste(af_data$sample, ":", af_data$pos)
	allvars$samplepos = paste(allvars$sample, ":", allvars$pos)
	af_data_genome <- merge(af_data,allvars,by = "samplepos")

	# optional filter to remove uncalled variants
	af_data_called <- filter(af_data,confidence!="uncalled")

	### COLOR PLOT BY CONFIDENCE FLAGS

	# plot of all variants
	p_all <- ggplot(af_data, aes(x=ont_AF, y=illumina_AF, color = homopolymer, shape=called)) + 
	  geom_point() +
	  annotate("rect", xmin=-Inf, xmax=0.25, ymin=-Inf, ymax=Inf, alpha=0.2, fill="gray") +
	  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0.25, alpha=0.2, fill="gray") +
	  geom_abline(slope=1, intercept=0, linetype = "dashed") +
	  scale_shape_manual(labels = c("Ambiguous", "Called", "Uncalled"), values=c(1, 19, 3))+
	  scale_color_manual(labels = c("Non-homopolymer", "Homopolymer"),values=c("skyblue2", "red2")) +
	  labs(x = "ONT allele frequency", y = "Illumina allele frequency") +
	  theme_classic() +
	  guides(color=FALSE) + # dont show color legend on this plot
	  theme(legend.position = c(0.75, 0.2), legend.title=element_blank(),legend.background = element_blank(),
		    panel.grid.major = element_line(colour="black", size = (0.2), linetype = "dotted"),
		    panel.grid.minor = element_line(colour="black", size = (0.1), linetype = "dotted"))

	# plot of called variants only
	p_called <- ggplot(af_data_called, aes(x=ont_AF, y=illumina_AF, color=confidence)) + 
	  geom_point() +
	  annotate("rect", xmin=-Inf, xmax=0.25, ymin=-Inf, ymax=Inf, alpha = 0.2, fill="gray") +
	  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0.25, alpha = 0.2, fill="gray") +
	  geom_abline(slope=1, intercept=0, linetype = "dashed") +
	  scale_color_manual(labels = c("variant", "low-confidence variant", "uncalled varaint"),values=c("blue", "skyblue2","lightgray")) +
	  labs(x = "ONT allele frequency", y = "Illumina allele frequency") +
	  theme_classic() +
	  theme(legend.position = c(0.75, 0.1),legend.title=element_blank(),legend.background = element_blank(),
		    panel.grid.major = element_line(colour="black", size = (0.2), linetype = "dotted"),
		    panel.grid.minor = element_line(colour="black", size = (0.1), linetype = "dotted"))

	af_data_called$absdiff = abs(af_data_called$freq_diff)
	importantDf <- af_data_called %>% group_by(pos) %>% top_n(1, absdiff) %>% filter(abs(freq_diff) > .20 | homopolymer == 'True')

	#importantPositions <- af_data_called %>% filter(abs(freq_diff) > .25 | homopolymer == 'True')
	#importantPositions <- unique(importantPositions$pos)
	#importantDf <- data.frame(importantPositions)
	#importantDf$homopolymer = "NA"
	## PLOT THE DIFFERENCE IN ALLELE FREQUENCIES BY POSITION

	p3 <- ggplot(af_data_called %>% arrange(desc(confidence)), aes(x=pos, y=freq_diff, color=homopolymer)) + 
	  #geom_vline(xintercept = importantDf$pos, colour="lightgray", size = (0.1), linetype = "longdash") +
	  geom_point() +
	  geom_hline(yintercept = 0) + 
	  #geom_text(aes(label=ifelse(homopolymer,pos,'')),hjust=0, vjust=0) +
	  scale_color_manual(labels = c("Non-homopolymer", "Homopolymer"),values=c("skyblue2", "red2")) +
	  labs(x = "Position", y = "ONT allele freq - Illumina allele freq") +
	  theme_classic() +
	  theme(legend.position = c(0.78,0.9), legend.title=element_blank(),legend.background = element_blank(),
		    panel.grid.major = element_line(colour="black", size = (0.2), linetype = "dotted"),
		    panel.grid.minor = element_line(colour="black", size = (0.05), linetype = "dotted"),
		    axis.line.x=element_blank()) +
	  geom_text(data = importantDf, color = "black", aes( x= pos - 220, y=freq_diff + .015 * absdiff/freq_diff, label=pos), size = 1.8, angle = 0)
	  #geom_text(data=importantDf, color = "black", mapping=aes(x=pos, label = pos, y= freq_diff + .06 * (absdiff / freq_diff - 1) + .02 * (absdiff / freq_diff + 1)), size=2, vjust=-0.4, hjust=0)

	plots <- list(p_all, p_called, p3)
	return(plots)
}
