#!/usr/bin/Rscript

# script: lendist.R
# Aim: plot contig cumulative length for one or several fasta assemblies
# input can be archived fasta files
# requires bioawk (https://github.com/lh3/bioawk) for rapid feature extraction
#
# Stéphane Plaisance - VIB-Nucleomics Core - 2017-05-11 v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# depends on
library("readr")
library("ggplot2")

# adapt the path to your bioawk executable
bioawk <- "/opt/biotools/bioawk/bioawk"

# provided arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one fasta must be supplied\n", call.=FALSE)
}

# allow blobbing several fasta files
file.list <- paste(args, collapse=" ")

# build command
cmd <- sprintf("cnt=0; cat /dev/null > /tmp/len_info.txt; for f in %s;  do cnt=$((cnt+1)); %s -c fastx -v id=$cnt -v asm=$(basename $f) \'{print id, asm, $name, length($seq)}\' $f >> /tmp/len_info.txt; done", file.list, bioawk)

cat("# collecting lengths from file(s)\n")

t1 <- try(system(cmd, intern = TRUE, wait = TRUE))

cat("# loading results and plotting\n")

data <- read_delim("/tmp/len_info.txt",
                   "\t",
                   escape_double = FALSE,
                   col_names = FALSE,
                   trim_ws = TRUE,
                   col_types = cols())

colnames(data) <- c("idx", "asm", "source","length")
data <- data[order(data$idx,data$length),]
# convert to Mb
data$length <- data$length/1000000

# add columns for plotting
data$cumulative.length <- unlist(tapply(data$length, data$idx, cumsum))
data$contig.index <- ave(data$idx, data$idx, FUN = seq_along)

plots <- pdf(file = "assembly_sizes.pdf", width=5, height=4, onefile = TRUE)

# display y-scale as integer
options(scipen=1000000)

# basis-plot
p1 <- ggplot(data=data, aes(x=contig.index, y=cumulative.length, group=asm)) +
  geom_point(aes(shape=asm), size=1, col="gray50") +
  scale_shape_manual(values = 0:length(unique(data$asm))) +
  theme(axis.text.x = element_text(colour="grey20",size=8,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=8,angle=0,hjust=1,vjust=0,face="plain"),
        axis.title.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=10,angle=90,hjust=.5,vjust=.5,face="plain"),
        #legend.justification = c(0,1),
        #legend.position = c(0,1),
        legend.title=element_blank(),
        legend.text=element_text(size=5),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.8, 'lines'),
        legend.background = element_rect(fill="transparent"),
        plot.title = element_text(margin=margin(b=0), size = 14))

# plot linear-scale
p1 +
  labs(title="Cumulative length plot\n", x="contig index", y="cumulative size (Mb)")

# plot log-scale
p1 +
  labs(title="Cumulative length plot (log-scale)\n", x="contig index (log10)", y="cumulative size (Mb)") +
  scale_x_log10()

# add boxplot with contig size distribution
data$loglength=log(data$length*1000,10)
formatBack <- function(x) 10^x
p2 <- ggplot(data=data, aes(asm, loglength), size=0.5)
p2 + geom_boxplot() + 
	scale_y_continuous(labels=formatBack) + 
	coord_flip() +
	labs(title="Contig size distribution", y="contig / chr size (kbps)", x="") + 
	theme(axis.text=element_text(size=8), axis.title=element_text(size=10,face="bold"))

close.plots <- dev.off()

cat("# results were plotted to assembly_sizes.pdf\n")
