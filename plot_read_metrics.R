#!/usr/bin/Rscript

# Plot from FastQ (ONT or PacBIO)
# usage: plot_read_metrics.R -i fastq_metrics.txt file
#
# Stephane Plaisance VIB-NC November-21-2023 v1.0
# depends on external command for read metrics extraction
# eg: zcat ${fastq}.gz | fastq2metrics > fastq_metrics.txt

# R libraries
suppressMessages(library("optparse")
suppressMessages(library("readr"))
suppressMessages(library("plyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggpubr"))

option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="read metrics text file", metavar="character"),
  make_option(c("-m", "--minlength"), type="numeric", default="0", 
              help="minimal read length [default= %default]", metavar="numeric"),
  make_option(c("-M", "--maxlength"), type="numeric", default="1000000", 
              help="maximal read length [default= %default]", metavar="numeric"),
  make_option(c("-q", "--minquality"), type="numeric", default="0", 
              help="minimal read average basecall quality [default= %default]", metavar="numeric"),
  make_option(c("-Q", "--maxquality"), type="numeric", default="99", 
              help="maximal read average basecall quality [default= %default]", metavar="character"),
  make_option(c("-f", "--format"), type="character", default="png", 
              help="output format (png or pdf) [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$infile)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# extract basename for output
basename_without_ext <- tools::file_path_sans_ext(basename(opt$infile))

Nvalue <- function(lim, x, na.rm = TRUE){
  # handle NA values
  if(isTRUE(na.rm)){
    x <- x[!is.na(x)]
  }
  cutval <- 100/lim
  # compute LXX and NXX
  sorted <- sort(x, decreasing = TRUE)
  SXX <- sum(x)/cutval
  csum <- cumsum(sorted)
  GTLXX <- as.vector(csum >= SXX)
  LXX=min(which(GTLXX == TRUE))
  NXX <- round(sorted[LXX], 1)
  # eg: get NXX with lst['NXX']
  NXX
}

full_data <- read_delim(opt$infile, delim = "\t", col_names = TRUE, show_col_types = FALSE)

# filter data based on user limits
slen=opt$minlength
elen=opt$maxlength
squal=opt$minquality
equal=opt$maxquality

read_data <- subset(full_data, length>=slen & length <= elen & meanq >= squal & meanq <= equal)

minl <- min(read_data$length)
maxl <- max(read_data$length)
n50reads <- Nvalue(50,read_data$length)

minQ <- min(c(40, min(read_data$meanq)))
maxQ <- min(c(99, max(read_data$meanq)))
medQ <- round(median(read_data$meanq),1)

# total yield in Mb
totbasecount <- round(sum(read_data$length/1E6),2)

# head(read_data)

# # plot lengths
#  scale_x_log10() +
p1 <- ggplot() + 
 geom_density(data=read_data, aes(length, colour="orange"), lwd=1.25, show.legend=FALSE) +
 xlim(c(minl, maxl)) +
 scale_x_log10() +
 labs(x = "", y = "density") +
 theme_linedraw() +
 theme(plot.title = element_text(margin=margin(b=0), size = 14)) +
 ggtitle(paste0("read length ([", minl, "..", maxl, "], N50=" ,n50reads, ")"))

# plot meanQ
p2 <- ggplot() + 
  geom_density(data=read_data, aes(meanq, colour="red"), lwd=1.25, show.legend=FALSE) +
  labs(x = "", y = "density") +
  xlim(c(minQ, maxQ)) +
  theme_linedraw() +
  theme(plot.title = element_text(margin=margin(b=0), size = 14)) +
  coord_flip() +
  ggtitle(paste0("average basecall quality ([", round(minQ,1), "..", round(maxQ,1), "], median=" , round(medQ,1), ")"))
   
# biplot meanQ x length
# coord_cartesian(xlim = c(minl, maxl), ylim = c(minQ, maxQ)) +
p3 <- ggplot(read_data, aes(x=length, y=meanq)) + 
 geom_point(pch=20, cex=0.75, col="grey60") +
 scale_x_log10() +
 coord_cartesian(xlim = c(minl, maxl), ylim = c(minQ, maxQ)) +
 labs(x = "read length", y = "average basecall quality") +
 stat_density_2d(aes(fill = after_stat(level)), geom="polygon") +
 scale_fill_gradient(low="blue", high="red") +
 geom_hline(aes(yintercept=20), linewidth=0.5, colour="green", lty=1) +
 geom_hline(aes(yintercept=30), linewidth=0.5, colour="blue", lty=2) +
 geom_vline(aes(xintercept=n50reads), linewidth=0.5, colour="grey25", lty=3) +
 theme_linedraw() +
 theme(plot.title = element_text(margin=margin(b=0), size = 14),
       legend.position = "none") +
 annotate(geom="text", 
          x=minl+2/5*(maxl-minl), 
          y=minQ+1/4*(maxQ-minQ), 
          label="green line:Q20\n dashed-blue line Q30\ndotted-grey line len N50",
          color="black",
          cex=4) +
 ggtitle(paste0("read-count: ", nrow(read_data), " - total base count (Mb): ", totbasecount))

p0 <- ggplot() +
  theme_void()  # Remove axes and labels


# Save plots in the chosen format (PDF or PNG)
if (opt$format == "pdf") {
  outfile <- paste0(basename_without_ext, ".pdf")
  pdf(outfile, width = 10, height = 10)
  } else {
    outfile <- paste0(basename_without_ext, ".png")
    png(outfile, width = 1600, height = 1600, res = 150)
    }

ggarrange(p1, p0, p3, p2, 
          labels = c("A", "", "B", "C"),
          ncol = 2, nrow = 2)

null <- dev.off()
