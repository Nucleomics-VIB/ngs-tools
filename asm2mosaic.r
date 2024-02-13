#!/usr/bin/Rscript

# scrip: asm2mosaic.r
# parse multifasta (argument)
# plot mosaic to PDF
# plot cumsize contig N-plot
# replace deprecated library fastaR with ShortRead
#
# Stephane Plaisance - VIB-Nucleomics Core - June-19-2023 v2.0
#
# visit our Git: https://github.com/Nucleomics-VIB

suppressMessages({
require("ShortRead")
require("treemap")
require("gridExtra")
require("scales")
require("ggrepel")
})

args = commandArgs(trailingOnly=TRUE)

sequences <- readFasta(fasta_file=args[1])
sizes <- width(sequences)

# create data.frame
seq_df <- data.frame(seqlen=as.numeric(sizes[order(sizes, decreasing=TRUE)]))
seq_df$cum <- cumsum(seq_df$seqlen)
seq_df$idx <- seq(nrow(seq_df))

tot <- max(seq_df$cum)
num <- max(seq_df$idx)

# mosaic plot
suppressWarnings({
pdf("mosaic_plot.pdf", width=6, height=3)
treemap(sizes,
        index="idx",
        vSize="seqlen",
        type="index",
        title="",
        fontsize.labels=8,
        lowerbound.cex.labels=0,
        border.lwds=0.5)
null <- dev.off()
})

# cumulative plot
suppressWarnings({
pdf("N_plot.pdf", width=6, height=6)
ggplot(data=seq_df, aes(x=idx, y=cum, group=1)) +
  geom_line(col="grey70")+
  geom_point()+
  geom_label_repel(aes(label = idx), max.overlaps=10, nudge_x=1, nudge_y=1)+
  theme_bw()+
  geom_hline(yintercept = tot*0.5, col="grey80", lwd=1, linetype='dotted')+
  geom_text(aes(x = num-1, y = tot*0.53, label = "L50"), col="grey80", size=8)+ 
  xlab("contig number")+
  ylab("assembly size")+
  scale_y_continuous(labels = unit_format(unit = "e+06", scale = 1 / 1e+06, digits = 2))
null <- dev.off()
})
