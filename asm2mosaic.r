#!/usr/bin/Rscript

# scrip: asm2mosaic.r
# parse multifasta (argument)
# plot mosaic to PDF
# plot cumsize contig N-plot
#
# Stephane Plaisance - VIB-Nucleomics Core - August-09-2022 v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

suppressMessages({
require("fastaR")
require("treemap")
require("gridExtra")
require("scales")
require("ggrepel")
})

args = commandArgs(trailingOnly=TRUE)

sizes <- fastaR::fa_size(fasta_file=args[1])

# sort decreasing and add extra columns
sizes <- sizes[order(sizes$Length, decreasing=TRUE),]
sizes$cum <- cumsum(sizes$Length)
sizes$idx <- seq(nrow(sizes))
tot <- max(sizes$cum)
num <- max(sizes$idx)

# mosaic plot
suppressWarnings({
pdf("mosaic_plot.pdf", width=6, height=3)
treemap(sizes,
        index="Seq_id",
        vSize="Length",
        type="index",
        title="",
        fontsize.labels=8,
        lowerbound.cex.labels=0,
        border.lwds=0.5)
null <- dev.off()
})

# cumulative plot
#suppressWarnings({
pdf("N_plot.pdf", width=6, height=6)
ggplot(data=sizes, aes(x=idx, y=cum, group=1)) +
  geom_line(col="grey70")+
  geom_point()+
  geom_label_repel(aes(label = Seq_id), max.overlaps=10, nudge_x=1, nudge_y=1)+
  theme_bw()+
  geom_hline(yintercept = tot*0.5, col="grey80", lwd=1, linetype='dotted')+
  geom_text(aes(x = num-1, y = tot*0.53, label = "L50"), col="grey80", size=8)+ 
  xlab("contig number")+
  ylab("assembly size")+
  scale_y_continuous(labels = unit_format(unit = "e+06", scale = 1 / 1e+06, digits = 2))
null <- dev.off()
#})
