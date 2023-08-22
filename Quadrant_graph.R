library(grid)
library(ggplot2) # ggplot2_3.4.0
library(ggrepel) # ggrepel_0.9.2      
library(RColorBrewer) # RColorBrewer_1.1-3
library(svglite) # svglite_2.1.0


### filter rows to common tested gene list, then merge.
int.genes = intersect(rownames(Visium.pseudobulk.hippocampus.SOR.for.quadrant), rownames(bulk.RNA.seq.hippocampus.SOR.for.quadrant))
res = cbind(Visium.pseudobulk.hippocampus.SOR.for.quadrantv2[int.genes,], bulk.RNA.seq.hippocampus.SOR.for.quadrant[int.genes, c("logFC", "FDR")])
colnames(res)[2:3] = paste0("visium.", colnames(res)[2:3])
colnames(res)[4:5] = paste0("bulk", colnames(res)[4:5])
# switch logFC sign for experiment1 so CTRL+SOR samples have higher expression of positive fold change genes in both lists.
res$visium.logFC = -1 * res$visium.logFC
res$bulk.logFC = -1 * res$bulk.logFC
# create signed log10(FDR) for visualization.
res$visium.slfdr = log10(res$visium.FDR)*sign(res$visium.logFC)
res$bulk.slfdr = log10(res$bulk.FDR)*sign(res$bulk.logFC)
# create max FDR between both analyses for mapping to color/size.
res$max.fdr = pmax(res$visium.FDR, res$bulk.FDR)

### set up color palette.
palette.rdylbu = (RColorBrewer::brewer.pal(11, "RdYlBu"))
max.fdr = 1 - res$max.fdr

# if signs are concordant, use RdYlBu palette. Otherwise, use gray.
sign.concord = sign(res$visium.slfdr) * (sign(res$bulk.slfdr) )
sign.direction = sign(res$visium.slfdr)
color.data = max.fdr * sign.direction

### function to map numeric to color vectors. Found here: https://stackoverflow.com/questions/15006211/how-do-i-generate-a-mapping-from-numbers-to-colors-in-r/18749392, but unsure of original source.
map2color = function (x, pal, limits = NULL) 
{
  if (is.null(limits)) 
    limits = range(x)
  pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal) + 
                            1), all.inside = TRUE)]
}

max.fdr.colors = map2color(x=color.data, pal=rev(colorRampPalette(palette.rdylbu)(101)), limits=c(-1, 1))
# change non-concordant genes to grey.
idx.discordant = sign.concord == -1
max.fdr.colors[idx.discordant] = scales::alpha(rep("#AAAAAA", sum(idx.discordant)), max.fdr[idx.discordant])

p = ggplot(res, aes(visium.slfdr, bulk.slfdr)) +
  geom_point(color = ifelse(res$visium.FDR <= 0.05 & res$bulk.FDR <= 0.05, "red", "black"))

p1 = p + geom_label_repel(
  aes(label=gene.name, size=0.1),
  #nudge_y       = 36 - subset(dat, mpg > 30)$mpg,
  # segment.size  = 0.1,
  segment.color = "grey50",
  force=.1,
  max.iter=10000,
  #direction     = "x",
  seed=777,
  xlim=c(4.65, 5.45),
  ylim=c(0.4, 2.95),
  box.padding=0.18,
  label.padding=0.05,
  point.padding=0.09,
  label.size=0.1,
  direction="y"
)

p2 = p + geom_label_repel(
  aes(label=gene.name, size=0.2),
  #nudge_y       = 36 - subset(dat, mpg > 30)$mpg,
  segment.size  = 0.1,
  segment.color = "grey50",
  force=0.4,
  max.iter=10000,
  # direction     = "x",
  seed=777,
  xlim=c(-2.45, log10(0.05)),
  ylim=c(-2.95, log10(0.05)),
  box.padding=0.7,
  label.padding=0.3,
  point.padding=0.1,
  label.size=0.5
) 

# Function: get the x and y positions of a single ggrepel label
get.xy.pos.labs <- function(n) {
  grb <- grid.get(n)
  data.frame(
    x = xrg[1]+diff(xrg)*convertX(grb$x, "native", valueOnly = TRUE),
    y = yrg[1]+diff(yrg)*convertY(grb$y, "native", valueOnly = TRUE)
  )
}

# Get x and y plot ranges 
xrg <- ggplot_build(p1)$layout$panel_params[[1]]$x.range
yrg <- ggplot_build(p1)$layout$panel_params[[1]]$y.range

# prep color data.
channels.min = 0.1
channel.coef = 1.2
palette.rdylbu = (RColorBrewer::brewer.pal(11, "RdBu"))[-6]
color.data = max.fdr * sign.direction
max.fdr.colors = map2color(x=color.data, pal=rev(colorRampPalette(palette.rdylbu)(101)), limits=c(log10(0.05), -log10(0.05)))

# change non-concordant genes to grey.
idx.discordant = sign.concord == -1
max.fdr.colors[idx.discordant] = scales::alpha(rep("#AAAAAA" ,  sum(idx.discordant)), pmax(max.fdr[idx.discordant], channels.min))

# svglite("figures/bulk_bulk_scatterplot.svg", height=9.5, width=10.75, system_fonts = list(sans = "Arial"))
par(mar=c(6,6,6,2), xpd=FALSE)
plot( 
  x=res$visium.slfdr,
  y=res$bulk.slfdr,
  pch=16, cex=0, cex.lab=1.2, cex.main=1.8,
  col=scales::alpha(max.fdr.colors, pmax(max.fdr, channels.min)*channel.coef),
  #ylim=c(-20,20), 
  #xlim=c(-10,20),
  ylab=expression("signed -log"[10]*"(FDR) :  bulk RNA-seq hippocampus"),
  xlab=expression("signed -log"[10]*"(FDR) :  pseudobulk visium hippocampal subregions"),
  main="bulk RNA-seq Vs excitatory neurons"
)

# add gridlines.
abline(h=0, v=0, lty=1, lwd=1.5, col=scales::alpha("#333333", 0.8))
abline(h=c(log10(0.05), -log10(0.05)), v=c(log10(0.05), -log10(0.05)), col=scales::alpha("#333333", 0.8), lty=3, lwd=1.2)

# add points.
points(
  x=res$visium.slfdr,
  y=res$bulk.slfdr,
  pch=16, cex=pmax(max.fdr, channels.min)*channel.coef,
  col=scales::alpha(max.fdr.colors, pmax(max.fdr, channels.min)*channel.coef)
)

# add labels.
text(res$visium.slfdr, res$bulk.slfdr, res$gene.name, cex = 0.5, col="#444444", font=2, adj=0)

# add legend.
legend(3.25, -1.27, expression(bold("FDR = 0.05")), lwd=1.3, lty=3, cex=0.9, bty="n")
legend(2.8, -1.7, c("discordant between experiments", "concordant: higher in both experiments", "concordant: lower in both experiments"), title=expression(bold("direction of effect")), fill=c("#AAAAAA", palette.rdylbu[c(3,8)]), cex=0.9)

dev.off()
