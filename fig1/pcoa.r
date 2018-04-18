# Fig.1 PCoA

# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory

# clean enviroment object
rm(list=ls()) 

# load related packages
library("ggplot2") 
library("vegan")

# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
main_theme = theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(size=.5, colour="black"),
                    axis.line.y=element_line(size=.5, colour="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(color="black", size=7),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    legend.text= element_text(size=7),
                    text=element_text(family="sans", size=7))

# Public file 1. "design.txt"  Design of experiment
design = read.table("../data/design.txt", header=T, row.names= 1, sep="\t") 

# setting subset design
if (TRUE){
	sub_design = subset(design,groupID %in% c("A50Cp0","A50Cp1","A50Cp10","A50Cp112","A50Cp119","A50Cp14","A50Cp2","A50Cp21","A50Cp28","A50Cp3","A50Cp35","A50Cp42","A50Cp49","A50Cp63","A50Cp7","A50Cp70","A50Cp77","A50Cp84","A50Cp91","A50Cp98","A50Sz0","A50Sz1","A50Sz10","A50Sz118","A50Sz13","A50Sz2","A50Sz20","A50Sz27","A50Sz3","A50Sz34","A50Sz41","A50Sz48","A50Sz5","A50Sz56","A50Sz62","A50Sz69","A50Sz7","A50Sz76","A50Sz83","A50Sz90","A50Sz97","HNCp112","HNCp119","HNSz118","IR24Cp0","IR24Cp1","IR24Cp10","IR24Cp112","IR24Cp119","IR24Cp14","IR24Cp2","IR24Cp21","IR24Cp28","IR24Cp3","IR24Cp35","IR24Cp42","IR24Cp49","IR24Cp63","IR24Cp7","IR24Cp70","IR24Cp77","IR24Cp84","IR24Cp91","IR24Cp98","IR24Sz0","IR24Sz1","IR24Sz10","IR24Sz118","IR24Sz13","IR24Sz2","IR24Sz20","IR24Sz27","IR24Sz3","IR24Sz34","IR24Sz41","IR24Sz48","IR24Sz5","IR24Sz56","IR24Sz62","IR24Sz69","IR24Sz7","IR24Sz76","IR24Sz83","IR24Sz90","IR24Sz97","soilCp0","soilCp42","soilCp49","soilCp57","soilCp63","soilCp77","soilCp84","soilCp91","soilCp98","soilSz0","soilSz41","soilSz48","soilSz54","soilSz62","soilSz76","soilSz83","soilSz90","soilSz97") ) # select group1
}else{
	sub_design = design
}
# if (TRUE){
# 	sub_design = subset(sub_design,site %in% c("Cp","Sz") ) # select group2
# }
# 
# # Set group style, single or combine
# if (FALSE){
# 	sub_design$group=paste(sub_design$groupID,sub_design$site,sep = "")
# }else{
# 	sub_design$group=sub_design$groupID
# }


print(paste("Number of group: ",length(unique(sub_design$group)),sep="")) # show group numbers

#  PCoA bray_curtis
bray_curtis = read.table("../data/bray_curtis_otu_table_css.txt", sep="\t", header=T, check.names=F)

# subset matrix and design
idx = rownames(sub_design) %in% colnames(bray_curtis) 
sub_design = sub_design[idx,]
bray_curtis = bray_curtis[rownames(sub_design), rownames(sub_design)] # subset and reorder distance matrix

# cmdscale {stats}, Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
pcoa = cmdscale(bray_curtis, k=4, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z","a") 
eig = pcoa$eig
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)), ])
#points$group = factor(points$group, levels=colors$group)

# plot PCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=day, shape = compartment))
p = p + geom_point(alpha=.7, size=2) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="bray_curtis PCoA") + main_theme

p
ggsave("beta_pcoa_day_bray_curtis_default.pdf", q, width = 4, height = 2.5)

# scale_color_gradientn 按多种颜色连续着色，如彩虹色
# topo.colors(), rainbow(), heat.colors(), terrain.colors(), cm.colors(), RColorBrewer::brewer.pal()
q= p + scale_color_gradientn(colours=rainbow(7))
q
ggsave("beta_pcoa_day_bray_curtis_rainbow.pdf", q, width = 4, height = 2.5)


# plot PCo 1 and 3
points$siteXcompt=paste(points$site,points$compartment,sep = "")

p = ggplot(points, aes(x=x, y=z, color=day, shape = siteXcompt))
p = p + geom_point(alpha=.7, size=2) +
  scale_color_gradientn(colours=rainbow(7)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep=""),
       title="bray_curtis PCoA") + main_theme
p
ggsave("beta_pcoa_day_bray_curtis3.pdf", p, width = 4, height = 2.5)

