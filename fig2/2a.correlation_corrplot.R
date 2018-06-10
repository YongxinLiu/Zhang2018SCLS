# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object

library("corrplot")
library("pheatmap")
library(ggcorrplot)

# Public file 1. "design.txt"  Design of experiment
design = read.table("../data/design.txt", header=T, row.names= 1, sep="\t") 

# Public file 2. "otu_table.txt"  raw reads count of each OTU in each sample
otu_table = read.delim("../data/otu_table.txt", row.names= 1,  header=T, sep="\t")

# setting subset design
if (TRUE){
	sub_design = subset(design,groupID %in% c("A50Cp0","A50Cp1","A50Cp2","A50Cp3","A50Cp7","A50Cp10","A50Cp14","A50Cp21","A50Cp28","A50Cp35","A50Cp42","A50Cp49","A50Cp63","A50Cp70","A50Cp77","A50Cp84","A50Cp91","A50Cp98","A50Cp112","A50Cp119") ) # select group1
}else{
	sub_design = design
}
sub_design$group=sub_design$groupID

# Set group order
if ("TRUE" == "TRUE") {
    sub_design$group  = factor(sub_design$group, levels=c("A50Cp0","A50Cp1","A50Cp2","A50Cp3","A50Cp7","A50Cp10","A50Cp14","A50Cp21","A50Cp28","A50Cp35","A50Cp42","A50Cp49","A50Cp63","A50Cp70","A50Cp77","A50Cp84","A50Cp91","A50Cp98","A50Cp112","A50Cp119"))
	}else{sub_design$group  = as.factor(sub_design$group)}

print(paste("Number of group: ",length(unique(sub_design$group)),sep="")) # show group numbers

# sub and reorder subdesign and otu_table
idx = rownames(sub_design) %in% colnames(otu_table)
sub_design = sub_design[idx,]
count = otu_table[, rownames(sub_design)]
norm = t(t(count)/colSums(count,na=T)) * 100 # normalization to total 100


# Pearson correlation among groups
sampFile = as.data.frame(sub_design$group,row.names = row.names(sub_design))
colnames(sampFile)[1] = "group"
mat = norm
mat_t = t(mat)

mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno

# 计算相关系数，并保留3位小数
sim=cor(mat_mean_final,method="pearson")
sim=round(sim,3)

pdf(file="ggcorrplot_pearson_A50Cp.pdf", height = 2.5, width = 2.5)
ggcorrplot(sim, type = "upper", colors = c("green", "yellow", "red")) # , method = "circle"
dev.off()

# 人为设置颜色度
col1 <- colorRampPalette(c("green", "green", "red"))

pdf(file="corplot_pie_pearson_A50Cp.pdf", height = 2.5, width = 2.5)
corrplot(sim, method="pie", type="lower", col=col1(100)) # , diag=F , na.label = "1"
dev.off()


## 生成图例 

col1 <- colorRampPalette(c("green", "red"))
corrplot(sim, method="pie", type="lower", col=col1(100)) # , diag=F , na.label = "1"

# 生成时间热图，分别为土和植物的
time1 = c(0,1,2,3,7,10,14,21,28,35,42,49,63,70,77,84,91,98,112,119)
time2 = c(0,41,48,54,62,77,84,90,97,119,0,0,0,0,0,0,0,0,0,0)
time=data.frame(time1,time2)
pheatmap(time, cluster_rows = F,  cluster_cols = F)
pheatmap(time, cluster_rows = F,  cluster_cols = F, filename = "corplot_pie_legend_time.pdf" ,width=2, height=4)

