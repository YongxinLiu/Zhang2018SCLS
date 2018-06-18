# Install related packages
if (FALSE){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggplot2","grid","scales","vegan","agricolae","ggrepel","dplyr"))
	install.packages("devtools", repo="http://cran.us.r-project.org")
	library(devtools)
	install_github("vqv/ggbiplot")
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
library("ggplot2") # load related packages
library("grid")
library("scales")
library("vegan")
library("agricolae")
library("ggbiplot")
library("dplyr")
library("ggrepel")

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
	design = subset(design,groupID %in% c("A50Cp0","A50Cp1","A50Cp10","A50Cp112","A50Cp119","A50Cp14","A50Cp2","A50Cp21","A50Cp28","A50Cp3","A50Cp35","A50Cp42","A50Cp49","A50Cp63","A50Cp7","A50Cp70","A50Cp77","A50Cp84","A50Cp91","A50Cp98","A50Sz0","A50Sz1","A50Sz10","A50Sz118","A50Sz13","A50Sz2","A50Sz20","A50Sz27","A50Sz3","A50Sz34","A50Sz41","A50Sz48","A50Sz5","A50Sz56","A50Sz62","A50Sz69","A50Sz7","A50Sz76","A50Sz83","A50Sz90","A50Sz97","IR24Cp0","IR24Cp1","IR24Cp10","IR24Cp112","IR24Cp119","IR24Cp14","IR24Cp2","IR24Cp21","IR24Cp28","IR24Cp3","IR24Cp35","IR24Cp42","IR24Cp49","IR24Cp63","IR24Cp7","IR24Cp70","IR24Cp77","IR24Cp84","IR24Cp91","IR24Cp98","IR24Sz0","IR24Sz1","IR24Sz10","IR24Sz118","IR24Sz13","IR24Sz2","IR24Sz20","IR24Sz27","IR24Sz3","IR24Sz34","IR24Sz41","IR24Sz48","IR24Sz5","IR24Sz56","IR24Sz62","IR24Sz69","IR24Sz7","IR24Sz76","IR24Sz83","IR24Sz90","IR24Sz97","soilCp0","soilCp42","soilCp49","soilCp57","soilCp63","soilCp77","soilCp83","soilCp91","soilCp98","soilCp118","soilSz0","soilSz41","soilSz48","soilSz54","soilSz62","soilSz76","soilSz84","soilSz90","soilSz97","soilSz119") ) # select group1
}
otu_table = read.delim("../data/otu_table.txt", row.names= 1,  header=T, sep="\t")

calc_dis = function(subset){
  
# 筛选名中包含部分关键字的组
idx = grepl(subset, design$groupID)
sub_design = design[idx,]
print(paste("Number of group: ",length(unique(sub_design$group)),sep="")) # show group numbers

# 筛选样品，计算距离矩阵
# 有些样品异常会剔除，需要交叉筛选
idx=rownames(sub_design) %in% colnames(otu_table)
sub_design=sub_design[idx,]
count = otu_table[, rownames(sub_design)]
norm = t(t(count)/colSums(count,na=T)) * 100 # normalization to total 100
norm=as.data.frame(norm)

# 筛选119天，作为终点，再计算各点到终点的距离，观察群落结构变化

final=rownames(sub_design[sub_design$day2 %in% 119,])
ck=norm[,final] # 筛选并计算均值
ck_mean= rowMeans(ck)
norm$ck_mean=ck_mean # 均值添加到OTU表

# 计算包括终点均值的所有样品bray距离
bray_curtis = vegdist(t(norm), method = "bray")
bray_curtis= as.matrix(bray_curtis)

# 计算各组内部差异程度
# 建立一个存储组内差异的数据框
dat = t(as.data.frame(c("sampleA","sampleB","0.15","group","genosite")))
colnames(dat) = c("sampleA","sampleB","distance","group","type")
rownames(dat) = c("test")

# 每个样品与final对应的距离
for (i in sort(unique(sub_design$day2))){
  group = rownames(sub_design[sub_design$day2 %in% i,])
  for (m in 1:(length(group)-1)) {
    x = c(group[m],"ck_mean",bray_curtis[group[m],"ck_mean"],i,subset)
    dat=rbind(dat,x)
  }
}
dat = as.data.frame(dat[-1,], stringsAsFactors=F) # 删除首行框架数据
# dat = dat[dat$distance != 0,] # 
# 距离转换为3位数值
dat$distance=round(as.numeric(as.character(dat$distance)), digits=3)
# 分组添加levels，分组有时间顺序
dat$group = factor(dat$group, levels=unique(dat$group))
return(dat)
}

# 计算第一大类A50品种在Cp地点下每个点到终点的距离
dat = calc_dis("A50Cp")
# 计算第一次结果赋给all，以后每类别修改品种和地点重新运行上面55 - 103行后，cbind追加dat至all中
all = dat

#分别按"A50","IR24", X "Cp","Sz"筛选数据
for (i in c("A50Sz", "IR24Cp", "IR24Sz", "soilCp", "soilSz")){
  dat = calc_dis(i)
  all = rbind(all, dat)
}

# 抽样部分检查
i=sample(dat$group,1)
head(dat[dat$group==i,])

all$group = factor(all$group, levels = sort(as.numeric(as.character(levels(all$group)))))


# 先用各组箱线图查看数据分布
p = ggplot(all, aes(x=group, y=distance, color=type)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="Groups", y="Bray-Curtis distance") + theme_classic()  
p=p+geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)
p=p+theme(axis.text=element_text(angle=45,vjust=1, hjust=1))
p
ggsave(paste("e.boxplot",".pdf", sep=""), p, width = 8, height =5 )

# 规律不方便，土壤变化小没有像根系那样连续取样，展示中也比较丑

# 开始尝试散点图拟合
require(splines)
require(MASS)
all$group = as.numeric(as.character(all$group))
p = ggplot(all, aes(x=group, y=distance, color=type)) +
  geom_point(alpha=1, size=0.7) + geom_jitter()+
  labs(x="Groups", y="Bray-Curtis distance") + 
  theme(axis.text=element_text(angle=45,vjust=1, hjust=1)) +
  # 拟合只有loess方法初始为曲线，lm方法默认为直线，但可以添加公式
  main_theme + geom_smooth(method = "lm", formula = y ~ poly(x,3)) # , formula = y ~ ns(x,3)
p
ggsave(paste("e.beta_scatter_bray_curtis_mean6",".pdf", sep=""), p, width = 4, height =2.5 )
# 为什么拟合曲线没有了


# 我们的另一个备选方案，拆线图也很好看
library(ggpubr)
p=ggline(all, x="group", y="distance", add = "mean_se",  color = "type",
       palette = "jco",legend = "right") +main_theme
p
ggsave(paste("e.line",".pdf", sep=""), p, width = 5, height = 3)
