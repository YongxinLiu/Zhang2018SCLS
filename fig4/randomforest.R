# 随机森林回归 Random Forest Regression


## 回归分析

# 我们主要有两个品种，种植在两个地点。这里先以A50建模，IR24验证的方案来演示。本实验较复杂 ，具体的方法会有多种组合。

# 读取实验设计、和物种分类文件
tc_map =read.table("../data/design.txt",header = T, row.names = 1)
# 物种分类文件，由qiime summarize_taxa.py生成，详见扩增子分析流程系列
# 本研究以纲水平进行训练，其实各层面都可以，具体那个层面最优，需要逐个测试寻找。推荐纲、科，不建议用OTU，差异过大
otu_table =read.table("../data/otu_table_tax_L3.txt",header = T, row.names = 1)
# 筛选品种作为训练集
sub_map = tc_map[tc_map$genotype %in% c("A50"),] # ,"IR24"
# 筛选OTU
idx = rownames(sub_map) %in% colnames(otu_table)
sub_map = sub_map[idx,]
sub_otu = otu_table[, rownames(sub_map)]   

## 随机森林回归
library(randomForest)
set.seed(315)
rf = randomForest(t(sub_otu), sub_map$day, importance=TRUE, proximity=TRUE, ntree = 1000)
print(rf)

## 交叉验证选择Features
set.seed(315) # 随机数据保证结果可重复，必须
# rfcv是随机森林交叉验证函数：Random Forest Cross Validation
result = rfcv(t(sub_otu), sub_map$day, cv.fold=10)
# 查看错误率表，23时错误率最低，为最佳模型
result$error.cv
# 绘制验证结果 
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

# 导出feature重要性
imp= as.data.frame(rf$importance)
imp = imp[order(imp[,1],decreasing = T),]
head(imp)
write.table(imp,file = "importance_class.txt",quote = F,sep = '\t', row.names = T, col.names = T)
# 简单可视化
varImpPlot(rf, main = "Top 23 - Feature OTU importance",n.var = 25, bg = par("bg"),
           color = par("fg"), gcolor = par("fg"), lcolor = "gray" )


## ggplot2美华feature贡献度柱状图

# 软件内部的varImpPlot可以快速可视化贡献度，简单全面，但发表还是要美美哒，美是需要代码的，就是花时间
# 基本思路同绘制Top 23 feature柱状图，按门着色，简化纲水平名字

# 读取所有feature贡献度
imp = read.table("importance_class.txt", header=T, row.names= 1, sep="\t") 
# 分析选择top23分组效果最好
imp = head(imp, n=23)
# 反向排序X轴，让柱状图从上往下画
imp=imp[order(1:23,decreasing = T),]

# imp物种名分解
# 去除公共部分
imp$temp = gsub("k__Bacteria;p__","",rownames(imp),perl=TRUE) 
# 提取门名称
imp$phylum = gsub(";[\\w-\\[\\]_]+","",imp$temp,perl=TRUE) # rowname unallowed same name
imp$phylum = gsub("[\\[\\]]+","",imp$phylum,perl=TRUE) 
# 提取纲名称
imp$class = gsub("[\\w\\[\\];_]+;c__","",imp$temp,perl=TRUE)  
imp$class = gsub("[\\[\\]]+","",imp$class,perl=TRUE)
# 添加纲level保持队形
imp$class=factor(imp$class,levels = imp$class)

# 图4.1. 绘制物种类型种重要性柱状图

p=ggplot(data = imp, mapping = aes(x=class,y=X.IncMSE,fill=phylum)) + 
  geom_bar(stat="identity")+coord_flip()+main_theme
p
ggsave(paste("rf_imp_feature",".pdf", sep=""), p, width = 4, height =4)



# 图4.2. 绘制时间序列热图

# 加载热图绘制包
library(pheatmap)

# 数据筛选23个feature展示
sub_abu = sub_otu[rownames(imp),]

# 简化名字
rownames(sub_abu)=imp[rownames(sub_abu),"class"]

# 直接自动聚类出图
pheatmap(sub_abu, scale = "row")
# 保存结果
pheatmap(sub_abu, scale = "row", filename = "heatmap_samples.pdf", width = 5, height = 5)


# 按时间为组合并均值
sampFile = as.data.frame(sub_map$day2,row.names = row.names(sub_map))
colnames(sampFile)[1] = "group"
mat_t = t(sub_abu)
mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
otu_norm_group = do.call(rbind, mat_mean)[-1,]
colnames(otu_norm_group) = mat_mean$group
pheatmap(otu_norm_group,scale="row",cluster_cols = F, cluster_rows = T)
pheatmap(otu_norm_group, scale="row",cluster_cols = F, cluster_rows = T, filename = "heatmap_groups.pdf", width = 5, height = 5)



## 求每组最大值
bak=otu_norm_group

otu_norm_group = otu_norm_group[as.character(imp$class),] # 按初始排序

for (i in 1:length(rownames(otu_norm_group))) {
#  i=1
  x=as.data.frame(sort(otu_norm_group[i,],decreasing = T))
  imp[i,"order"]=rownames(x)[1]
}
library(dplyr)
imp$order2 =  as.numeric(gsub("A50Sz","",imp$order,perl=TRUE) )
taxonomy = arrange(imp, desc(order2), class)

otu_norm_group1 = otu_norm_group[match(taxonomy$class,rownames(otu_norm_group)),] # 按初始排序

# 按初始排序
pheatmap(otu_norm_group1,scale="row",cluster_cols = F, cluster_rows = F)

pheatmap(otu_norm_group1,scale="row",cluster_cols = F, cluster_rows = F,filename ="fig4/pheatmap_order_all.pdf",width=8, height=4)

# 再用新顺序画样品
# X轴按从小到大排序
sub_abu1 = sub_abu[match(taxonomy$class,rownames(sub_abu)),] # 按初始排序
# Y轴按时间排序
sub_design2 = arrange(sub_design, day2, groupID)

sub_abu2 = sub_abu1[,match(sub_design2$Description,colnames(sub_abu1))] # 按初始排序

pheatmap(sub_abu2,scale="row",cluster_cols = F, cluster_rows = F)
# 重绘所有样品，效果也没有均值好
pheatmap(sub_abu2,scale="row",cluster_cols = F, cluster_rows = F,filename ="fig4/pheatmap_order_all_sample.pdf",width=8, height=4)
