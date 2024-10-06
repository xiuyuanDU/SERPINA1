setwd('F:/scRNA/workpath/5.cellchat')
thyroid.combined=readRDS("cellchat.rds")
Idents(thyroid.combined) <- thyroid.combined@meta.data[["stim"]]
af<-thyroid.combined
cellchat = createCellChat(object = af , group.by = "cellType")
levels(cellchat@idents) #展示以下现在的细胞分组
groupSize <- as.numeric(table(cellchat@idents)) #细胞亚群各组数量
groupSize
#设置配体受体交互数据库 
CellChatDB = CellChatDB.human 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
###选择使用的细胞通讯集合
#使用所有集合进行细胞通讯分析
CellChatDB.use = CellChatDB
#在对象中设置需要使用的细胞通讯集合 
cellchat@DB = CellChatDB.use

#预处理表达数据以进行细胞间通讯分析 
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 检查是否存在没有过表达基因的细胞群，如果有，则进一步处理
if(any(is.na(cellchat@netP$prob))) {
  cellchat <- removeCellGroups(cellchat, groups = which(is.na(cellchat@netP$prob)))
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
}
# 确保所有的因子水平是正确的，并删除未使用的因子水平
cellchat@idents <- droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents),unique(cellchat@idents)))
#细胞通信网络的推断 
cellchat <- computeCommunProb(object=cellchat,raw.use = TRUE)
#在信号通路级别推断细胞-细胞通信 
cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞通信网络 
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents)) #细胞亚群各组数量
groupSize
#汇总结果贝壳图 
pdf(file = "汇总结果贝壳图.pdf", width =5 , height = 10)
par(mfrow = c(1,2), xpd=TRUE)     #设置画板 
#交互的数量 
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
#交互的权重 
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#单一贝壳图 
mat <- cellchat@net$weight
par(mfrow = c(1,1), xpd=TRUE)        #设置画板，两行三列 
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#cellchat@netP$pathways 
#head(cellchat@LR$LRsig) 
#显示重要通信的信号都可以通过cellchat@netP$pathways. 
pathways.show <- cellchat@netP$pathways 
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()
#显示从某些细胞组（由“sources.use”定义）到其他细胞组（由“targets.use”定义）的所有重要相互作用 L-R 对
table(cellchat@idents) #简单查看细胞群序 
pdf(file = "Moderately.pdf", width =10 , height = 20)
netVisual_bubble(cellchat, sources.use = 16, targets.use = c(1:33), remove.isolate = FALSE) 
dev.off()
#########################
table(cellchat@idents) #简单查看细胞群序 
pdf(file = "Poorly.pdf", width =10 , height = 20)
netVisual_bubble(cellchat, sources.use = 21, targets.use = c(1:33), remove.isolate = FALSE)  
dev.off()
#########################
pdf(file = "state1.pdf", width =10 , height = 20)
table(cellchat@idents) #简单查看细胞群序 
netVisual_bubble(cellchat, sources.use = 22, targets.use = c(1:33), remove.isolate = FALSE) 
dev.off()
#########################
pdf(file = "state2.pdf", width =10 , height = 20)
table(cellchat@idents) #简单查看细胞群序 
netVisual_bubble(cellchat, sources.use = 23, targets.use = c(1:33), remove.isolate = FALSE) 
dev.off()
#########################
pdf(file = "state3.pdf", width =10 , height = 20)
table(cellchat@idents) #简单查看细胞群序 
netVisual_bubble(cellchat, sources.use = 24, targets.use = c(1:33), remove.isolate = FALSE)、
dev.off()
save.image("细胞通讯.RData")
load("细胞通讯.RData")
