library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(AUCell)
library(clusterProfiler)
library(org.Mm.eg.db)
library(devtools)
library(RColorBrewer)
library(tidydr)
library(reshape2)
library(paletteer) 
library(usethis)
library(ggstatsplot)
options(Seurat.object.assay.version = "v4")
#导入数据
setwd("E:/EXP document/single_cell_data")
PBS_basal.data <- Read10X(data.dir = c("PBS_basal"))
Py_basal.data <- Read10X(data.dir = c("Py_basal"))
PBS_clp.data <- Read10X(data.dir = c("PBS_clp"))
Py_clp.data <- Read10X(data.dir = c("Py_clp"))

#创建seurat对象
PBS_basal.data <- CreateSeuratObject(PBS_basal.data,project = "PBS_basal.data",assay = "RNA", min.cells = 3, min.features = 200)
Py_basal.data <- CreateSeuratObject(Py_basal.data ,project = "Py_basal.data",assay = "RNA", min.cells = 3, min.features = 200)
PBS_clp.data <- CreateSeuratObject(PBS_clp.data,project = "PBS_clp.data",assay = "RNA", min.cells = 3, min.features = 200)
Py_clp.data <- CreateSeuratObject(Py_clp.data,project = "Py_clp.data",assay = "RNA", min.cells = 3, min.features = 200)

#合并seurat对象
data <- merge(PBS_clp.data, add.cell.ids = c("PBS_clp"))
data[['RNA']] = as(object = data[['RNA']], Class = "Assay")
data <- NormalizeData(data)

#颜色设置
mycols <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
            '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
            '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
            '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
            '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
            '#968175')

#数据质量控制
data[["percent.mt"]] <- PercentageFeatureSet(data,pattern = "^mt-")
VlnPlot(data,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

#标准化
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(data), 10)
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
Top10_VaiableFeature <- plot2
Top10_VaiableFeature

#降维
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
data <- RunPCA(data, features = VariableFeatures(object = data),nfeatures.print = 10)
print(data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca")
DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(data)

#细胞分群
data <- data %>% 
  RunUMAP(reduction = "pca", dims = 1:6) %>% 
  FindNeighbors(reduction = "pca", dims = 1:6) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
DimPlot(data, label = T, pt.size = 1) + 
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 0,hjust = 0.01)) +NoLegend()

DimPlot(data, label = F,pt.size = 1,split.by = "orig.ident",label.size = 4,ncol = 2) + 
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 0,hjust = 0.01)) +theme(legend.position = "right")
#找富集基因
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25)
top5 <- data.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(data, features = top5$gene) + NoLegend()
#细胞分群

Idents(data) <- factor(Idents(Myeloid), levels= My_levels)
B_cell = DotPlot(data,features = c("Cd19", "Cd79a","Ms4a1","Ly6d"))
T_cell = DotPlot(data,features = c("Cd3g","Cd3e","Cd8a","Cd8b1","Cd4","Ctla4","Klrb1","Nkg7"))
T_cell
B_cell
Neutrophils = DotPlot(data,features = c("S100a8","S100a9","Mmp9","Ly6g","Elane","Cxcr2","Il1b","Csf3r","Il6"))
Neutrophils
Monocytes = DotPlot(data,features = c("Ccr2","Itgam","Ly6c2","Crip1","Lgals3","Cd14","S100a10","Il1b","Csf2","Csf3","Csf1"))
Monocytes
Basophil = DotPlot(data,features = c("Fcer1a","Cd63","Cd49b","Singlecf","Epx"))
Basophil

data$cell_type <- factor(data$cell_type,levels = c("B_cells","Pre_B_cells","Cd4+_T_cell","Cd4+_T_cell","Nk_T_cells","Nk_cells",
                                               "Neutrophils","Immature_neutrophils","Ly6c_hi_monocytes","Macrophage","Unknow"))

My_levels <- c("B_cells","Pro_B_cells","Cd4+_T_cell","Cd8+_T_cell","Nk_T_cells","Nk_cells","Neutrophils","Ly6c_hi_monocytes","Ly6c_low_monocytes", "Macrophage","Un1","Un2") 
Idents(data) <- factor(Idents(data), levels= My_levels)
data$cell_type <- factor(data$cell_type,labels = c("B_cells","Pro_B_cells","Cd4+_T_cell","Cd8+_T_cell","Nk_T_cells","Nk_cells","Neutrophils","Ly6c_hi_monocytes","Ly6c_low_monocytes", "Macrophage","Un1","Un2"))
DotPlot(data,features = c("Cd19", "Cd79a","Ms4a1","Ly6d",
                                         "Cd3g","Cd3e","Cd8b1","Cd4","Ctla4","Klrb1","Nkg7","Klrc1",
                                         "S100a8","S100a9","Ly6g",
                                         "Ccr2","Ly6c2","Crip1", "Cd14","S100a10",
                                        "S100a4","Ear2","Cst3",
                                        "Cx3cr1"))+  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签
  labs(x="cell_type",y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  viridis::scale_color_viridis(option="D")

DotPlot(data,features = c("Il6"))
#定群
new.cluster.ids <- c("B_cells","B_cells","Neutrophils","T_cell","T_cell","Neutrophils","T_cell","nk","B_cells","monocytes",
                      "T_cell","T_cell","UN1","monocyte","UN2","un3","UN4")
names(new.cluster.ids) <- levels(data)
data<- RenameIdents(data, new.cluster.ids)
data[['cell_type']] = unname(new.cluster.ids[data@meta.data$seurat_clusters])
data@meta.data$seurat_clusters = data@meta.data$cell_type  
data$seurat_clusters <- data@active.ident
#颜色选择
display.brewer.all()
mycolor2<-brewer.pal(14, "Paired")
mycolor3<- brewer.pal(11,"PRGn")
mycolor4<- brewer.pal(12,"Set3")

p1 = DimPlot(data, label = T, reduction = 'umap',pt.size = 0.5) + 
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 0,hjust = 0.01)) +NoLegend()

p2 = DimPlot(data, label = F,reduction = 'umap',pt.size = 0.5,split.by = "orig.ident",label.size = 4,cols = mycolor2,ncol = 2) + 
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 0,hjust = 0.01)) +theme(legend.position = "right")
p1+p2
#细胞比例
table(data$orig.ident)
prop.table(table(Idents(data)))
table(Idents(data), data$orig.ident)

Cellratio <- prop.table(table(Idents(data), data$orig.ident), margin = 2)

Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = "black")+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")) + scale_fill_manual(values = mycols)
 
#细胞因子
VlnPlot(data,features = c("Il1b","Tnf"),group.by = "cell_type",ncol = 1,pt.size =0.1,split.by = "orig.ident")+
                        theme(legend.position = "right")
VlnPlot(data,features = c("Il6","Csf1"),group.by = "cell_type",ncol = 1,pt.size =0.1,split.by = "orig.ident")+
  theme(legend.position = "right")
VlnPlot(data,features = c("Csf2","Csf3"),group.by = "cell_type",ncol = 1,pt.size =0.1,split.by = "orig.ident")+
  theme(legend.position = "right")
VlnPlot(data,features = c("Ccl2","Ccl4"),group.by = "cell_type",ncol = 1,pt.size =0.1,split.by = "orig.ident")+
  theme(legend.position = "right")
VlnPlot(data,features = c("Tnf","Ccl3"),group.by = "cell_type",ncol = 1,pt.size =0.1,split.by = "orig.ident")+
  theme(legend.position = "right")
FeaturePlot(data,features = 
              c("Ifng","Gzmb","Tnf"),split.by = "orig.ident")
#基因富集分析

BM.marker <- FindAllMarkers(data)
jjVolcano(diffData = BM.marker,
          tile.col =mycols,
          size  = 3.5,
          fontface = 'italic',
          topGeneN = 10)

data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25)
top5<- data.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(data, features = top5$gene,
          group.by = "cell_type",
          assay = "RNA",
          group.colors = mycols) +
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+
  theme(legend.position = "right")


#Go分析total

cluster <- subset(data,ident = "Neutrophils")
  cluster1 <- FindMarkers(cluster,min.pct = 0.1,
                              logfc.threshold = 0.25,
                              group.by = "orig.ident",
                              ident.1 = "PBS_clp.data",
                              ident.2 = "PBS_basal.data")
  up <-rownames(cluster1[intersect(which(cluster1 [,1]<0.05),which(cluster1 [,2]>=0.25)),])
 # down <-rownames(diff_Ly6c_hi_monocytes_cluster[intersect(which(diff_Ly6c_hi_monocytes_cluster [,1]<0.05),which(diff_Ly6c_hi_monocytes_cluster [,2]<=(-0.25))),])
  gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  head(gs)
  ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
  head(ego.bp)
  write.csv(ego.bp, file ="PBS Neutrophils_go.csv")
  pdf("PBS Neutrophils", height =12, width = 8)
  print(dotplot(ego.bp, showCategory=20,title="PBS Neutrophils GoTerm"))
  dev.off()
  

Neutrophil_cluster <- subset(data,ident = "Neutrophils")
diff_Neutrophil_cluster <- FindMarkers(Neutrophil_cluster,min.pct = 0.1,
                                              logfc.threshold = 0.25,
                                              group.by = "orig.ident",
                                              ident.1 = "PBS_clp.data",
                                              ident.2 = "PBS_basal.data")
up <-rownames(diff_Neutrophil_cluster[intersect(which(diff_Neutrophil_cluster [,1]<0.05),which(diff_Neutrophil_cluster [,2]>=0.25)),])
down <-rownames(diff_Neutrophil_cluster[intersect(which(diff_Neutrophil_cluster [,1]<0.05),which(diff_Neutrophil_cluster [,2]<=(-0.25))),])
gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(gs)
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
head(ego.bp)
write.csv(ego.bp, file ="PBSNeutrophils_up_go.csv")
pdf("diff_PBneutrophil_clusterf", height =5, width = 10)
print(dotplot(ego.bp, showCategory=10,title="PBS_clpVsPBSbasal_Neutrophil_gene GoTerm"))
dev.off()
getwd()

diff$Cluster <- factor(diff$Cluster,levels = c("B_cells","Pre_B_cells","Cd4+_T_cell","Cd8+_T_cell","Nk_T_cells","Nk_cells",
                                                   "Neutrophils","Immature_neutrophils","Ly6c_hi_monocytes","Macrophage","Unknow"))
ggplot(diff,aes(x = Cluster,y = Description,colour = pvalue))+
  geom_point(aes(color = -log10(pvalue),
                 size = Count))+
  viridis::scale_color_viridis(option="D")+
  xlab("cell_type")+
  theme(axis.text = element_text( size = 10),
        axis.text.x = element_text(size = 10,angle =45,hjust = 1))+
  theme(panel.grid = element_blank())+
  theme(panel.background = element_rect(
    fill = "white", color = "black", size = 0.5))


#DEGs
Fibroblast <- subset(scedata, celltype=="Fibroblast")
diff_Fibroblast <- FindMarkers(Fibroblast, min.pct = 0.25, 
                    logfc.threshold = 0.25,
                    group.by = "group",
                    ident.1 ="GM",
                    ident.2="BM")

 EnhancedVolcano(diff_Fibroblast,
                lab = rownames(diff_Fibroblast),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'diff_Fibroblast')

                
#GO分析



cluster0.markers1 <- FindMarkers(CLP_monocyte, ident.1 = "0", ident.2 = "1", min.pct = 0.25)
up <-rownames(cluster0.markers1[intersect(which(cluster0.markers [,1]<0.05),which(cluster0.markers[,2]>=0.25)),])
down <-rownames(cluster0.markers1[intersect(which(cluster0.markers [,1]<0.05),which(cluster0.markers [,2]<=(-0.25))),])
gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(gs)
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
head(ego.bp)
write.csv(ego.bp, file ="cluster0.markers1_up_go.csv")
pdf("cluster0.markers1_up_go.pdf", height =10, width = 8)
print(dotplot(ego.bp, showCategory=20,title="Cluster0 Vs.Cluster2 up gene GoTerm"))
dev.off()

table(data$orig.ident)
prop.table(table(Idents(data)))
table(Idents(data), data$orig.ident)

Cellratio <- prop.table(table(Idents(data), data$orig.ident), margin = 2)

Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = "black")+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

Ly6c_hi_monocytes <- subset(data,ident="Ly6c_hi_monocytes",)
Ly6c_hi_monocytes <- NormalizeData(Ly6c_hi_monocytes, normalization.method = "LogNormalize", scale.factor = 10000)
Ly6c_hi_monocytes<- FindVariableFeatures(Ly6c_hi_monocytes, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Ly6c_hi_monocytes), 10)
plot1 <- VariableFeaturePlot(Ly6c_hi_monocytes)
plot2 <- LabelPoints(plot = plot1, points = top10)
Top10_VaiableFeature <- plot2
Top10_VaiableFeature
plot1+plot2

#CLP单独分析
CLP_data <- merge(PBS_clp.data,c(Py_clp.data ), add.cell.ids = c("PBS_clp","Py_clp"))
CLP_data[['RNA']] = as(object = CLP_data[['RNA']], Class = "Assay")
CLP_data <- NormalizeData(CLP_data)

CLP_data[["percent.mt"]] <- PercentageFeatureSet(CLP_data,pattern = "^mt-")
VlnPlot(CLP_data,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,cols = mycols)
plot1 <- FeatureScatter(CLP_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CLP_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

CLP_data <- subset(CLP_data , subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)
CLP_data  <- NormalizeData(CLP_data , normalization.method = "LogNormalize", scale.factor = 10000)
CLP_data  <- FindVariableFeatures(CLP_data , selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CLP_data)
CLP_data<- ScaleData(CLP_data, features = all.genes)
CLP_data <- RunPCA(CLP_data, features = VariableFeatures(object = CLP_data),nfeatures.print = 10)
print(CLP_data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(CLP_data, dims = 1:2, reduction = "pca")
DimPlot(CLP_data, reduction = "pca")
DimHeatmap(CLP_data, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(CLP_data)

CLP_data <- CLP_data %>% 
  RunUMAP(reduction = "pca", dims = 1:6) %>% 
  FindNeighbors(reduction = "pca", dims = 1:6) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
DimPlot(CLP_data, label = T, pt.size = 1,cols = mycols) + 
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 0,hjust = 0.01)) +theme(legend.position = "right")

DimPlot(CLP_data, label = T,pt.size = 1,split.by = "orig.ident",label.size = 4,ncol = 2,cols = mycols) + 
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 0,hjust = 0.01)) +theme(legend.position = "right")

data.markers <- FindAllMarkers(CLP_data, only.pos = TRUE, min.pct = 0.25)
top5 <- data.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(CLP_data, features = top5$gene) + NoLegend()

Clp_Total_marker = DotPlot(CLP_data,features = c("Cd19", "Cd79a","Ms4a1","Ly6d",
                                         "Cd3g","Cd3e","Cd8a","Cd8b1","Cd4","Ctla4","Klrb1","Nkg7",
                                         "S100a8","S100a9","Mmp9","Ly6g",
                                         "Ccr2","Ly6c2","Crip1", "Cd14","S100a10",
                                         "S100a4","Ear2","Adgre180","Cst3",
                                         "Cd11c","H2","Cd74","Fscn1","Ccl5","Irf8","Fcgr3","Cx3cr1"))
Clp_Total_marker
new.cluster.ids <- c("B_cells","Cd4+_T_cell","Inflammatory_neutrophils","Basal_neutrophils","Inflammatory_neutrophils","B_cells","Cd8+_T_cell",
                     "Inflammatory_neutrophils", "Inflammatory_monocytes","Inflammatory_neutrophils","Nk_T_cell","Nk_cell","Un1", "Un2","Pro_B_cell","Macrophage")
names(new.cluster.ids) <- levels(CLP_data)
CLP_data <- RenameIdents(CLP_data, new.cluster.ids)
CLPcelltype <- c("0" = "B_cells",
                      "1" = "Cd4+_T_cell",
                      "2" ="Inflammatory_neutrophils",
                      "3" = "Basal_neutrophils",
                      "4"= "Inflammatory_neutrophils", 
                      "5"= "B_cells",
                      "6"= "Cd8+_T_cell", 
                      "7"= "Inflammatory_neutrophils", 
                      "8"= "Inflammatory_monocytes",
                      "9"= "Inflammatory_neutrophils",
                      "10"= "Nk_T_cell",
                      "11"= "Nk_cell",
                      "12"= "Un1",
                      "13"= "Un2",
                      "14"= "Pro_B_cell",
                      "15"= "Macrophage")
CLP_data[['cell_type']] = unname(CLPcelltype[CLP_data@meta.data$seurat_clusters])
CLP_data@meta.data$seurat_clusters = CLP_data@meta.data$cell_type  

DimPlot(CLP_data, reduction = 'umap',pt.size = 1,cols = mycols,label = T) + 
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 0,hjust = 0.01)) +theme(legend.position = "right")

DimPlot(CLP_data, label = T,reduction = 'umap',pt.size = 1,split.by = "orig.ident",label.size = 4,ncol = 2,cols = mycols) + 
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 0,hjust = 0.01)) +theme(legend.position = "right")

VlnPlot(CLP_data,features = c("Il1b","Il6","Ccl5","Csf1","Csf2","Csf3","Ccl2","Ccl3","Tnf"),group.by = "cell_type",ncol = 3,pt.size = 0.1,split.by = "orig.ident")   
monocytes <- subset(PBMC, ident= "Ly6c_hi_monocytes")
VlnPlot(monocytes,features = "Ccr5",group.by = "seurat_clusters",split.by = "orig.ident",pt.size = 0)

table(CLP_data$orig.ident)
prop.table(table(Idents(CLP_data)))
table(Idents(CLP_data), CLP_data$orig.ident)

Cellratio <- prop.table(table(Idents(CLP_data), CLP_data$orig.ident), margin = 2)

Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
p3 <- ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = "black")+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
p3+scale_fill_manual(values = mycols)

#Monocyte分析

CLP_monocyte <-subset(CLP_data,ident="Inflammatory_monocytes")
CLP_monocyte  <- NormalizeData(CLP_monocyte , normalization.method = "LogNormalize", scale.factor = 10000)
CLP_monocyte <- FindVariableFeatures(CLP_monocyte , selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(CLP_monocyte ), 10)
plot1 <- VariableFeaturePlot(CLP_monocyte )
plot2 <- LabelPoints(plot = plot1, points = top10)
Top10_VaiableFeature <- plot2
Top10_VaiableFeature
plot1+plot2

all.genes <- rownames(CLP_monocyte)
CLP_monocyte <- ScaleData(CLP_monocyte, features = all.genes)
CLP_monocyte <- RunPCA(CLP_monocyte, features = VariableFeatures(object = CLP_monocyte),nfeatures.print = 10)
print(CLP_monocyte[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(CLP_monocyte, dims = 1:2, reduction = "pca")
DimPlot(CLP_monocyte, reduction = "pca")
DimHeatmap(CLP_monocyte, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(CLP_monocyte)

CLP_monocyte <- CLP_monocyte %>% 
  RunUMAP(reduction = "pca", dims = 1:6) %>% 
  FindNeighbors(reduction = "pca", dims = 1:6) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
DimPlot(CLP_monocyte, label = T, pt.size = 2) + 
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inche=812s"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 0,hjust = 0.01)) +theme(legend.position = "right")

DimPlot(CLP_monocyte, label = T,pt.size = 2,split.by = "orig.ident",label.size = 4,ncol = 2) + 
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 0,hjust = 0.01)) +theme(legend.position = "right")

FeaturePlot(CLP_monocyte,features =c("Il1b","Il6","Tnf","Csf1","Csf3","Ccl2") ,split.by = "orig.ident")

CLP_data.markers <- FindAllMarkers(CLP_monocyte, only.pos = TRUE, min.pct = 0.25)
top10 <- CLP_data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(CLP_monocyte, features = top10$gene) + NoLegend()

cluster0.markers1 <- FindMarkers(CLP_monocyte, ident.1 = "0", ident.2 = "1", min.pct = 0.25)
up <-rownames(cluster0.markers1[intersect(which(cluster0.markers [,1]<0.05),which(cluster0.markers[,2]>=0.25)),])
down <-rownames(cluster0.markers1[intersect(which(cluster0.markers [,1]<0.05),which(cluster0.markers [,2]<=(-0.25))),])
gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(gs)
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
head(ego.bp)
write.csv(ego.bp, file ="cluster0.markers1_up_go.csv")
pdf("cluster0.markers1_up_go.pdf", height =10, width = 8)
print(dotplot(ego.bp, showCategory=20,title="Cluster0 Vs.Cluster2 up gene GoTerm"))
dev.off()

FeaturePlot(CLP_monocyte,features = c("Ly6c2","Ccr2","Il1b","Fcgr1","Ccl2","Nfkb1","Relb","Ccl3","Cxcr2","Cxcr3"))


test <- read.csv("clpmonocytetop40.csv")
test_2 <-test[order(test[,4],decreasing = T),]
test_2$Description <- factor(test_2$Description,levels = c(rev(test_2$Description)))
max <- max(test_2 $GeneRatio)
min <- min(test_2 $GeneRatio)
p = ggplot(test_2,aes(x =GeneRatio, y = Description))+
  geom_point(aes(size = Count,color = -log10(p.adjust)))+
  theme_bw()+
  scale_colour_gradient(low ="#443B84FF", high = "#B9DE28FF")+
  xlim(min,max)+
  scale_y_discrete(labels = function(x) str_wrap(x,width = 40))+
  labs(x = "GeneRatio",y = "",title = "Dotplot",
       color = expression(-log10(p.adjust)),size = "Count")+
  theme(axis.title = element_text(size = 12),axis.text = element_text(size = 10),
        plot.title = element_text(size = 1,hjust = 0.5,face = "bold"),
        legend.title = element_text(size = 9),legend.text = element_text(size = 10))

#小提琴堆叠图
Myeloid$cell_type <- factor(Myeloid$cell_type,labels = c(
                                                   "Neutrophils","Immature_neutrophils","Ly6c_hi_monocytes","Macrophage"))

Idents(Myeloid) <- factor(Idents(Myeloid), levels= My_levels)

Myeloid <- subset(data,idents = c("Macrophage","Neutrophils","Immature_neutrophils","Ly6c_hi_monocytes"))
My_levels <- c("Neutrophils","Immature_neutrophils","Ly6c_hi_monocytes","Macrophage")                    
Myeloid <- factor(Myeloid$cell_type, ordered  = My_levels)
?factor()
VlnPlot(Myeloid, features = c("Il6","Il2","Ccl3","Ccl2","Ccl4","Tnf","Il1b","Il18","Il1a","Il4","Il10","Ifng","Il17f","Crp","Il27"), stack = TRUE, sort = TRUE,flip = TRUE,split.by  = "orig.ident",cols =  c("#57C3F3","#E95C59")) +
  theme(legend.position = "top") +scale_color_manual(values = c("#57C3F3","#E95C59"))
 
VlnPlot(Myeloid, features = c("Trim30a","Lbp","Trim5","Clec4e","Gbp3","Ifi35","Ncf1","Il1r2","Syk","Zbp1","Il4ra","Ctr9","Gab1","Jak1","Gpr35"), stack = TRUE, sort = TRUE,flip = TRUE,split.by  = "orig.ident",cols =  c("#57C3F3","#E95C59")) +
  theme(legend.position = "top") +scale_color_manual(values = c("#57C3F3","#E95C59"))


 
#GGPLOT2多细胞信号通路  
ggplot(top10,aes(x = Cluster,y = Description,colour = pvalue))+
  geom_point(aes(color = -log10(pvalue),
                 size = Count))+
  viridis::scale_color_viridis(option="E")+
  xlab("Cluster")+
  theme(axis.text = element_text( size = 10),
        axis.text.x = element_text(size = 10,angle = 45,hjust = 1))+
  theme(panel.grid = element_blank())+
  theme(panel.background = element_rect(
    fill = "white", color = "black", size = 0.5))



# Neutrophils
Neutrophis <- subset(data,ident = "Neutrophils")
all.genes <- rownames(Neutrophis)
Neutrophis <- ScaleData(Neutrophis, features = all.genes)
Neutrophis <- RunPCA(Neutrophis, features = VariableFeatures(object = Neutrophis),nfeatures.print = 10)
print(Neutrophis[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Neutrophis, dims = 1:2, reduction = "pca")
DimPlot(Neutrophis, reduction = "pca")
DimHeatmap(Neutrophis, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(Neutrophis)

#细胞分群
Neutrophis <- Neutrophis %>% 
  RunUMAP(reduction = "pca", dims = 1:6) %>% 
  FindNeighbors(reduction = "pca", dims = 1:6) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
new.cluster.ids <- c("Neu4","Neu1","Neu8","Neu6","Neu5","Neu3","Neu7","Neu4")
names(new.cluster.ids) <- levels(Neutrophis)
Neutrophis<- RenameIdents(Neutrophis, new.cluster.ids)
Neutrophis[['cell_type']] = unname(new.cluster.ids[Neutrophis@meta.data$seurat_clusters])
Neutrophis@meta.Neutrophis$seurat_clusters = Neutrophis@meta.data$cell_type  
Neutrophis$seurat_clusters <- Neutrophis@active.ident

p1 = DimPlot(Neutrophis, label = T, reduction = 'umap',pt.size = 0.5) + 
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 0,hjust = 0.01)) +NoLegend()

p2 = DimPlot(Neutrophis, label = F,reduction = 'umap',pt.size = 0.5,split.by = "orig.ident",label.size = 4,ncol = 2) + 
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 0,hjust = 0.01)) +theme(legend.position = "right")
p1+p2
#AUCell

mmu <- read.gmt("mh.all.v2023.2.Mm.symbols.gmt")
cells_rankings <- AUCell_buildRankings(Neutrophis@assays$RNA@data)
unique(mmu$term)
geneSets <- lapply(unique(mmu$term), function(x){mmu$gene[mmu$term == x]})
names(geneSets) <- unique(mmu$term)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
length(rownames(cells_AUC@assays@data$AUC))
grep("HALLMARK_INFLAMMATORY_RESPONSE",rownames(cells_AUC@assays@data$AUC),value = T)

geneSet <- "HALLMARK_INFLAMMATORY_RESPONSE"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
Neutrophis$AUC <- aucs
df<- data.frame(Neutrophis@meta.data,Neutrophis@reductions$umap@cell.embeddings)
head(df)

class_avg <- df %>%
  group_by(cell_type) %>%
  summarise(
    umap_1 = median(umap_1),
   umap_2 = median(umap_2)
  )
ggplot(df, aes(umap_1, umap_2))+
  geom_point(aes(colour = AUC)) + 
  viridis::scale_color_viridis(option="A") +
  ggrepel::geom_label_repel(aes(label = cell_type),
                            data = class_avg,
                            size = 4,
                            label.size = 0,
                            segment.color = NA)+
  theme(legend.position = "none") + 
  theme_bw()

Idents(data) <- factor(Idents(Myeloid), levels= My_levels)
FeaturePlot(data,features = "Xkr4")

Neutrophis <- subset(PBMC, ident= "Nk_cells")
exprSet <- Neutrophis@assays[["RNA"]]@data  
exprSet<-as.data.frame(t(exprSet))

#基因相关性分析
y <- as.numeric(exprSet[,"Klrg1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)
for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")
View(cor_data_df)

cor_data_sig_pos <- cor_data_df %>%
  filter(pvalue < 0.01) %>% filter(correlation > 0)%>%
  arrange(desc(correlation))
cor_data_sig_neg <- cor_data_df %>%
  filter(pvalue < 0.01) %>% filter(correlation < 0)%>%
  arrange(desc(abs(correlation)))


ggscatterstats(data = exprSet,
               y = Itgam,
               x = Gzma,
               centrality.para = "mean",
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               marginal.type = "density" # #类型可以换成density,boxplot,violin,densigram
               )
