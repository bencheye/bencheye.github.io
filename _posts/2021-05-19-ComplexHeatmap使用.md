layout: post
title: "complexHeatmap热图包使用"
date: 2021-05-19
description: ""

tag: 
---   

```R
library(ComplexHeatmap)
# data
set.seed(123)
nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3
mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
)
mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
rownames(mat) = paste0("row", seq_len(nr))
colnames(mat) = paste0("column", seq_len(nc))
# 两个热图组合在一起
f1 <- colorRamp2(seq(min(mat), max(mat), length=3), c("blue","#EEEEEE", "red"))
f2 <- colorRamp2(seq(min(mat), max(mat), length=3), c("blue","#EEEEEE", "red"), space = "RGB")
H1 <- Heatmap(mat, col = f1, column_title = "LAB color space")
H2 <- Heatmap(mat, col = f2, column_title = "RGB color space")
H1+H2
# 设置颜色
Heatmap(mat2, col = colorRamp2(c(-3,0,3), c("green","white","red")), cluster_rows = FALSE, cluster_columns = FALSE)
Heatmap(mat, col = rev(rainbow(10)))
# 根据特定分组来分割热图
Heatmap(mat, name = "mat", row_split = yy, column_split = factor(rep(c("C", "D"), 12)))
Heatmap(mat, name = "mat", row_km = 3, row_gap = unit(2, "mm"))
# 突出重要基因
# 由于基因很多直接展示出来，根本看不清，我们可以强调几个标记基因。用到两个函数是rowAnnotation和anno_mark
mark_gene <- c("IL7R","CCR7","IL7R","S100A4","CD14","LYZ","MS4A1","CD8A","FCGR3A","MS4A7","GNLY","NKG7","FCER1A", "CST3","PPBP")
gene_pos <- which(rownames(mat) %in% mark_gene)
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, 
                                                 labels = mark_gene))
Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        top_annotation = top_anno,
        right_annotation = row_anno,
        column_title = NULL)
# 为了增加聚类注释，我们需要用到HeatmapAnnotation函数，它对细胞的列进行注释，
# 而rowAnnotation函数可以对行进行注释。这两个函数能够增加各种类型的注释，
# 包括条形图，点图，折线图，箱线图，密度图等等，这些函数的特征是anno_xxx，例如anno_block就用来绘制区块图。
top_anno <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = col), # 设置填充色
                       labels = levels(cluster_info), 
                       labels_gp = gpar(cex = 0.5, col = "white"))) # 设置字体
# 其中anno_block中的gp参数用于设置各类图形参数，labels设置标签，labels_gp设置和标签相关的图形参数。可以用?gp来了解有哪些图形参数。
Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        top_annotation = top_anno, # 在热图上边增加注释
        column_title = NULL ) # 不需要列标题
# 修改聚类树不同分组的颜色
library(dendextend)
row_dend = as.dendrogram(hclust(dist(mat)))
row_dend = color_branches(row_dend, k = 2) # `color_branches()` returns a dendrogram object
Heatmap(mat, name = "mat", cluster_rows = row_dend)
```

