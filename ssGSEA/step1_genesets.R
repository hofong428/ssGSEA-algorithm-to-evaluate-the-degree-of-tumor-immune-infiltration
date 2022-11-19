rm(list=ls())
###矩阵信息提取
a <- read.table('GSE112996_merged_fpkm_table.txt.gz',
                header = T,
                row.names=1)
raw_data<- a[,-1]
###表型信息提取
pheno <- read.csv(file = 'GSE112996_series_matrix.txt')
pheno <- data.frame(num1 = strsplit(as.character(pheno[42,]),split='\t')[[1]][-1],
                    num2 = gsub('patient: No.','P',strsplit(as.character(pheno[51,]),split='\t')[[1]][-1]))
table(pheno$num2)
####数据过滤
data<- a[!apply(raw_data,1,sum)==0,] 
####去除重复基因名的行，归一化
uni_matrix<- data[!duplicated(data$GeneName),]
rownames(uni_matrix) <- uni_matrix[,1]
uni_matrix <- uni_matrix[,-1]
uni_matrix <- log2(uni_matrix+1)
colnames(uni_matrix)<- gsub('X','',gsub('\\.','\\-',colnames(uni_matrix)))
uni_matrix<- uni_matrix[,order(colnames(uni_matrix))]
uni_matrix[1:4,1:4]
save(uni_matrix,file = 'uni_matrix.Rdata')

library(genefilter)
library(GSVA)
library(Biobase)
load('gene_set.Rdata')
gsva_matrix<- gsva(as.matrix(uni_matrix), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
library(pheatmap)
gsva_matrix1<- t(scale(t(gsva_matrix)))
gsva_matrix1[gsva_matrix1< -2] <- -2
gsva_matrix1[gsva_matrix1>2] <- 2
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}
nor_gsva_matrix1 <- normalization(gsva_matrix1)
annotation_col = data.frame(patient=pheno$num2)
rownames(annotation_col)<-colnames(uni_matrix)
bk = unique(c(seq(0,1, length=100)))
pheatmap(nor_gsva_matrix1,show_colnames = F,
         cluster_rows = F,cluster_cols = F,annotation_col = annotation_col,
         breaks=bk,cellwidth=5,cellheight=5,fontsize=5,filename = 'ssgsea.png')
save(gsva_matrix,gsva_matrix1,pheno,file = 'score.Rdata')
