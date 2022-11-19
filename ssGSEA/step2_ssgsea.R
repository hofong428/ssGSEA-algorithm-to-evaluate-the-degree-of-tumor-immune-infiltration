rm(list=ls())
library(pdftools)
options(stringsAsFactors = F)

b <- pdf_text('SupplementaryTables.pdf')

tmp = unlist(lapply(20:36, function(i){
  trimws(strsplit(b[[i]],split = '\n')[[1]])
}))
tmp=tmp[-c(1,2)]
library(stringr)
tmp=do.call(rbind,lapply(str_split(tmp,' '), function(x){
  x=x[-length(x)]
  gene_name<- x[1]
  cell_type<- trimws(paste(x[-1],collapse = ' '))
  return(c(gene_name,cell_type))
}))
immune_list <- split(tmp[,1],tmp[,2])

geneset_substract<- function(tmp){split_to_line<- gsub('\r','',strsplit(tmp,split = '\n')[[1]])
              gene_name<- apply(data.frame(split_to_line),1,function(x){ line<- strsplit(x,split=' ')[[1]]
                                                   pos<- grep('[A-Za-z\\d+]|\\d+',line)
                                                   res <- line[pos[1]]})
              cell_type<- apply(data.frame(split_to_line),1,function(x){ line<- strsplit(x,split=' ')[[1]]
                                                           pos<- grep('[A-Za-z\\d+]|\\d+',line)
                                                           res <- line[pos]
                                                           res <- res[c(-1,-length(res))]
                                                           s <- ''
                                                           for (i in 1:length(res)){
                                                             s<- paste(s,res[i])}
                                                           return(s)})
              result<- data.frame(gene_name,cell_type)
              return(result)
              }
gene_set <- data.frame()
for(i in 20:36){
  gene_set<- rbind(gene_set,geneset_substract(b[i]))
}
gene_set <- gene_set[c(-1,-2),]
list <- list()
for(i in 1:length(unique(gene_set$cell_type))){
  list[[i]] <- gene_set$gene_name[gene_set$cell_type== (unique(gene_set$cell_type)[i])]
}
names(list)<- unique(gene_set$cell_type)
save(list,file='gene_set.Rdata')
