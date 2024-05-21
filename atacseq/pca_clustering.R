# last modified 2024-05-22
# Author: Feng Yan

## plot PCA and Hiearchical clustering using dba object (from DiffBind)

pca_clustering <- function(db, sampleid, prefix, top=nrow(db$binding), peakset=NULL){
  
  ## hcluster or pca
  if (is.null(peakset)) {
    idx=rowSds(db$binding[,3+sampleid]) %>% order(decreasing = T) %>% head(n=top)
  } else if (class(peakset)=="GRanges") {
    idx=findOverlaps(peakset, db$peaks[[1]] %>% makeGRangesFromDataFrame()) %>%
      subjectHits()
  }
  
  library(amap) ## load Dist function
  pdf(paste0("plot/hc_dend_atac_",prefix,".pdf"), width = 5, height = 5)
  
  db$binding[idx,3+sampleid] %>% t() %>% log2() %>%
    Dist(method = "pearson") %>% 
    hclust(method = "complete") %>% as.dendrogram() %>% 
    dendextend::set("labels_col", value = as.numeric(db$samples$Treatment)[sampleid], order_value=T) %>% 
    #dendextend::set("branches_col", value = as.numeric(db$samples$Treatment)[sampleid], order_value=T) %>%
    set("branches_lwd", 2) %>%
    place_labels(labels = db$samples %>% slice(sampleid) %>%
                   unite(label, Condition, Treatment) %>% pull(label)) %>%
    plot 
  dev.off()
  
  autoplot(db$binding[idx,3+sampleid] %>% t() %>% log2()  %>% prcomp(scale.=T) ,
           data=db$samples[sampleid,], shape = "Condition",
           colour="Treatment", size=2 , label=T, label.label="SampleID",label.repel=T
  )
  ggsave(paste0("plot/pca_atac_",prefix,".pdf"), width = 4, height = 3, units = "in")
}
