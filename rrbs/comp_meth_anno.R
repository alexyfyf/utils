# last modified 2024-01-25
# Author: Feng Yan

## compute annotation of genomic ranges using annotatr (for CGI) and ChIPseeker (for genomic location)

comp_meth_anno <- function(gr, cpgs_info = cpgs_info, column = "mean.meth.diff", 
                           simplify = F, # don't split up and down
                           plot = "bar",
                           genome = "mm10"){
  library(annotatr)
  library(ChIPseeker)
  library(RColorBrewer)
  #library(pheatmap)
  library(cowplot)
  library(tidyverse)
  library(rlang)
  library(plyranges)
  
  seqlevelsStyle(gr) <- "UCSC"
  
  if (is.null(cpgs_info)) {
    library("annotatr")
    cpgs_info <- build_annotations(genome = genome, annotations = paste0(genome, "_cpgs"))
  }
  
  if (genome == "mm10") {
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    txdb=TxDb.Mmusculus.UCSC.mm10.knownGene
    org="org.Mm.eg.db"
    print(paste0('use txdb ', packageVersion('TxDb.Mmusculus.UCSC.mm10.knownGene'), 
                 'and orgdb ', packageVersion('org.Mm.eg.db')))
  } else if (genome == "hg38") {
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    txdb=TxDb.Hsapiens.UCSC.hg38.knownGene
    org="org.Hs.eg.db"
    print(paste0('use txdb ', packageVersion('TxDb.Hsapiens.UCSC.hg38.knownGene'), 
                 ' and orgdb ', packageVersion('org.Hs.eg.db')))
  } else if (genome == "hg19") {
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
    txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
    org="org.Hs.eg.db"
    print(paste0('use txdb ', packageVersion('TxDb.Hsapiens.UCSC.hg19.knownGene'), 
                 'and orgdb ', packageVersion('org.Hs.eg.db')))
  } else {
    print("Only hg38, hg19 and mm10 are currently supported.")
  }
  
  cpg_anno <- annotate_regions(regions = gr,
                               annotations = cpgs_info,
                               ignore.strand = TRUE,
                               quiet = FALSE)
  tss_anno <- annotatePeak(gr, TxDb = txdb, #version 3.10
                           annoDb = org, overlap = "TSS") ## version 3.16
  
  cpg_anno_sum <- cpg_anno %>% summarize_annotations()
  
  tss_anno <- tss_anno@anno %>% data.frame() %>% 
    mutate(annotation2=str_replace_all(annotation, 
                                       c("(Intron).*"="\\1",
                                         "(Exon).*"="\\1" ,
                                         "(Downstream).*"="\\1"
                                       )))
  tss_anno_sum <- tss_anno %>% dplyr::count(annotation2)  
  
  if (simplify) {
    
    colnames(cpg_anno_sum) <- c("Anno","total")
    colnames(tss_anno_sum) <- c("Anno","total")
    
    cpgpeaknum <- cpg_anno_sum[,-1] %>% as.matrix() %>% colSums()
    tsspeaknum <- tss_anno_sum[,-1] %>% as.matrix() %>% colSums()
    
    cpgtitle <- paste("total", cpgpeaknum, sep = "=")
    tsstitle <- paste("total", tsspeaknum, sep = "=")
    
  } else if (!is.null(column)) {
    
    cpg_anno_sum <- cpg_anno_sum %>% 
      left_join(cpg_anno %>% filter(!!sym(column) > 0) %>% summarize_annotations(),by=c("annot.type")) %>% 
      left_join(cpg_anno %>% filter(!!sym(column) < 0) %>% summarize_annotations(),by=c("annot.type")) %>%
      replace_na(list(n.y=0, n=0)) ## fix missing annotation type
    
    tss_anno_sum <- tss_anno_sum %>% 
      left_join(tss_anno %>% filter(!!sym(column) > 0) %>% dplyr::count(annotation2),by=c("annotation2")) %>% 
      left_join(tss_anno %>% filter(!!sym(column) < 0) %>% dplyr::count(annotation2),by=c("annotation2")) %>%
      replace_na(list(n.y=0, n=0))
    
    colnames(cpg_anno_sum) <- c("Anno","total","up","down")
    colnames(tss_anno_sum) <- c("Anno","total","up","down")
    
    cpgpeaknum <- cpg_anno_sum[,-1] %>% as.matrix() %>% colSums()
    tsspeaknum <- tss_anno_sum[,-1] %>% as.matrix() %>% colSums()
    
    cpgtitle <- paste(c("up","down","total"), cpgpeaknum[c(2,3,1)], sep = "=") %>% paste(., collapse = " ")
    tsstitle <- paste(c("up","down","total"), tsspeaknum[c(2,3,1)], sep = "=") %>% paste(., collapse = " ")
  }
  
  if (plot=="bar") {
    # The palette with grey:
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    p1 <- cpg_anno_sum %>% gather(key = "group", value = "number", 2:ncol(cpg_anno_sum)) %>%
      mutate(CGIanno=factor(Anno, 
                            levels = paste(genome, c("cpg_inter","cpg_shelves","cpg_shores", "cpg_islands"), sep = "_"),
                            labels = c("OpenSea","Shelves","Shores", "Islands"))) %>%
      mutate(group=factor(group, levels = c("up","down","total"))) %>% 
      ggplot(aes(x = group, y = number, fill = CGIanno)) + 
      geom_bar(stat = "identity", position = "fill") + 
      scale_fill_manual(values=cbPalette[c(1,2,8,7)], drop = F) +
      coord_flip() + 
      geom_text(aes(label=number),  vjust = 0.5, position = position_fill(vjust = 0.5)) +
      theme(legend.position="bottom") +
      ggtitle(cpgtitle)
    
    p2 <- tss_anno_sum %>% gather(key = "group", value = "number", 2:ncol(tss_anno_sum)) %>%
      mutate(TSSanno=factor(Anno, 
                            levels = c("3' UTR","5' UTR","Distal Intergenic", "Downstream","Exon","Intron",
                                       "Promoter (<=1kb)", "Promoter (1-2kb)","Promoter (2-3kb)"))
      ) %>%
      mutate(group=factor(group, levels = c("up","down","total"))) %>% 
      ggplot(aes(x = group, y = number, fill = TSSanno)) + 
      geom_bar(stat = "identity", position = "fill") + 
      scale_fill_manual(values=brewer.pal(9, 'Set1'), drop = F) +
      coord_flip() + 
      geom_text(aes(label=number),  vjust = 0.5, position = position_fill(vjust = 0.5)) +
      theme(legend.position="bottom") +
      ggtitle(tsstitle)
    
    p3 <- plot_grid(p1,p2, nrow = 2)
    
  } else if (plot=="pie") {
    
    p1 <- cpg_anno_sum %>% gather(key = "group", value = "number", 2:ncol(cpg_anno_sum)) %>%
      mutate(CGIanno=factor(Anno, 
                            levels = paste(genome, c("cpg_islands","cpg_inter","cpg_shelves","cpg_shores"), sep = "_"),
                            labels = c("Islands", "OpenSea","Shelves","Shores"))) %>%
      mutate(group=factor(group, levels = c("up","down","total"))) %>% 
      ggplot(aes(x = "", y = number, fill = CGIanno)) + 
      geom_bar(stat="identity", position="fill") + ## without position="fill", will generate absolute value not pct
      coord_polar(theta="y") + scale_fill_discrete(drop=F) + 
      facet_wrap(~group) + theme(legend.position="none") +
      ggtitle(cpgtitle)
    #geom_text(aes(label=number),  vjust = 1, position = position_fill(vjust = 0.5))
    
    p2 <- tss_anno_sum %>% gather(key = "group", value = "number", 2:ncol(tss_anno_sum)) %>%
      mutate(TSSanno=factor(Anno, 
                            levels = c("3' UTR","5' UTR","Distal Intergenic", "Downstream","Exon","Intron",
                                       "Promoter (<=1kb)", "Promoter (1-2kb)","Promoter (2-3kb)"))
      ) %>% ## in case drop of factors with 0 counts, and make color consistent
      mutate(group=factor(group, levels = c("up","down","total"))) %>% 
      ggplot(aes(x = "", y = number, fill = TSSanno)) + 
      geom_bar(stat="identity", position="fill") + ## without position="fill", will generate absolute value not pct
      coord_polar(theta="y") + scale_fill_discrete(drop=F) + 
      facet_wrap(~group) + theme(legend.position="none") +
      ggtitle(tsstitle)
    #geom_text(aes(label=number),  vjust = 2, position = position_fill(vjust = 0.5)) 
    
    p3 <- plot_grid(p1,p2, nrow = 2)
  }
  
  hm <- lapply(list(cpg_anno_sum, tss_anno_sum), function(x){
    if (ncol(x)==4) {
      x %>% mutate(total_pct=total/sum(total),
                   up_pct=up/sum(up),
                   down_pct=down/sum(down))
    } else if (ncol(x)==2) {
      x %>% mutate(total_pct=total/sum(total))
    }}) %>% 
    Reduce(rbind,.) 
  
  hmat <- as.matrix(hm[,-1])
  rownames(hmat) <- hm$Anno
  
  print(p3)
  
  return(list(cgi=cpg_anno, tss=tss_anno %>% makeGRangesFromDataFrame(keep.extra.columns = T),
              summary=hmat, plot=p3))
}

