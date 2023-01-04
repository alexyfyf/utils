# last modified 2022-12-16
# Author: Feng Yan

## plot correlation between logFC of DEGs and meth.diff of DMC or DMR

plot_meth_exp <- function(anno, deg, plot=TRUE, pval.DEG = 0.05){
  ## merge annoation df

  deg <- topTags(deg, n=Inf, p.value = pval.DEG)$table
  df <- anno %>% data.frame() %>% 
    mutate(annotation2=str_replace_all(annotation, 
                                       c("(Intron).*"="\\1",
                                         "(Exon).*"="\\1" ,
                                         "(Downstream).*"="\\1"))) %>%
    inner_join(deg, by=c("SYMBOL"="Gene.Name")) %>% 
    group_by(SYMBOL, geneId, annotation2) %>% 
    summarise(mean.meth.diff = mean(mean.meth.diff),
              logFC = mean(logFC), FDR = mean(FDR))
  ## if multiple regions match same gene, this will calculate mean mean.meth.diff across the same annotated region fot the gene

  
  ## filter on significant DEGs don't change the correlation
  
  ## get count in each segments
  df2<- df %>% mutate(logFC_cut=cut(logFC, breaks = c(-Inf,-1,1, Inf), 
                                    labels = c(-5,0,5)),
                      meanmethdiff_cut=cut(mean.meth.diff, breaks = c(-Inf,0, Inf), 
                                           labels = c(-90,90))) 
  df3 <- df2 %>% ungroup() %>%
    dplyr::count(logFC_cut, meanmethdiff_cut, annotation2) %>%
    complete(logFC_cut, meanmethdiff_cut, annotation2, fill=list(n=0)) %>%
    mutate(logFC_cut=as.numeric(as.character(logFC_cut)),
           meanmethdiff_cut=as.numeric(as.character(meanmethdiff_cut)))
  
  p <- df %>%
    ggplot(aes(x=mean.meth.diff,y=logFC))+
    geom_point(aes(col=(FDR>0.05), 
                   #alpha=(FDR<0.05)
    ), 
    show.legend = F)+
    geom_smooth(method = "lm", se = F)+
    stat_cor()+
    geom_vline(xintercept = c(-20,20))+
    geom_hline(yintercept = c(-1,1))+
    facet_wrap(~annotation2) +
    geom_text(data=df3, aes(x=meanmethdiff_cut, y=logFC_cut, label=n))
  print(p)
  
  return(df2)
}
