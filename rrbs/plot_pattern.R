# last modified 2022-12-16
# Author: Feng Yan

## plot epiallele patterns in samples ordered by order
## loci in the format of for example "chr4_+_55531465:55531477:55531489:55531502"
## order must have a column sample

plot_pattern <- function(loci, grlist=epi.gr, order=order){
  tmp=loci %>% str_split(pattern = "_|:", simplify = T) 
  
  pt=lapply(grlist, function(x){
    pt=data.frame(x) %>% 
      filter(seqnames==tmp[1] & strand==tmp[2] & start==tmp[3] & end==tmp[6])
    df=data.frame(pt_pct=t(pt)[13:28] %>% as.numeric())
    df
  }) %>% Reduce(cbind, .)
  
  colnames(pt) <- names(grlist)
  
  pattern=(values(grlist[[1]]) %>% colnames)[8:23] 
  pattern=str_split(pattern,":",n=2, simplify = T)[,2]
  
  pt2=pt %>% mutate(pattern=pattern) %>% 
    gather(key = sample, value = pct, 1:length(grlist)) %>% 
    left_join(order)
  
  if (!is.null(order)) {
    pt2$sample <- factor(pt2$sample, levels = as.character(order$sample))
  }
  
  # mycolors=colorRampPalette(c("blue","white","red"))(16)
  mycolors=colorRampPalette(brewer.pal(name="Dark2", n = 8))(16)
  p=ggplot(pt2, aes(x=sample, y=pct, fill=pattern)) + 
    geom_bar(stat="identity") + ggtitle(loci) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values = mycolors)
  p
}
