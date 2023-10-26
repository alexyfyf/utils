# last modified 2022-12-16
# Author: Feng Yan

# This code plot Sankey plot using list of elements

plot_sankey <- function(gr_list, column = "meth.diff", fulljoin = T){
  require(ggalluvial)
  require(tidyverse)
  require(plyranges)
  # source("http://118.138.241.167/data/Lmo2_ATAC/atac_diff/code/join_overlap_full.R")
  
  subset <- lapply(gr_list, function(x) {
    x <- x %>% plyranges::mutate(status=ifelse(!!sym(column)>0,"up","down")) %>% 
      plyranges::select(status)
  })
  
  z <- Reduce(join_overlap_full, subset) ## use custome full join function
  colnames(mcols(z)) <- rev(names(gr_list))
  
  p <- z %>% mcols() %>%
    data.frame() %>% 
    mutate_all(function(fac) {
      factor(fac, levels = levels(addNA(fac)), labels = c("down","up","nochange"), 
             exclude = NULL)
    }) %>% 
    dplyr::select(length(gr_list):1) %>% ## to reverse the merged data.frame
    group_by_all() %>%
    dplyr::summarise(count = dplyr::n()) %>%
    to_lodes_form(axes=1:length(gr_list), id="id") %>%
    ggplot(aes(x=x, y=count, stratum=stratum, alluvium=id,
               fill=stratum, label=stratum)) +
    geom_flow(stat = "alluvium", color="darkgray")+
    geom_stratum(alpha=0.5) + 
    scale_fill_manual(values = c("#F8766D","#00BA38","#619CFF"))
  print(p)
  return(list(df = z, plot = p))
}
