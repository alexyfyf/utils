# last modified 2022-12-16
# Author: Feng Yan

# This code performs full join on GRanges (under development, especially to address the 'within' and 'directed' suffix)
# see this issue https://github.com/sa-lee/plyranges/issues/68

join_overlap_full <- function(x,y){
  options(stringsAsFactors = F)
  hit <- findOverlaps(x,y)
  colnames(mcols(x)) <- paste0(colnames(mcols(x)),".x")
  colnames(mcols(y)) <- paste0(colnames(mcols(y)),".y")
  cns.x <- mcols(x) %>% colnames()
  cns.y <- mcols(y) %>% colnames()
  
  z <- x[queryHits(hit),]
  hitz <- cbind(mcols(y)[subjectHits(hit),],
                mcols(x)[queryHits(hit),])
  colnames(hitz) <- c(cns.y,cns.x)
  mcols(z) <- hitz
  uniqex <- x[-queryHits(hit),]
  uniqey <- y[-subjectHits(hit),]
  mcols(uniqex)[, (length(cns.x)+1):(length(cns.x)+length(cns.y))] <- NA
  colnames(mcols(uniqex)) <- c(cns.x,cns.y)
  mcols(uniqex) <- mcols(uniqex)[,c((length(cns.x)+1):(length(cns.x)+length(cns.y)),1:(length(cns.x)))]
  mcols(uniqey)[, (length(cns.y)+1):(length(cns.x)+length(cns.y))] <- NA
  colnames(mcols(uniqey)) <- c(cns.y,cns.x)
  res <- c(z, uniqex,uniqey)
  mcols(res) <- mcols(res) %>% data.frame() %>% mutate_all(as.character())
  return(res)
}