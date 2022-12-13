# last modified 2022-12-12

parse_gct <- function(x){
  # x is a DGElist
  library(edgeR)
  gctfile <- cbind(as.character(x$genes$Gene.Name) %>% toupper(), as.character(x$genes$Gene.ID), 
                   cpm(x, log = T))
  colnames(gctfile)[1:2] <- c("SYMBOL","ENSEMBL")
  write("#1.2", "lcpm.gct")
  write(dim(x), "lcpm.gct", append = T, sep = "\t")
  write.table(gctfile, "lcpm.gct", quote = F, row.names = F, append = T, sep = "\t")
}

parse_cls <- function(groupfct){
  # group is a factor, usually x$sample$group
  ## need to reorder because GSEA recognize groups by order, not by matching
  library(glue)
  group_cha <- as.character(groupfct)
  group_cha <- group_cha[!duplicated(group_cha)]
  group <- factor(as.character(groupfct), levels = group_cha)
  clsfile <- glue('{nsample} {ngroup} 1 
                  # {name}
                  {group}',
                  nsample = length(group),
                  ngroup = length(levels(group)),
                  name = glue_collapse(levels(group), sep = " "),
                  group = glue_collapse(as.character(group), sep = " ")) ## use names rather than numbers
  write(clsfile, "pheno.cls", sep = "\t")
}