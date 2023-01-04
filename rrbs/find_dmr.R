# last modified 2022-12-16
# Author: Feng Yan

# This code runs DMC and DMR identification together using methylKit and edmr
# the DMC and DMR cutoff used here is relative liberal compared to the default

find_dmr <- function(obj, condition, min.per.group = 3L, cores = 10) {
  library(methylKit)
  library(edmr)
  ## input methylBase or methylRawlist
  if (class(obj) == "methylBase") {
    ## subset directly
    meth <- obj 
    meth_sub <- methylKit::reorganize(meth,
                                      sample.ids = obj@sample.ids[obj@treatment %in% condition],
                                      treatment = obj@treatment[obj@treatment %in% condition])
    
  } else if (class(obj) == "methylRawList") {
    ## first subset, then unite 
    cov_temp <- methylKit::reorganize(obj,
                                      sample.ids = getSampleID(obj)[obj@treatment %in% condition],
                                      treatment = obj@treatment[obj@treatment %in% condition])
    meth_sub <- methylKit::unite(cov_temp, destrand = FALSE, min.per.group = min.per.group)
  } else 
    print("Input is not a methylBase or methylRawList object")
  
  ## this is to avoid seqname inconsitent issue

  ## DMC 
  meth.diff <- methylKit::calculateDiffMeth(meth_sub, mc.cores = cores)
  cpg.sig <- getMethylDiff(meth.diff, difference = 25, qvalue = 0.01) ## default difference and qvalue
  cpg.sig <- as(cpg.sig, "GRanges")
  cpg.sig <- keepStandardChromosomes(cpg.sig, pruning.mode = "coarse")
  
  ## DMR, use liberal cutoff to output unfiltered DMRs
  meth.dmr <- edmr::edmr(meth.diff, 
                         DMC.methdiff = 20, ## default is 25
                         DMR.methdiff = 10, ## default is 20
                         DMC.qvalue = 0.01, num.DMCs = 1, num.CpGs = 3)
  dmr.sig <- filter.dmr(meth.dmr, mean.meth.diff = 20,
                        DMR.qvalue = 0.001, num.DMCs = 3, 
                        num.CpGs = 3)  ## default is 5
  dmr.sig <- keepStandardChromosomes(dmr.sig, pruning.mode = "coarse")
  return(list(cpg.sig = cpg.sig, dmr.sig = dmr.sig, dmr = meth.dmr))
}
