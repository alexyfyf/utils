# last modified 2022-12-16
# Author: Feng Yan

## plot profile meta line plot using genomation package

plot_signal <- function(region, flank, filename, cores=12, select=1:4,
                        meth.gr = meth_grlist, atac.files = bwfiles3, 
                        width=15, height=8) {
  require(scales)
  require(GenomicRanges)
  require(genomation)
  #sml <- ScoreMatrixList(bwfiles, region+flank, cores = 12) ## HINT-ATAC already normalize and correct signal
  #sml2 <- ScoreMatrixList(bwfiles2, region+flank, cores = 12)
  sml3 <- ScoreMatrixList(atac.files[select], region+flank, cores = cores
                          #, rpm = T ##reads normalized to library size
                          ) 
  sml4 <- ScoreMatrixList(meth.gr[select], #bwfiles4[c(5,17,1,13)],
                          region+flank, cores = cores,
                          ## not sure here
                          is.noCovNA = T, weight.col = "score")
  
  #pdf(filename, width = 15, height = 15)
  #par(mfrow=c(2,2))

  # browser()
  pdf(filename, width = width, height = height)
  par(mfrow=c(1,2))
  xlim=flank+ floor(width(region)[1]/2)
  # plotMeta(sml, profile.names = c("LSC","WL","preLSC","WE"),
  #          xcoords = -xlim:xlim, ylab = "HINT-corrected signal",
  #          line.col = hue_pal()(4)[c(1,2,4,3)])
  # plotMeta(sml2, profile.names = c("LSC","WL","preLSC","WE"),
  #          xcoords = -xlim:xlim, ylab = "TOBIAS-corrected signal",
  #          line.col = hue_pal()(4)[c(1,2,4,3)])
  colors <- c(hue_pal()(4)[c(1,2,4,3)], "#808080")
  colors <- colors[select]
  names <- c("LSC","WL","preLSC","WE","TALL")[select]
  plotMeta(sml3, profile.names = names,
           xcoords = -xlim:xlim, ylab = "Normalized read counts",
           line.col = colors)
  plotMeta(sml4, profile.names = names, 
           #smoothfun=function(x) stats::lowess(x, f = 1/5),
           xcoords = -xlim:xlim, ylab = "Methylation %",
           line.col = colors)
  dev.off()
}


plot_signal_bin <- function(region, flank, filename, cores = 12, bin.num = 50,
                            select=1:4,
                            meth.gr = meth_grlist, atac.files = bwfiles3,
                            width=15, height=8) {

  require(scales)
  require(GenomicRanges)
  require(genomation)
  #sml <- ScoreMatrixList(bwfiles, region+flank, cores = 12) ## HINT-ATAC already normalize and correct signal
  #sml2 <- ScoreMatrixList(bwfiles2, region+flank, cores = 12)
  sml3 <- ScoreMatrixList(atac.files[select], region+flank, cores = cores, 
                          bin.num = bin.num
                          #, rpm = T ##reads normalized to library size
  ) 
  sml4 <- ScoreMatrixList(meth.gr[select], #bwfiles4[c(5,17,1,13)], 
                          region+flank, cores = cores, bin.num = bin.num,
                          ## not sure here
                          is.noCovNA = T, weight.col = "score")
  
  # pdf(filename, width = 15, height = 15)
  # par(mfrow=c(2,2))
  
  # browser()
  pdf(filename, width = width, height = height)
  par(mfrow=c(1,2))
  # xlim=flank+ floor(width(region)[1]/2)
  # plotMeta(sml, profile.names = c("LSC","WL","preLSC","WE"),
  #          xcoords = -xlim:xlim, ylab = "HINT-corrected signal",
  #          line.col = hue_pal()(4)[c(1,2,4,3)])
  # plotMeta(sml2, profile.names = c("LSC","WL","preLSC","WE"),
  #          xcoords = -xlim:xlim, ylab = "TOBIAS-corrected signal",
  #          line.col = hue_pal()(4)[c(1,2,4,3)])
  colors <- c(hue_pal()(4)[c(1,2,4,3)], "#808080")
  colors <- colors[select]
  names <- c("LSC","WL","preLSC","WE","TALL")[select]
  plotMeta(sml3, profile.names = names,
           #xcoords = -xlim:xlim, 
           ylab = "Normalized read counts",
           line.col = colors)
  plotMeta(sml4, profile.names = names, 
           #smoothfun=function(x) stats::lowess(x, f = 1/5),
           #xcoords = -xlim:xlim, 
           ylab = "Methylation %",
           line.col = colors)
  dev.off()
}

