## calculate heatmap size for plotting in pdf format
## https://jokergoo.github.io/2020/05/11/set-cell-width/height-in-the-heatmap/

calc_ht_size = function(ht, unit = "inch") {
  pdf(NULL)
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()
  
  c(w, h)
}

## usage

# mat2 = mat[1:5, 1:5]
# ht2 = Heatmap(mat2, name = "mat2", 
#         width = ncol(mat2)*unit(5, "mm"), 
#         height = nrow(mat2)*unit(5, "mm"))
# size2 = calc_ht_size(ht2)

# pdf(..., width = size2[1], height = size2[2])
# ht2 # or draw(ht2)