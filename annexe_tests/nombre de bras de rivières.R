library(rgdal)

setwd("D:/WKSP/150916_SBEP_JGL_COURS_DEAU/")

riv = readOGR("IN/rivieres.shp", "rivieres")

nBras <- sapply (1:nrow(riv), function(x) {
  
  m = mergeThenSplitLine(riv[x, ]  )
  return(length(m))
  
})
ids = riv$id[which(nBras > 1)]