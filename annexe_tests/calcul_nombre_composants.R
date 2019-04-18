library(rgdal)
library(igraph)
setwd("D:/WKSP/150916_SBEP_JGL_COURS_DEAU/")
source("R/functions.R")


ce = readOGR("OUT/bdtopo_ce_appariee.shp", "bdtopo_ce_appariee")

compterComposants = function(riv) {

  out=list()
  for (i in 1:nrow(riv)) {
    
    f= riv[i, ]
    f = mergeThenSplitLine(f)
    g = construireReseau(f)
    dec = decompose.graph(g)
    ncomps=length(dec)
    out[[i]]=ncomps
    
  }
  
  res = unlist(out)
  return(res)

}

compterComposants(riv)

ce$nLignes = compterComposants(ce)

writeOGR(ce, "OUT", "bdtopo_ce_appariee_nLignes", "ESRI Shapefile")