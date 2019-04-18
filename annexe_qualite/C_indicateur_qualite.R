library(rgdal)
library(rgeos)

setwd("D:/WKSP/150916_SBEP_JGL_COURS_DEAU")

coucheref <- readOGR("IN", "rivieres")
bdtopo <- readOGR("OUT", "bdtopo_riviere_appariee")

# LES IDS à CALCULER
ids = bdtopo$id

pts1 = spsample(coucheref, n=100000, type="regular")
ids1 = over(gBuffer(pts1, width=1, byid=T), coucheref)$id
coords1 = coordinates(pts1)

pts2 = spsample(bdtopo, n=100000, type="regular")
ids2 = over(gBuffer(pts2, width=1, byid=T), bdtopo)$id
coords2 = coordinates(pts2)


# DISTANCES PAR POINTS
out=list()
dists = get.knnx(data=coords1, query=coords2, k=1)$nn.dist
pts = SpatialPointsDataFrame(pts2, data=data.frame(id=ids2, distance=dists))

pts$classe_distance = 0
pts$classe_distance[which((dists - mean(dists)) > 2*sd(dists))] = 1
pts$classe_distance[which((dists - mean(dists)) > 3*sd(dists))] = 2
pts$classe_distance[which((dists - mean(dists)) > 4*sd(dists))] = 3
pts$classe_distance[which((dists - mean(dists)) > 5*sd(dists))] = 4
pts$classe_distance[which((dists - mean(dists)) > 6*sd(dists))] = 5

writeOGR(pts, "OUT_QUALITE", "pts_riviere", "ESRI Shapefile")


# DISTANCES MOYENNES PAR LIGNES
out=list()

for (i in 1:length(ids)) {
  
  id=ids[i]
  
  coucheref.sample = coucheref[coucheref$id==id, ]
  bdtopo.sample = bdtopo[bdtopo$id==id,  ]
  
  pts.ref = spsample(coucheref.sample, n=gLength(coucheref.sample)/10, type="regular") # tous les 10 m
  pts.bdtopo = spsample(bdtopo.sample, n=gLength(bdtopo.sample)/10, type="regular")
  
  coords.ref = coordinates(pts.ref)
  coords.bdtopo = coordinates(pts.bdtopo)
  
  dists = get.knnx(data=coords.ref, query=coords.bdtopo, k=1)$nn.dist
  meandist = mean(dists)
    
  out[[i]] = meandist
  
}

meandists = unlist(out)

df = data.frame(ids, meandist=meandists, classe_distance=0)

df$classe_distance[which((meandists-mean(meandists)) > 1*sd(meandists))] = 1
df$classe_distance[which((meandists-mean(meandists)) > 2*sd(meandists))] = 2
df$classe_distance[which((meandists-mean(meandists)) > 3*sd(meandists))] = 3
df$classe_distance[which((meandists-mean(meandists)) > 4*sd(meandists))] = 4
df$classe_distance[which((meandists-mean(meandists)) > 5*sd(meandists))] = 5

write.csv(df, "OUT_QUALITE/meandists_riviere.csv")


# exportshp(pts1, "pts1")
# 
# 
# # coucheref <- coucheref[coucheref$id == 494, ]
# pts1 <- readOGR("OUT/BDTOPO_coucherefERE//300915/FUSION/pts_bdtopo.shp", "pts_bdtopo")
# pts2 <- readOGR("OUT/BDTOPO_coucherefERE//300915/FUSION/pts_coucherefere.shp", "pts_coucherefere")
# 
# distanceEntreLignes <- function(ligne1, ligne2, every=100) {
# 
#  lg <- gLength(ligne1)
# 
#  nPts <- lg/every
# 
#  ligne1.pts <- spsample(ligne1, n=nPts, type="regular")
#  ligne2.pts <- spsample(ligne2, n=nPts, type="regular")
# 
#  d <- gDistance(ligne1.pts, ligne2.pts, byid=T)
#  mind <- apply(d, 2, min)
#  meand <- mean(mind)
#  sdd <- sd(mind)
#  outlier <- (mind - meand) > sdd
#  df <- data.frame(dist=mind, diff=(mind-meand), outlier=outlier)
#  ligne1.pts.df <- SpatialPointsDataFrame(ligne1.pts, data=df)
# 
#  return(ligne1.pts.df)
# 	
# }
# 
# distanceEntrePoints <- function(ligne1.pts, ligne2.pts, every=100) {
#   
#   d <- gDistance(ligne1.pts, ligne2.pts, byid=T)
#   mind <- apply(d, 2, min)
#   meand <- mean(mind)
#   sdd <- sd(mind)
#   outlier <- (mind - meand) > sdd
#   df <- data.frame(dist=mind, diff=(mind-meand), outlier=outlier)
#   ligne1.pts.df <- SpatialPointsDataFrame(ligne1.pts, data=df)
#   
#   return(ligne1.pts.df)
#   
# }
# 
# pts <- distanceEntrePoints(pts1, pts2)
# 
# pts <- distanceEntreLignes(river, coucheref, every=10)
# pts@data
# 
# writeOGR(pts , "OUT_QUALITE", "494_5", "ESRI Shapefile")
# 
# 
# ligne1.pts <- spsample(river, n=nPts, type="regular")
# 
# ligne1.pts <- spsample(river, n=nPts, type="regular")
# ligne2.pts <- spsample(coucheref, n=nPts, type="regular")
# 
# d <- gDistance(ligne1.pts, ligne2.pts, byid=T)
# mind <- apply(d, 2, min)
# 
# df <- data.frame(dist=mind)
# ligne1.pts.df <- SpatialPointsDataFrame(ligne1.pts, data=df)
# 
# meand <- mean(mind)
# sdd <- sd(mind)
# 
# 
# 
# distanceEntreLignes <- function(ligne1, ligne2, n) {
# 
# 	ligne1.pts <- spsample(ligne1, n=n, type="regular")
# 	ligne2.pts <- spsample(ligne2, n=n, type="regular")
# 
# 	d <- gDistance(ligne1.pts, ligne2.pts, byid=T)
# 	d <- apply(d, 2, min)
# 
# 	meanDist <- mean(d)
# 
# 	return(meanDist)
# }
# 
# distanceEntreLignes(river, coucheref, 100)



