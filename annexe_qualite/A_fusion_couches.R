library(rgdal)

setwd("D:/WKSP/160208_DDTM83_BCAE")

rootFolder = "OUT/APPAR_BCAE/"

l <- list.files(rootFolder, "*.shp$")

layerFileNames = unlist(l)
layerNames = sapply(l, function(x) gsub(".shp", "", x))

out=list()
ncol=list()

for (i in 1:length(l)) {
  print(i)
  spLine =  readOGR(file.path(rootFolder, layerFileNames[i]), layerNames[i])
  spLineM = gLineMerge(spLine) # FUSION
  spLineM = spChFIDs(spLineM, as.character(i)) # CHANGEMENT de L'ID
  df = spLine@data[1, ]; row.names(df)=i # RECUP DES ATTRIBUTS
  spLine = SpatialLinesDataFrame(spLineM, data=df) # CREATION DE SPDF
  out[[i]] <- spLine 
}

save(out, file="OUT/out.rda")

# RESULTATS NULS
which(sapply(out, is.null))

# NOMBRE DE COLONNES
ncols = sapply(out, function(x) ncol(x@data))

out1 = out[ncols == 4]
out2 = out[ncols == 5]

res1 = do.call(rbind, out1)
res2 = do.call(rbind, out2)

res2$id_pce=NULL
names(res1)=names(res2)

res = spRbind(res1, res2)

names(res) = names(out[[1]])

for (i in 1:length(out)) {
  names(out[[i]]) <- c("idLigne", "id", "id_pce", "nLignes", "propLong") 
}

# FUSION DES DONNEES
res= do.call(rbind, out)
res = res[order(res$id), ]

writeOGR(res, "OUT", "bdtopo_bcae_83_appariee", "ESRI Shapefile")
