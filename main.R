# détection des segments peu recouvrants
# amélioration pour la détection des multiples
# http://www.rpubs.com/geospacedman/routing

# AM° : stop function, arcs multiples, snapping des extrêmités selon un seuil d'éloignement (si pas de croisement)
setwd("D:/WKSP/160208_DDTM83_BCAE")

source("R/functions.R")
# source("R/script_old.R")
source("R/script.R")


## CHARGEMENT DES DONNEES ##
# bcae <- readOGR("IN", "TRONCON_HYDROGRAPHIQUE__V3_modif")
# save(bcae, file="IN/bcae.rda")
load("IN/bcae.rda")
# bdtopo_bcae <- readOGR("IN", "bdtopo_bcae_83")
# save(bdtopo_bcae, file="IN/bdtopo_bcae.rda")
load("IN/bdtopo_bcae.rda")

lg <- gLength(bcae, byid=T)
o <- order(lg)

for (i in 1:length(o)) {
  
  if (!(i %in% c(271, 437, 438))) {
  
  # EXTRACTION DES DONNEES
  id <- bcae$id[o][i] 
  
  # INFOS
  print(paste(">>", i, "sur", length(o), ": bcae" , id, "|", round(lg[o][i]/1000), sep=" "))
  
  # SELECTION DU COURS D'EAU
  couche2 <- bcae[bcae$id==id, ]  
  
  # SELECTION DE LA BDTOPO PRE-IDENTIFIEE DANS LA ZONE DE LA FRAYERE
  regle <- paste("^", id, "[^(0-9)]+", "|", # id en début de ligne
               "[^(0-9)]+", id, "$", "|", # id en fin de ligne
               "[^(0-9)]+",id, "[^(0-9)]+", "|", # id entre deux autres ids
               "^",id, "$", sep="") # seul
#   rivers <- bdtopo_frayere[grep(id, bdtopo_frayere$id_frayere), ]  
#   bdtopo_frayere$id_frayere[grep(reg, bdtopo_frayere$id_frayere)] 
  couche1 <- bdtopo_bcae[grep(regle, bdtopo_bcae$id_bcae_83), ] 
  
  # CONFIGURATION DE L'EXPORT
  outputDir = "OUT/APPAR_BCAE"
  outputName = paste("bdtopo_bcae", id, sep="_")  
  fullPath = file.path(outputDir, paste(outputName, "shp", sep="."))
  
  # LAUNCH
  if(!(file.exists(fullPath))) {
    
    if (gIntersects(couche1, couche2)) {
      
      trouverSegments(couche1=couche1,
                      couche2=couche2, 
                      id=id, outputDir=outputDir, outputName=outputName)
    } else {
      line=paste("bca", id, "n'intersecte pas")
      write(line, file=file.path(outputDir, "log.txt"), append=TRUE)         
    }  
  
  } else {
    print(paste(fullPath, "existe"))
  } 
  }
}