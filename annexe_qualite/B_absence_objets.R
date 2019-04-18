library(rgdal)

setwd("D:/WKSP/160208_DDTM83_BCAE")

appar = readOGR("OUT", "bdtopo_bcae_83_appariee")

bcae = readOGR("IN", "TRONCON_HYDROGRAPHIQUE__V3_modif")

m=match(bcae$id, appar$id)

ids = bcae$id[which(is.na(m))]

ids[order(ids)]

# [1] 863  28 266 280
