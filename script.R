# remplacer labels garder et supprimer par booléen keep=T ou F
# mettre noms plus parlants

trouverSegments <- function(couche1, couche2, id, outputDir, outputName) {

time1 <- Sys.time()
   
##>> GRAPHE ##
  
print("construction du graphe...")

couche1.fus = mergeThenSplitLine(couche1)
fusion <- croiserLignes(couche1.fus, couche2)

couche1_buff_100m <- gBuffer(couche1, width=100)
couche2_buff_100m <- gBuffer(couche2, width=100)

g <- construireReseau(fusion)
E(g)$name <- attributes(E(g))$vnames


##>> LIBELLER LES NOEUDS SELON CROISEMENT OU RIVIERE

print("étiquetage des noeuds...")

V(g)$label = NA
# INTERSECTIONS
# !! des fois, pas d'intersection trouvée
cross <- croisements(g, couche1.fus, couche2)
V(g)$label[cross] <- "croisement"

# RIVIERES
delacouche <- deLaCouche(g, couche1.fus)
V(g)$label[delacouche] <- "riviere"
# FRAYERES
V(g)[which(is.na(V(g)$label))]$label = "frayere"


# ARCS JOIGNANT DEUX CROISEMENTS ET APPARTENANT AUX FRAYERES
vertexes <- V(g)[which(V(g)$label=="croisement")]
pasdelapremierecouche <- pasDansLaCouche(g, fusion, couche1.fus, vertexes) # fonction calculant si les arcs joignant deux vertexes appartiennent ou pas à une couche

# LABEL DES ARCS
E(g)$garder = TRUE
E(g)[inc(V(g)[which(V(g)$label == "frayere")])]$garder = FALSE
E(g)[pasdelapremierecouche]$garder <- FALSE

# SOUS-GRAPHE SANS LES ARCS RELATIFS AUX FRAYERES
g.premierecouche <- subgraph.edges(g, E(g)[E(g)$garder])
premierecouche <- fusion[E(g)$garder, ]

# DEBUT ET FIN
vertexes = V(g.premierecouche)[which(V(g.premierecouche)$label == 'riviere')]
extremitesCouche2 <- procheExtremitesCouche(g.premierecouche, vertexes, couche2, min=10)
V(g.premierecouche)[extremitesCouche2]$label <- "début et fin"

##>> RESEAU ÉPURÉ
print("épure du graphe...")

pivots <- V(g.premierecouche)[which(V(g.premierecouche)$label %in% c("croisement", "début et fin"))]
nPivots <- length(pivots)

out <- list()
for (i in 1:nPivots) { 
  
  # FROM
  vfrom <- pivots[[i]]  
  # TO
  tonodes <- pivots
  tonodes <- tonodes[as.numeric(tonodes) >= as.numeric(vfrom)]
  
  # TROUVER LES PROCHAINS NOEUDS
  res <- findNextNodes(g.premierecouche, vfrom, tonodes)
  epaths <- res$epaths
  ends <- res$ends
  
  # FILTRAGE PAR PROPORTION DE RECOUVREMENT
  props <- sapply(epaths, function(p) {
    premierecouche.lg <- gLength(premierecouche[p, ])
    if(gIntersects(premierecouche[p, ], couche2_buff_100m)) {
      overlapped.lg <- gLength(gIntersection(premierecouche[p, ], couche2_buff_100m))
      prop <- overlapped.lg / premierecouche.lg
    } else {
      prop <- 0
    }
    return(prop)  
  })  
  p <- which(props > 0.1)
  out[[i]] <- unlist(epaths[p])
}

if(length(out) > 0) {
  
selec <- unique(unlist(out)); selec <- selec[order(selec)]

g.epure <- subgraph.edges(g.premierecouche, E(g.premierecouche)[selec])
coucheepure <- premierecouche[selec, ]

# save(g.epure, file="OUT/g.epure_ce863.rda")
# save(coucheepure, file="OUT/couceepure_ce863.rda")

##>> DETECTION DU MEILLEUR CHEMIN

print("détection du meilleur chemin...")
comps <- decompose.graph(g.epure, min.vertices=2)

# ON FILTRE LES COMPOSANTS QUI TOUCHENT LA COUCHE REFERENTIELLE
ok <- which(sapply(1:length(comps), function(i) {
  g.comp <- comps[[i]]
  sp.comp <- coucheepure[match(E(g.comp)$name, E(g.epure)$name), ]
  return(gIntersects(sp.comp, couche2))
}))

# POUR CHAQUE COMPOSANT, CALCUL DU MEILLEUR CHEMIN
out=list()
for (i in ok) {
#   print(i)
  g.comp <- comps[[i]]
  sp.comp <- coucheepure[match(E(g.comp)$name, E(g.epure)$name), ]
  couchecompare = gBuffer(couche2, width=25)
  out[[i]] <- trouverMeilleurChemin(g.comp, sp.comp, couchecompare)
} 

if(length(out) > 0) {
  best = SpatialLines(list(Lines(unlist(out), ID="id")))
} else {
  #  !! quand couche epure, compiler tous les objets
  best = coucheepure  
}

     
##>> PLOTTING ##

plot(couche2, lwd=5, col="lightblue"); plot(best, add=T, col="red")
title(id)


##>> EXPORT

# MAIN 
df <- data.frame(idLigne=1:length(best), 
                 id = id, 
                 id_pce = NA, 
                 nLignes = length(comps),
                 propLong = gLength(best) / gLength(couche2) * 100,
                 row.names=row.names(best))

# EXPORT
best.df <- SpatialLinesDataFrame(best, data=df)

# EXPORT
writeOGR(best.df, outputDir, outputName, "ESRI Shapefile")

}# if(length(e) > 0)
else {
  print("pas de plus court chemin trouvé")  
}
# INFORMATION
time2 <- Sys.time()
print(difftime(time2, time1))

}# function

