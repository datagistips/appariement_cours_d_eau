require(igraph)
require(sp)
require(FNN)
require(rgeos)
require(rgdal)
require(maptools)
require(R.utils)

findStartEndNodes <- function(spLines) {
  coords <- spLines@lines[[1]]@Lines[[1]]@coords
  startEndNodes <- rbind(head(coords, 1), tail(coords, 1))
  rownames(startEndNodes) <- apply(startEndNodes, 1, paste, collapse=" ")
  return(startEndNodes)
}

construireReseau <- function(spLines) {  
  # EDGES
  edges = do.call(rbind, lapply(spLines@lines, function(ls) {
    cbind(head(ls@Lines[[1]]@coords, 1), 
          tail(ls@Lines[[1]]@coords, 1))}))
  
  # EDGE IDS
  # il y a mieux ?
  IDS <- sapply(1:nrow(edges), function(x) paste(edges[x, ][order(edges[x, ])], collapse=" "))
  
  # LONGUEUR DES SEGMENTS
  lengths = gLength(spLines, byid=T)
  
  # GRAPHE
  froms = paste(edges[, 1], edges[, 2])
  tos = paste(edges[, 3], edges[, 4])
  
  graph = graph.edgelist(cbind(froms, tos), directed = FALSE)
  
  E(graph)$weight = lengths
  
  # COORDONNEES DES VERTICES
  xy = do.call(rbind, strsplit(V(graph)$name, " "))
  
  V(graph)$x = as.numeric(xy[, 1])
  V(graph)$y = as.numeric(xy[, 2])
  
  return(graph)
}

routePoints <- function(graph, from, to) {
  require(FNN)
  xyg = cbind(V(graph)$x, V(graph)$y)
  
  ifrom = get.knnx(xyg, from, 1)$nn.index[1, 1]
  ito = get.knnx(xyg, to, 1)$nn.index[1, 1]
  
#   p = get.shortest.paths(graph, ifrom, ito, output = "both")
  p = get.shortest.paths(graph, ifrom, ito)
  p[[1]]
}

croiserLignes_old <- function(couche1, couche2) {  
  proj4string(couche1) <- proj4string(couche2)
  
  couche1 <- spChFIDs(couche1, paste("1", 1:length(couche1), sep="_"))
  couche2 <- spChFIDs(couche2, paste("2", 1:length(couche2), sep="_"))
  # MERGE
  fus0 <- spRbind(as(couche1, "SpatialLines"), as(couche2, "SpatialLines"))
  # sans le TRUE c'est plus rapide
  fus1 <- gIntersection(fus0, fus0)
  #   s <- sapply(row.names(fus1), function(x) {
  #     split <- strsplit(x, " ")
  #     return(split[[1]][1] == split[[1]][2])})
  #   df <- data.frame(identicalSegs = s, pairs = row.names(fus1), row.names=row.names(fus1))
  #   fus <- SpatialLinesDataFrame(fus1, data=df)
  # EXPLODE
  nLines <- length(fus1@lines[[1]]@Lines); nLines
  l <- lapply(1:nLines, function(x) Lines(list(fus1@lines[[1]]@Lines[[x]]), ID=as.character(x)))
  fus <- SpatialLines(l)
  return(fus)
}

croiserLignes <- function(couche1, couche2) {  
  proj4string(couche1) <- proj4string(couche2)
  
  couche1 <- spChFIDs(couche1, paste("1", 1:length(couche1), sep="_"))
  couche2 <- spChFIDs(couche2, paste("2", 1:length(couche2), sep="_"))
  # MERGE
  fus0 <- spRbind(as(couche1, "SpatialLines"), as(couche2, "SpatialLines"))
  # sans le TRUE c'est plus rapide
  fus1 <- gIntersection(fus0, fus0)
  #   s <- sapply(row.names(fus1), function(x) {
  #     split <- strsplit(x, " ")
  #     return(split[[1]][1] == split[[1]][2])})
  #   df <- data.frame(identicalSegs = s, pairs = row.names(fus1), row.names=row.names(fus1))
  #   fus <- SpatialLinesDataFrame(fus1, data=df)
  # EXPLODE
#   nLines <- length(fus1@lines[[1]]@Lines); nLines
#   l <- lapply(1:nLines, function(x) Lines(list(fus1@lines[[1]]@Lines[[x]]), ID=as.character(x)))
#   fus <- SpatialLines(l)
  
  fus <- mergeThenSplitLine(fus1)
  
  return(fus)
}

mergeThenSplitLine <- function(spLine) {
  m <- gLineMerge(spLine)
  nLines <- length(m@lines[[1]]@Lines)
  l <- lapply(1:nLines, function(x) Lines(list(m@lines[[1]]@Lines[[x]]), ID=x))
  merged <- SpatialLinesDataFrame(SpatialLines(l), data=data.frame(1:length(l)))  
  return(merged)
}

getNodes <- function(spL) {
  l <- lapply(1:length(spL), function(x) spL@lines[[x]]@Lines[[1]]@coords)
  nodes <- do.call(rbind, l)
  rownames(nodes) <- apply(nodes, 1, paste, collapse=" ")
  return(nodes)
}

SpatialLinesFromGraph <- function(graph) {
  e <- get.edges(graph, E(graph))
  l <- lapply(1:nrow(e), function(i) {
    first <- cbind(V(graph)[[e[i,1]]]$x, V(graph)[[e[i,1]]]$y)
    second <- cbind(V(graph)[[e[i,2]]]$x, V(graph)[[e[i,2]]]$y)
    coords <- rbind(as.numeric(first), as.numeric(second))
    return(SpatialLines(list(Lines(list(Line(coords)), ID=i))))
  })
  return(do.call(rbind, l))
}
  
makeLineFromCoords <- function(coords, i) {
    Sl1 = Line(coords)
    S1 = Lines(list(Sl1), ID=as.character(i))
    Sl = SpatialLines(list(S1))  
    return(Sl)
  }

SpatialPointsFromVertexes <- function(vertexes) {
  coords <- cbind(vertexes$x, vertexes$y)
#   df <- data.frame(id=as.numeric(vertexes), label=ifelse(is.null(vertexes$label), "-rien-", vertexes$label))
  df <- data.frame(id=as.numeric(vertexes), label=ifelse(is.null(vertexes$label), "-rien-", vertexes$label))
  pts <- SpatialPointsDataFrame(SpatialPoints(coords), data=df)
  return(pts)
}

edgeFromVertices <- function(graph, vpath) {
  vsSuite <- as.vector(sapply(1:(length(vpath)-1), function(i) c(vpath[[i]], vpath[[i+1]])))
  return(get.edge.ids(graph, vsSuite))  
}

findNextNodes <- function(graph, node, nodes) {  
  w <- setdiff(nodes, node)
  p <- get.shortest.paths(graph, node, V(graph)[w])$vpath  
  p <- p[which(sapply(p, length) > 0)]
  
  vpaths <- unique(Map(function(p) as.vector(p[1:which(!is.na(match(p, nodes)))[2]]), p))
  ends <- lapply(vpaths, tail, 1)  
  epaths <- lapply(vpaths, function(x) edgeFromVertices(graph, x))
  
  return(list(vpaths=vpaths, ends=ends, epaths=epaths))
}

exportshp <- function(spatfile, name, folder="OUT_IGRAPH") {
  if(class(spatfile) %in% c("SpatialLinesDataFrame", "SpatialPointsDataFrame")) {
    spatfile.df <- spatfile
  }
  if(class(spatfile) == "SpatialLines") {
    df <- data.frame(id=1:length(spatfile), row.names=row.names(spatfile))
    spatfile.df <- SpatialLinesDataFrame(spatfile, data=df)
  }
  if(class(spatfile) == "SpatialPoints") {
    df <- data.frame(id=1:length(spatfile), row.names=row.names(spatfile))
    spatfile.df <- SpatialPointsDataFrame(spatfile, data=df)
  }  
  writeOGR(spatfile.df, folder, name, "ESRI Shapefile")
}

distanceEntreLignes <- function(ligne1, ligne2, n) {
  
  ligne1.pts <- spsample(ligne1, n=n, type="regular")
  ligne2.pts <- spsample(ligne2, n=n, type="regular")
  
  d <- gDistance(ligne1.pts, ligne2.pts, byid=T)
  d <- apply(d, 2, min)
  
  meanDist <- mean(d)
  
  return(meanDist)
}

cutPoints <- function(g.comp, sp.comp) {
    
  nCompIni <- no.clusters(g.comp)
  out <- list()
  coords <- list()
  for (i in 1:ecount(g.comp)) {
    nClusters <- no.clusters(delete.edges(g.comp, i))
    #     print(nClusters)
    if(nClusters == nCompIni) {
      out[[i]]=FALSE # PAS PONT
    } else {
      out[[i]]=TRUE # PONT
    }
  }
     
  pont <- unlist(out)
  
  if(any(pont==FALSE) & any(pont)==TRUE) {
        
    bridges <- sp.comp[which(pont), ]    
    notBridges <- sp.comp[which(!pont), ]
        
    int <- gIntersection(bridges, notBridges, byid=T)
    coordsint <- apply(coordinates(int), 1, paste, collapse="-") 
    coordsv <- apply(cbind(V(g.comp)$x,V(g.comp)$y), 1, paste, collapse="-")  
    w <- which(!is.na(match(coordsv, coordsint)))        
  } else {
    w = NULL # tous les arcs sont des ponts
  } 
    
  return(w)  
}

trouverMeilleurChemin <- function(g.comp, sp.comp, couchecompare) {
  
  # FUSION DU GRAPHE
  sp.compM <- mergeThenSplitLine(sp.comp)
  g.compM <- construireReseau(sp.compM)
  
  # START et END NODES
  V(g.compM)[get_diameter(g.compM)]$diametre = TRUE
  V(g.compM)[farthest.nodes(g.compM)$vertices]$label="début ou fin"
  
  # NOEUDS DE DEGRE UN
  deg1 <- which(degree(g.compM)==1 & (is.na(V(g.compM)$label)))
  
  # GARDER OU PAS
  V(g.compM)$garder=TRUE
  V(g.compM)[deg1]$garder = FALSE
  
  # SOUS GRAPHE
  E(g.compM)$garder = TRUE
  E(g.compM)[inc(V(g.compM)[which(!V(g.compM)$garder)])]$garder = FALSE
  
  g.nodeg1 <- subgraph.edges(g.compM, E(g.compM)[E(g.compM)$garder])
  sp.nodeg1 <- sp.compM[E(g.compM)$garder, ]
      
  # SUPPRESSION DES BOUCLES EXTERIEURES QUI NE TOUCHENT PAS LE DIAMETRE
  bc = biconnected.components(g.nodeg1)
  diamnodes = which(V(g.nodeg1)$diametre)
  
  w <- which(sapply(bc$component_edges, length) > 2)
   
  if (length(w) > 0) {
    # DIAMETRES
    dansdiam <- sapply(w, function(x) any(!is.na(match(diamnodes, bc$components[[x]])))) # composant sur le diamètre
    toremove <- unlist(bc$component_edges[w[which(!dansdiam)]])
    if(!is.null(toremove)) {
      sp.nobc <- sp.nodeg1[-toremove, ]
      g.nobc <- subgraph.edges(g.nodeg1, E(g.nodeg1)[-toremove])
    } else {
      sp.nobc <- sp.nodeg1
      g.nobc <- g.nodeg1
    }    
  } else {
    toremove = NULL
    sp.nobc <- sp.nodeg1
    g.nobc <- g.nodeg1
  }  
  
  
  # TRAITEMENT DES ARCS MULTIPLES
  wmult <- which(is.multiple(g.nobc))
  
  if(length(wmult) > 0) {
    vpairs <- get.edges(g.nobc, E(g.nobc)[wmult])
    l <- lapply(1:nrow(vpairs), function(i) {
      es <- E(g.nobc)[from(vpairs[i, 1])][to(vpairs[i, 2])]
      
      props = sapply(1:length(es), function(j) {     
        e = es[j]      
        if(gIntersects(sp.nobc[e, ], couchecompare)) {
          lg = gLength(sp.nobc[e, ])
          overlapped.lg = gLength(gIntersection(sp.nobc[e, ], couchecompare))
          prop <- overlapped.lg / lg
        } else {
          prop = 0
        }
        return(prop)
      })      
      sel = which.min(props)
      
      return(es[[sel]])   
    })
    
    toremove = unlist(l)
    sp.nomult = sp.nobc[-toremove, ]
    g.nomult <- subgraph.edges(g.nobc, E(g.nobc)[-toremove])
  } else {
    sp.nomult = sp.nobc
    g.nomult <- g.nobc
  }  
  
  # FUSION DES COUCHES
  sp.nomultM <- mergeThenSplitLine(sp.nomult)
  g.nomultM <- construireReseau(sp.nomultM)
  
  # DETERMINATION DES PIVOTS
  cutpoints <- cutPoints(g.nomultM, sp.nomultM)
  
  # START et END NODES
  if(length(cutpoints) > 0) { 
    pivots = c(cutpoints, as.numeric(farthest.nodes(g.nomultM)$vertices))
  } else {
    pivots = as.numeric(farthest.nodes(g.nomultM)$vertices) 
  }  
    
  # DETERMINATION DU MEILLEUR CHEMIN  
  out = list()
  for (i in 1:length(pivots)) {
    print(paste(i, "/", length(pivots)))

    vfrom = pivots[[i]]
    tos = findNextNodes(g.nomultM, vfrom, pivots)$ends
    tos = tos[which(tos > as.numeric(vfrom))]
    
    if(length(tos) > 0) {
      for (j in 1:length(tos)) {
        vto <- V(g.nomultM)[tos[[j]]]
        
        # TOUS LES CHEMINS OU LE PLUS COURT SI TIME OUT
        timea <- Sys.time()
        
        setTimeLimit(cpu=30) # on augmente la taille 
        
        alls = tryCatch(
          expr = {
            evalWithTimeout({all_simple_paths(g.nomultM, vfrom, vto)}, timeout = 15)
          }, TimeoutException = function(ex) {
            cat("timeout\n") 
            return(get.shortest.paths(g.nomultM, vfrom, vto)[[1]])
          }, error = function(e) {
            cat("timeout CPU\n")
            return(get.shortest.paths(g.nomultM, vfrom, vto)[[1]])
          }
        )
        
        timeb <- Sys.time()
        print(difftime(timeb, timea))
                
        # CALCUL DU MEILLEUR CHEMIN
        if(length(alls) == 1) {
          edges <- edgeFromVertices(g.nomultM, alls[[1]])              
        } else {          
#           dists <- sapply(alls, function(p) {
#             e <- edgeFromVertices(g.nomultM, p)
#             d <- mean(distanceEntreLignes(sp.nomulT[e, ], couche2)$dist)
#             return(d)
#           })
#           mindist <- which.min(dists)
#           edges <- edgeFromVertices(g.nomultM, alls[[mindist]])                     
          props <- sapply(alls, function(p) {
            e <- edgeFromVertices(g.nomultM, p)
            lg <- gLength(sp.nomultM[e, ])
            if (gIntersects(sp.nomultM[e, ], couchecompare)) {
              overlapped.lg <- gLength(gIntersection(sp.nomultM[e, ], couchecompare))
              prop <- overlapped.lg / lg
            }else{
              prop=0  
            }
            return(prop)
          })
          
          maxoverlap <- which.max(props)
          edges <- edgeFromVertices(g.nomultM, alls[[maxoverlap]])   
          #         out[[i]] <- edges
        } # if
        out <- c(out, list(edges))      
      } # for j 
    } # if
  } # for i
  
  ok <- as.numeric(unique(unlist(out))); ok <- ok[order(ok)]
    
  g.best <- subgraph.edges(g.nomultM, E(g.nomultM)[ok])
  sp.best <- gLineMerge(sp.nomultM[ok, ])
  
#   # NETTOYAGE FINAL
#   deg1 = which(degree(g.best)==1 & is.na(V(g.best)$label))
#   
#   if(length(deg1) > 0) {    
#     toremove <- E(g.best)[inc(V(g.best)[which(V(g.best) %in% deg1)])][[1]]
#     sp.bestnodeg1 <- gLineMerge(sp.best[-toremove, ]) 
#   } else {
#     sp.bestnodeg1 <- gLineMerge(sp.best)
#   }
    
  res = sp.best@lines[[1]]@Lines
  
  return(res)
}

find_all_paths <- function(graph, start, end, mypath=vector()) {
  mypath = append(mypath, start)
  
  if (start == end) {
    return(mypath)
  }
  
  paths = list()
  
  for (node in graph[[start]][[1]]) {
    if (!(node %in% mypath)){
      newpaths <- find_all_paths(graph, node, end, mypath)
      for (newpath in newpaths){
        paths <- append(paths, newpath)
      }
    }
  }
  return(paths)
}

croisements <- function(g, couche1.fus, couche2) {
  int = gIntersection(couche1.fus, couche2)
  if(class(int)=="SpatialCollections") int=int@pointobj
  coords <- coordinates(int)
  m <- match(V(g)$name, apply(coords, 1, paste, collapse=" "))  
  res <- which(!is.na(m))
  return(res)
}

deLaCouche <- function(g, couche1.fus) {
  coords <- getNodes(couche1.fus)
  m <- match(V(g)$name, rownames(coords)) 
  res <- which(!is.na(m))
  return(res)
}

procheExtremitesCouche <- function(g.premierecouche, vertexes, couche2, min=10) {
  
  # EXTREMITES DE LA COUCHE
  startEnd <- findStartEndNodes(couche2)
  
  # COORDONNEES DES NOEUDS
  coords <- cbind(vertexes$x, vertexes$y)
  
  # NOEUDS LES PLUS PROCHES DES EXTREMES
  starts = get.knnx(coords, cbind(startEnd[1,1], startEnd[1,2]), length(vertexes))$nn.index
  ends = get.knnx(coords, cbind(startEnd[2,1], startEnd[2,2]), length(vertexes))$nn.index
  
  # STARTS
  m = match(vertexes[starts], V(g.premierecouche))
  nodes <- sapply(m, function(n) {
    if(length(subcomponent(g.premierecouche, n)) >= min) return(n)
    })
  start=nodes[1]
  
  if(is.null(start[[1]])) {
    start = starts[[1]]    
  }
    
  # ENDS
  m = match(vertexes[ends], V(g.premierecouche))
  nodes <- sapply(m, function(n) {
    if(length(subcomponent(g.premierecouche, n)) >= min) return(n)
  })
  end=nodes[1]
  
  if(is.null(end[[1]])) {
    end = ends[[1]]    
  }
  
  res = as.numeric(c(start, end))
  return(res)
  
}

pasDansLaCouche <- function(g, fusion, couche1, vertexes) {
  
  couche1_buff_1m = gBuffer(couche1, width=1)
  
  l <- lapply(1:length(vertexes), function(i) {    
    v = vertexes[[i]]
    neighs <- neighbors(g, v)
    
    s <- sapply(1:length(neighs), function(j) {
      neigh = neighs[[j]]
      if(neigh$label == "croisement") {
        eid <- get.edge.ids(g, c(v, neigh))
        if(!gContains(couche1_buff_1m, fusion[eid, ])) {
          return(eid)          
        }
      }
    })      
    return(s)
  })
  
  return(unique(unlist(l)))
}