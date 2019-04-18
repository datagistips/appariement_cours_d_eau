library(rgdal)

setwd("D:/WKSP/150916_SBEP_JGL_COURS_DEAU/")

f = readOGR("OUT", "bdtopo_frayere_appariee")

# UTILISATION DES VARIABLES propLong, distance
# http://www.statmethods.net/advstats/cluster.html
# https://onlinecourses.science.psu.edu/stat857/node/125

mydata = f@data[, c("proplong", "nlignes", "distance")]

wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))

# for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
#                                      centers=i)$withinss)
# 
# plot(1:15, wss, type="b", xlab="Number of Clusters",
#      ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(mydata, 5) # 5 cluster solution
# get cluster means 
aggregate(mydata,by=list(fit$cluster), FUN=mean)
# PLOT
# plot(pc.comp1, pc.comp2,col=cl$cluster)
# append cluster assignment
df = data.frame(id=f$id, classe=fit$cluster)

write.csv(df, "OUT_QUALITE/cluster_frayere.csv")