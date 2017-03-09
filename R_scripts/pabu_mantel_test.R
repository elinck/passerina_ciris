#pabu IBD analysis
setwd("~/Dropbox/Passerina_ciris/")
library(sp);library(ggplot2);library(plyr);library(magrittr);library(ape)

loc <- read.csv("specimen_data/map_data_kmeans.csv")
loc.brd <- subset(loc,season=="breeding")
loc.wnt <- subset(loc,season=="nonbreeding")
pts <- SpatialPoints(data.frame(loc.brd$Longitude,loc.brd$Latitude),proj4string = crs(range))
geo_dist <- spDists(pts,pts)

phy <- read.dna("pyrad/outfiles/pabu_c60h60_noCOA.unlinked_snps",format = "sequential")
phy.brd <- read.dna("pyrad/outfiles/pabu_c60h60_breeding.unlinked_snps",format="sequential")
phy.wnt <- read.dna("pyrad/outfiles/pabu_c60h60_nonbreeding.unlinked_snps",format="sequential")
gen_dist <- dist.dna(phy.brd,"K80",as.matrix=T,pairwise.deletion = T)

mantel.test(gen_dist,geo_dist)

df <- data.frame(c(geo_dist),c(gen_dist))
fit <- lm(df$c.gen_dist.~df$c.geo_dist.)
ggplot(data=df,aes(x=c.geo_dist.,y=c.gen_dist.))+theme_minimal()+
  xlab("Geographic Distance (km)")+
  ylab("Genetic Distance")+
  geom_point(size=1.5,shape=1,col="grey60")+stat_smooth(col="black",method = "lm")
  #annotate("text",x=10,y=0,label=paste0("p < 2e-16, R^2 = ",round(summary(fit)$adj.r.squared,3)),cex=3)


