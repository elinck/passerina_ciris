#passerina PCA + plots
library(adegenet);library(plyr);library(ggplot2);library(viridis);library(grid);library(magrittr);library(raster)
setwd("~/Dropbox/Passerina_ciris/")

########### load data #########################
#load geographic/specimen data
map <- map_data("world")
state <- map_data("state")
pabu.palette <- c("#2B59A8","#E94F24","#C9DA2B")
bunting <- read.csv("specimen_data/Pciris_map_data.csv")
bunting <- subset(bunting,sampleID!="COA_CWT097") #remove outgroup

#load range shapefile
range <- shapefile("~/Dropbox/BirdlifeIntnl_Range_Maps/Passerina_ciris_22723957.shp")
range.west.df <- crop(range,c(-113,-84,12,38)) %>% fortify() %>% subset(id != 0)
range.west.df$season <- gsub("1","Nonbreeding",range.west.df$id) %>% gsub("2","Breeding",.) %>% gsub("3","Migration",.)
range.df <- fortify(range)
range.df <- subset(range.df,id != 0)
range.df$season <- gsub("1","Nonbreeding",range.df$id) %>% gsub("2","Breeding",.) %>% gsub("3","Migration",.)
  
#read in SNP data (.str or .gen formats)
seq <- read.structure("./pyrad/outfiles/pabu_c60h60.str",n.ind=96,n.loc=4552,onerowperind = F,col.lab=1,col.pop=0,row.marknames=0,NA.char="-9",ask=F)
seq <- seppop(seq,pop=(rownames(seq@tab) %in% c("COA_CWT097")))[[1]] #remove outgroup

############# cluster and merge locations for better plotting ######### ####################
# 
# pts <- SpatialPoints(data.frame(bunting$Longitude,bunting$Latitude),proj4string = crs(range))
# dist <- spDists(pts,pts) %>% as.dist()
# tmp <- hclust(dist,method="complete")
# tmp <- cutree(tmp,h=10)
# bunting$locality <- tmp
# bunting <- arrange(bunting,locality)
# maxlocality <- 1
# tmp <- bunting[1,]
# for(i in 1:nrow(bunting)){
#   row <- bunting[i,]
#   if(row$locality>maxlocality){
#     row$locality.long <- row$Longitude
#     row$locality.lat <- row$Latitude
#     maxlocality <- row$locality
#   } else if(row$locality==maxlocality){
#     row$locality.long <- subset(tmp,locality==maxlocality)$locality.long[1]
#     row$locality.lat <- subset(tmp,locality==maxlocality)$locality.long[1]
#   } else if(row$locality<maxlocality){
#     row$locality.long <- subset(tmp,locality==maxlocality)$locality.long[1]
#     row$locality.lat <- subset(tmp,locality==maxlocality)$locality.lat[1]
#   }
#   tmp <- rbind(tmp,row)
# }
# bunting <- tmp
#######
###########find k-means clusters for k=2-4, merge with specimen locations
clust.k2 <- find.clusters(seq,n.pca=95,n.clust = 2,choose.n.clust = F)
clust.k2 <- data.frame(sampleID=names(clust.k2$grp),clust.k2=clust.k2$grp)
clust.k3 <-  find.clusters(seq,n.pca=95,n.clust = 3,choose.n.clust = F)
clust.k3 <- data.frame(sampleID=names(clust.k3$grp),clust.k3=clust.k3$grp)
clust.k4 <-  find.clusters(seq,n.pca=95,n.clust = 4,choose.n.clust = F)
clust.k4 <- data.frame(sampleID=names(clust.k4$grp),clust.k4=clust.k4$grp)
clust <- merge(clust.k2,clust.k3,by="sampleID")
clust <- merge(clust,clust.k4,by="sampleID")
bunting <- merge(clust,bunting,by="sampleID")

#subset by geography (E/W) or season (breeding/nonbreeding) 
seq.west <- seppop(seq,pop=rownames(seq@tab) %in% subset(bunting,clust.k2==1)$sampleID)[[1]]
seq.east <- seppop(seq,pop=rownames(seq@tab) %in% subset(bunting,clust.k2==1)$sampleID)[[2]]
seq.brd <- seppop(seq,pop=(rownames(seq@tab) %in% subset(bunting,season=="breeding")$sampleID))[[2]]
seq.wnt <- seppop(seq,pop=(rownames(seq@tab) %in% subset(bunting,season=="breeding")$sampleID))[[1]]

#kmeans on western population only
clust.west <- find.clusters(seq.west,n.pca=72,n.clust = 2,choose.n.clust = F)
clust.west <- data.frame(sampleID=names(clust.west$grp),clust.west=clust.west$grp)
bunting <- merge(bunting,clust.west,by="sampleID",all.x=T)
clust.east <- find.clusters(seq.east,n.pca=72,n.clust=2,choose.n.clust=F)
clust.east <- data.frame(sampleID=names(clust.east$grp),clust.east=clust.east$grp)

#center and replace missing data with mean allele frequency
seq.scaled <- scaleGen(seq,NA.method=c("mean"),center=T,scale=F)
west.scaled <- scaleGen(seq.west,NA.method=c("mean"),center=T,scale=F)
east.scaled <- scaleGen(seq.east,NA.method=c("mean"),center=T,scale=F)
brd.scaled <- scaleGen(seq.brd,NA.method=c("mean"),center=T,scale=F)
wnt.scaled <- scaleGen(seq.wnt,NA.method=c("mean"),center=T,scale=F)

#run PCA, merge output with specimen locations
pca <- prcomp(west.scaled,center=F,scale=F)
screeplot(pca)
pc <- data.frame(pca$x[,1:3])
pc$sampleID <- rownames(pc)
pc <- merge(pc,bunting,by="sampleID")

#plot PCA with sample ID's
ggplot(data=pc,aes(x=PC1,y=PC2,col=clust.k3))+geom_text(aes(label=sampleID))

####################
######## k3 clusters on a map with inset PCA plot
clust_per_locality <- ddply(pc,.(clust.k3,locality.long,locality.lat),summarize,n=length(PC1))
clust_per_locality$n <- as.numeric(clust_per_locality$n)
range.crop <- crop(range,c(-118,-73,7,39))
range.df <- fortify(range.crop)
range.df <- subset(range.df,id != 0)
range.df$season <- gsub("1","Nonbreeding",range.df$id) %>% gsub("2","Breeding",.) %>% gsub("3","Migration",.)

inset <- ggplot(pc,aes(x=PC1,y=PC2,col=clust.k3))+geom_point(size=0.5)+xlab("PC1")+ylab("PC2")+
  scale_color_manual(values = pabu.palette)+
  theme(legend.position = "none",axis.ticks = element_blank(),axis.text=element_blank(),
        axis.title = element_text(size=6),
        plot.background = element_rect(fill = "transparent", colour = NA))

mapPlot <- ggplot()+coord_map()+theme_bw()+ylim(7,39)+xlim(-118,-73.5)+
  scale_size(breaks=c(1,4,7,10),name="Samples",range=c(1,10))+
  scale_color_manual(values = pabu.palette,name="Cluster")+
  scale_shape_discrete(solid = F)+
  scale_fill_manual(name="Season",values = c("grey45","grey85"))+
  theme(panel.grid=element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
        axis.text = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  geom_polygon(data=subset(range.df,season %in% c("Breeding","Nonbreeding")),
               aes(x=long,y=lat,group=group,fill=season),col=NA)+
  geom_path(data=state,aes(x=long,y=lat,group=group),col="grey",lwd=0.25)+
  geom_path(data=map,aes(x=long,y=lat,group=group))+
  geom_point(data=clust_per_locality,aes(x=locality.long,y=locality.lat,size=n,col=clust.k3),alpha=0.85,stroke=0)
  
vp <- viewport(height=0.3,width=0.27,x=0.15,y=0.23)
print(mapPlot)
print(inset,vp=vp)
####

######black and white version#######
clust_per_locality <- ddply(pc,.(clust.k3,locality.long,locality.lat),summarize,n=length(PC1))
clust_per_locality$n <- as.numeric(clust_per_locality$n)
range.crop <- crop(range,c(-118,-73,7,39))
range.df <- fortify(range.crop)
range.df <- subset(range.df,id != 0)
range.df$season <- gsub("1","Nonbreeding",range.df$id) %>% gsub("2","Breeding",.) %>% gsub("3","Migration",.)

inset <- ggplot(pc,aes(x=PC1,y=PC2,col=clust.k3))+geom_point(size=0.5)+xlab("PC1")+ylab("PC2")+
  scale_color_manual(values = c("gray20","gray50","gray80"))+
  theme(legend.position = "none",axis.ticks = element_blank(),axis.text=element_blank(),
        axis.title = element_text(size=6),
        plot.background = element_rect(fill = "transparent", colour = NA))

mapPlot <- ggplot()+coord_map()+theme_bw()+ylim(7,39)+xlim(-118,-73.5)+
  scale_size(breaks=c(1,4,7,10),name="Samples",range=c(1,10))+
  scale_color_manual(values = c("gray20","gray50","gray80"))+
  scale_shape_discrete(solid = F)+
  scale_fill_manual(name="Season",values = c("grey40","grey88"))+
  theme(panel.grid=element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
        axis.text = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  geom_polygon(data=subset(range.df,season %in% c("Breeding","Nonbreeding")),
               aes(x=long,y=lat,group=group,fill=season),col=NA)+
  geom_path(data=state,aes(x=long,y=lat,group=group),col="black",lwd=0.2)+
  geom_path(data=map,aes(x=long,y=lat,group=group),col="black",lwd=0.4)+
  geom_point(data=clust_per_locality,aes(x=locality.long,y=locality.lat,size=n,col=clust.k3),alpha=0.85,stroke=0)

vp <- viewport(height=0.3,width=0.27,x=0.15,y=0.23)
print(mapPlot)
print(inset,vp=vp)
####

########regional maps############
#western breeding population, k2 split, PCA + kmeans clusters on a map
clust_per_locality <- ddply(pc,.(clust.west,locality.long,locality.lat),summarize,n=length(PC1))
inset <- ggplot(pc,aes(x=PC1,y=PC2,col=clust.west))+geom_point(size=0.5)+xlab("PC1")+ylab("PC2")+
  scale_color_manual(values = pabu.palette)+
  theme(legend.position = "none",axis.ticks = element_blank(),axis.text=element_blank(),
        axis.title = element_text(size=6),
        plot.background = element_rect(fill = "transparent", colour = NA))
mapPlot <- ggplot()+coord_map()+theme_bw()+ylim(5,39)+xlim(-115,-85)+
  scale_size(breaks=c(1,4,7,10),name="Samples",range=c(1,10))+
  scale_color_manual(values = pabu.palette,name="Genotype\nCluster")+
  scale_shape_discrete(solid = F)+
  scale_fill_manual(name="Season",values = c("grey45","grey85"))+
  theme(panel.grid=element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
        axis.text = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  geom_polygon(data=subset(range.df,season %in% c("Breeding","Nonbreeding")),
               aes(x=long,y=lat,group=group,fill=season),col=NA)+
  geom_path(data=state,aes(x=long,y=lat,group=group),col="grey",lwd=0.25)+
  geom_path(data=map,aes(x=long,y=lat,group=group))+
  geom_point(data=subset(clust_per_locality,clust.west==1),aes(x=locality.long,y=locality.lat,size=n,col=clust.west),alpha=0.7)+
  geom_point(data=subset(clust_per_locality,clust.west==2),aes(x=locality.long,y=locality.lat,size=n,col=clust.west),alpha=0.7)+
  geom_point(data=subset(clust_per_locality,clust.west==3),aes(x=locality.long,y=locality.lat,size=n,col=clust.west),alpha=0.7)
vp <- viewport(height=0.32,width=0.34,x=0.19,y=0.2)
print(mapPlot)
print(inset,vp=vp)

############################################
###### Plot of color-scaled PC1 by locality
meanPC1 <- ddply(pc,.(locality.long,locality.lat,season),summarize,meanPC1=mean(PC1),n=length(PC1))
ggplot()+coord_map()+theme_bw()+theme(panel.grid=element_blank())+ylim(12,38)+xlim(-110,-84)+
  theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),
        legend.spacing = unit(0,"cm"))+
  #scale_color_gradient2(high="lawngreen",mid = "orangered",low="slateblue4",name="PC1")+
  scale_color_gradientn(colors = c("slateblue4","orangered","yellow","lawngreen"))+
  scale_size_continuous(breaks=c(1,4,7,10))+
  geom_polygon(data=subset(range.west.df,season!="Migration"),
               aes(x=long,y=lat,group=group,fill=season),alpha=0.2)+
  geom_path(data=map,aes(x=long,y=lat,group=group))+
  geom_point(data=meanPC1,aes(x=locality.long,y=locality.lat,col=meanPC1,size=n))


#map with elevation raster
alt <- raster("~/Documents/worldclim/alt_2-5m_bil/alt.bil") %>% crop(c(-120,-65,10,40))
alt <- rasterToPoints(alt)
alt <- data.frame(alt)
mapplot <- ggplot()+coord_cartesian()+theme_bw()+ylim(12,38)+xlim(-116,-75)+
  theme(panel.grid=element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
        axis.text = element_blank(),legend.title=element_text(size=9),legend.text=element_text(size=7))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  scale_fill_gradient(low="grey99",high="grey1",name="Elevation (M)")+
  scale_color_manual(values = pabu.palette,name="Cluster")+
  scale_size_continuous(breaks=c(2,4,6,8),range=c(1,10),name="Samples")+
  geom_raster(data=alt,aes(x=x,y=y,fill=alt))+
  geom_path(data=map,aes(x=long,y=lat,group=group))+
  geom_point(data=clust_per_locality,aes(x=locality.long,y=locality.lat,size=n,col=clust.k3),stroke=0,alpha=0.9)
inset <- ggplot(pc,aes(x=PC1,y=PC2,col=clust.k3))+geom_point(size=0.5)+xlab("PC1")+ylab("PC2")+
  scale_color_manual(values = pabu.palette)+
  theme(legend.position = "none",axis.ticks = element_blank(),axis.text=element_blank(),
        axis.title = element_text(size=6),plot.background = element_rect(fill = "transparent", colour = NA))
vp <- viewport(x=0.145,y=0.155,width=.26,height=.28)
print(mapplot)
print(inset,vp=vp)




  