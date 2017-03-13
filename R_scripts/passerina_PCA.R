#passerina PCA + plots
library(adegenet);library(plyr);library(ggplot2);library(viridis);library(grid);library(magrittr);library(raster)
library(foreach);library(doMC);library(stringr)
registerDoMC(cores=8)
setwd("~/Dropbox/Passerina_ciris/")

###############################################
########### load data #########################
###############################################
map <- map_data("world")
state <- map_data("state")
pabu.palette <- c("#2B59A8","#E94F24","#C9DA2B",viridis(5))
range <- shapefile("~/Dropbox/BirdlifeIntnl_Range_Maps/Passerina_ciris_22723957.shp")

#SNP data
seq <- read.structure("./structure_runs/pabu_c48_full.str",n.ind=95,n.loc=3615,onerowperind = F,col.lab=1,col.pop=0,row.marknames=0,NA.char="-9",ask=F)

#specimen data
bunting <- read.csv("specimen_data/Pciris_specimen_data.csv")

#######################################################################
#################### Data Cleaning & Prep #############################
#######################################################################
#find k-means clusters for k=2-6, merge with specimen locations for plotting
clust.k2 <- find.clusters(seq,n.pca=95,n.clust = 2,choose.n.clust = F)
clust.k3 <- find.clusters(seq,n.pca=95,n.clust = 3,choose.n.clust = F)
clust.k4 <- find.clusters(seq,n.pca=95,n.clust = 4,choose.n.clust = F)
clust.k5 <- find.clusters(seq,n.pca=95,n.clust = 5,choose.n.clust = F)
clust.k6 <- find.clusters(seq,n.pca=95,n.clust = 6,choose.n.clust = F)
clust <- cbind(sampleID=rownames(seq@tab),clust.k2=unname(clust.k2$grp),clust.k3=unname(clust.k3$grp),
               clust.k4=unname(clust.k4$grp),clust.k5=unname(clust.k5$grp),clust.k6=unname(clust.k6$grp)) %>% data.frame()
bunting <- merge(clust,bunting,by="sampleID")

#subset by geography (E/W) or season (breeding/nonbreeding) 
seq.west <- seppop(seq,pop=rownames(seq@tab) %in% subset(bunting,k2pop=="W")$sampleID)[[2]] #[[2]] selects condition==T
seq.east <- seppop(seq,pop=rownames(seq@tab) %in% subset(bunting,k2pop=="E")$sampleID,drop=T)[[2]] #NOTE drop=T required for c48. Loci will not share #/order with seq.west
seq.brd <- seppop(seq,pop=(rownames(seq@tab) %in% subset(bunting,season=="breeding")$sampleID))[[2]]
seq.wnt <- seppop(seq,pop=(rownames(seq@tab) %in% subset(bunting,season=="nonbreeding")$sampleID))[[2]]
seq.west.brd <- seppop(seq,pop=rownames(seq@tab) %in% subset(bunting,k2pop=="W" & season=="breeding")$sampleID)[[2]]
seq.west.wnt <- seppop(seq,pop=rownames(seq@tab) %in% subset(bunting,k2pop=="W" & season=="nonbreeding")$sampleID)[[2]]

#kmeans on e/w populations only
clust.west <- find.clusters(seq.west,n.pca=72,n.clust = 2,choose.n.clust = F)
clust.west <- data.frame(sampleID=names(clust.west$grp),clust.west=clust.west$grp)
bunting <- merge(bunting,clust.west,by="sampleID",all.x=T)
clust.east <- find.clusters(seq.east,n.pca=72,n.clust=2,choose.n.clust=F)
clust.east <- data.frame(sampleID=names(clust.east$grp),clust.east=clust.east$grp)
bunting <- merge(bunting,clust.east,by="sampleID",all.x=T)

#center and replace missing data with mean allele frequency
seq.scaled <- scaleGen(seq,NA.method=c("mean"),center=T,scale=F)
west.scaled <- scaleGen(seq.west,NA.method=c("mean"),center=T,scale=F)
east.scaled <- scaleGen(seq.east,NA.method=c("mean"),center=T,scale=F)
brd.scaled <- scaleGen(seq.brd,NA.method=c("mean"),center=T,scale=F)
wnt.scaled <- scaleGen(seq.wnt,NA.method=c("mean"),center=T,scale=F)

####################################################################
########################### PCA & DAPC #############################
####################################################################
#run PCA, merge output with specimen locations
pca <- prcomp(seq.scaled,center=F,scale=F)
screeplot(pca)
pc <- data.frame(pca$x[,1:3])
pc$sampleID <- rownames(pc)
pc <- merge(pc,bunting,by="sampleID")

#plot PCA with sample ID's
ggplot(data=pc,aes(x=PC1,y=PC2,col=clust.k3))+geom_text(aes(label=sampleID))

#linear regression of PC1 against longitude
tmp <- subset(pc,!is.na(clust.west))
model <- lm(tmp$PC1~tmp$Longitude)
summary(model)
ggplot(data=tmp,aes(y=PC1,x=Longitude))+theme_minimal()+geom_point()+stat_smooth(method="lm")

#check for correlation bw missing data and PC1 score (at c60: no for PC1, p=0.41; yes for PC2, p<1.9e-5 R^2=0.17)
a <- fread("~/Dropbox/Passerina_ciris/pyrad/outfiles/pabu_c48.unlinked_snps",header=T) #unlinked SNP's from same pyrad step 7 run as seq file. 
nchar <- as.numeric(names(a)[2])
colnames(a) <- c("sample","seq")
b <- lapply(a$seq,FUN=function(e) str_count(e,"N"))
c <- as.numeric(b)/nchar
md <- data.frame(sample=a$sample,md=c) %>% cbind(pc=pc$PC1)
lm(abs(md$pc)~md$md) %>% summary()
plot(y=abs(md$pc),x=md$md)

#dapc
dapc <- dapc(seq,pop=bunting$clust.k3)

#DAPC x-validation
dapc.xval <- function(){
  train <- dlply(data.frame(sampleID=names(clust.k2$grp),clust.k2=clust.k2$grp),
                 .(clust.k2),function(e) sample(e$sampleID,10)) %>% unlist()
  predict <- subset(data.frame(sampleID=names(clust.k2$grp),clust.k2=clust.k2$grp)
                    ,sampleID %in% train==F)$sampleID
seq.train <- seppop(seq,pop=rownames(seq@tab) %in% train)[[2]] 
  seq.predict <- seppop(seq,pop=rownames(seq@tab) %in% predict)[[2]]
  dapc.obs <- dapc(seq.train,
                   pop=subset(bunting,sampleID %in% train)$clust.k2,
                   n.pca=length(subset(bunting,sampleID %in% train)$clust.k2)/3-1,
                   n.da=2)
  dapc.pred <- predict(dapc.obs,seq.predict)
  comp <- data.frame(sampleID=rownames(dapc.pred$ind.scores),grp=dapc.pred$assign)
  comp <- merge(comp,bunting,by="sampleID",all.x=T,all.y=F)
  1-nrow(subset(comp,grp!=clust.k2))/length(predict)
}
pc_assign_match <-foreach(i=1:1000,.combine="c") %dopar% dapc.xval()
mean(pc_assign_match) # % assigned to correct cluster in DAPC

#xval results summary
xval_c48 <- data.frame(k=2:6,match=c(0.989375,0.9356042),md=round(48/95,2))
xval_c60 <- data.frame(k=2:6,match=c(0.9812,0.9224,0.7103,0.6871,0.5664),md=round(60/95,2))
xval_c85 <- data.frame(k=2:6,match=c(0.9545208,0.8667292,0.6827143,0.6005102,0.60729),md=round(85/95,2))
xval_c91 <-  data.frame(k=2:6,match=c(0.6926,0.650375,0.50733,0.3817,0.397),md=round(91/95,2)) #redo this one
xval_c95 <- data.frame(k=2:6,match=c(0.902,0.7993542,0.6817083,0.5253673,0.4982041),md=round(95/95,2))
xval <- rbind(xval_c60,xval_c85,xval_c91,xval_c95)
ggplot(data=xval,aes(x=k,y=match*100,col=factor(md)))+
  theme_bw()+
  scale_color_manual(values = brewer.pal(name="Set1",n=4),name="Minimum\nSample\nCoverage (%)")+
  geom_point()+geom_line()

#DAPC on breeding range, predict wintering samples
dapc.brd <- dapc(seq.brd,
                 pop=subset(bunting,season=="breeding")$clust.k2,
                 n.pca=length(subset(bunting,season=="breeding")$clust.k2)/3-1,
                 n.da=2)
pred.wnt <- predict(dapc.brd,seq.wnt)
ld1.brd <- data.frame(dapc.brd$ind.coord,grp=dapc.brd$grp,season="breeding")
ld1.wnt <- data.frame(pred.wnt$ind.scores,grp=pred.wnt$assign,season="nonbreeding")
ld1 <- rbind(ld1.brd,ld1.wnt)
ld1$sampleID <- rownames(ld1)

assigncomp <- merge(ld1,bunting,by="sampleID")[,c("sampleID","grp","clust.k2","season.y")]
assigncomp$grp <- as.numeric(as.character(assigncomp$grp))
assigncomp$clust.k2 <- as.numeric(as.character(assigncomp$clust.k2))
assigncomp$diff <- assigncomp$clust.k2 - assigncomp$grp
1-sum(assigncomp$diff!=0)/nrow(pred.wnt$ind.scores)# number of conflicting assignments

ggplot(data=pc,aes(x=LD1.y,fill=xpred_grp))+theme_minimal()+facet_grid(~season)+
  scale_fill_manual(values = pabu.palette)+
  geom_density(alpha=0.7)+
  geom_point(aes(y=0),shape=124,size=4)

#train PCA on breeding range & project wintering range (western samples only)
# pca.brd <- prcomp(seq.scaled[rownames(seq.scaled)  %in% intersect(rownames(brd.scaled),rownames(west.scaled)),],center=F,scale=F)
# pc.wnt <- predict(pca.brd,seq.scaled[rownames(seq.scaled)  %in% intersect(rownames(wnt.scaled),rownames(west.scaled)),])
# pc1 <- data.frame(pca.brd$x[,1:3],pcatype="trained")
# pc2 <- data.frame(pc.wnt[,1:3],pcatype="predicted")
# pc <- rbind(pc1,pc2)
# pc$sampleID <- rownames(pc)
# pc <- merge(pc,bunting,by="sampleID")
# ggplot(data=pc)+theme_bw()+facet_grid(~pcatype)+
#   #geom_label(data=pc,aes(x=PC1,y=PC2,col=clust.k3,label=sampleID))
#   geom_point(data=pc,aes(x=PC1,y=PC2,col=clust.k3))

###################################################################
############################ Plots ################################
###################################################################
#k3 clusters on a map with inset PCA plot
clust_per_locality <- ddply(pc,.(clust.k3,locality.long,locality.lat),summarize,n=length(PC1))
clust_per_locality$n <- as.numeric(clust_per_locality$n)
range.crop <- crop(range,c(-118,-73,7,39))
range.df <- fortify(range.crop)
range.df <- subset(range.df,id != 0)
range.df$season <- gsub("1","Nonbreeding",range.df$id) %>% gsub("2","Breeding",.) %>% gsub("3","Migration",.)
inset <- ggplot(pc,aes(x=PC1,y=PC2,col=clust.k3))+xlab("PC1")+ylab("PC2")+
  scale_color_manual(values = viridis(3))+
  theme(panel.border = element_rect(color="black",fill=NA),
        legend.position="none",
        panel.background = element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_text(size=7))+
  geom_vline(aes(xintercept=0),lwd=0.25)+
  geom_hline(aes(yintercept=0),lwd=0.25)+
  geom_point(size=.75)
 
mapPlot <- ggplot()+coord_map()+theme_bw()+ylim(7,39)+xlim(-118,-73.5)+
  scale_size(breaks=c(1,4,7,10),name="Samples",range=c(1.5,10))+
  scale_color_manual(values = viridis(3),name="Cluster")+
  scale_shape_discrete(solid = F)+
  scale_fill_manual(name="Season",values = c("grey50","grey85"))+
  theme(panel.grid=element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
        axis.text = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  geom_polygon(data=subset(range.df,season %in% c("Breeding","Nonbreeding")),
               aes(x=long,y=lat,group=group,fill=season),col=NA)+
  geom_path(data=state,aes(x=long,y=lat,group=group),col="grey",lwd=0.25)+
  geom_path(data=map,aes(x=long,y=lat,group=group))+
  geom_point(data=clust_per_locality,aes(x=locality.long,y=locality.lat,size=n,col=clust.k3),alpha=0.7,stroke=0)
vp <- viewport(height=0.3,width=0.27,x=0.17,y=0.2)
print(mapPlot)
print(inset,vp=vp)

#Map with localities color-scaled by PC1
meanPC1 <- ddply(pc,.(locality.long,locality.lat),summarize,meanPC1=mean(PC1),n=length(PC1))
range.crop <- crop(range,c(-118,-73,7,39))
range.df <- fortify(range.crop)
range.df <- subset(range.df,id != 0)
range.df$season <- gsub("1","Nonbreeding",range.df$id) %>% gsub("2","Breeding",.) %>% gsub("3","Migration",.)
inset <- ggplot(pc,aes(x=PC1,y=PC2,col=PC1))+geom_point(size=0.5)+xlab("PC1")+ylab("PC2")+
  scale_color_gradient2(low=pabu.palette[1],mid=pabu.palette[2],high=pabu.palette[3],name="Cluster")+
  theme(legend.position = "none",axis.ticks = element_blank(),axis.text=element_blank(),
        axis.title = element_text(size=6),
        plot.background = element_rect(fill = "transparent", colour = NA))
mapPlot <- ggplot()+coord_map()+theme_bw()+ylim(7,39)+xlim(-118,-73.5)+
  scale_size(breaks=c(1,4,7,10),name="Samples",range=c(1,10))+
  scale_color_gradient2(low=pabu.palette[1],mid=pabu.palette[2],high=pabu.palette[3],name="Cluster")+
  scale_shape_discrete(solid = F)+
  scale_fill_manual(name="Season",values = c("grey45","grey85"))+
  theme(panel.grid=element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
        axis.text = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  geom_polygon(data=subset(range.df,season %in% c("Breeding","Nonbreeding")),
               aes(x=long,y=lat,group=group,fill=season),col=NA)+
  geom_path(data=state,aes(x=long,y=lat,group=group),col="grey",lwd=0.25)+
  geom_path(data=map,aes(x=long,y=lat,group=group))+
  geom_point(data=meanPC1,aes(x=locality.long,y=locality.lat,size=n,col=meanPC1),alpha=0.85,stroke=0)
vp <- viewport(height=0.3,width=0.27,x=0.15,y=0.23)
print(mapPlot)
print(inset,vp=vp)

#DAPC cross-prediction map
clust_per_locality <- ddply(pc,.(xpred_grp,locality.long,locality.lat),summarize,n=length(PC1))
clust_per_locality$n <- as.numeric(clust_per_locality$n)
mapPlot <- ggplot()+coord_map()+theme_bw()+ylim(7,39)+xlim(-118,-73.5)+
  scale_size(breaks=c(1,4,7,10),name="Samples",range=c(1.5,10))+
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
  geom_point(data=clust_per_locality,aes(x=locality.long,y=locality.lat,size=n,col=xpred_grp),alpha=0.85,stroke=0)
inset <- ggplot(data=ld1,aes(x=LD1,fill=grp))+theme_blank()+facet_grid(~season)+
  scale_fill_manual(values = pabu.palette)+
  geom_density(alpha=0.7)+
  geom_point(aes(y=0),shape=124,size=2)
vp <- viewport(height=0.27,width=0.32,x=0.2,y=0.17)
print(mapPlot)
print(inset,vp=vp)


#black and white version k3 clusters map
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

#western breeding population, k2 split, PCA + kmeans clusters on a map
range.crop <- crop(range,c(-115,-84,12,38))
range.df <- fortify(range.crop)
range.df <- subset(range.df,id != 0)
range.df$season <- gsub("1","Nonbreeding",range.df$id) %>% gsub("2","Breeding",.) %>% gsub("3","Migration",.)
clust_per_locality <- ddply(pc,.(locality.long,locality.lat,clust.west),summarize,n=length(PC1),meanPC1=mean(PC1))
inset <- ggplot(pc,aes(x=PC1,y=PC2,col=clust.west))+geom_point(size=0.5)+xlab("PC1")+ylab("PC2")+
  #scale_color_gradient2(low=pabu.palette[1],mid=pabu.palette[2],high=pabu.palette[3],name="Mean PC1")+
  scale_color_manual(values = pabu.palette)+
  theme(legend.position = "none",axis.ticks = element_blank(),axis.text=element_blank(),
        axis.title = element_text(size=6),
        plot.background = element_rect(fill = "transparent", colour = NA))
mapPlot <- ggplot()+coord_map()+theme_bw()+ylim(12,38)+xlim(-115,-84)+
  scale_size(breaks=c(1,4,7,10),name="Samples",range=c(1.5,10))+
  #scale_color_gradient2(low=pabu.palette[1],mid=pabu.palette[2],high=pabu.palette[3],name="Mean PC1")+
  scale_color_manual(values = pabu.palette)+
  scale_shape_discrete(solid = F)+
  scale_fill_manual(name="Season",values = c("grey45","grey85"))+
  theme(panel.grid=element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
        axis.text = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  geom_polygon(data=subset(range.df,season %in% c("Breeding","Nonbreeding")),
               aes(x=long,y=lat,group=group,fill=season),col=NA)+
  geom_path(data=state,aes(x=long,y=lat,group=group),col="grey",lwd=0.25)+
  geom_path(data=map,aes(x=long,y=lat,group=group))+
  geom_point(data=clust_per_locality,aes(x=locality.long,y=locality.lat,size=n,col=clust.west),stroke=0,alpha=0.9)
vp <- viewport(height=0.25,width=0.23,x=0.13,y=0.205)
print(mapPlot)
print(inset,vp=vp)

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




  