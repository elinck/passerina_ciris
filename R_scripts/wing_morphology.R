#wing length analysis
setwd("~/Dropbox/Passerina_ciris/morpho/")
library(data.table)
source("/R/ggthemes.R")
# wing <- read.csv("uwbm_pabu_morphology.csv")
# uwbm <- read.csv("../specimen_data/extra/uwbm_west_brd.csv")
# wing <- merge(wing,uwbm,by="Mus_Number",all.x=T,all.y=F)
# write.csv(wing,"wing_uwbmDB.csv",row.names = F)
# wing <- read.csv("wing_uwbmDB.csv")
# bunting <- read.csv("../specimen_data/Pciris_specimen_data.csv")
# wing <- merge(wing,bunting,by="Mus_Number",all.x=T,all.y=F)

wing <- read.csv("wing_uwbmDB.csv")
wing <- subset(wing,wing$X=="" & !is.na(wing.chord|tarsus.length)) #drop measurements with notes about poor specimen quality (n=10)
str <- fread("../structure_runs/structure_k3_final.txt") #read in structure results
colnames(str) <- c("sampleID","q1","q2","q3")
tmp <- merge(pc,str,by="sampleID")
wing <- merge(wing,tmp,by="Mus_Number",all.x=T,all.y=T) #run passerina_PCA to get object "pc"
wing <- wing[,c("Mus_Number","wing.chord","PC1","q1","Longitude.x","Longitude.y","Latitude.x","Latitude.y","Weight","tarsus.length","MO.x")]
for(i in 1:nrow(wing)){                               #merge coordinates for genotype & morphology samples
  if(is.na(wing[i,"Longitude.x"])){
    wing[i,"Longitude.x"] <- wing[i,"Longitude.y"]
    wing[i,"Latitude.x"] <- wing[i,"Latitude.y"]
  }
}
wing$Longitude <- wing$Longitude.x
wing$Latitude <- wing$Latitude.x
wing <- wing[,-c(5:8)]
wing <- subset(wing,MO.x %in% c(4,5,6) & Longitude>-100)

pcdata <- wing[,c("Mus_Number","wing.chord","tarsus.length")] %>% na.omit()
wingpca <- prcomp(pcdata[,2:3],center=T,scale=T) %>% .$x %>% .[,1]
wingpca <- data.frame(Mus_Number=pcdata$Mus_Number,wing_tarsus_PC1=-wingpca)
wing <- merge(wing,wingpca,by="Mus_Number")

melt.wing <- melt(wing,id.vars = c("Longitude","Latitude"))
melt.wing <- subset(melt.wing,variable %in% c("wing.chord","tarsus.length","wing_tarsus_PC1"))
melt.wing$value <- as.numeric(melt.wing$value)
wing.sum <- ddply(wing,.(Longitude,Latitude),summarize,n=length(wing.chord),wing_tarsus_PC1=mean(na.omit(wing_tarsus_PC1)),wing=mean(na.omit(wing.chord)),tarsus=mean(na.omit(tarsus.length)))

lm(wing$wing.chord~wing$q1) %>% summary
ggplot(data=wing,aes(x=q1,y=wing_tarsus_PC1))+theme_minimal()+
  geom_point()+
  geom_smooth(method="lm",col="grey",fill=NA,lwd=0.5)

lm(wing$wing.chord~wing$Longitude) %>% summary()
ggplot(data=melt.wing,aes(x=Longitude,y=value))+
  theme_minimal()+
  facet_grid(variable~.,scales="free_y")+
  geom_smooth(col="grey",fill=NA,lwd=0.35,size=0.5,method="lm")+
  geom_point(size=.9)

map <- map_data("state")
map2 <- map_data("world")
wing.sum <- subset(wing.sum,!is.na(wing))
ggplot()+coord_map()+theme_minimal()+ylim(25,38)+xlim(-102,-88)+
  theme(legend.position = "right")+
  scale_size_continuous(range=c(1,6),breaks=c(1,4,7,10),labels=c(1,4,7,10))+
  scale_color_gradient(low="grey95",high="grey20")+
  geom_path(data=map2,aes(x=long,y=lat,group=group),col="grey",lwd=0.5)+
  geom_path(data=map,aes(x=long,y=lat,group=group),col="grey",lwd=0.5)+
  geom_point(data=wing.sum,aes(x=Longitude,y=Latitude,size=n,col=wing_tarsus_PC1))


w2 <- subset(wing,!is.na(wing.chord)&!is.na(PC1))
ggplot(w2,aes(x=PC1,y=wing.chord))+geom_point()
  
  