#bunting maps
library(ggplot2);library(magrittr);library(sp);library(raster);library(rgeos)
setwd("~/Dropbox/Passerina_ciris/")
map <- map_data("world")
bunting <- read.csv("Pciris_map_data.csv")
range <- shapefile("~/Dropbox/BirdlifeIntnl_Range_Maps/Passerina_ciris_22723957.shp")
range.df <- fortify(range)
range.df <- subset(range.df,id != 0)
range.df$season <- gsub("1","Nonbreeding",range.df$id) %>% gsub("2","Breeding",.) %>% gsub("3","Migration",.)

count_per_locality <- ddply(bunting,.(Longitude,Latitude),summarize,n=length(Country))

ggplot()+coord_map()+theme_bw()+theme(panel.grid=element_blank())+ylim(7,40)+xlim(-115,-75)+
  scale_fill_manual(values=c("orangered","yellow","cornflowerblue"),name="Seasonal Status")+
  geom_polygon(data=range.df,aes(x=long,y=lat,group=group,fill=season),alpha=0.5)+
  geom_path(data=map,aes(x=long,y=lat,group=group))+
  geom_point(data=count_per_locality,aes(x=Longitude,y=Latitude,size=n),alpha=0.7)

  