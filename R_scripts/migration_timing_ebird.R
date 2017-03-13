#pabu eBird data analysis:
#q1: Do P. ciris in Louisiana/Mississippi follow trans-gulf or circum-gulf routes from the Yucatán? 
#method1: compare spring arrival on Gulf Coast vs. S. Texas - if circum-gulf, should see S. Texas first; if trans-, Louisiana.
library(ggplot2);library(data.table);library(plyr);library(foreach);library(doMC);library(raster);library(magrittr)
registerDoMC(cores=8)
setwd("~/Dropbox/Passerina_ciris/gbif")
source("/R/ggthemes.R")

#load sampling regions
stex <- shapefile("~/Dropbox/Passerina_ciris/range_shp/S-Texas-pabu.shp") %>% buffer(.,width=0.1)
slou <- shapefile("~/Dropbox/Passerina_ciris/range_shp/S-Lou-pabu.dbf") %>% buffer(.,width=0.1)

#load and trim pabu occurrence reports
gbif <- fread("pabu_gbif_28Feb2017.csv") %>% data.frame()
gbif <- gbif[,c("countrycode","locality","decimallongitude","decimallatitude","day","month","year","institutioncode","basisofrecord","recordnumber")]
gbif$date <- paste(gbif$day,gbif$month,gbif$year,sep="/") %>% as.POSIXlt(.,format="%d/%m/%Y")
gbif$yday <- gbif$date$yday
gbif <- subset(gbif,!is.na(yday) & !is.na(decimallatitude))
gbif <- gbif[,names(gbif)[names(gbif)!="date"]] #drop date format for plyr

#load eBird effort 
effort <- fread("~/Dropbox/ebird_effort/ebd_sampling_latlongdate_Feb2015.txt")
effort$yday <- as.POSIXlt(effort$`OBSERVATION DATE`,format="%Y-%m-%d")$yday
effort$year <- substr(effort$`OBSERVATION DATE`,1,4) %>% as.numeric()
effort.pts <- SpatialPointsDataFrame(data.frame(effort$LONGITUDE,effort$LATITUDE),data=effort,proj4string = crs(stex))
effort.stex <- intersect(effort.pts,stex)
effort.stex <- effort.stex@data
effort.stex$bin <- "Texas"
effort.slou <- intersect(effort.pts,slou)
effort.slou <- effort.slou@data
effort.slou$bin <- "Louisiana"
effort <- rbind(effort.slou,effort.stex)
effort <- subset(effort,year>2004&year<2016)

#find pabu reports in S. Lousiana and S. Texas
gbif.pts <- SpatialPointsDataFrame(data.frame(gbif$decimallongitude,gbif$decimallatitude),proj4string = crs(stex),data=gbif)
gbif.pts.stex <- intersect(gbif.pts,stex)
gbif.pts.slou <- intersect(gbif.pts,slou)
tex <- gbif.pts.stex@data
tex$bin <- "Texas"
lou <- gbif.pts.slou@data
lou$bin <- "Louisiana"
df <- rbind(tex,lou)
df <- subset(df,year>2004&year<2016)

#calculate frequency of pabu reports by julian day
pabuCount <- ddply(df,.(yday,bin),summarize,n.pabu=length(year))
effortCount <- ddply(effort,.(yday,bin),summarize,n.reports=length(year))
freq <- merge(pabuCount,effortCount,by=c("yday","bin"),all.x=T,all.y=F)
freq$freq <- freq$n.pabu/freq$n.reports

#plot frequency histogram
ggplot(data=freq,aes(x=yday,y=freq))+theme_bw()+facet_grid(bin~.)+geom_bar(stat="identity",width=1)+xlim(60,121)
ggplot(data=freq,aes(x=yday,y=freq))+facet_grid(bin~.)+geom_point()+xlim(60,121)+stat_smooth()

#stats on report frequencies
aov(freq$freq~freq$bin+freq$yday) %>% summary() #Is this appropriate? must be better analyses...

###analysis of raw count data
w <- wilcox.test(subset(df,bin=="Louisiana")$yday,subset(df,bin=="Texas")$yday)

maxdays <- ddply(df,.(bin,yday),summarize,n=length(day)) %>% 
            ddply(.(bin),function(e) e$yday[e$n==max(e$n)][1]) 

ggplot(data=df,aes(x=yday,fill=bin))+theme_minimal()+
  facet_grid(bin~.,scales="free_y")+
  xlim(60,152)+
  ggtitle("Passerina ciris\neBird Report Count")+
  xlab("Day of the Year")+ylab("")+
  scale_fill_manual(values = c("grey45","grey80"))+
  geom_histogram(binwidth=1)+
  geom_vline(data=maxdays,aes(xintercept=V1))

#map sub-figure (combine in illustrator)
map <- map_data("world")
states <- map_data("state")
ggplot()+
  theme_blank()+coord_map()+
  xlim(-102,-80)+ylim(18,34)+
  scale_fill_manual(values = c("grey45","grey80"))+
  geom_polygon(data=fortify(slou),aes(x=long,y=lat,group=group),fill="grey50")+
  geom_polygon(data=fortify(stex),aes(x=long,y=lat,group=group),fill="grey75")+
  geom_path(data=states,aes(x=long,y=lat,group=group),col="grey")+
  geom_path(data=map,aes(x=long,y=lat,group=group))
