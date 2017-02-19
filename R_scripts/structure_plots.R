### script to quickly visualize fastStructure results 

#setwd("~/Desktop/")
#library(ggplot2)

#standard cjbattey colorscheme
colors <- c("gold","forestgreen","magenta3","orangered","cornflowerblue", "orange","sienna","dodgerblue4")

#plot K2
k2 <- read.table("pabu_fs_K2.2.meanQ")
barplot(t(as.matrix(k2[,2:3])),axes=FALSE,col=colors,border=NA,names=k2[,1],las=2,cex.names=0.6,space=0.05,xpd=FALSE)

#plot K3
k3 <- read.table("pabu_fs_K3.3.meanQ")
barplot(t(as.matrix(k3[,3:4])),axes=FALSE,col=colors,border=NA,names=k3[,1],las=2,cex.names=0.6,space=0.05,xpd=FALSE)
