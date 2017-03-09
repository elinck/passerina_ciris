#pabu wing length plots, data via Storer 1951
library(rpsychi)
wing <- read.csv("~/Dropbox/Passerina_ciris/lit/storer_wing_data.csv")
wing$Number <- as.numeric(wing$Number)
anova <- ind.oneway.second(m=wing$wing.mean,sd=wing$wing.sd,n=wing$Number,unbiased=F)
pf(anova$anova.table$F[1],11,368)

wing <- arrange(wing,wing.mean)
wing <- mutate(wing, Locality = factor(Locality,Locality))

ggplot(data=wing)+facet_grid(~season,scales = "free")+theme_bw()+ylab("Wing Length")+
  theme(axis.text=element_text(angle = 45,hjust=1,vjust=1))+
  geom_point(data=wing,aes(x=Locality,y=wing.mean))+
  geom_errorbar(aes(x=Locality,ymin=wing.mean-(wing.sd/sqrt(Number)),ymax=wing.mean+(wing.sd/sqrt(Number))))

