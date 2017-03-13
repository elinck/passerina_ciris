### perform Evanno test on STRUCTURE replicates

setwd("~/Dropbox/Passerina_ciris/structure_runs/final/raw")
# install devtools package from CRAN
install.packages('devtools',dependencies=T)
library(devtools)

# install pophelper package from GitHub
install_github('royfrancis/pophelper')

# load library for use
library(pophelper)

sfiles <- list.files("~/Dropbox/Passerina_ciris/structure_runs/final/raw")
slist <- readQ(sfiles,filetype="structure")
tabq <- tabulateQ(qlist = slist)
sumq <- summariseQ(tabq)
evanno <- evannoMethodStructure(sumq) #same issue w/ SDs
