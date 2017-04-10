#warning: written for pyRAD vcf's, partial test on STACKS vcf, no guarantees otherwise. 
#
#vcf2dadi converts a vcf file to the file format required by dadi (basically, allele counts by population)
#
#required parameters: 
#vcf=path to vcf file (tested with vcf's from pyrad; limited tests with stacks)
#populations=path to a two-column text file with first column listing sample ID's (must exactly match vcf) 
#            and second listing population assignment as an integer (same format as STACKS popmap files). Any sep should work. 
#
#optional parameters:
#outgroup=text string giving sample ID of a single outgroup sample (default NULL)
#min.maf=minimum minor allele frequency (default NULL). Not recommended for use with dadi. 
#outfile=file name for output (default NULL)
#minsamples=numeric vector of minimum individuals per population to retain a locus. Currently only works for 2 populations. 
#bootstrap=produce bootstrap datasets? (default F). NOTE: bootstraps over SNP's, so only independent if sampling 1 SNP per locus. 
#n.bootstrap=number of bootstrapped datasets to write (default NULL)
#bootstrap.dir=directory path to output bootstrapped datasets (default NULL)
#biallelic=remove non-biallelic loci? (default True)
#cores=number of processing cores to use (default 1)
#oneSNPperLocus=sample a single SNP per locus? (default F). NOTE: bootstraps over SNP's, so only independent if sampling 1 SNP per locus.

#utility functions to parse vcf files
makeDadiRow_k1 <- function(i,vcf,outgroup,outcol,populations,min.maf,allele.freq){
  row <- vcf[i,]
  if(is.null(outgroup)){
    outseq <- "---"
    refseq <- "---"
  } else {
    outseqcol <- outcol[i] %>% strsplit(.,"|") %>% unlist() %>% .[1] %>% as.numeric()+4
    outseq <- paste0("-",row[outseqcol],"-")
    refseq <- paste0("-",row$REF,"-")
  }
  allele1 <- row$REF
  allele2 <- row$ALT
  genotypes <- row[10:ncol(row)] %>% lapply(.,function(e) strsplit(e,":")[[1]][1])
  if(allele.freq==F){
    genotype.counts <- ddply(populations,.(pop),function(e){
      pop <- genotypes[names(genotypes) %in% e$id]
      n.allele1 <- sum(str_count(pop,"0")) 
      n.allele2 <- sum(str_count(pop,"1"))
      c(n.allele1=n.allele1,n.allele2=n.allele2)
    })
  } else{
    genotype.counts <- ddply(populations,.(pop),function(e){
      pop <- genotypes[names(genotypes) %in% e$id]
      n.allele1 <- sum(str_count(pop,"0")) 
      n.allele2 <- sum(str_count(pop,"1"))
      c(n.allele1=n.allele1/(n.allele1+n.allele2),n.allele2=n.allele1/(n.allele1+n.allele2))
    })
  }
  #concatenate output for one SNP
  dadiRow <- c(outseq,
               refseq, 
               allele1,genotype.counts$n.allele1[1],
               allele2,genotype.counts$n.allele2[1],
               row$X.CHROM,row$POS)
  
  if(max(genotype.counts[2:3])!=0){ #remove loci that are all missing data after outgroup filter
    if(!is.null(min.maf)){ #optional minor allele filter
      maf <- min(sum(genotype.counts$n.allele1),sum(genotype.counts$n.allele2))/(sum(genotype.counts$n.allele1)+sum(genotype.counts$n.allele2))
      if(maf>=min.maf){
        return(dadiRow)
      }
    } else if(is.null(min.maf)){
      return(dadiRow)
    }
  }
}
makeDadiRow_k2 <- function(i,vcf,outgroup,outcol,populations,min.maf,allele.freq){
  row <- vcf[i,]
  if(is.null(outgroup)){
    outseq <- "---"
    refseq <- "---"
  } else {
    outseqcol <- outcol[i] %>% strsplit(.,"|") %>% unlist() %>% .[1] %>% as.numeric()+4
    outseq <- paste0("-",row[outseqcol],"-")
    refseq <- paste0("-",row$REF,"-")
  }
  allele1 <- row$REF
  allele2 <- row$ALT
  genotypes <- row[10:ncol(row)] %>% lapply(.,function(e) strsplit(e,":")[[1]][1]) #get genotypes, remove comments/qual info
  if(allele.freq==F){
    genotype.counts <- ddply(populations,.(pop),function(e){
      pop <- genotypes[names(genotypes) %in% e$id] #e$id is the vector of sample ID's for a population
      n.allele1 <- sum(str_count(pop,"0")) 
      n.allele2 <- sum(str_count(pop,"1"))
      c(n.allele1=n.allele1,n.allele2=n.allele2)
    })
  } else{
    genotype.counts <- ddply(populations,.(pop),function(e){ #not sure this will work with downstream filters
      pop <- genotypes[names(genotypes) %in% e$id]
      n.allele1 <- sum(str_count(pop,"0")) 
      n.allele2 <- sum(str_count(pop,"1"))
      c(n.allele1=n.allele1/(n.allele1+n.allele2),n.allele2=n.allele1/(n.allele1+n.allele2))
    })
  }
  #concatenate output for one SNP
  dadiRow <- c(outseq,
               refseq, 
               allele1,genotype.counts$n.allele1[1],genotype.counts$n.allele1[2],
               allele2,genotype.counts$n.allele2[1],genotype.counts$n.allele2[2],
               row$X.CHROM,row$POS)
  
  if(max(genotype.counts[2:3])!=0){ #remove loci that are all missing data after outgroup filter
    if(!is.null(min.maf)){ #optional minor allele filter
      maf <- min(sum(genotype.counts$n.allele1),sum(genotype.counts$n.allele2))/(sum(genotype.counts$n.allele1)+sum(genotype.counts$n.allele2))
      if(maf>=min.maf){
        return(dadiRow)
      }
    } else if(is.null(min.maf)){
      return(dadiRow)
    }
  }
}
makeDadiRow_k3 <- function(i,vcf,outgroup,outcol,populations,min.maf,allele.freq){
  row <- vcf[i,]
  if(is.null(outgroup)){
    outseq <- "---"
    refseq <- "---"
  } else {
    outseqcol <- outcol[i] %>% strsplit(.,"|") %>% unlist() %>% .[1] %>% as.numeric()+4
    outseq <- paste0("-",row[outseqcol],"-")
    refseq <- paste0("-",row$REF,"-")
  }
  allele1 <- row$REF
  allele2 <- row$ALT
  genotypes <- row[10:ncol(row)] %>% lapply(.,function(e) strsplit(e,":")[[1]][1])
  if(allele.freq==F){
    genotype.counts <- ddply(populations,.(pop),function(e){
      pop <- genotypes[names(genotypes) %in% e$id]
      n.allele1 <- sum(str_count(pop,"0")) 
      n.allele2 <- sum(str_count(pop,"1"))
      c(n.allele1=n.allele1,n.allele2=n.allele2)
    })
  } else{
    genotype.counts <- ddply(populations,.(pop),function(e){
      pop <- genotypes[names(genotypes) %in% e$id]
      n.allele1 <- sum(str_count(pop,"0")) 
      n.allele2 <- sum(str_count(pop,"1"))
      c(n.allele1=n.allele1/(n.allele1+n.allele2),n.allele2=n.allele1/(n.allele1+n.allele2))
    })
  }
  #concatenate output for one SNP
  dadiRow <- c(outseq,
               refseq, 
               allele1,genotype.counts$n.allele1[1],genotype.counts$n.allele1[2],genotype.counts$n.allele1[3],
               allele2,genotype.counts$n.allele2[1],genotype.counts$n.allele2[2],genotype.counts$n.allele2[3],
               row$X.CHROM,row$POS)
  
  if(max(genotype.counts[2:3])!=0){ #remove loci that are all missing data after outgroup filter
    if(!is.null(min.maf)){ #optional minor allele filter
      maf <- min(sum(genotype.counts$n.allele1),sum(genotype.counts$n.allele2))/(sum(genotype.counts$n.allele1)+sum(genotype.counts$n.allele2))
      if(maf>=min.maf){
        return(dadiRow)
      }
    } else if(is.null(min.maf)){
      return(dadiRow)
    }
  }
}

#main function
vcf2dadi <- function(vcf,populations,allele.freq=F,outgroup=NULL,min.maf=NULL,outfile=NULL,minsamples=NULL,
                     bootstrap=F,n.bootstraps=NULL,bootstrap.dir=NULL,biallelic=T,cores=1,oneSNPperLocus=F){
  library(data.table);library(magrittr);library(stringr);library(plyr);library(R.utils);
  library(foreach);library(doMC)
  registerDoMC(cores=cores)
  
  print("reading in data")
  #read in data
  vcf <- data.frame(fread(vcf))
  nSNP <- nrow(vcf)
  nLoci <- nlevels(factor(vcf$X.CHROM))
  
  #remove sites w>2 alleles
  if(biallelic==T){
    print("removing non-biallelic sites")
    vcf <- vcf[-grep(",",vcf$ALT),] 
  }
  
  #filter for only sites with outgroup present
  if(!is.null(outgroup)){
    print("removing sites not sequenced in outgroup")
    outcolnum <- grep(outgroup,names(vcf))
    outcol <- vcf[,outcolnum]
    no_outgroup_sites <- grep("\\.",outcol)
    vcf <- vcf[-no_outgroup_sites,-outcolnum]
    outcol <- outcol[-no_outgroup_sites]
  }
  
  #read in populations file
  populations <- read.delim(populations,header=F)
  colnames(populations) <- c("id","pop")
  nPop <- nlevels(factor(populations$pop))

  print("processing vcf")
  #convert vcf to dadi
  if (nPop==1){
    dadiOut <- foreach(i=1:nrow(vcf),.combine = rbind) %dopar% makeDadiRow_k1(i,vcf,outgroup,outcol,populations,min.maf,allele.freq)
    colnames(dadiOut) <- c("OutgroupSeq","IngroupSeq","Allele1","pop1","Allele2","pop1","Locus","Site")
  } else if (nPop==2){
    dadiOut <- foreach(i=1:nrow(vcf),.combine = rbind) %dopar% makeDadiRow_k2(i,vcf,outgroup,outcol,populations,min.maf,allele.freq)
    colnames(dadiOut) <- c("OutgroupSeq","IngroupSeq","Allele1","pop1","pop2","Allele2","pop1","pop2","Locus","Site")
  } else if (nPop==3){
    dadiOut <- foreach(i=1:nrow(vcf),.combine = rbind) %dopar% makeDadiRow_k3(i,vcf,outgroup,outcol,populations,min.maf,allele.freq)
    colnames(dadiOut) <- c("OutgroupSeq","IngroupSeq","Allele1","pop1","pop2","pop3","Allele2","pop1","pop2","pop3","Locus","Site")
    
  }
  out <- dadiOut
  
  #filter for minimum sample coverage (currently only works at k=2)
  if(!is.null(minsamples)){
    print(paste("removing sites not present in at least ",minsamples," samples in", levels(factor(populations$pop))))
    nloci_pre_minsamples <- nlevels(factor(out[,ncol(out)-1]))
    sampleCoverage <- apply(dadiOut,1,function(e) list((as.numeric(e[4])+as.numeric(e[7]))/2,(as.numeric(e[5])+as.numeric(e[8]))/2))
    keep <- lapply(sampleCoverage,function(e) !any(unlist(e) < minsamples)) %>% unlist()
    out <- out[keep,]
  }
  
  #take a random snp from each locus
  if(oneSNPperLocus==T){ 
    print("sampling one SNP per locus")
    pb <- ProgressBar(max=nlevels(factor(out[,ncol(out)-1])),ticks=10,stepLength=100/nlevels(factor(out[,ncol(out)-1])),newlineWhenDone=F)
    tmp <- out[0,]
    for(i in levels(factor(out[,ncol(out)-1]))){
      locus <- out[out[,ncol(out)-1]==i,]
      if(class(locus)=="matrix"){
        tmp <- rbind(tmp,locus[sample(1:length(locus[1]),1),])
      } else {
        tmp <- rbind(tmp,locus)
      }
      increase(pb)
    }
    out <- tmp
  }
  
  #bootstrap datasets (note this bootstraps over snp's, which are only independent if oneSNPperLocus=T)
  if(bootstrap==T){
    print("writing bootstrapped datasets")
    pb <- ProgressBar(max=n.bootstraps,ticks=10,stepLength=100/n.bootstraps,newlineWhenDone=F)
    pwd <- getwd()
    setwd(bootstrap.dir)
    for(i in 1:n.bootstraps){
      boot <- out[sample(1:nrow(out),nrow(out),replace = T),] 
      boot[,ncol(boot)-1] <- 1:nrow(boot) #change to unique Locus+Site bc dadi removes duplicate SNP's based on these columns
      boot[,ncol(boot)] <- 1:nrow(boot)
      write.table(boot,paste0("bootstrap_",i,".dadi"),sep="\t",quote = F,row.names = F)
      increase(pb)
    }
    setwd(pwd)
  }
  
  #final message
  print(paste("printed",nrow(out)," (out of",nSNP,"total)"," snp's from ",nLoci,"loci"))
  
  #write to file
  if(!is.null(outfile)){
    write.table(out,outfile,sep="\t",quote = F,row.names = F)
  }
  
  return(out)
  
}


