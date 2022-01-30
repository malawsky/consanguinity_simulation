require(data.table)
require(tidyr)
require(genio)
require(stringr)

### individuals with unrelated parents
consang_outbred <- function(map,ped,nn,chr,output){
  
  geneticmap <- read.table(map)
  ped1 <- as.data.frame(read.table(ped,sep=" ",colClasses = "character"), as.factor = F)
  odd <- seq(1,dim(ped1)[2]*2+1,2)
  even <- seq(2,dim(ped1)[2]*2+2,2)
  relist2 <- list()
  num <- dim(ped1)[2]/2
  for (i in 1:num) {
    relist2[[i]] <- ped1[,c(odd[i],even[i])]
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("rbind.data.frame", consang)
  
  small_ped <- consangdf[-c(1:3),1:nn]
  
  colsplit <- lapply(1:nn, function(i) as.data.frame(small_ped[,i]))
  
  splitgeno <- lapply(1:dim(small_ped)[2], function(i) str_split_fixed(colsplit[[i]][,1], " ", 2))
  
  recombine <- function(i,geno){
    data1 <- geno[1:i,]
    data1 <- cbind.data.frame(as.character(data1[,1]),as.character(data1[,2]))
    colnames(data1) <- c("V1","V2")
    data2 <- geno[(i+1):dim(small_ped)[1],]
    fliped <- cbind.data.frame(as.character(data2[,2]),as.character(data2[,1]))
    colnames(fliped) <- c("V1","V2")
    final <-  rbind.data.frame(data1,fliped)
    return(final)
  }
  
  multirec <- function(geno){
    difference <- diff(geneticmap$V3)*.01
    flips <- rbinom(length(difference),1,difference)
    vectorflips <- match(c(1),flips)
    data <- geno
    data <- cbind.data.frame(as.character(data[,1]),as.character(data[,2]))
    colnames(data) <- c("V1","V2")
    if (!is.na(vectorflips)) {
      for (i in vectorflips) {
        data <- recombine(i,data)
      }
    }
    
    return(data)
  }
  
  randomchrom <- function(df) {
    return(as.data.frame(df[, sample(2, 1)]))
  }
  
  #### make siblings and unrelated
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m1 <- cbind.data.frame(df_1, df_2)   
  
  df1 <- m1[, c(matrix(1:ncol(m1), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-4,4), times =1 ,  each = 4) + rep(c(1,3,2,4),nn/2)
  
  firstbreed <- df1[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  
  consang <- lapply(1:length(relist), function(i) apply(relist[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("cbind.data.frame", consang)
  colnames(consangdf) <- paste0('sample',1:nn)
  
  data <- consangdf
  data <- data.frame(lapply(data, function(x) {gsub("G C", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G G", "0", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("C C", "2", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G T", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T G", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T T", "2", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("G A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A G", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A A", "0", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C G", "1", x)})) 
  
  data[] <- lapply(data, function (i) as.numeric(as.character(i)))
  
  colnames(geneticmap) <- c('chr','id','posg','pos','ref','alt')
  fam <- cbind.data.frame(1,paste0('sibmate',1:nn),0,0,0,0)
  colnames(fam) <- c('fam','id','pat','mat','sex','pheno')
  write_plink(paste0(output,chr,"_outbred"),bim=geneticmap,as.matrix(data) )
}

### individuals with parents that are second cousins

consang_secondcous <- function(map,ped,nn,chr,output){
  
  geneticmap <- read.table(map)
  ped1 <- as.data.frame(read.table(ped,sep=" ",colClasses = "character"), as.factor = F)
  odd <- seq(1,dim(ped1)[2]*2+1,2)
  even <- seq(2,dim(ped1)[2]*2+2,2)
  relist2 <- list()
  num <- dim(ped1)[2]/2
  for (i in 1:num) {
    relist2[[i]] <- ped1[,c(odd[i],even[i])]
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("rbind.data.frame", consang)
  
  small_ped <- consangdf[-c(1:3),1:nn]
  colsplit <- lapply(1:nn, function(i) as.data.frame(small_ped[,i]))
  
  splitgeno <- lapply(1:dim(small_ped)[2], function(i) str_split_fixed(colsplit[[i]][,1], " ", 2))
  
  recombine <- function(i,geno){
    data1 <- geno[1:i,]
    data1 <- cbind.data.frame(as.character(data1[,1]),as.character(data1[,2]))
    colnames(data1) <- c("V1","V2")
    data2 <- geno[(i+1):dim(small_ped)[1],]
    fliped <- cbind.data.frame(as.character(data2[,2]),as.character(data2[,1]))
    colnames(fliped) <- c("V1","V2")
    final <-  rbind.data.frame(data1,fliped)
    return(final)
  }
  
  multirec <- function(geno){
    difference <- diff(geneticmap$V3)*.01
    flips <- rbinom(length(difference),1,difference)
    vectorflips <- match(c(1),flips)
    data <- geno
    data <- cbind.data.frame(as.character(data[,1]),as.character(data[,2]))
    colnames(data) <- c("V1","V2")
    if (!is.na(vectorflips)) {
      for (i in vectorflips) {
        data <- recombine(i,data)
      }
    }
    return(data)
  }
  
  randomchrom <- function(df) {
    return(as.data.frame(df[, sample(2, 1)]))
  }
  
  #### make siblings and unrelated
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m1 <- cbind.data.frame(df_1, df_2)   
  
  df1 <- m1[, c(matrix(1:ncol(m1), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-12,12), times =1 ,  each = 12) + rep(c(1,3,2,4,5,6,7,8,9,10,11,12),nn/6)
  
  firstbreed <- df1[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  ### make cousins
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-12,12), times =1 ,  each = 12) + rep(c(1,5,2,6,3,7,4,8,9,10,11,12),nn/6)
  
  firstbreed <- df[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-12,12), times =1 ,  each = 12) + rep(c(1,9,2,10,3,9,5,11,6,12,7,11),nn/6)
  
  cousinsconsang <- df[,reorder]
  
  relist2 <- list()
  for (i in 1:nn) {
    relist2[[i]] <- as.data.frame(cousinsconsang[,c(odd[i],even[i])])
  }
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist2[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist2[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-12,12), times =1 ,  each = 12) + rep(c(1,7,2,8,3,9,4,10,5,11,6,12),nn/6)
  
  cousinsconsang <- df[,reorder]
  
  relist2 <- list()
  for (i in 1:nn) {
    relist2[[i]] <- as.data.frame(cousinsconsang[,c(odd[i],even[i])])
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("cbind.data.frame", consang)
  colnames(consangdf) <- paste0('sample',1:nn)
  
  data <- consangdf
  data <- data.frame(lapply(data, function(x) {gsub("G C", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G G", "0", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("C C", "2", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G T", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T G", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T T", "2", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("G A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A G", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A A", "0", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C G", "1", x)})) 
  
  data[] <- lapply(data, function (i) as.numeric(as.character(i)))
  colnames(geneticmap) <- c('chr','id','posg','pos','ref','alt')
  fam <- cbind.data.frame(1,paste0('sibmate',1:nn),0,0,0,0)
  colnames(fam) <- c('fam','id','pat','mat','sex','pheno')
  write_plink(paste0(output,chr,"_secondcousins"),bim=geneticmap,as.matrix(data) )
}

### individuals with parents that are first cousins once removed 

consang_cousremoved <- function(map,ped,nn,chr,output){
  
  geneticmap <- read.table(map)
  ped1 <- as.data.frame(read.table(ped,sep=" ",colClasses = "character"), as.factor = F)
  odd <- seq(1,dim(ped1)[2]*2+1,2)
  even <- seq(2,dim(ped1)[2]*2+2,2)
  relist2 <- list()
  num <- dim(ped1)[2]/2
  for (i in 1:num) {
    relist2[[i]] <- ped1[,c(odd[i],even[i])]
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("rbind.data.frame", consang)
  
  small_ped <- consangdf[-c(1:3),1:nn]
  colsplit <- lapply(1:nn, function(i) as.data.frame(small_ped[,i]))
  
  splitgeno <- lapply(1:dim(small_ped)[2], function(i) str_split_fixed(colsplit[[i]][,1], " ", 2))
  
  recombine <- function(i,geno){
    data1 <- geno[1:i,]
    data1 <- cbind.data.frame(as.character(data1[,1]),as.character(data1[,2]))
    colnames(data1) <- c("V1","V2")
    data2 <- geno[(i+1):dim(small_ped)[1],]
    fliped <- cbind.data.frame(as.character(data2[,2]),as.character(data2[,1]))
    colnames(fliped) <- c("V1","V2")
    final <-  rbind.data.frame(data1,fliped)
    return(final)
  }
  
  multirec <- function(geno){
    difference <- diff(geneticmap$V3)*.01
    flips <- rbinom(length(difference),1,difference)
    vectorflips <- match(c(1),flips)
    data <- geno
    data <- cbind.data.frame(as.character(data[,1]),as.character(data[,2]))
    colnames(data) <- c("V1","V2")
    if (!is.na(vectorflips)) {
      for (i in vectorflips) {
        data <- recombine(i,data)
      }
    }
    return(data)
  }
  
  randomchrom <- function(df) {
    return(as.data.frame(df[, sample(2, 1)]))
  }
  
  #### make siblings and unrelated
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m1 <- cbind.data.frame(df_1, df_2)   
  
  df1 <- m1[, c(matrix(1:ncol(m1), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-10,10), times =1 ,  each = 10) + rep(c(1,3,2,4,5,6,7,8,9,10),nn/5)
  
  firstbreed <- df1[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  ### make cousins
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-10,10), times =1 ,  each = 10) + rep(c(1,5,2,6,3,7,4,8,9,10),nn/5)
  
  firstbreed <- df[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  for (i in 1:nn) {
    if (i %% 10 == 5 || i %% 10 == 6 || i %% 10 == 7 || i %% 10 == 8 ){
      df[,i] = firstbreed[,i]
    }
  }
  
  reorder <- rep(seq(0,nn*2-10,10), times =1 ,  each = 10) + rep(c(1,9,2,10,3,9,5,6,7,8),nn/5)
  
  
  cousinsconsang <- df[,reorder]
  
  
  relist2 <- list()
  for (i in 1:nn) {
    relist2[[i]] <- as.data.frame(cousinsconsang[,c(odd[i],even[i])])
  }
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist2[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist2[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-10,10), times =1 ,  each = 10) + rep(c(1,7,2,8,3,9,4,10,5,7),nn/5)
  
  cousinsconsang <- df[,reorder]
  
  relist2 <- list()
  for (i in 1:nn) {
    relist2[[i]] <- as.data.frame(cousinsconsang[,c(odd[i],even[i])])
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("cbind.data.frame", consang)
  colnames(consangdf) <- paste0('sample',1:nn)
  
  data <- consangdf
  data <- data.frame(lapply(data, function(x) {gsub("G C", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G G", "0", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("C C", "2", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G T", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T G", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T T", "2", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("G A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A G", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A A", "0", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C G", "1", x)})) 
  
  data[] <- lapply(data, function (i) as.numeric(as.character(i)))
  colnames(geneticmap) <- c('chr','id','posg','pos','ref','alt')
  fam <- cbind.data.frame(1,paste0('sibmate',1:nn),0,0,0,0)
  colnames(fam) <- c('fam','id','pat','mat','sex','pheno')
  write_plink(paste0(output,chr,"_removedcousins"),bim=geneticmap,as.matrix(data) )
}

### individuals with parents that are avuncular for three generations

consang_avan3rddgen <- function(map,ped,nn,chr,output){
  
  geneticmap <- read.table(map)
  ped1 <- as.data.frame(read.table(ped,sep=" ",colClasses = "character"), as.factor = F)
  odd <- seq(1,dim(ped1)[2]*2+1,2)
  even <- seq(2,dim(ped1)[2]*2+2,2)
  relist2 <- list()
  num <- dim(ped1)[2]/2
  for (i in 1:num) {
    relist2[[i]] <- ped1[,c(odd[i],even[i])]
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("rbind.data.frame", consang)
  
  small_ped <- consangdf[-c(1:3),1:nn]
  colsplit <- lapply(1:nn, function(i) as.data.frame(small_ped[,i]))
  
  splitgeno <- lapply(1:dim(small_ped)[2], function(i) str_split_fixed(colsplit[[i]][,1], " ", 2))
  
  recombine <- function(i,geno){
    data1 <- geno[1:i,]
    data1 <- cbind.data.frame(as.character(data1[,1]),as.character(data1[,2]))
    colnames(data1) <- c("V1","V2")
    data2 <- geno[(i+1):dim(small_ped)[1],]
    fliped <- cbind.data.frame(as.character(data2[,2]),as.character(data2[,1]))
    colnames(fliped) <- c("V1","V2")
    final <-  rbind.data.frame(data1,fliped)
    return(final)
  }
  
  multirec <- function(geno){
    difference <- diff(geneticmap$V3)*.01
    flips <- rbinom(length(difference),1,difference)
    vectorflips <- match(c(1),flips)
    data <- geno
    data <- cbind.data.frame(as.character(data[,1]),as.character(data[,2]))
    colnames(data) <- c("V1","V2")
    if (!is.na(vectorflips)) {
      for (i in vectorflips) {
        data <- recombine(i,data)
      }
    }
    return(data)
  }
  
  randomchrom <- function(df) {
    return(as.data.frame(df[, sample(2, 1)]))
  }
  
  #### make siblings and unrelated
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m1 <- cbind.data.frame(df_1, df_2)   
  
  df1 <- m1[, c(matrix(1:ncol(m1), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-10,10), times =1 ,  each = 10) + rep(c(1,3,2,4,5,6,7,8,9,10),nn/5)
  
  firstbreed1 <- df1[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed1[,c(odd[i],even[i])])
  }
  
  ### make cousins
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  
  reorder <- rep(seq(0,nn*2-10,10), times =1 ,  each = 10) + rep(c(1,2,3,5,4,6,7,8,9,10),nn/5)
  
  
  firstbreed <- df[,reorder]
  
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  
  reorder <- rep(seq(0,nn*2-10,10), times =1 ,  each = 10) + rep(c(1,3,2,4,2,3,7,8,9,10),nn/5)
  
  
  firstbreed <- df[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  
  reorder <- rep(seq(0,nn*2-10,10), times =1 ,  each = 10) + rep(c(1,2,3,4,5,7,6,8,9,10),nn/5)
  
  
  firstbreed <- df[,reorder]
  
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-10,10), times =1 ,  each = 10) + rep(c(5,6,7,9,8,10,7,10,8,9),nn/5)
  
  
  firstbreed <- df[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  for (i in 1:nn) {
    if ( i %% 10 == 9 || i %% 10 == 0){
      df[,i] = firstbreed1[,i]
    }
  }
  
  reorder <- rep(seq(0,nn*2-10,10), times =1 ,  each = 10) + rep(c(1,3,2,4,1,5,2,6,1,7),nn/5)
  
  cousinsconsang <- df[,reorder]
  
  relist2 <- list()
  for (i in 1:nn) {
    relist2[[i]] <- as.data.frame(cousinsconsang[,c(odd[i],even[i])])
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("cbind.data.frame", consang)
  colnames(consangdf) <- paste0('sample',1:nn)
  
  data <- consangdf
  data <- data.frame(lapply(data, function(x) {gsub("G C", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G G", "0", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("C C", "2", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G T", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T G", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T T", "2", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("G A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A G", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A A", "0", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C G", "1", x)})) 
  
  data[] <- lapply(data, function (i) as.numeric(as.character(i)))
  colnames(geneticmap) <- c('chr','id','posg','pos','ref','alt')
  fam <- cbind.data.frame(1,paste0('sibmate',1:nn),0,0,0,0)
  colnames(fam) <- c('fam','id','pat','mat','sex','pheno')
  write_plink(paste0(output,chr,"_3_gen_avuncular"),bim=geneticmap,as.matrix(data) )
}

### individuals with parents that are avuncular for two generations

consang_avan2ndgen <- function(map,ped,nn,chr,output){
  
  geneticmap <- read.table(map)
  ped1 <- as.data.frame(read.table(ped,sep=" ",colClasses = "character"), as.factor = F)
  odd <- seq(1,dim(ped1)[2]*2+1,2)
  even <- seq(2,dim(ped1)[2]*2+2,2)
  relist2 <- list()
  num <- dim(ped1)[2]/2
  for (i in 1:num) {
    relist2[[i]] <- ped1[,c(odd[i],even[i])]
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("rbind.data.frame", consang)
  
  small_ped <- consangdf[-c(1:3),1:nn]
  colsplit <- lapply(1:nn, function(i) as.data.frame(small_ped[,i]))
  
  splitgeno <- lapply(1:dim(small_ped)[2], function(i) str_split_fixed(colsplit[[i]][,1], " ", 2))
  
  recombine <- function(i,geno){
    data1 <- geno[1:i,]
    data1 <- cbind.data.frame(as.character(data1[,1]),as.character(data1[,2]))
    colnames(data1) <- c("V1","V2")
    data2 <- geno[(i+1):dim(small_ped)[1],]
    fliped <- cbind.data.frame(as.character(data2[,2]),as.character(data2[,1]))
    colnames(fliped) <- c("V1","V2")
    final <-  rbind.data.frame(data1,fliped)
    return(final)
  }
  
  multirec <- function(geno){
    difference <- diff(geneticmap$V3)*.01
    flips <- rbinom(length(difference),1,difference)
    vectorflips <- match(c(1),flips)
    data <- geno
    data <- cbind.data.frame(as.character(data[,1]),as.character(data[,2]))
    colnames(data) <- c("V1","V2")
    if (!is.na(vectorflips)) {
      for (i in vectorflips) {
        data <- recombine(i,data)
      }
    }
    return(data)
  }
  
  randomchrom <- function(df) {
    return(as.data.frame(df[, sample(2, 1)]))
  }
  
  #### make siblings and unrelated
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m1 <- cbind.data.frame(df_1, df_2)   
  
  df1 <- m1[, c(matrix(1:ncol(m1), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-8,8), times =1 ,  each = 8) + rep(c(1,3,2,4,5,6,7,8),nn/4)
  
  firstbreed1 <- df1[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed1[,c(odd[i],even[i])])
  }
  
  ### make cousins
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  
  reorder <- rep(seq(0,nn*2-8,8), times =1 ,  each = 8) + rep(c(1,2,3,5,4,6,7,8),nn/4)
  
  
  firstbreed <- df[,reorder]
  
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  
  reorder <- rep(seq(0,nn*2-8,8), times =1 ,  each = 8) + rep(c(1,3,1,4,2,3,7,8),nn/4)
  
  
  firstbreed <- df[,reorder]
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  
  reorder <- rep(seq(0,nn*2-8,8), times =1 ,  each = 8) + rep(c(1,2,3,4,5,7,6,8),nn/4)
  
  
  firstbreed <- df[,reorder]
  
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-8,8), times =1 ,  each = 8) + rep(c(1,5,2,6,3,7,4,8),nn/4)
  
  cousinsconsang <- df[,reorder]
  
  relist2 <- list()
  for (i in 1:nn) {
    relist2[[i]] <- as.data.frame(cousinsconsang[,c(odd[i],even[i])])
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("cbind.data.frame", consang)
  colnames(consangdf) <- paste0('sample',1:nn)
  
  data <- consangdf
  data <- data.frame(lapply(data, function(x) {gsub("G C", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G G", "0", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("C C", "2", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G T", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T G", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T T", "2", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("G A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A G", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A A", "0", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C G", "1", x)})) 
  
  data[] <- lapply(data, function (i) as.numeric(as.character(i)))
  colnames(geneticmap) <- c('chr','id','posg','pos','ref','alt')
  fam <- cbind.data.frame(1,paste0('sibmate',1:nn),0,0,0,0)
  colnames(fam) <- c('fam','id','pat','mat','sex','pheno')
  write_plink(paste0(output,chr,"_2_gen_avuncular"),bim=geneticmap,as.matrix(data) )
}

### individuals with parents that are avuncular for one generation

consang_avan <- function(map,ped,nn,chr,output){
  
  geneticmap <- read.table(map)
  ped1 <- as.data.frame(read.table(ped,sep=" ",colClasses = "character"), as.factor = F)
  odd <- seq(1,dim(ped1)[2]*2+1,2)
  even <- seq(2,dim(ped1)[2]*2+2,2)
  relist2 <- list()
  num <- dim(ped1)[2]/2
  for (i in 1:num) {
    relist2[[i]] <- ped1[,c(odd[i],even[i])]
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("rbind.data.frame", consang)
  
  small_ped <- consangdf[-c(1:3),1:nn]
  colsplit <- lapply(1:nn, function(i) as.data.frame(small_ped[,i]))
  
  splitgeno <- lapply(1:dim(small_ped)[2], function(i) str_split_fixed(colsplit[[i]][,1], " ", 2))
  
  recombine <- function(i,geno){
    data1 <- geno[1:i,]
    data1 <- cbind.data.frame(as.character(data1[,1]),as.character(data1[,2]))
    colnames(data1) <- c("V1","V2")
    data2 <- geno[(i+1):dim(small_ped)[1],]
    fliped <- cbind.data.frame(as.character(data2[,2]),as.character(data2[,1]))
    colnames(fliped) <- c("V1","V2")
    final <-  rbind.data.frame(data1,fliped)
    return(final)
  }
  
  multirec <- function(geno){
    difference <- diff(geneticmap$V3)*.01
    flips <- rbinom(length(difference),1,difference)
    vectorflips <- match(c(1),flips)
    data <- geno
    data <- cbind.data.frame(as.character(data[,1]),as.character(data[,2]))
    colnames(data) <- c("V1","V2")
    if (!is.na(vectorflips)) {
      for (i in vectorflips) {
        data <- recombine(i,data)
      }
    }
    return(data)
  }
  
  randomchrom <- function(df) {
    return(as.data.frame(df[, sample(2, 1)]))
  }
  
  #### make siblings and unrelated
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m1 <- cbind.data.frame(df_1, df_2)   
  
  df1 <- m1[, c(matrix(1:ncol(m1), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-6,6), times =1 ,  each = 6) + rep(c(1,3,2,4,5,6),nn/3)
  
  firstbreed1 <- df1[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed1[,c(odd[i],even[i])])
  }
  
  ### make cousins
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  for (i in 1:nn) {
    if (i %% 6 == 1 || i %% 6 == 2){
      df[,i] = firstbreed1[,i]
    }
  }
  
  reorder <- rep(seq(0,nn*2-6,6), times =1 ,  each = 6) + rep(c(1,2,3,5,4,6),nn/3)
  
  
  firstbreed <- df[,reorder]
  
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-6,6), times =1 ,  each = 6) + rep(c(1,3,2,4,1,5),nn/3)
  
  cousinsconsang <- df[,reorder]
  
  relist2 <- list()
  for (i in 1:nn) {
    relist2[[i]] <- as.data.frame(cousinsconsang[,c(odd[i],even[i])])
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("cbind.data.frame", consang)
  colnames(consangdf) <- paste0('sample',1:nn)
  
  data <- consangdf
  data <- data.frame(lapply(data, function(x) {gsub("G C", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G G", "0", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("C C", "2", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G T", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T G", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T T", "2", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("G A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A G", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A A", "0", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C G", "1", x)})) 
  
  data[] <- lapply(data, function (i) as.numeric(as.character(i)))
  colnames(geneticmap) <- c('chr','id','posg','pos','ref','alt')
  fam <- cbind.data.frame(1,paste0('sibmate',1:nn),0,0,0,0)
  colnames(fam) <- c('fam','id','pat','mat','sex','pheno')
  write_plink(paste0(output,chr,"_avunc"),bim=geneticmap,as.matrix(data) )
}


### individuals with parents that are first cousins for three generations

consang_cous_3gen <- function(map,ped,nn,chr,output){
  geneticmap <- read.table(map)
  ped1 <- as.data.frame(read.table(ped,sep=" ",colClasses = "character"), as.factor = F)
  odd <- seq(1,dim(ped1)[2]*2+1,2)
  even <- seq(2,dim(ped1)[2]*2+2,2)
  relist2 <- list()
  num <- dim(ped1)[2]/2
  for (i in 1:num) {
    relist2[[i]] <- ped1[,c(odd[i],even[i])]
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("rbind.data.frame", consang)
  
  small_ped <- consangdf[-c(1:3),1:nn]
  
  colsplit <- lapply(1:nn, function(i) as.data.frame(small_ped[,i]))
  
  splitgeno <- lapply(1:dim(small_ped)[2], function(i) str_split_fixed(colsplit[[i]][,1], " ", 2))
  
  recombine <- function(i,geno){
    data1 <- geno[1:i,]
    data1 <- cbind.data.frame(as.character(data1[,1]),as.character(data1[,2]))
    colnames(data1) <- c("V1","V2")
    data2 <- geno[(i+1):dim(small_ped)[1],]
    fliped <- cbind.data.frame(as.character(data2[,2]),as.character(data2[,1]))
    colnames(fliped) <- c("V1","V2")
    final <-  rbind.data.frame(data1,fliped)
    return(final)
  }
  
  multirec <- function(geno){
    difference <- diff(geneticmap$V3)*.01
    flips <- rbinom(length(difference),1,difference)
    vectorflips <- match(c(1),flips)
    data <- geno
    data <- cbind.data.frame(as.character(data[,1]),as.character(data[,2]))
    colnames(data) <- c("V1","V2")
    if (!is.na(vectorflips)) {
      for (i in vectorflips) {
        data <- recombine(i,data)
      }
    }
    return(data)
  }
  
  randomchrom <- function(df) {
    return(as.data.frame(df[, sample(2, 1)]))
  }
  
  #### make siblings and unrelated
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m1 <- cbind.data.frame(df_1, df_2)   
  
  df1 <- m1[, c(matrix(1:ncol(m1), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-16,16), times =1 ,  each = 16) + rep(c(1,2,3,5,4,6,7,9,8,10,11,13,12,14,15,16),nn/8)
  
  firstbreed <- df1[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  ### make cousins
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-16,16), times =1 ,  each = 16) + rep(c(1,3,2,4,5,7,6,8,9,11,10,12,13,15,14,16),nn/8)
  
  firstbreed <- df[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-16,16), times =1 ,  each = 16) + rep(c(1,2,3,5,4,6,7,9,8,10,11,13,12,14,15,16),nn/8)
  
  firstbreed <- df[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-16,16), times =1 ,  each = 16) + rep(c(1,2,3,4,5,7,6,8,9,11,10,12,13,14,15,16),nn/8)
  
  firstbreed <- df[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  
  
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-16,16), times =1 ,  each = 16) + rep(c(5,9,7,11,5,10,7,12,6,9,8,11,6,10,8,12),nn/8)
  
  cousinsconsang <- df[,reorder]
  
  relist2 <- list()
  for (i in 1:nn) {
    relist2[[i]] <- as.data.frame(cousinsconsang[,c(odd[i],even[i])])
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("cbind.data.frame", consang)
  colnames(consangdf) <- paste0('sample',1:nn)
  
  data <- consangdf
  data <- data.frame(lapply(data, function(x) {gsub("G C", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G G", "0", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("C C", "2", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G T", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T G", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T T", "2", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("G A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A G", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A A", "0", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C G", "1", x)})) 
  
  data[] <- lapply(data, function (i) as.numeric(as.character(i)))
  colnames(geneticmap) <- c('chr','id','posg','pos','ref','alt')
  fam <- cbind.data.frame(1,paste0('sibmate',1:nn),0,0,0,0)
  colnames(fam) <- c('fam','id','pat','mat','sex','pheno')
  write_plink(paste0(output,chr,"_3_gen_cousins"),bim=geneticmap,as.matrix(data) )
}

### individuals with parents that are first cousins for 2 generations

consang_cous_2gen <- function(map,ped,nn,chr,output){
  geneticmap <- read.table(map)
  ped1 <- as.data.frame(read.table(ped,sep=" ",colClasses = "character"), as.factor = F)
  odd <- seq(1,dim(ped1)[2]*2+1,2)
  even <- seq(2,dim(ped1)[2]*2+2,2)
  relist2 <- list()
  num <- dim(ped1)[2]/2
  for (i in 1:num) {
    relist2[[i]] <- ped1[,c(odd[i],even[i])]
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("rbind.data.frame", consang)
  
  small_ped <- consangdf[-c(1:3),1:nn]
  
  colsplit <- lapply(1:nn, function(i) as.data.frame(small_ped[,i]))
  
  splitgeno <- lapply(1:dim(small_ped)[2], function(i) str_split_fixed(colsplit[[i]][,1], " ", 2))
  
  recombine <- function(i,geno){
    data1 <- geno[1:i,]
    data1 <- cbind.data.frame(as.character(data1[,1]),as.character(data1[,2]))
    colnames(data1) <- c("V1","V2")
    data2 <- geno[(i+1):dim(small_ped)[1],]
    fliped <- cbind.data.frame(as.character(data2[,2]),as.character(data2[,1]))
    colnames(fliped) <- c("V1","V2")
    final <-  rbind.data.frame(data1,fliped)
    return(final)
  }
  
  multirec <- function(geno){
    difference <- diff(geneticmap$V3)*.01
    flips <- rbinom(length(difference),1,difference)
    vectorflips <- match(c(1),flips)
    data <- geno
    data <- cbind.data.frame(as.character(data[,1]),as.character(data[,2]))
    colnames(data) <- c("V1","V2")
    if (!is.na(vectorflips)) {
      for (i in vectorflips) {
        data <- recombine(i,data)
      }
    }
    return(data)
  }
  
  randomchrom <- function(df) {
    return(as.data.frame(df[, sample(2, 1)]))
  }
  
  #### make siblings and unrelated
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m1 <- cbind.data.frame(df_1, df_2)   
  
  df1 <- m1[, c(matrix(1:ncol(m1), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-12,12), times =1 ,  each = 12) + rep(c(1,2,3,5,4,6,7,9,8,10,11,12),nn/6)
  
  firstbreed <- df1[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  ### make cousins
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-12,12), times =1 ,  each = 12) + rep(c(1,3,2,4,5,7,6,8,9,11,10,12),nn/6)
  
  firstbreed <- df[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-12,12), times =1 ,  each = 12) + rep(c(1,2,3,5,4,6,7,9,8,10,11,12),nn/6)
  
  firstbreed <- df[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-12,12), times =1 ,  each = 12) + rep(c(5,9,3,9,5,7,6,8,4,10,3,7),nn/6)
  
  cousinsconsang <- df[,reorder]
  
  relist2 <- list()
  for (i in 1:nn) {
    relist2[[i]] <- as.data.frame(cousinsconsang[,c(odd[i],even[i])])
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("cbind.data.frame", consang)
  colnames(consangdf) <- paste0('sample',1:nn)
  
  data <- consangdf
  data <- data.frame(lapply(data, function(x) {gsub("G C", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G G", "0", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("C C", "2", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G T", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T G", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T T", "2", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("G A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A G", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A A", "0", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C G", "1", x)})) 
  
  data[] <- lapply(data, function (i) as.numeric(as.character(i)))
  colnames(geneticmap) <- c('chr','id','posg','pos','ref','alt')
  fam <- cbind.data.frame(1,paste0('sibmate',1:nn),0,0,0,0)
  colnames(fam) <- c('fam','id','pat','mat','sex','pheno')
  write_plink(paste0(output,chr,"_2_gen_cousins"),bim=geneticmap,as.matrix(data) )
}

### individuals with parents that are first cousins for one generation

consang_cous <- function(map,ped,nn,chr,output){
  
  geneticmap <- fread(map)
  ped1 <- as.data.frame(fread(ped,sep=" ",colClasses = "character"), as.factor = F)
  odd <- seq(1,dim(ped1)[2]*2+1,2)
  even <- seq(2,dim(ped1)[2]*2+2,2)
  relist2 <- list()
  num <- dim(ped1)[2]/2
  for (i in 1:num) {
    relist2[[i]] <- ped1[,c(odd[i],even[i])]
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("rbind.data.frame", consang)
  
  small_ped <- consangdf[-c(1:3),1:nn]
  colsplit <- lapply(1:nn, function(i) as.data.frame(small_ped[,i]))
  
  splitgeno <- lapply(1:dim(small_ped)[2], function(i) str_split_fixed(colsplit[[i]][,1], " ", 2))
  
  recombine <- function(i,geno){
    data1 <- geno[1:i,]
    data1 <- cbind.data.frame(as.character(data1[,1]),as.character(data1[,2]))
    colnames(data1) <- c("V1","V2")
    data2 <- geno[(i+1):dim(small_ped)[1],]
    fliped <- cbind.data.frame(as.character(data2[,2]),as.character(data2[,1]))
    colnames(fliped) <- c("V1","V2")
    final <-  rbind.data.frame(data1,fliped)
    return(final)
  }
  
  multirec <- function(geno){
    difference <- diff(geneticmap$V3)*.01
    flips <- rbinom(length(difference),1,difference)
    vectorflips <- match(c(1),flips)
    data <- geno
    data <- cbind.data.frame(as.character(data[,1]),as.character(data[,2]))
    colnames(data) <- c("V1","V2")
    if (!is.na(vectorflips)) {
      for (i in vectorflips) {
        data <- recombine(i,data)
      }
    }
    return(data)
  }
  
  randomchrom <- function(df) {
    return(as.data.frame(df[, sample(2, 1)]))
  }
  
  #### make siblings and unrelated
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m1 <- cbind.data.frame(df_1, df_2)   
  
  df1 <- m1[, c(matrix(1:ncol(m1), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-8,8), times =1 ,  each = 8) + rep(c(1,3,2,4,5,6,7,8),nn/4)
  
  firstbreed <- df1[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  ### make cousins
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-8,8), times =1 ,  each = 8) + rep(c(1,5,2,6,3,7,4,8),nn/4)
  
  firstbreed <- df[,reorder]
  
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  reorder <- rep(seq(0,nn*2-8,8), times =1 ,  each = 8) + rep(c(1,5,2,6,3,7,4,8),nn/4)
  
  cousinsconsang <- df[,reorder]
  
  relist2 <- list()
  for (i in 1:nn) {
    relist2[[i]] <- as.data.frame(cousinsconsang[,c(odd[i],even[i])])
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("cbind.data.frame", consang)
  colnames(consangdf) <- paste0('sample',1:nn)
  
  data <- consangdf
  data <- data.frame(lapply(data, function(x) {gsub("G C", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G G", "0", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("C C", "2", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G T", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T G", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T T", "2", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("G A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A G", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A A", "0", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C G", "1", x)})) 
  
  data[] <- lapply(data, function (i) as.numeric(as.character(i)))
  colnames(geneticmap) <- c('chr','id','posg','pos','ref','alt')
  fam <- cbind.data.frame(1,paste0('sibmate',1:nn),0,0,0,0)
  colnames(fam) <- c('fam','id','pat','mat','sex','pheno')
  write_plink(paste0(output,chr,"_1_gen_cousins"),bim=geneticmap,as.matrix(data) )
}

### individuals with parents that are siblings 

consang_sib <- function(map,ped,nn,chr,output){
  
  geneticmap <- read.table(map)
  ped1 <- as.data.frame(read.table(ped,sep=" ",colClasses = "character"), as.factor = F)
  odd <- seq(1,dim(ped1)[2]*2+1,2)
  even <- seq(2,dim(ped1)[2]*2+2,2)
  relist2 <- list()
  num <- dim(ped1)[2]/2
  for (i in 1:num) {
    relist2[[i]] <- ped1[,c(odd[i],even[i])]
  }
  
  consang <- lapply(1:length(relist2), function(i) apply(relist2[[i]], 1, paste, collapse=" "))
  consangdf <- do.call("rbind.data.frame", consang)
  
  small_ped <- consangdf[-c(1:3),1:nn]
  
  ## turn sample dataframe into list of individual samples, columns look like e.g. "A G"
  colsplit <- lapply(1:nn, function(i) as.data.frame(small_ped[,i]))
  
  ## split each individual into a dataframe with two columns, one for each haplotype
  splitgeno <- lapply(1:dim(small_ped)[2], function(i) str_split_fixed(colsplit[[i]][,1], " ", 2))
  
  ## recombination function, given an individual geno and position i the columns will be flipped after position i
  recombine <- function(i,geno){
    data1 <- geno[1:i,]
    data1 <- cbind.data.frame(as.character(data1[,1]),as.character(data1[,2]))
    colnames(data1) <- c("V1","V2")
    data2 <- geno[(i+1):dim(small_ped)[1],]
    fliped <- cbind.data.frame(as.character(data2[,2]),as.character(data2[,1]))
    colnames(fliped) <- c("V1","V2")
    final <-  rbind.data.frame(data1,fliped)
    return(final)
  }
  
  ## modelling multiple recombinations, calculate probability of flip between genetic positions using genetic map, model recombination loci 
  ##using binomial distribution, create vector with positions i to flip and then use recombination function to recombine at a given position
  multirec <- function(geno){
    difference <- diff(geneticmap$V3)*.01
    flips <- rbinom(length(difference),1,difference)
    vectorflips <- match(c(1),flips)
    data <- geno
    data <- cbind.data.frame(as.character(data[,1]),as.character(data[,2]))
    colnames(data) <- c("V1","V2")
    if (!is.na(vectorflips)) {
      for (i in vectorflips) {
        data <- recombine(i,data)
      }
    }
    return(data)
  }
  
  ## randomly select a single gamete from a recombination event
  randomchrom <- function(df) {
    return(as.data.frame(df[, sample(2, 1)]))
  }
  
  ## recombine individuals and randomly choose a gamete 
  
  recomb1_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb1_1 <- lapply(recomb1_2, function(i) randomchrom(i))
  
  ## do it again
  recomb2_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(splitgeno[[i]]))
  
  recomb2_1 <- lapply(recomb2_2, function(i) randomchrom(i))
  
  df_1 <- do.call("cbind.data.frame", recomb1_1)
  df_2 <- do.call("cbind.data.frame", recomb2_1)
  
  ## combine into a dataframe where the order is ind1_gamete1 ind1_gamete2 ind2_gamete1....
  m <- cbind.data.frame(df_1, df_2)   
  
  df <- m[, c(matrix(1:ncol(m), nrow = 2, byrow = T))] 
  
  ## make sibilings by reordering gametes systematically
  reorder <- rep(seq(0,nn*2-4,4), times =1 ,  each = 4) + rep(c(1,3,2,4),nn/2)
  
  firstbreed <- df[,reorder]
  
  ## make it into list of individuals again
  relist <- list()
  odd <- seq(1,nn*2-1,2)
  even <- seq(2,nn*2,2)
  for (i in 1:nn) {
    relist[[i]] <- as.data.frame(firstbreed[,c(odd[i],even[i])])
  }
  
  ###make two gamaetes per individual again
  recomb3_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb3_1 <- lapply(recomb3_2, function(i) randomchrom(i))
  
  recomb4_2 <- lapply(1:dim(small_ped)[2], function(i) multirec(relist[[i]]))
  
  recomb4_1 <- lapply(recomb4_2, function(i) randomchrom(i))
  
  df_3 <- do.call("cbind.data.frame", recomb3_1)
  df_4 <- do.call("cbind.data.frame", recomb4_1)
  
  m1 <- cbind.data.frame(df_3, df_4)   
  
  df1 <- m1[, c(matrix(1:ncol(m1), nrow = 2, byrow = T))] 
  
  ## mate siblings
  reorder <- rep(seq(0,nn*2-4,4), times =1 ,  each = 4) + rep(c(1,3,2,4),nn/2)
  
  consangsib <- df1[,reorder]
  
  relist1 <- list()
  for (i in 1:nn) {
    relist1[[i]] <- as.data.frame(consangsib[,c(odd[i],even[i])])
  }
  
  ### combine genotypes into single column and make dataframe
  breed2 <- lapply(1:length(relist1), function(i) apply(relist1[[i]], 1, paste, collapse=" "))
  df11 <- do.call("cbind.data.frame", breed2)
  colnames(df11) <- paste0('sample',1:nn)
  
  ### recode as homosygous 0,2 or hetero 1
  data <- df11
  data <- data.frame(lapply(data, function(x) {gsub("G C", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G G", "0", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("C C", "2", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("G T", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T G", "1", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("T T", "2", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("T C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("G A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A C", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A G", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("A A", "0", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C A", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C T", "1", x)})) 
  data <- data.frame(lapply(data, function(x) {gsub("C G", "1", x)})) 
  
  ### write out as plink file of simualted individuals
  data[] <- lapply(data, function (i) as.numeric(as.character(i)))
  colnames(geneticmap) <- c('chr','id','posg','pos','ref','alt')
  fam <- cbind.data.frame(1,paste0('sibmate',1:nn),0,0,0,0)
  colnames(fam) <- c('fam','id','pat','mat','sex','pheno')
  write_plink(paste0(output,chr,"_sib"),bim=geneticmap,as.matrix(data) )
}
