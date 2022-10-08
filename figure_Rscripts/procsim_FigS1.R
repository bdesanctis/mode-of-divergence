filedir <- 'complete_sims/simsFigS1/'


maxd <- 500

lf <- list.files(filedir, pattern='res',recursive=T)
# lf <- lf[grep(lf, pattern='scenario[1-2]/', invert=T)]
# lf <- lf[grep(lf, pattern='variable_dominance/')]
# #lf <- lf[grep(lf, pattern='E2', invert=T)]
# #lf <- lf[grep(lf, pattern='E1', invert=T)]
# lf <- lf[grep(lf, pattern='E0')]
# lf <- lf[grep(lf, pattern='k2')]
lf1 <- lf[grep(lf,pattern='rep1.txt')]

getparams <- function(fn){
  p <- unlist(strsplit(unlist(strsplit(fn,'res'))[2],"_"))[2:11]
  params <- sapply(p, FUN=substring, 2)
  names(params) <- sapply(p, FUN=substr, start=1,stop=1)
  return(params)
}

process_sim <- function(f1, f2, maxd){

  fc <- rbind(cbind(f1,P=rep(1,dim(f1)[1])),
              cbind(f2,P=rep(2,dim(f2)[1])))
  fc <- fc[order(fc$fix_time),]
  fc <- fc[1:maxd,]
  
  mcols <- grep('mut_effect',colnames(fc))
  dcols <- grep('dom',colnames(fc))
  
  #PROCESS SIMS
  #1. scale by sqrt(1/2) to deal with annoying alpha
  fc[,mcols] <- fc[,mcols]*sqrt(1/2)
  
  #Compute d_ij
  fc[which(fc$P==1),dcols] <- fc[which(fc$P==1),mcols]*(fc[which(fc$P==1),dcols] - 1/2)
  fc[which(fc$P==2),dcols] <- fc[which(fc$P==2),mcols]*(fc[which(fc$P==2),dcols] - 1/2)
  
  #Compute a_ij
  fc[which(fc$P==2),mcols] <- fc[which(fc$P==2),mcols]/2
  #For subs that have fixed in P1 (if both populations evolve) we need to flip them, to go from P1 to P2.
  fc[which(fc$P==1),mcols] <- -fc[which(fc$P==1),mcols]/2
  
  return(fc)
}

getstats <- function(fc,mcols,dcols,D=25){
  n <- length(mcols)
  
  Ma <- cumsum(rowSums(fc[,mcols]^2)) #sum across traits, then subs
  ma <- rowSums(apply(fc[,mcols],2,cumsum)^2) #sum across subs, then traits
  fa <- (ma-Ma)

  mat <- cbind(Ma=Ma,ma=ma,fa=fa)

  return(mat[D,])
}


res <- matrix(NA, ncol=14,nrow=0)

for(i in 1:length(lf1)){
  if(i %% 10 == 0) cat(i,'/',length(lf1),'\n')
  
  fn1 <- lf1[i]
  fn2 <- gsub(fn1, pattern='rep1', replacement = 'rep2')

  params <- getparams(fn1)
  
  f1 <- read.table(paste0(filedir,fn1),header=T)
  f2 <- read.table(paste0(filedir,fn2), header=T)

  fc <- process_sim(f1,f2, maxd=maxd)
  mcols <- grep('mut_effect',colnames(fc))
  dcols <- grep('dom',colnames(fc))
  
  stats <- getstats(fc, mcols,dcols,maxd)
  res <- rbind(res, c(params,stats,maxd))
}
colnames(res) <- c(names(params), names(stats), 'd')

df <- as.data.frame(res, stringsAsFactors = T)
for(i in 11:14){
  df[,i] <- as.numeric(res[,i])
}


save(df, file='procSim_FigS1.RData')


