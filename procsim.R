getstats <- function(fc,mcols,dcols,overD,D=25){
  n <- length(mcols)

  Ma <- cumsum(rowSums(fc[,mcols]^2)) #sum across traits, then subs
  Md <- cumsum(rowSums(fc[,dcols]^2))
  Mad <- cumsum(rowSums(fc[,mcols]*fc[,dcols]))

  ma <- rowSums(apply(fc[,mcols],2,cumsum)^2) #sum across subs, then traits
  md <- rowSums(apply(fc[,dcols],2,cumsum)^2)
  mad <- rowSums(apply(fc[,mcols],2,cumsum)*apply(fc[,dcols],2,cumsum))

  fa <- (ma-Ma)
  fd <- (md-Md)
  fad <- (mad-Mad)

  mat <- cbind(Ma=Ma,ma=ma,fa=fa,
        Md=Md,md=md,fd=fd,
        Mad=Mad,mad=mad,fad=fad)
  if(overD){
    return(mat)
  }else{
    return(mat[D,])
  }
}

process_sim <- function(rep1, rep2, dir, params){
  maxd <- params$D
  anc <- params$P1ancestral
  scenario <- params$sc
  
  fstem <- paste0('res_N',params$N,'_L',params$L,'_U',params$U,'_k',params$k,
                  '_s',formatC(params$s),'_n', params$n,'_F',params$f,
                  '_q',params$q,'_M',params$M,'_E',params$E,'_rep')
  f1 <- read.table(paste0(dir,'/',fstem,rep1,'.txt'),header=T)
  f2 <- read.table(paste0(dir,'/',fstem,rep2,'.txt'),header=T)

  if(anc){
    if(scenario=='drift'){
      fc <- cbind(f1[1:(maxd),],P=rep(1,maxd))
    }else{
      fc <- rbind(cbind(f1[1:(maxd/2),],P=rep(1,maxd/2)),
                  cbind(f2[1:(maxd/2),],P=rep(2,maxd/2)))
     }

  }else{
    fc <- rbind(cbind(f1,P=rep(1,dim(f1)[1])),
                cbind(f2,P=rep(2,dim(f2)[1])))
    fc <- fc[order(fc$fix_time),]
    fc <- fc[1:maxd,]
  }

  mcols <- grep('mut_effect',colnames(fc))
  dcols <- grep('dom',colnames(fc))

  #PROCESS SIMS
  #1. scale by sqrt(1/2) to deal with annoying alpha
  fc[,mcols] <- fc[,mcols]*sqrt(1/2)

  #2. Deal with scenarios
  if(scenario == 'return' & anc){
      fc[which(fc$P==1),mcols] <- -fc[which(fc$P==1),mcols]
  }else if(scenario == 'orth'){
      fc[which(fc$P==1),mcols[c(1,2)]] <- fc[which(fc$P==1),mcols[c(2,1)]]
      fc[which(fc$P==1),dcols[c(1,2)]] <- fc[which(fc$P==1),dcols[c(2,1)]]
  }else if(scenario == 'straight' & !anc){
      fc[which(fc$P==1),mcols] <- -fc[which(fc$P==1),mcols]
  }

  if(anc){
    fc$P <- 2
  }

  #Compute d_ij
  fc[which(fc$P==1),dcols] <- fc[which(fc$P==1),mcols]*(fc[which(fc$P==1),dcols] - 1/2)
  fc[which(fc$P==2),dcols] <- fc[which(fc$P==2),mcols]*(fc[which(fc$P==2),dcols] - 1/2)

  #Compute a_ij
  fc[which(fc$P==2),mcols] <- fc[which(fc$P==2),mcols]/2
  #For subs that have fixed in P1 (if both populations evolve) we need to flip them, to go from P1 to P2.
  fc[which(fc$P==1),mcols] <- -fc[which(fc$P==1),mcols]/2

  return(fc)
}


#PARAMETER definitions
# L=-1: free recombination among all loci
# k=2:  quadratic decline in fitness with distance from optimum
# q=0.5 dominance coefficients neither phenotypically recessive nor dominant ON AVERAGE
# M=1: homozygous mutation effects follow multivariate normal distribution
# f:  variance in dominance coefficients of new mutations. 
# f=-1 indicates uniform distribution. f=123 indicates 'realistic' distribution. f=-999 indicates per-mutation dominance coefficient
# U=0.01 mutation rate per haploid genome
# s=0.01 average selection coefficient of new mutations in optimal background
# n number of traits under selection / dimensions
# D: total number of substitutions
# overD: record stats at the end, or throughout
# anc: if true, all subs have occurred in P2, and P1 stays in the ancestral state
# mutmodel: independent dominance coefficients per trait, or a single dominance coefficient for all traits

filedir <- 'complete_sims'

parammat <- read.csv(file='sims_param_combinations.csv')
parammat <- parammat[13:16,]
#COMPILE SIMULATIONS
a <- list()
for(i in 1:nrow(parammat)){

    params <- parammat[i,]
    pairing <- sample.int(200,200) #randomize pairing
    np <- 100
    
    cat(i,'/',nrow(parammat),' ',as.character(params),'\n')

    if(params$overD){
      res <- list()
    }else{
      res <- as.data.frame(matrix(NA,ncol=9,nrow=np))
    }

    for(j in 1:np)
    {
      if((j %% 10) == 1) cat('|')
      fc <- process_sim(rep1=pairing[j], rep2=pairing[j+np], dir=filedir, params=params)
      mcols <- grep('mut_effect',colnames(fc))
      dcols <- grep('dom',colnames(fc))

      if(params$overD){
        stats <- getstats(fc,mcols,dcols,overD=T)
        res[[length(res)+1]] <- stats
      }else{
        stats <- getstats(fc,mcols,dcols,overD=F,D=params$D)
        res[j,] <- stats
        colnames(res) <- names(stats)
      }
      rm(fc)
    }
    a[[length(a)+1]] <- list(params=params,res=res)
    cat('\n')
  }

save(a, file="procSim_Fig4.RData")
