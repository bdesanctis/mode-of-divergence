library(Hmisc)
library(RColorBrewer)
library(scales)

f <- 123 #this flag indicates that dominance coefficients are drawn from a 'realistic' distribution where large effect mutations have more extreme coefficients
mut <- 'perTrait'
overD <- TRUE
M <- 1
D <- 1000
n <- 20
scenarios <- 'drift'

paths <- c('M','m','f')
stats <- c('a','d','ad')

load("procSim.RData")
params <- as.data.frame(t(sapply(a, '[[','params')))

#initiate colour palette
palpaired <- brewer.pal(12,name='Paired')
palpaired[11] <- brewer.pal(11,name="BrBG")[4]
cols <- matrix(palpaired,ncol=2,byrow=T)


getstat <- function(x, diff=F, term, path,sub=T,maxd=D){
  if(sub){
    lines <- c(1:maxd)[which((c(1:maxd) %% 20) == 1 )]
  }else{
    lines <- c(1:maxd)
  }
                  
  if(path=='f'){
    if(diff){
      return(x[lines,paste0('m',term)]-x[lines,paste0('M',term)])
    }else{
      return(x[lines,paste0('M',term)]/x[lines,paste0('m',term)])
    }
  }else{
    return(x[lines,paste0(path,term)])
  }
  
}
getylab <- function(path, stat,diff=F){
  if(path=='m'){
    if(stat=='a') return(expression(paste('(B) ',italic(m),'(',bold(A),',',bold(A),')')))
    if(stat=='d') return(expression(paste('(E) ',italic(m),'(',bold(Delta),',',bold(Delta),')')))
    if(stat =='ad') return(expression(paste('(H) ',italic(m),'(',bold(A),',',bold(Delta),')')))
  }else if(path=='M'){
    if(stat=='a') return(expression(paste('(A) ',italic(M),'(',bold(A),',',bold(A),')')))
    if(stat=='d') return(expression(paste('(D) ',italic(M),'(',bold(Delta),',',bold(Delta),')')))
    if(stat =='ad') return(expression(paste('(G) ',italic(M),'(',bold(A),',',bold(Delta),')')))
  }else if(path=='f' & diff==T){
    if(stat=='a') return(expression(paste('(C) ',italic(m),'(',bold(A),',',bold(A),') - ',italic(M),'(',bold(A),',',bold(A),')')))
    if(stat=='d') return(expression(paste('(F) ',italic(m),'(',bold(Delta),',',bold(Delta),') - ',italic(M),'(',bold(Delta),',',bold(Delta),')')))
    if(stat =='ad') return(expression(paste('(I) ',italic(m),'(',bold(A),',',bold(Delta),') - ',italic(M),'(',bold(A),',',bold(Delta),')')))
  }else if(path=='f' & diff==F){
    if(stat=='a') return(expression(paste('(C) ',italic(M),'(',bold(A),',',bold(A),') / ',italic(m),'(',bold(A),',',bold(A),')')))
    if(stat=='d') return(expression(paste('(F) ',italic(M),'(',bold(Delta),',',bold(Delta),') / ',italic(m),'(',bold(Delta),',',bold(Delta),')')))
    if(stat =='ad') return(expression(paste('(I) ',italic(M),'(',bold(A),',',bold(Delta),') / ',italic(m),'(',bold(A),',',bold(Delta),')')))
  }
}

dplotall <- function(a, cols,sc=0,shaded=F,lm=F, sub=F,D,n, paths,stats){
  if(sub){
    lines <- c(1:D)[which((c(1:D) %% 20) == 1 )]
  }else{
    lines <- c(1:D)
  }
  
  for(i in 1:length(paths)){
    for(j in 1:length(stats)){
      
      #get y limits
      dat <- lapply(a, FUN=function(x) {sapply(x$res, FUN=getstat, term=stats[j], path=paths[i], diff=diff,sub=sub,maxd=D)})
      means <- sapply(dat, FUN=rowMeans)
      if (stats[j]=='ad') means <- -means #needs to be flipped if P2 ancestral instead of P1.
      sds <- sapply(dat, FUN=function(x) {apply(x, 1, sd)*2})
      ymin <- min(means-sds); ymax <- max(means+sds)
      
      
      for(w in 1:length(a)){
        mycol <- cols[1+(a[[w]]$params$P1ancestral==TRUE),2]
        dat <- sapply(a[[w]]$res, FUN=getstat, term=stats[j], path=paths[i], diff=diff,sub=sub,maxd=D)
        means <- rowMeans(dat)
        if (stats[j]=='ad') means <- -means
        sds <- apply(dat, 1, sd)*2
        
        if(w==1){
          plot(0,type='n',
             ylim=c(ymin,ymax),
             xlim=c(1,D),
             xlab='Divergence (D)',
             ylab='',
             main='')
          mtext(side=3,adj=0, getylab(paths[i],stats[j],diff=diff),cex=1,padj=-.5)
          if(stats[j]=='ad' | paths[i]=='f' & stats[j]=='d') abline(h=0, lty=3, lwd=2)
          if(i==1 & j==1) {
            legend('topleft',
                 fill=cols[1:2,2],
                 c(expression('Similar N'['e']),expression('P1 much lower N'['e'])),
                 bty='n',cex=1.2,border=F)
          }
        }
        
        if(shaded) polygon(c(lines,rev(lines)),c(means+sds,rev(means-sds)),col=alpha(mycol,.1),border=F)
        lines(lines,means, col=mycol,
              lty=1,
              lwd=2)
      }
    }
  }
}

quartz(width=10,height=8)
layout(matrix(1:9,ncol=3))

w <- which(params$n==n & 
             params$sc %in% scenarios & 
             params$mutmodel==mut & 
             params$M==M &
             params$overD == overD &
             params$f == f)
dplotall(a=a[rev(w)], cols=cols[c(1,3),],shaded=T,sub=T,D=D,paths=paths,stats=stats)

dev.copy2pdf(file='FigS2.pdf',
             width=10, height=8)
