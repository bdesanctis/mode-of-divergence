library(Hmisc)
library(RColorBrewer)
library(scales)

fs <- c(-1,-999) # uniform distribution of dominance coefficients, either per-trait or per-mutation
n <- 20
muts <- c('perTrait','perMut')
overD <- TRUE
D <- 1e4
M <- 1
scenarios <- 'drift'

#which stats to plot
paths <- c('M','m','f')
stats <- c('a','d','ad')

load('procSim_FigS3.RData')
params <- as.data.frame(t(sapply(a, '[[','params')))

#set up colour palette
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

dplotall <- function(a, cols,sc=0,shaded=F,lm=F, sub=F,D,paths,stats, diff=T){
  for(i in 1:length(paths)){
    for(j in 1:length(stats)){

      for(w in 1:length(a)){
        mycol <- cols[1+(a[[w]]$params$P1ancestral==TRUE),2]
        dat <- sapply(a[[w]]$res, FUN=getstat, term=stats[j], path=paths[i], diff=diff,sub=sub)
        means <- rowMeans(dat)
        if (stats[j]=='ad') means <- -means #needs to be flipped if P2 ancestral instead of P1.
        sds <- apply(dat, 1, sd)
        
        ymin <- min(c(0,means-sds)); ymax <- max(means+sds)
        if(stats[j]=='ad') {
          ymin <- min(means-sds/4)
          ymax <- max(means+sds/4)
        } else if(paths[i]=='f' & stats[j]=='d'){
          ymin <- min(c(0,means-sds/4))
          ymax <- max(means+sds/4)#1.4
        }
        
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
                 fill=c(cols[c(1,2),2]),
                 c(expression('Similar N'['e']),expression('P1 much lower N'['e'])),
                 bty='n',cex=1,border=F, seg.len = 1)
            legend('bottomright', 
                   lwd=c(1,2),
                   c('per-trait dominance','per-mutation dominance'),
                   cex=1, bty='n')
          }
        }
        
        if(sub){
          Dvals <- c(1:D)[which((c(1:D) %% 20) == 1 )]
        }else{
          Dvals <- 1:D
        }
        if(shaded) polygon(c(Dvals,rev(Dvals)),c(means+sds,rev(means-sds)),col=alpha(cols[1+(a[[w]]$params$P1ancestral==TRUE),1],.1),border=F)
        lines(Dvals,means, col=mycol,
              lwd=c(1,2)[1+(a[[w]]$params$f==-999)])
      }
    }
  }
}

#set up plotting device
quartz(width=10,height=8)
layout(matrix(1:9,ncol=3))

#make plot
w <- which(params$n==n & 
             params$sc %in% scenarios & 
             params$mutmodel %in% muts & 
             params$M==M &
             params$overD == overD &
             params$f %in% fs 
             )

dplotall(a=a[rev(w)], cols=cols[c(1,3),],shaded=T,sub=T,D=D,paths=paths,stats=stats)

dev.copy2pdf(file='FigS3.pdf',
             width=10, height=8)
