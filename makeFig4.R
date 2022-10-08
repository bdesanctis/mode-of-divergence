library(scales)
library(RColorBrewer)
library(Hmisc)

n <- 20
overD <- TRUE
D <- 1000
M <- 1
load('procSim.RData')
params <- as.data.frame(t(sapply(a,'[[','params')))
cex <- 1

terms <- c('a','d','ad') #which term to compute stats for, additive dominance or additive by dominance
diff <- T #compute difference or ratio
shaded=FALSE # include 1 sd around the mean

#set up colour palette
palpaired <- brewer.pal(12,name='Paired')
palpaired[11] <- brewer.pal(11,name="BrBG")[4]
cols <- matrix(palpaired,ncol=2,byrow=T)

getstat <- function(x, diff=F, term){
  if(diff){
    return(x[,paste0('m',term)]-x[,paste0('M',term)])
  }else{
    return(x[,paste0('M',term)]/x[,paste0('m',term)])
  }
}

getylab <- function( diff=F, term){
  if(diff){
    if(term=='a') return(expression(paste(italic(m),'(',bold(A),',',bold(A),') - ',italic(M),'(',bold(A),',',bold(A),')')))
    if(term=='d') return(expression(paste(italic(m),'(',bold(Delta),',',bold(Delta),') - ',italic(M),'(',bold(Delta),',',bold(Delta),')')))
    if(term =='ad') return(expression(paste(italic(m),'(',bold(A),',',bold(Delta),') - ',italic(M),'(',bold(A),',',bold(Delta),')')))
  }else{
    if(term=='a') return(expression(paste(italic(M),'(',bold(A),',',bold(A),') / ',italic(m),'(',bold(A),',',bold(A),')')))
    if(term=='d') return(expression(paste(italic(M),'(',bold(Delta),',',bold(Delta),') / ',italic(m),'(',bold(Delta),',',bold(Delta),')')))
    if(term =='ad') return(expression(paste(italic(M),'(',bold(A),',',bold(Delta),') / ',italic(m),'(',bold(A),',',bold(Delta),')')))
  }
}

dplot <- function(a,scenario,n=20, terms=c('a','d','ad'), diff=T, anc=T,D,overD,invert=F,
                  log=F,shaded=T,colid=c(1,3,5), M=1, f=-1,mutmodel='perTrait',dplot=T,lty=1,ann=F,...){
  
  params <- as.data.frame(t(sapply(a, '[[','params')))
  w <- which(params$n==n & 
               params$sc == scenario & 
               params$P1ancestral==anc & 
               params$mutmodel==mutmodel & 
               params$M==M &
               params$overD == overD &
               params$f == f)
  print(w)
  xmax <- ifelse(log, log10(D),D)
  
  if(dplot){
    plot(0,type='n',
         ylim=c(-.5,.2),
         xlim=c(0,xmax),
         xlab='Divergence (D)',
         ylab='',
         main='',...)
    abline(h=0, col='light gray')
    if(log) axis(side=1,at=log10(c(1,2,5,10,20,50,100)),labels=c(1,2,5,10,20,50,100))
  }
  
  
  for(term in terms){
    if(term == 'd') lty <- 2
    if(term == 'ad') lty <- 3
    dat <- sapply(a[[w]]$res, FUN=getstat, term=term, diff=diff)
    means <- rowMeans(dat)
    if(invert & term=='ad') means <- -means
    sds <- apply(dat, 1, sd)
    
    mycol <- cols[colid[which(terms==term)],]
    
    if(log) {
      if(shaded) polygon(log10(c(1:D,rev(1:D))),c(means+sds,rev(means-sds)),col=alpha(mycol[2],.08),border=F)
      lines(log10(c(1:D)),means, col=mycol[2],lty=lty,lwd=2)
      #if(P1inverse & term =='ad')  lines(log10(c(1:D)),-means, col=mycol[2],lty=lty,lwd=2)
    }else{
      if(shaded) polygon(c(1:D,rev(1:D)),c(means+sds,rev(means-sds)),col=alpha(mycol[2],.08),border=F)
      lines(1:D,means, col=mycol[2],lty=lty,lwd=2)
      #if(P1inverse & term =='ad')  lines(1:D,-means, col=mycol[2],lty=lty,lwd=2)
    }
    
    if(dplot & ann){
      if(term=='a') {
        segments(x0=930,x1=1010, y0=-.2,y1=-.2, col=mycol[2], lty=lty)
        text(900,-.2,getylab(diff=diff,term=term), col=mycol[2],adj=1)
      }else if(term=='d'){
        segments(x0=930,x1=1010, y0=-.3,y1=-.3, col=mycol[2], lty=lty)
        text(900,-.3,getylab(diff=diff,term=term), col=mycol[2],adj=1)
      }else if(term=='ad'){
        segments(x0=930,x1=1010, y0=-.4,y1=-.4, col=mycol[2], lty=lty)
        text(900,-.4,getylab(diff=diff,term=term), col=mycol[2],adj=1)
      }
    }
  }
}


#set up plotting device
graphics.off()
quartz(width=7,height=6)
layout(matrix(1:4,byrow=T,ncol=2))

par(mai=c(0,1,1,0))
dplot(a, scenario='drift',D=D,M=M,f=0.1667,n=n,anc=F,overD=overD,diff=diff, xaxt='n', yaxt='n')
axis(side=2,at=c(-.4,-.2,0,.2))
mtext(side=3, '(A)',adj=.01,cex=.8, padj=2)

mtext(side=3, expression('Similar N'['e']),adj=.5,cex=cex, padj=-1.5)
mtext(side=2, 'No mutational bias',adj=0.5,cex=cex, padj=-4)

par(mai=c(0,0,1,1))
dplot(a, scenario='drift',D=D,M=M,f=0.1667,n=n,anc=T,overD=overD,diff=diff, lty=1, xaxt='n', yaxt='n', invert=T)
mtext(side=3, '(B)',adj=.01,cex=.8, padj=2)

mtext(side=3, expression('P1 much lower N'['e']),adj=.5,cex=cex, padj=-1.5)


par(mai=c(1,1,0,0))
dplot(a, scenario='drift',D=D,M=M,f=123,n=n,anc=F,overD=overD,diff=diff, yaxt='n')
axis(side=2,at=c(-.4,-.2,0,.2))
mtext(side=3, '(C)',adj=.01,cex=.8, padj=2)

mtext(side=2, 'Strong bias for',adj=.5,cex=cex,padj=-5)
mtext(side=2, 'recessive mutations',adj=.5,cex=cex,padj=-3.5)

par(mai=c(1,0,0,1))
dplot(a, scenario='drift',D=D,M=M,f=123,n=n,anc=T,overD=overD,diff=diff, lty=1, yaxt='n',invert=T, ann=T)
mtext(side=3, '(D)',adj=.01,cex=.8, padj=2)

dev.copy2pdf(file='Fig4.pdf')
