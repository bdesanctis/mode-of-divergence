library(scales)
library(RColorBrewer)
library(Hmisc)
load('procSim.RData')
params <- as.data.frame(t(sapply(a,'[[','params')))

ns <- c(2,20)
anc <- FALSE
mut <- 'perTrait'
f <- 0.1667
overD <- TRUE
scenarios <- c('orth','straight')


ltys <- c(1,5)
D <- 100


#select colour palette
palpaired <- brewer.pal(12,name='Paired')
palpaired[11] <- brewer.pal(11,name="BrBG")[4]
cols <- matrix(palpaired,ncol=2,byrow=T)[c(2,3),]

graphics.off()
quartz(width=5,height=6)

layout(matrix(1:2,ncol=1), heights=c(1,2))

par(mai=c(0,.8,.5,0.1))
plot(0,type='n',
     ylim=c(0,.5),
     xlim=c(0,log10(D)),
     xaxt='n',yaxt='n',
     xlab='',
     ylab=expression(paste(italic(m),'(',bold(A),',',bold(A),')')), 
     main='')
mtext(side=3, '(A)',adj=.01,padj=-.1)

axis(side=2, at=c(0,.25,.5),las=T)
for(n in ns){
  for(s in 1:length(scenarios)){
    w <- which(params$n==n & 
                 params$sc == scenarios[s] & 
                 params$P1ancestral==anc & 
                 params$mutmodel==mut &
                 params$f==f & 
                 params$overD == overD &
                 params$D ==D)
    dat <- sapply(a[[w]]$res, FUN=function(x) x[,'ma'])
    means <- rowMeans(dat)
    sds <- apply(dat, 1, sd)
    polygon(log10(c(1:D,rev(1:D))),c(means+sds,rev(means-sds)),col=alpha(cols[s,1],.5),border=F)
    lines(log10(c(1:D)),means, col=cols[s,2],lty=ltys[which(ns==n)],lwd=2)
  }
}
text(log10(1),.48,'II:  Different traits',col=cols[1,2],adj=0)
text(log10(1),.42,'III: Opposite directions',col=cols[2,2],adj=0)
legend(x=log10(0.8), y=.4, lty=c(1,2), legend=c('n=2', 'n=20'), bty='n')

par(mai=c(1,.8,0,0.1))
plot(0,type='n',
     ylim=c(0,1.6),
     xlim=c(0,log10(D)),
     xaxt='n',
     xlab='Divergence (D)',
     ylab=expression(paste(italic(M),'(',bold(A),',',bold(A),') / ',italic(m),'(',bold(A),',',bold(A),')')), 
     main='')
mtext(side=3, '(B)',adj=.01,padj=1.6)
lines(log10(c(1:D)),1/(1:D),lwd=2,col='gray')

axis(side=1,at=log10(c(1,2,5,10,20,50,100)),labels=c(1,2,5,10,20,50,100))
for(n in ns){
  for(s in 1:length(scenarios)){
    w <- which(params$n==n & 
                 params$sc == scenarios[s] & 
                 params$P1ancestral==anc & 
                 params$mutmodel==mut &
                 params$f==f & 
                 params$overD == overD &
                 params$D ==D)
    m <- sapply(a[[w]]$res, FUN=function(x) x[,'ma'])
    M <- sapply(a[[w]]$res, FUN=function(x) x[,'Ma'])
    means <- rowMeans(M)/rowMeans(m)
    lines(log10(c(1:D)),means, col=cols[s,2],lty=which(ns==n),lwd=2) # ratio of mean M/m
  }
}

dy <- .08

#plot 'direct trajectory' label
lab <- unlist(strsplit('direct trajectory',split=''))
subD <- seq(log10(1.5), log10(10), length.out=length(lab))
points(subD,1/(10^subD)-dy,pch=lab,srt=-45,cex=.7, col='gray') 


dev.copy2pdf(file='Fig3.pdf')

