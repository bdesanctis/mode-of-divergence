library(scales)
library(RColorBrewer)
greycol <- 'darkgrey'
cex=.8

f <- 0.1667 #variance in dominance coefficients
n <- 20 # number of traits
ancs <- c(FALSE,TRUE) #whether both evolve (FALSE), or only P2 (TRUE)
D <- 50
mut <- 'perTrait'
overD <- FALSE
scenarios <- c('return','orth','straight') #select only divergence scenarios involving directional selection
nscen <- length(scenarios)

load("procSim.RData")
params <- as.data.frame(t(sapply(a, '[[','params')))

botax <- 3; topax <- 3

#set up colour palette
palpaired <- brewer.pal(nscen*2,name='Paired')
if(nscen*2>10) palpaired[11] <- brewer.pal(11,name="BrBG")[4]
cols <- matrix(palpaired,ncol=2,byrow=T)


#manually set range of y-axis
mins <- matrix(
  c(0,  0,-.01,
    0,  0,-.1,
    -.35,-.05,-.1),
  byrow=T,ncol=3)
maxs <- matrix(
  c(.2,.1/2,.02,
    .7,.15,.2,
    .7,.1,.2),
  byrow=T,ncol=3)

#combinations of stats to plot
paths <- c('a','d','ad')
stats <- c('M','m','f')

#titles
mains <- c('Additive effects', 'Dominance effects','Additive-by-dominance')

subtitls <- matrix(c(
  expression(paste('(A) Total amount, ',italic(M),'(',bold(A),',',bold(A),')')),
  expression(paste('(D) ',italic(M),'(',bold(Delta),',',bold(Delta),')')),
  expression(paste('(G) ',italic(M),'(',bold(A),',',bold(Delta),')')),
  
  expression(paste('(B) Net effect, ',italic(m),'(',bold(A),',',bold(A),')')),
  expression(paste('(E) ',italic(m),'(',bold(Delta),',',bold(Delta),')')),
  expression(paste('(H) ', italic(m),'(',bold(A),',',bold(Delta),')')),
  
  expression(paste('(C) ',italic(m),'(',bold(A),',',bold(A),') - ',italic(M),'(',bold(A),',',bold(A),')')),
  expression(paste('(F) ',italic(m),'(',bold(Delta),',',bold(Delta),') - ',italic(M),'(',bold(Delta),',',bold(Delta),')')),
  expression(paste('(I) ',italic(m),'(',bold(A),',',bold(Delta),') - ',italic(M),'(',bold(A),',',bold(Delta),')'))
  
), nrow=length(stats), ncol=length(paths),byrow=T)


dplot_path_cartoon_short <- function(scenario,cols,myy,anc){
  if(anc){
    mycol <- cols[2]
  }else{
    mycol <- cols[c(1,2)]
  }
  
  if(scenario=='return'){
    a.xs <- c(-1,-.1,-.1,-1)
    a.ys <- c(.05,.05,-.05,-.05)
    dy1 <- .2; dx1 <- .1
    dy2 <- -.2; dx2 <- .1

    if(anc){
      d.xs <- c(.2,-.2)
      d.ys <- c(.1,-.1)
      MRCAx <- a.xs[1]; MRCAy <- a.ys[1]
    }else{
      d.xs <- c(-.2,-.2)
      d.ys <- c(.1,-.1)
      MRCAx <- 0; MRCAy <- 0
      text(x=MRCAx+.15, y=MRCAy+.07+myy,'ancestral',adj=0,cex=.8)
      text(x=MRCAx+.15, y=MRCAy-.07+myy,'phenotype',adj=0,cex=.8)
    }
    
  }else if(scenario=='orth'){
    a.xs <- c(0,0,.1,1)
    a.ys <- c(-.7,-.1,0,0)
    if(anc){
      d.xs <- c(.1,.2)
      d.ys <- c(.2,.1)
      MRCAx <- a.xs[1]; MRCAy <- a.ys[1]
    }else{
      d.xs <- c(.1,.2)
      d.ys <- c(-.2,.1)
      MRCAx <- 0; MRCAy <- 0
    }
    dx1 <- -.3; dx2 <- 0
    dy1 <- 0; dy2 <- -.2
    
    
  }else if(scenario =='straight'){
    a.xs <- c(-1,-.1,.1,1)
    a.ys <- c(0,0,0,0)
    if(anc){
      d.xs <- c(.2,.2)
      d.ys <- c(.1,-.1)
      MRCAx <- a.xs[1]; MRCAy <- a.ys[1]
    }else{
      d.xs <- c(-.2,.2)
      d.ys <- c(.1,-.1)
      MRCAx <- 0; MRCAy <- 0
    }
    dx1 <- 0; dx2 <- 0
    dy1 <- .2; dy2 <- .2
  }

  points(x=a.xs[1],y=a.ys[1]+myy,col=cols[1], pch=16)
  points(x=a.xs[4], y=a.ys[4]+myy, col=cols[2], pch=16)
  points(x=MRCAx, y=MRCAy+myy, pch=1)
  
  text(x=a.xs[1]+dx1, y=a.ys[1]+dy1+myy, col=cols[1],'P1')
  text(x=a.xs[4]+dx2, y=a.ys[4]+dy2+myy, col=cols[2],'P2')
  
  sc <- .9
  arrows(x0=a.xs[1]*sc,
         x1=a.xs[2]*sc,
         y0=a.ys[1]*sc+myy,
         y1=a.ys[2]*sc+myy,
         length = .06,
         col=alpha(mycol[1],1),lwd=2)
  arrows(x0=a.xs[3]*sc,
         x1=a.xs[4]*sc,
         y0=a.ys[3]*sc+myy,
         y1=a.ys[4]*sc+myy,
         length = .06,
         col=alpha(mycol[length(mycol)],1),lwd=2)
  
  arrows(x0=(a.xs[1]/2+a.xs[2]/2)*sc,
         x1=(a.xs[1]/2+a.xs[2]/2)*sc+d.xs[1],
         y0=(a.ys[1]/2+a.ys[2]/2)*sc+myy,
         y1=(a.ys[1]/2+a.ys[2]/2)*sc+d.ys[1]+myy,
         length = .04,
         col=alpha(mycol[1],.7))
  arrows(x0=(a.xs[3]/2+a.xs[4]/2)*sc,
         x1=(a.xs[3]/2+a.xs[4]/2)*sc+d.xs[2],
         y0=(a.ys[3]/2+a.ys[4]/2)*sc+myy,
         y1=(a.ys[3]/2+a.ys[4]/2)*sc+d.ys[2]+myy,
         length = .04,
         col=alpha(mycol[length(mycol)],.7))
  points(x=MRCAx, y=MRCAy+myy, pch=1)
  
}

#initiate plotting device
graphics.off()
quartz(width=9,height=6)
layout(matrix(c(rep(1,3),rep(2,3),3,4,5,6,7,8,9,10,11),3,5),widths=c(0.6,0.6,1,1,1))
par(mar=c(botax+4,1,topax+4,0.1),xaxs='i',yaxs='i',xpd=NA)

#MSET UP GRID FOR LANDSCAPE CARTOONS
for(anc in ancs){
  plot.new()
  plot.window(xlim=c(-1.3,1.3),ylim=c(0,nscen*2))
  lines(c(-1.3,-1.3),c(0,nscen*2))
  lines(c(1.3,1.3),c(0,nscen*2))
  for(i in seq(0,nscen*2,by=2))
    lines(c(-1.3,1.3),c(i,i))
  for(i in seq(0,nscen*2-2,by=2))
    lines(c(0,0),c(i,i+1.75),lty='dotted',lwd=.5,col='gray')
  for(i in seq(1,nscen*2-1,by=2))
    lines(c(-1.3,1.3),c(i,i),lty='dotted',lwd=.5,col='gray')
  
  #MAKE LANDSCAPE CARTOONS
  for(i in 1:length(unique(scenarios))){
    scenario <- unique(scenarios)[i]
    dplot_path_cartoon_short(scenario, cols[i,], myy=nscen*2+1-2*i,anc=anc)
  }
  
  #ANNOTATE LANSCAPE CARTOONS WITH SCENARIO CLASSIFICATIONS
 if(anc){
  titls <- c(
    'VI: P2 distant opt.','V: Different traits','IV: P2 moving opt.')
  }else{
    titls <- c(
    'III: Opposite dir.','II: Different traits','I: Same direction')
  }
  
  eps <- 0.2
  for(i in 1:(nscen))
    text(y=(i-1)*2+2-eps,x=-1.2,titls[i],col=cols[nscen-i+1,2],adj=0)
  
  if(anc){
    text(x=0, y=nscen*2+.1, 'Only P2 evolves', col=greycol)
  }
  if(!anc){
    text(x=1.5,y=nscen*2+.4,'Divergence Scenarios',cex=1.2)
    text(x=0, y=nscen*2+.1, 'Both evolve', col=greycol)
  }
}

#PREPARE DATA FOR BOXPLOTS
par(mar=c(botax,3.5,topax,1),xaxs='i',yaxs='i')
w <- which(params$n==n & 
             params$sc %in% scenarios &
             params$mutmodel==mut &
             params$f==f & 
             params$overD == overD &
             params$D ==D)
boxdata <- do.call(rbind,lapply(a[w], FUN=function(x) {cbind(x$res,x$params['sc'],x$params['n'],x$params['P1ancestral'])}))
boxdata <- as.data.frame(boxdata)
colnames(boxdata) <- c(colnames(a[[w[1]]]$res),'sc','n','P1ancestral')

boxdata$sc2 <- boxdata$sc
boxdata$sc2[which(boxdata$sc2=='straight' & boxdata$P1ancestral==FALSE)] <- 'III'
boxdata$sc2[which(boxdata$sc2=='orth' & boxdata$P1ancestral==FALSE)] <- 'II'
boxdata$sc2[which(boxdata$sc2=='return' & boxdata$P1ancestral==FALSE)] <- 'I'
boxdata$sc2[which(boxdata$sc2=='straight' & boxdata$P1ancestral==TRUE)] <- 'VI'
boxdata$sc2[which(boxdata$sc2=='orth' & boxdata$P1ancestral==TRUE)] <- 'V'
boxdata$sc2[which(boxdata$sc2=='return' & boxdata$P1ancestral==TRUE)] <- 'IV'



#MAKE BOXPLOTS
for(i in 1:length(stats)) 
{
  for(j in 1:length(paths))
  {
    plot(0,type='n',axes=F,
         ylim=c(mins[i,j],maxs[i,j]),
         xlim=c(0.5,nscen*2+.5),xlab='',ylab='')
    if(mins[i,j]!=0) lines(c(0.5,nscen*2+.5),c(0,0),lty='dotted')
    boxplot(boxdata[,paste0(stats[i],paths[j])]~boxdata[,'sc2'],
            col=cols[,1],
            ann=F,
            axes=F,
            border=cols[,2],
            outline=F,add=T,range=0)
    if(stats[i]=='m' & paths[j]=='a'){
      axis(side=2,
           at=seq(0,maxs[i,j],by=.25),
           mgp=c(3, .5, 0),
           labels=seq(0,maxs[i,j],by=.25),
           padj=0)
    }else{
      axis(side=2,
           at=c(mins[i,j],0,maxs[i,j]),
           mgp=c(3, .5, 0),
           labels=c(as.character(mins[i,j]),0,as.character(maxs[i,j])),
           padj=0)
    }
    mtext(side=3,subtitls[i,j],adj=0,cex=0.9)
    
    if(i==1){
      mtext(side=3,mains[j],adj=0,padj=-1.5,cex=1.1)
    }

    ts <- c('I','II','III','IV','V', 'VI')
    mtext(side=1,at=1:length(ts),text=ts,col=cols[,2],padj=.3,cex=.9)
    lines(c(.75,3.25),rep(mins[i,j]-(maxs[i,j]-mins[i,j])/4.5,2),lwd=1.5,col=greycol)
    lines(c(3.75,6.25),rep(mins[i,j]-(maxs[i,j]-mins[i,j])/4.5,2),lwd=1.5,col=greycol)
    mtext(side=1,at=(.75+3.25)/2,text='Both evolve',col=greycol,padj=2,cex=0.6)
    mtext(side=1,at=(3.75+6.25)/2,text='Only P2 evolves',col=greycol,padj=2,cex=0.6)
    
    box()
  }
  
}

dev.copy2pdf(file=paste0("Fig2.pdf"))

