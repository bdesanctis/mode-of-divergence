library(e1071)
library(scales)
library(plotrix)
library(RColorBrewer)

palpaired <- brewer.pal(12,name='Paired')
palpaired[11] <- brewer.pal(11,name="BrBG")[4]
cols <- matrix(palpaired,ncol=2,byrow=T)[c(1,3,5),]

draw_triangle <- function(mm,d,dom,arrowl,col,col2,draw=T){
  lastpos <- mm[1,]

  for(i in 2:d){
    newpos <- c(lastpos[1]+(mm[i,1]-mm[i-1,1])*dom[i,1],lastpos[2]+(mm[i,2]-mm[i-1,2])*dom[i,2])
    if(draw) arrows(x0=lastpos[1],x1=newpos[1],y0=lastpos[2],y1=newpos[2],col=col,length=arrowl)
    lastpos <- newpos
  }
  
  f1pos <- lastpos
  if(draw){
    for(i in 2:d){
    newpos <- c(lastpos[1]+(mm[i,1]-mm[i-1,1])*(1-dom[i,1]),lastpos[2]+(mm[i,2]-mm[i-1,2])*(1-dom[i,2]))
    arrows(x0=lastpos[1],x1=newpos[1],y0=lastpos[2],y1=newpos[2],col=col2,length=arrowl)
    lastpos <- newpos
    }
  }
  return(f1pos)
}

dplotA <- function(P1.pos, P2.pos, MRCA.pos,mm,dom,main='',xlab='trait 1', ylab='trait 2',arrowl=.15,doA=T,cloud=T){
  
  #draw basic landscape
  plot(0,0, pch=15, 
       xlim=c(-.8,.8), ylim=c(-.8,.8),
       xlab='',
       ylab='', axes=F, frame.plot=T,
       main=main)
  mtext(side=1,xlab,.4)
  mtext(side=2,ylab)
  
  for(i in 1:10){
    rad <- seq(1,.1, by=-.1)[i]
    bgcol <- colorRampPalette(c("grey50", "grey100"))(11)[11-i]
    symbols(x=0, y=0, circles=rad, inches=F, add=T, bg=bgcol, fg=bgcol)
  }
  box()
  abline(v=0, lty=1, col='dark grey')
  abline(h=0, lty=1, col='dark grey')
  
  
  
  d <- nrow(mm)
  
  
  #subs
  arrows(x0=mm[1:(d-1),1],x1=mm[2:d,1],y0=mm[1:(d-1),2],y1=mm[2:d,2],length=arrowl,lwd=2)#,col='red')

  # #dom subs
  F1.pos <- draw_triangle(mm=mm,d=d,dom=dom,arrowl/2,col=cols[2,2],col2=cols[1,2])
  #F1
  points(x=F1.pos[1], y=F1.pos[2], col='black', pch=16,cex=1)
  text(x=F1.pos[1]-.05, y=F1.pos[2]+.07, expression(bold(bar(z)['F1'])), col='black')

  
  
  #P1
  if(cloud){
    points(x=rnorm(n=30,mean=P1.pos[1],sd=.1), 
         y=rnorm(n=30,mean=P1.pos[2],sd=.1), 
         col='black', pch=16,cex=1)
  }
  points(x=P1.pos[1], y=P1.pos[2], col='black', pch=16,cex=1)
  text(x=P1.pos[1], y=P1.pos[2]-.08, expression(bold(bar(z)['P1'])), col='black')
  
  #P2
  if(cloud){
    points(x=rnorm(n=30,mean=P2.pos[1],sd=.1), 
           y=rnorm(n=30,mean=P2.pos[2],sd=.1), 
           col='black', pch=16,cex=1)
  }
  points(x=P2.pos[1], y=P2.pos[2], col='black', pch=16,cex=1)
  text(x=P2.pos[1], y=P2.pos[2]+.08, expression(bold(bar(z)['P2'])), col='black')
  
  
  #OPT
  points(0,0, pch=4, cex=2, lwd=3)
  text(x=.075,y=-.075, 'opt')
  
  
  
  #shortest paths
  segments(x0=P1.pos[1],x1=P2.pos[1],y0=P1.pos[2],y1=P2.pos[2],lty='dashed') #add lines
  segments(x0=P1.pos[1],x1=F1.pos[1],y0=P1.pos[2],y1=F1.pos[2],lty='dashed',cols[2,2]) #add lines
  segments(x0=F1.pos[1],x1=P2.pos[1],y0=F1.pos[2],y1=P2.pos[2],lty='dashed',cols[1,2]) #add lines

  legend('bottomright', 
         legend=c(expression(paste(italic(M),'(2',bold(A),',2',bold(A),')')),
                  expression(paste(italic(m),'(2',bold(A),',2',bold(A),')')),
                  expression(paste(italic(M),'(',bold(A),'+',bold(Delta),',',bold(A),'+',bold(Delta),')')),
                  expression(paste(italic(m),'(',bold(A),'+',bold(Delta),',',bold(A),'+',bold(Delta),')')),
                  expression(paste(italic(M),'(',bold(A),'-',bold(Delta),',',bold(A),'-',bold(Delta),')')),
                  expression(paste(italic(m),'(',bold(A),'-',bold(Delta),',',bold(A),'-',bold(Delta),')'))),
         lty=rep(c(1,2),3),
         col=rep(c('black',
                   cols[2,2],
                   cols[1,2]),each=2),
         text.col=rep(c('black',
                        cols[2,2],
                        cols[1,2]),each=2),
         bg='white'
  )
}

dplotB  <- function(P1.pos, P2.pos, MRCA.pos,mm,dom,main='',xlab='trait 1', ylab='trait 2',arrowl=.15,doA=T){
  
  MP.pos <- c(mean(c(P1.pos[1],P2.pos[1])), mean(c(P1.pos[2], P2.pos[2])))
  #draw basic landscape
  plot(0,0, pch=15, 
       xlim=c(-.8,.8), ylim=c(-.8,.8),
       xlab='',
       ylab='', axes=F, frame.plot=T,
       main=main)
  mtext(side=1,xlab,.4)
  mtext(side=2,ylab)
  
  for(i in 1:10){
    rad <- seq(1,.1, by=-.1)[i]
    bgcol <- colorRampPalette(c("grey50", "grey100"))(11)[11-i]
    symbols(x=0, y=0, circles=rad, inches=F, add=T, bg=bgcol, fg=bgcol)
  }
  box()
  abline(v=0, lty=1, col='dark grey')
  abline(h=0, lty=1, col='dark grey')
  
  
  
  d <- nrow(mm)
  
  
  #subs
  startpos <- P1.pos
  for(i in 2:d){
    newpos <- c(startpos[1]+(mm[i,1]-mm[i-1,1])/2, startpos[2]+(mm[i,2]-mm[i-1,2])/2)
    arrows(x0=startpos[1], x1=newpos[1],y0=startpos[2], y1=newpos[2], length=arrowl)
    startpos <- newpos
  }
  
  startpos <- MP.pos
  delta <- dom - 1/2
  for(i in 2:d){
    newpos <- c(startpos[1]+delta[i,1]*(mm[i,1]-mm[i-1,1]), startpos[2]+delta[i,2]*(mm[i,2]-mm[i-1,2]))
    arrows(x0=startpos[1], x1=newpos[1],y0=startpos[2], y1=newpos[2], col=cols[3,2],length=arrowl/2)
    startpos <- newpos
  }
  
  arrows(.3,.3,.6,.3, length=arrowl)
  arrows(.6,.3,.8,.5, length=arrowl)
  segments(.6,.3,.8,.3,col='gray')
  draw.arc(.6,.3,angle1=0,
           angle2=pi/4,radius=.12,lwd=.7)
  text(.75,.37,expression(theta), cex=1.1)
  
  F1.pos <- draw_triangle(mm=mm,d=d,dom=dom,arrowl=arrowl,
                          col=cols[2,2],
                          col2=cols[1,2],
                          draw=F)
  
  
  #F1
  points(x=F1.pos[1], y=F1.pos[2], col='black', pch=16,cex=1)
  text(x=F1.pos[1]-.08, y=F1.pos[2], expression(bold(bar(z)['F1'])), col='black')
  
  #P1
  points(x=P1.pos[1], y=P1.pos[2], col='black', pch=16,cex=1)
  text(x=P1.pos[1], y=P1.pos[2]-.08, expression(bold(bar(z)['P1'])), col='black')
  
   #MP
  points(x=MP.pos[1], y=MP.pos[2], col='black', pch=16,cex=1)
  text(x=MP.pos[1]+.1, y=MP.pos[2], expression(bold(bar(z)['mp'])), col='black')
  
  
  #OPT
  points(0,0, pch=4, cex=2, lwd=3)
  text(x=.075,y=-.075, 'opt')
  
  segments(x0=P1.pos[1],x1=MP.pos[1],y0=P1.pos[2],y1=MP.pos[2],lty='dashed') #add lines
  segments(x0=MP.pos[1],x1=F1.pos[1],y0=MP.pos[2],y1=F1.pos[2],lty='dashed',cols[3,2]) #add lines
  
  legend('bottomright', 
         legend=c(expression(paste(italic(M),'(',bold(A),',',bold(A),')')),
                  expression(paste(italic(m),'(',bold(A),',',bold(A),')')),
                  expression(paste(italic(M),'(',bold(Delta),',',bold(Delta),')')),
                  expression(paste(italic(m),'(',bold(Delta),',',bold(Delta),')'))),
         
         lty=rep(c(1,2),3),
         col=rep(c('black',
                   cols[3,2]),each=2),
         text.col=rep(c('black',
                        cols[3,2]),each=2),
         bg='white'
  )
}


#set position of P1, P2, and MRCA
P1.pos <- c(-.7,-0.45)
P2.pos <- c(.3,.7)
MRCA.pos <- c(-.25,.22)
pos <- rbind(P1.pos,P2.pos,MRCA.pos)

#positions of substitutions
mm <- matrix(c(P1.pos,
                -0.3,-0.05,
                MRCA.pos,
                0.1,.12,
                .05,.45,
                P2.pos),byrow=T,ncol=2)


dom <- matrix(c(NA,NA,
                0.3, 0.8,
                0.05,.9,
                0.28,.8,
                0.9, 0.4,
                .5,.2,
                .3,.8),ncol=2,byrow=T)


#make figure
graphics.off()
quartz(width=10*.9,height=6*.9)

layout(matrix(1:2,ncol=2))
dplotA(P1.pos, P2.pos, MRCA.pos,mm,dom,arrowl=.1,cloud=F)
mtext(side=3,adj=0,'(A)')

dplotB(P1.pos, P2.pos, MRCA.pos,mm,dom,arrowl=.1)
mtext(side=3,adj=0,'(B)')

dev.copy2pdf(file='Fig1.pdf')
