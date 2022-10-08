library(RColorBrewer)
library(scales)

load('procSim_FigS1.RData')


palpaired <- brewer.pal(12,name='Paired')
palpaired[11] <- brewer.pal(11,name="BrBG")[4]
cols <- matrix(palpaired,ncol=2,byrow=T)[c(1,3),]


pchs <- c(1:2)

unmat=as.data.frame(cbind(s=rep(c('0.0001','0.01'),each=2,times=2),N=rep(c('10','1000'),each=4),n=rep(c('2','20'),4)),stringsAsFactors = T)
df$col <- apply(df[,c('s','N','n')],1, FUN=function(x) {which(unmat$s ==x['s'] & unmat$N == x['N'] & unmat$n==x['n'])})

mylabels <- c(expression(italic(N)),10,expression(10^3),
              expression(italic(bar(s))['mut']),expression(10^-2),expression(10^-4),
              expression(italic(n)),'2','20')

graphics.off()
quartz(width=10,height=6)
layout(matrix(1:2,ncol=2))

for (MVN in c(1,0)){
  df2 <- df[which(df$k==k & df$M == MVN),]
  plot(x=df2$ma, 
       y=df2$Ma,
       type='n',
       log='xy',
       ylim=c(.002,1.5),
       xlim=c(min(df$ma),1),
       xlab=expression(paste(italic(m),'(',bold(A),',',bold(A), ')  Net effect of evolutionary change')),
       ylab=expression(paste(italic(M),'(',bold(A),',',bold(A), ')  Total amount of evolutionary change')),
       xaxt='n')

  axis(1,at=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^0),c(expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1),expression(10^0)))
  
  df2.n2 <- df2[which(df2$n==2),]
  df2.n20 <- df2[which(df2$n==20),]
  
  points(x=df2.n2$ma,
         y=df2.n2$Ma,
         col=cols[df2.n2$N,1],
         pch=pchs[df2.n2$n]+15*c(0,1)[df2.n2$s])
  points(x=df2.n20$ma,
         y=df2.n20$Ma,
         col=cols[df2.n20$N,2],
         pch=pchs[df2.n20$n]+15*c(0,1)[df2.n20$s])

  predload <- outer(as.numeric(levels(df$N)),as.numeric(levels(df$n)), FUN=function(y,x) {x/(4*y)/2})
  segments(x0=predload,x1=predload,y0=predload,y1=20, col=cols,lty=3)
  
  mtext(side=3,ifelse(MVN, '(A) Multivariate normal distribution','(B) Exponential distribution'),adj=0)
  if(MVN==1){
    legend('bottomright',legend=mylabels,
           col=c('white',cols[1,2],cols[2,2],
                 'white','black',alpha('black',1),
                 'white','black','black'),
           pch=c(NA,15,15,
                 NA,15,0,
                 NA,pchs),
           cex=.8,ncol=3,
           text.width = .4,y.intersp=1)
  }
  abline(a=0,b=1,col='gray')
  text(x=.07,y=.1,expression(paste(italic(m),'(',bold(A),',',bold(A), ') = ',
                                   italic(M),'(',bold(A),',',bold(A), ')')),
       col='gray',srt=62,cex=.7)
  
}

mtext(side=3,expression(bold('Properties of fixed differences under stabilizing selection')),adj=4,padj=-2,cex=1)

dev.copy2pdf(file='FigS1.pdf')
