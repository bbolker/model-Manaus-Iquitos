dir <- getwd()
dir <- substr(dir,nchar(dir)-6,nchar(dir))

vacc<-read.csv('vacc20.csv')
vacc[,1]<-as.numeric(as.Date(vacc[,1])-as.Date('2019-12-31'))/366+2020
file <- paste('~/fig3','_',dir,'.pdf',sep="")
  pdf(file,width=10,height=5)
###
par(mar=c(3,4,2,4),mfrow=c(1,2),las=1,xaxs='i',yaxs='i')
  source('spline.R')
  require('pomp.orig')
  d1 <- as.numeric(as.Date(c('2020-1-31','2021-12-5'))-as.Date("2019-12-31"))/365.25+2020

#  par(mfrow=c(3,4),las=1,xaxs='i',yaxs='i',mar=c(2,4,2,4))
##################################################################################################

#cases<-read.csv('sari_city2.csv')
  load('mle_ncov001.rda')
  nd<-ncol(mle[[1]]@data)
  x <- read.csv('best_ncov.csv',row=1)
mm<-read.csv('model_parameters.csv',r=1)
ad<-sum(mm$type=='par')-length(grep('log.beta',rownames(mm)))
  y <- x
  np<-x$nm+ad+4
#nd<-nrow(cases)
print(nd)
  x$bic<- -2*x$loglik + np*log(nd)
  x$aicc<- -2*x$loglik+2*np*nd/(nd-np-1)
  y$bic<-x$bic
  y$aicc<-x$aicc
#  j <-  which.min(x$bic[1:10])
 j<-c(which.max(x$loglik[1:10]),10+which.max(x$loglik[11:20]))
 print(j)
jj<-j
#jj[1]<-1
#jj[2]<-11
jj<-c(1,12)
p<-0
Ro2<-as.null()
for(j in c(jj)){
p<-p+1
#if(p>3)par(mar=c(4,4,1,4))
if(p==4)plot(c(0,1),c(0,1),type='n',ann=F,axe=F)
#if(p>2)par(mar=c(4,4,1,4))
	print(j)
#country<-gsub(" ",".",read.csv('pop20.csv')[ifelse(j %% 20 ==0, 20, j%%20),1])
country<-'MANAUS'

	for(i in 1){
	y[i,] <- x[j[i],]
  load('allregion.rda')
        all[[i]] <- all[[j[i]]]
	}
  dyn.load('ncovmodel.so')
  source('ncovmodel.R')
  mle <- all[[1]]
  m<-mle@pred.mean
#  ar<-max(0.99-m[1,]/sum(m[,1]))
  params <- as.numeric(y[1,])
  names(params) <- colnames(y)
  mle@coef <- par.trans(params[match(names(coef(mle)),names(params))])
# mle@coef['p0']<- logit(1e-6)
  n <- 1000
  s <- simulate(mle,nsim=n,seed=100)
  print("succ")

  times <- mle@data[1,]
  death <- mle@data[3,]

  simul2 <- apply(sapply(s,function(x)x@data[3,]),1,median)
  s3 <- apply(sapply(s,function(x)x@data[3,]),1,function(x)quantile(x,0.025))
  s4 <- apply(sapply(s,function(x)x@data[3,]),1,function(x)quantile(x,0.975))

    matplot(times,sqrt(cbind(death,simul2,s3,s4)/y$bpop[1]*1e6),col='grey',lty=1,type='l',xlim=d1,ylim=c(0,35),
  ann=F,axe=F,frame=T,xlab='',ylab='')
#  mtext(side=3,line=0,paste('IAR',round(ar,2)))
  mtext(side=3,line=0,cex=0.8,bquote(psi==.(round(y$psi[1],2))))
  polygon(c(times,times[length(times):1]),sqrt(c(s3,s4[length(s3):1])/y$bpop[1]*1e6),col='lightgrey',border=NA)
#  lines(times,sqrt(death/y$bpop[1]*1e6),pch=20,col='red',cex=0.5)
  points(times,sqrt(death/y$bpop[1]*1e6),pch=21,bg='white',lty=1,col='red',lwd=1.5,cex=0.75)
  lines(times,sqrt(simul2/y$bpop[1]*1e6),col='black',lwd=2)
 month <- as.numeric(as.Date(c(paste(2020,seq(4,12,by=3),1,sep="-"),paste(2021,seq(1,12,by=3),1,sep="-")))-as.Date('2019-12-31'))/365.25+2020
  axis(1,at=month,lab=month.abb[c(seq(4,12,by=3),seq(1,12,by=3))],cex.axis=0.85,tck=-0.01,padj=-1.2)
  axis(1,at=month[seq(1,7,by=4)],lab=2020:2021,padj=0.2)
  axis(2,hadj=0.6,tck=-0.02)
  mtext(side=2,line=1.8,text=expression(sqrt(weekly~deaths~per~mil)),cex=0.75,las=0,col='black')

  mtext(side=3,line=0.5,adj=-0.1,text=letters[ifelse(p>3,1,0)+p],font=2,cex=1)

  #lines(times,sqrt(simul2/y$bpop[1]*1e6),col='black',lwd=2)
        col<-c('red','black','blue','brown','green')
  if(p==1)
        legend(times[38],28,bty='n',lty=c(1,1,1,1,1),lwd=1.5,pch=c(21,NA,3,NA,2),pt.bg=c('white',NA,NA,NA,NA),pt.cex=c(0.75,1,1,1,1),
	       col=col,text.col=col, cex=0.8,legend=c('reported deaths','simulation median','transmission rate','currently immunized','cumulative infected'),seg.len=3)
#text(times[30],4.5,expression(italic(R)[0](t)),col='black')
#text(times[30],4.0,expression(italic(S)(t)),col='black')
  a <- as.numeric(y[1,match('log.beta',names(y))+c(1:y$nm[1])-1])
 # slope <- a[length(a)]
 # a[length(a)]<-a[1]

  ms <- spline(seq(0,1,len=y$nm[1]/2),a[1:(y$nm[1]/2)],y$nm[1]/2)
  cut <- 2021-(3.36-y$eta[1])/33.6
  time2 <- seq(min(times),cut,by=1/366)
  my <- predict(seq(0,1,len=y$nm[1]/2),a[1:(y$nm[1]/2)],ms,y$nm[1]/2,(time2-min(time2))/(max(time2)-min(time2)))
  ms <- spline(seq(0,1,len=y$nm[1]/2),a[(y$nm[1]/2+1):y$nm[1]],y$nm[1]/2)
  time2 <- seq(cut,max(times),by=1/366)
  my <-c(my,predict(seq(0,1,len=y$nm[1]/2),a[(y$nm[1]/2+1):y$nm[1]],ms,y$nm[1]/2,(time2-min(time2))/(max(time2)-min(time2))))
  time2<- c(seq(min(times),cut,by=1/366),seq(cut,max(times),by=1/366))
  par(new=T)
  R0<-exp(my)*(1/y$gammah[1])
  plot(time2,R0, xlim=d1,type='l',lwd=0.5,cex=0.8,ylim=c(0,8),col='blue',ann=F,axe=F,frame=F,xlab='',ylab='')
  points(time2[seq(1,length(R0),by=15)],R0[seq(1,length(R0),by=15)],pch=3,lwd=0.5,col='blue')
  abline(h=1,lty=2,col='red')
  axis(4,tck=-0.02,hadj=0.2,at=0:5,col.lab='blue',col.ticks='blue')
  mtext(side=4,line=2.5,adj=0.3,text=expression(italic(R)[0](t)),cex=0.75,las=0,col='blue')
#  mtext(side=4,line=2.5,adj=0.95,text='vacc coverage',cex=0.75,las=0,col='brown')
  mtext(side=4,line=2,adj=0.9,text='immunized',cex=0.75,las=0,col='brown')
  mtext(side=4,line=2.7,adj=0.9,text='infected',cex=0.75,las=0,col='green')
#points(vacc[,1],vacc[,match(country,names(vacc))]/40+5.5,type='l',col='brown',lwd=2)
abline(h=5,col='lightgrey')
axis(4,tck=-0.02,hadj=0.2,at=seq(5,8,len=6),lab=c('',seq(0,1,len=6)[-1]),col.axis='brown') 
#load(sprintf('mle_ncov%3.3d.rda',jj))
ml<-as.null()
for(n in 0:1){
load(sprintf('mle_ncov%3.3d.rda',j+n*20))
ml<-c(ml,mle[[1]]@loglik)
}
n<-which.max(ml)-1
load(sprintf('mle_ncov%3.3d.rda',j+n*20))
#par(new=T)
points(mle[[1]]@data[1,],5+3*(0.99-mle[[1]]@pred.mean[1,]/y$bpop[1]),type='l',lwd=3,col='brown')
#points(mle[[1]]@data[1,],5+3*(mle[[1]]@pred.mean['CC',]/(y$bpop[1])),type='l',lwd=3,col='blue')
points(mle[[1]]@data[1,seq(1,89,by=4)],5+3*(mle[[1]]@pred.mean['BC',seq(1,89,by=4)]/(y$bpop[1])),type='b',lwd=1,pch=2,col='green')
lines(mle[[1]]@data[1,seq(1,89,by=4)],5+3*(mle[[1]]@pred.mean['BC',seq(1,89,by=4)]/(y$bpop[1])),lwd=1,col='green')
#mtext(side=3,line=-1,adj=0.1,cex=0.8,bquote(MLL==.(round(mle[[1]]@loglik,3))))
abline(h=0.7*3+5,col='red',lty=2)
abline(h=0.6*3+5,col='red',lty=2)
abline(h=0.4*3+5,col='red',lty=2)
abline(h=0.2*3+5,col='red',lty=2)
abline(v=(as.Date("2020-10-1")-as.Date("2019-12-31"))/366+2020,lty=2)
abline(v=(as.Date("2020-7-1")-as.Date("2019-12-31"))/366+2020,lty=2)

#plot(times,mle[[1]]@pred.mean[1,]/6543000, xlim=d1,type='l',lty=4,lwd=2,cex=0.85,ylim=c(0,0.7),col='darkgreen',ann=F,axe=F,frame=F,xlab='',ylab='')

# points(times2,temp-min(temp),type='b',pch=23,bg='yellow',lwd=1,cex=0.5)
# points(times2,rain/5,type='b',pch=22,bg='green',lwd=1,cex=0.75)
#  axis(4,tck=-0.02,hadj=0.2,at=seq(0,0.7,by=0.2),lab=ifelse(j==2,T,T),col='blue',col.ticks='blue')
#  mtext(side=4,line=2,adj=0.5,text=expression(italic(R)[0](t)),cex=1.2,las=0,col='blue')
#  mtext(side=4,line=2.2,adj=0.5,text=expression(italic(S)(t)),cex=1.2,las=0,col='darkgreen')
#par(new=T)
#  plot(times2,sqrt(R2), xlim=d1,type='l',lty=2,lwd=1,cex=0.85,ylim=c(0,20),col='blue',ann=F,axe=F,frame=F,xlab='',ylab='')


#for(k in 1:length(ts))
#lines(c(ts[k],ts[k]),c(0,w[j,k]/200),lwd=3,col='green')
print(y[1,])
if(p %in% c(1,2)) mtext(side=3,line=-1,adj=0.95,ifelse(p<2,'Manaus','Iquitos'))
if(FALSE){
for(p in c(3:4)){
if(p==3) par(mar=c(3,4,2,4),new=T,fig=c(0.75,1,0.5,1))
if(p==4) par(mar=c(4,4,1,4),new=T,fig=c(0.75,1,0,0.5))
x<-read.csv('best_ncov.csv',r=1)[c(1:10)+10*(p-3),]
x<-x[sort(x$psi,dec=T,index=T)$ix,]
#print(x)
plot(x$psi,x$loglik,type='b',log='',ylab='log likelihood',xlab=expression(psi),ylim=c(ifelse(p==3,-8,-7),1)+round(max(x$loglik)),lwd=2)
abline(h=max(x$loglik)-0.5*qchisq(0.95,1),lty=2)
mtext(side=1,line=3,adj=0.7,'')
par(new=T)
plot(x$psi,x$theta*x$theta*100,type='b',col='red',log='',pch=22,ann=F,axe=F,frame=T,ylim=c(0.25, 1.5))
points(x$psi,x$alpha*x$theta*x$theta*100,type='b',col='blue',pch=23)
axis(4)
mtext(side=4,line=3,las=0,'IFR (%)')
  mtext(side=3,line=0.5,adj=-0.15,text=letters[ifelse(p==3,4,8)],font=2,cex=1)
if(p==3)
legend("bottomleft",bty='n',pch=21:23,col=c('black','red','blue'),legend=c('log likelihood','IFR pre','IFR post'),lty=1,lwd=c(2,1,1))

mtext(side=3,line=-1,adj=0.05,ifelse(p==3,'Manaus','Iquitos'))
#abline(v=x$psi[6],lty=2,col='lightgrey')
}
}
}
dev.off()

