x<-read.csv('best_ncov.csv',r=1)
x<-x[,-match('phi',names(x))]
for(i in 1:ncol(x)){
if(abs(x[1,i])>1)x[,i]<-round(x[,i],3)
else x[,i]<-round(x[,i],5)
}
for(i in 1:7)
write.csv(x[,i*6+c(-5:0)],paste('out',i,'.csv',sep=""),quote=F)
