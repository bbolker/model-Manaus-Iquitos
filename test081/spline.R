predict<-function(xa, ya, y2a,  n, x)
	unlist(sapply(x, function(x)splint(xa, ya, y2a,  n,  x)))
splint<-function(xa, ya, y2a,  n,  x) {
klo=1; 
khi=n;
while (khi-klo > 1) {
k=(khi+klo)/2;
if(xa[k] > x) khi=k
else klo=k;
}                           
h=xa[khi]-xa[klo];
a=(xa[khi]-x)/h;                                          
b=(x-xa[klo])/h;                 
a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


spline<-function( x,  y, n,  yp1=0,  ypn=0)
{
u<-rep(0,n); 
y2<-rep(0,n); 
y2[1] = -0.5;
u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
qn=0.5;
un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
for (i in 2:(n-1)) { 
sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
p=sig*y2[i-1]+2.0;
y2[i]=(sig-1.0)/p;
u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
}
y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
for (k in (n-1):1)     
y2[k]=y2[k]*y2[k+1]+u[k];   
y2
}


