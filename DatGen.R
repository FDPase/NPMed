#tt<-3
#vvv<-0 		# also change vv in code.R
#ie<-0
#de<-0

g.0=0.2*0	
a.0=0.3*0; a.x=0.5; a.z=0.5
b.0=0.4*0; b.x=0.5; b.z=1.5; b.m=2.0

a.z=a.z*ie
b.z=b.z*de
b.m=b.m*(1-(1-ie)*(1-vvv))

genM<-function(a0=a.0, ax=a.x, az=a.z, Z, rZ, X, n, vv=vvv, type=tt, the=0){
	if (the==1) X<-g.0
	if (type==1) M<-a0+ax*X+az*Z+rnorm(n, sd=1+Z*vv)
	if (type==2) M<-a0+ax*X+az*Z+rt(n, df=7)*(2-(1-Z*vv))	# var=df/(df-2)
	if (type==3) M<-a0+ax*X+az*Z+ifelse(rZ==1, rnorm(n, sd=2), rnorm(n, sd=1))
	if (type==4) M<-a0+ax*X+az*Z+ifelse(rZ==1, 2*rt(n, df=7), rt(n, df=7))
	return(M)
}

genY<-function(b0=b.0, bx=b.x, bz=b.z, bm=b.m, Z, rZ, M, X, n, vv=vvv, type=tt, the=0){
	if (the==1) X<-g.0
	if (type==1) Y<-b0+bx*X+bz*Z+bm*M+rnorm(n, sd=1+Z*vv)
	if (type==2) Y<-b0+bx*X+bz*Z+bm*M+rt(n, df=7)*(2-(1-Z*vv))
	if (type==3) Y<-b0+bx*X+bz*Z+bm*M+ifelse(rZ==1, rnorm(n, sd=3), rnorm(n, sd=1))
	if (type==4) Y<-b0+bx*X+bz*Z+bm*M+ifelse(rZ==1, 2*rt(n, df=7), rt(n, df=7))
	return(Y)
}

GenDat<-function(
	seed=37, 
	n=1000, 
	type=tt,
	vv=vvv,
#	de=DE,
#	ie=IE,
	g0=g.0, gx=0.5,
	a0=a.0, ax=a.x, az=a.z,
	b0=b.0, bx=b.x, bz=b.z, bm=b.m
){

	set.seed(seed)
#n<-1000
	X<-rnorm(n)
#type<-1
#vv<-1
#de<-1
#ie<-1
	
#g0=0.2; gx=0.5
#a0=0.3; ax=0.5; az=0.5*ie
#b0=0.4; bx=0.5; bz=1.5*de; bm=2.0*(1-(1-ie)*(1-vv))
	
#	az=az*ie
#	bz=bz*de
#	bm=bm*(1-(1-ie)*(1-vv))
	
	eta<-g0+gx*X
	Z<-rbinom(size=1, n=n, prob=exp(eta)/(1+exp(eta)))
	if (vv==1) rZ<-rbinom(size=1, n=n, prob=ifelse(Z>0, 0.8, 0.2))
	if (vv==0) rZ<-rbinom(size=1, n=n, prob=0.3)
	
	M<-genM(a0, ax, az, Z, rZ, X, n, vv, type)
	Y<-genY(b0, bx, bz, bm, Z, rZ, M, X, n, vv, type)
	
	out<-list(M=M, Y=Y, Z=Z, X=X)
	return(out)
}
