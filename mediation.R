mediation.kern<-function(Yr, Mr, Zr){

	n<-length(Yr)
	h<-n^(-1/3)
	z1<-1
	z0<-0
	
	### Fm
	
	kbi0<-dnorm((z0-Zr)/h)/h
	Fm0.temp<-cumsum(kbi0[order(Mr)])/sum(kbi0)
	dFm0.temp<-Fm0.temp-c(0, Fm0.temp[1:(length(Fm0.temp)-1)])
	dFm0<-dFm0.temp[rank(Mr)]
	
	kbi1<-dnorm((z1-Zr)/h)/h
	Fm1.temp<-cumsum(kbi1[order(Mr)])/sum(kbi1)
	dFm1.temp<-Fm1.temp-c(0, Fm1.temp[1:(length(Fm1.temp)-1)])
	dFm1<-dFm1.temp[rank(Mr)]
	
	### Fy
	
	Mrdiff.temp<-lapply(Mr, function(x)x-Mr)
	Mrdiff.mat<-matrix(unlist(Mrdiff.temp), nc=n)
	Kh.2<-dnorm(-Mrdiff.mat/h)
	
	Kh0.1.temp<-dnorm((z0-Zr)/h)
	Kh0.1<-matrix(rep(Kh0.1.temp, n), nc=n)
	Kh1.1.temp<-dnorm((z1-Zr)/h)
	Kh1.1<-matrix(rep(Kh1.1.temp, n), nc=n)
	
	Fy0.den.temp<-Kh.2*Kh0.1/h^2
	Fy1.den.temp<-Kh.2*Kh1.1/h^2
	Fy0.den<-apply(Fy0.den.temp, 2, sum)
	Fy1.den<-apply(Fy1.den.temp, 2, sum)
	
	
	#####################################################################
		
	Fy0.num.yy<-apply(Fy0.den.temp[order(Yr),], 2, cumsum)
	Fy0.yy<-t(Fy0.num.yy)/Fy0.den
	Fy1.num.yy<-apply(Fy1.den.temp[order(Yr),], 2, cumsum)
	Fy1.yy<-t(Fy1.num.yy)/Fy1.den
	
	F00<-apply(Fy0.yy*dFm0, 2, sum)
	F10<-apply(Fy1.yy*dFm0, 2, sum)
	F11<-apply(Fy1.yy*dFm1, 2, sum)
	FDE<-apply((Fy1.yy-Fy0.yy)*dFm0, 2, sum)
	FIE<-apply(Fy1.yy*(dFm1-dFm0), 2, sum)
	dFIE<-FIE-c(0, FIE[1:length(FIE)-1])
	dFDE<-FDE-c(0, FDE[1:length(FDE)-1])
	Q1.IE<-sum(sort(Yr)*dFIE)
	Q1.DE<-sum(sort(Yr)*dFDE)
	
	Q2.IE<-sum(sort(Yr)^2*dFIE)
	Q2.DE<-sum(sort(Yr)^2*dFDE)
	
	Q3.IE<-sum(sort(Yr)^3*dFIE)
	Q3.DE<-sum(sort(Yr)^3*dFDE)
	Q.DE<-c(Q1.DE, Q2.DE, Q3.DE)
	Q.IE<-c(Q1.IE, Q2.IE, Q3.IE)
	
	out<-list(FDE=FDE, FIE=FIE, y=sort(Yr), Q.DE=Q.DE, Q.IE=Q.IE, F00=F00, F10=F10, F11=F11)
	return(out)	
}

mediation.boot<-function(Yr, Mr, Zr, B){

	orig<-mediation.kern(Yr, Mr, Zr)
	
	FDE.b<-orig$FDE
	FIE.b<-orig$FIE
	QDE.b<-orig$Q.DE
	QIE.b<-orig$Q.IE
	n<-length(Yr)

	if (B>0){
		for (b in 1:B){
			set.seed(b*17)
			select<-sample(n, re=TRUE)
			Yboot<-Yr[select]
			Mboot<-Mr[select]
			Zboot<-Zr[select]
			
			boot<-mediation.kern(Yboot, Mboot, Zboot)
			
			boot.y.mat<-t(matrix(rep(boot$y, length(orig$y)), nr=length(boot$y)))
			ind<-apply(boot.y.mat<=orig$y, 1, sum)
			FDE.boot<-c(rep(0, sum(ind==0)), boot$FDE[ind])
			FIE.boot<-c(rep(0, sum(ind==0)), boot$FIE[ind])
			FDE.b<-rbind(FDE.b, FDE.boot)
			FIE.b<-rbind(FIE.b, FIE.boot)
			QDE.b<-cbind(QDE.b, boot$Q.DE)
			QIE.b<-cbind(QIE.b, boot$Q.IE)
			print(b); flush.console()
		}
	}
	out<-list(FDE.b=FDE.b, FIE.b=FIE.b, QDE.b=QDE.b, QIE.b=QIE.b, y=orig$y)
	return(out)
}
