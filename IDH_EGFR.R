source("mediation.R")

load(file="df_GBMLGG_IDH1_EGFR.rda")

dat<-as.data.frame(df_GBMLGG_IDH1_EGFR)

M<-as.numeric(dat[, "EGFR_cg18071865"])
Z<-as.numeric(dat[, "IDH1"])
Y<-as.numeric(dat[, "EGFR_RPPA"])

##

par(mfrow=c(2,1))
hist(Y[Z==0], xlim=c(min(Y), max(Y)), br="Scott")
hist(Y[Z==1], xlim=c(min(Y), max(Y)), br="Scott")

par(mfrow=c(2,1))
Mstar<-Mr
hist(Mstar[Z==0], xlim=c(min(Mstar), max(Mstar)), br="Scott")
hist(Mstar[Z==1], xlim=c(min(Mstar), max(Mstar)), br="Scott")


##

Zr<-residuals(lm(Z~years_to_birth+years_to_birth_squared+histological_type+gender+radiation_therapy+race, data=dat))+mean(Z)
Mr<-residuals(lm(M~years_to_birth+years_to_birth_squared+histological_type+gender+radiation_therapy+race, data=dat))
Yr<-residuals(lm(Y~years_to_birth+years_to_birth_squared+histological_type+gender+radiation_therapy+race, data=dat))

med<-mediation.boot(Yr, Mr, Zr, B=500)
de<-med$QDE.b
ie<-med$QIE.b
rbind(de[,1], apply(de, 1, function(x) quantile(x, probs=c(0.025, 0.975))))
rbind(ie[,1], apply(ie, 1, function(x) quantile(x, probs=c(0.025, 0.975))))

ci<-apply(med$FIE.b[-1,], 2, function(x) quantile(x, probs=c(0.025, 0.975)))

med0<-mediation.kern(Yr, Mr, Zr)
FF<-c(med0$F00, med0$F10, med0$F11)

##
f10<-NULL; f11<-NULL
gp<-20
yy<-seq(min(med0$y), max(med0$y), by=(max(med0$y)-min(med0$y))/gp)
yy1<-yy
yy2<-yy+0.5*(max(med0$y)-min(med0$y))/gp
for (i in 1:(gp-1)) f10<-c(f10, max(med0$F10[med0$y<=yy1[i+1]])-max(med0$F10[med0$y<=yy1[i]]))
for (i in 1:(gp-1)) f11<-c(f11, max(med0$F11[med0$y<=yy2[i+1]])-max(med0$F11[med0$y<=yy2[i]]))
ff<-c(f10, f11)
f102<-rep(f10, each=2)
f112<-rep(f11, each=2)
yy11<-c(yy1[1], rep(yy1[2:(gp-1)], each=2), yy1[gp])
yy22<-c(yy2[1], rep(yy2[2:(gp-1)], each=2), yy2[gp])

#png(file="TCGA_GBM_IDH_EGFR_dF.png", pointsize=16, width=480, height=480)
plot(y=f10, x=yy1[2:gp]+mean(Y), type="S", ylim=c(min(ff), max(ff)), col=rgb(0.9, 0.5, 0.1),
	main="", bty="n", 
	ylab=NA, xlab="Protein expression of EGFR")
legend(x=0.25, y=0.30, fill=c(rgb(0.9, 0.1, 0.1, 0.2), rgb(0.6, 0.1, 0.9, 0.2)), 
	border=c("red", "purple"),
	legend=c(expression(italic(d*hat(F)[Y(list(1,M(0)))])), 
			expression(italic(d*hat(F)[Y(list(1,M(1)))]))
		),
	bty="n")
polygon(y=c(0, f102, 0, 0), x=c(yy11[1], yy11, yy11[(gp-1)*2], yy11[1])+mean(Y), col=rgb(0.9, 0.1, 0.1, 0.2), border=2)
polygon(y=c(0, f112, 0, 0), x=c(yy22[1], yy22, yy22[(gp-1)*2], yy22[1])+mean(Y), col=rgb(0.6, 0.1, 0.9, 0.2), border="purple")

#dev.off()

##

#png(file="TCGA_GBM_IDH_EGFR.png", pointsize=16, width=480, height=600)

layout(mat=matrix(c(1,2), nc=1), heights=c(2,3))
par(mar=c(3,5,1,1))
plot(y=med$FIE.b[1,], x=med$y+mean(Y), type="l", bty="n", ylim=c(-0.2, 0.2), 
	ylab=expression(italic(hat(F)[IE])
	==italic(hat(F)[Y(list(1,M(1)))])
	-italic(hat(F)[Y(list(1,M(0)))])))
points(y=ci[1,], x=med$y+mean(Y), type="l", col=rgb(0.5, 0.5, 0.5))
points(y=ci[2,], x=med$y+mean(Y), type="l", col=rgb(0.5, 0.5, 0.5))
abline(h=0, col=1, lty=2)

par(mar=c(5,5,1,1))
plot(y=med0$F00, x=med0$y+mean(Y), type="l", ylim=c(min(FF), max(FF)), col=rgb(0.9, 0.5, 0.1),
	main="", bty="n", 
	ylab=NA, xlab="Protein expression of EGFR")
legend(x=0.25, y=0.45, col=c(rgb(0.6, 0.1, 0.9), rgb(0.9, 0.1, 0.1), rgb(0.9, 0.5, 0.1)), lty=1, 
	legend=c(expression(italic(hat(F)[Y(list(1,M(1)))])), 
			expression(italic(hat(F)[Y(list(1,M(0)))])),
			expression(italic(hat(F)[Y(list(0,M(0)))]))
		),
	bty="n")
points(y=med0$F10, x=med0$y+mean(Y), type="l", col=rgb(0.9, 0.1, 0.1))
points(y=med0$F11, x=med0$y+mean(Y), type="l", col=rgb(0.6, 0.1, 0.9))

dev.off()

########################################################################

#### Parametric

fit.m<-lm(M~Z+years_to_birth+years_to_birth_squared+histological_type+gender+radiation_therapy+race, data=dat)
fit.y<-lm(Y~M+Z+years_to_birth+years_to_birth_squared+histological_type+gender+radiation_therapy+race, data=dat)	
a<-fit.m$coef["Z"]
b<-fit.y$coef["M"]
c<-fit.y$coef["Z"]
sd.a<-summary(fit.m)$coef["Z", "Std. Error"]
sd.b<-summary(fit.y)$coef["M", "Std. Error"]
sd.c<-summary(fit.y)$coef["Z", "Std. Error"]
ie.par<-a*b
de.par<-c
v.IE<-a^2*sd.b^2+b^2*sd.a^2

low.qde<-de.par-1.96*sd.c
upp.qde<-de.par+1.96*sd.c

low.qie<-ie.par-1.96*sqrt(v.IE)
upp.qie<-ie.par+1.96*sqrt(v.IE)

c(de.par, low.qde, upp.qde)
c(ie.par, low.qie, upp.qie)

