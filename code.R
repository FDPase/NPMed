adjust<-1

for (tt in 1:4){

for (deievv in 0:3){

if (deievv==0) {de=0; ie=0; vvv=0}
if (deievv==1) {de=1; ie=0; vvv=0}
if (deievv==2) {de=0; ie=1; vvv=0}
if (deievv==3) {de=0; ie=0; vvv=1}

source("DatGen.R")
source("mediation.R")

S<-1000
FF<-NULL; QQ<-NULL
for (s in 1:S){

	dat<-GenDat(seed=37*s)
	X<-dat$X
	Z<-dat$Z
	M<-dat$M
	Y<-dat$Y
	
	if (adjust==1){ 
		Zr<-residuals(lm(Z~X))+mean(Z)
		Mr<-residuals(lm(M~X))+mean(M)
		Yr<-residuals(lm(Y~X))+mean(Y)
	}
	if (adjust==0) Zr<-Z; Mr<-M; Yr<-Y
	
	med<-mediation.boot(Yr, Mr, Zr, B=0)
	
	FF<-rbind(FF, med$y, med$FDE.b, med$FIE.b)
	QQ<-rbind(QQ, c(med$QDE.b, med$QIE.b))
	print(s); flush.console()
}

if (adjust==1) save(FF, QQ, file=paste("DE", de, "IE", ie, "V", vvv, "T", tt, ".RData", sep=""))
if (adjust==0) save(FF, QQ, file=paste("DE", de, "IE", ie, "V", vvv, "T", tt, "_unadj.RData", sep=""))

}}

##################################################################

## Theoretical values
setwd("C:/Users/ythuang/Desktop/nonparametric_mediation/code")

for (tt in 1:4){
for (deievv in 0:3){

if (deievv==0) {de=0; ie=0; vvv=0}
if (deievv==1) {de=1; ie=0; vvv=0}
if (deievv==2) {de=0; ie=1; vvv=0}
if (deievv==3) {de=0; ie=0; vvv=1}

source("DatGen.R")

z1<-1
z0<-0
vv<-vvv

m<-1000
if (vv==1){
	rZ0<-rep(c(rep(0, 400), rep(1, 100)), m)
	rZ1<-rep(c(rep(0, 100), rep(1, 400)), m)
} else {
	rZ0<-rep(c(rep(0, 350), rep(1, 150)), m)
	rZ1<-rep(c(rep(0, 350), rep(1, 150)), m)
}
M0<-genM(Z=0, rZ=rZ0, n=ifelse(tt<=2, m*50, m*500), the=1)
M1<-genM(Z=1, rZ=rZ1, n=ifelse(tt<=2, m*50, m*500), the=1)
Y00<-genY(Z=0, rZ=rZ0, M=M0, n=ifelse(tt<=2, m*m*50, m*m*500), the=1)
Y10<-genY(Z=1, rZ=rZ1, M=M0, n=ifelse(tt<=2, m*m*50, m*m*500), the=1)
Y11<-genY(Z=1, rZ=rZ1, M=M1, n=ifelse(tt<=2, m*m*50, m*m*500), the=1)

F00<-ecdf(Y00)
F10<-ecdf(Y10)
F11<-ecdf(Y11)

yy<-seq(-20, 20, 0.01)
FF.t<-rbind(yy, 
	F10(yy)-F00(yy),
	F11(yy)-F10(yy))
QQ.t<-c(mean(Y10-Y00), mean(Y10^2-Y00^2), mean(Y10^3-Y00^3), mean(Y11-Y10), mean(Y11^2-Y10^2), mean(Y11^3-Y10^3))

save(FF.t, QQ.t, file=paste("True_DE", de, "IE", ie, "V", vvv, "T", tt, ".RData", sep=""))

}}


