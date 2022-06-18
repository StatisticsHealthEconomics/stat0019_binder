setwd("C:/Users/epnjw/nickyw/ConDuCT-II/VOICourseFlorence/Intro EVPI")
library(BCEA)

#Read in 150,000 simulated values of Qwi, Qnwi, Cwi, Cnwi, Cdrug, phi, psi
par<-read.table("wounds.txt",header=TRUE)
Nsim<-nrow(par); Nt<-2; wtp=20000	#Nt=no. trts

#Net Benefit based on current information
effects<-matrix(rep(0,Nsim*Nt),Nsim,Nt)
costs<- matrix(rep(0,Nsim*Nt),Nsim,Nt)
p1<-exp(par$psi)/(1+exp(par$psi))
p2<-exp(par$psi+par$phi)/(1+exp(par$psi+par$phi))

effects[,1]<- p1*par$Qwi+(1-p1)*par$Qnwi
effects[,2]<- p2*par$Qwi+(1-p2)*par$Qnwi

costs[,1]<- p1*par$Cwi+(1-p1)*par$Cnwi
costs[,2]<- p2*par$Cwi+(1-p2)*par$Cnwi + par$Cdrug

### EVPI in BCEA ###

bcea.out <- bcea (e=effects, c=costs, ref=2, wtp=seq(1000,40000, 1000), 
interventions=c("no antibiotics", "antibiotics"),plot=F)
summary(bcea.out)

bcea.out$eib[20]
bcea.out$ceac[20]
bcea.out$evi[20]

ceac.plot(bcea.out)
evi.plot(bcea.out)

PopEVPI.bcea<- bcea.out$evi*7.7217*166081

PopEVPI.bcea[20]

### EVPI without BCEA###

NB<-effects*wtp-costs

#Expected Net Benefit
ENB<-apply(NB,2,mean)

#Optimal intervention based on current information
tstar<- which.max(ENB)

#Prob tstar is cost-effective
CE<-ifelse(NB[,2]>NB[,1],1,0)
probCE<-mean(CE)

#Find maximum NB gain for each simulation
max.NBgain<-apply(NB,1,max) - NB[,tstar]

#Compute EVPI (per 166081 women per year for 10 year technology horizon, discounted)
PopEVPI<- 166081*7.7217*mean(max.NBgain)

EVPI
ENB
tstar
probCE

#Running mean to assess convergence
PopEVPI.run<-c(rep(0,150))
for (i in 1:150){
	PopEVPI.run[i]<-166081*7.7217*mean(max.NBgain[1:(i*1000)])
}

#Convergence check plot
plot(seq(1000,150000,1000), PopEVPI.run, type="l",lty=1, xlab="Simulation", ylab="EVPI")
plot(seq(101000,150000,1000), PopEVPI.run[101:150], type="l",lty=1, xlab="Simulation", ylab="EVPI")




