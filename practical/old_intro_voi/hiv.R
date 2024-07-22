library(BCEA)

setwd(here::here("11_intro_voi"))

#Read in 150,000 simulated values of M, a, b, e, h
par<-read.table("hiv150.txt",header=TRUE)
N<-105000; T<-3; Nsim<-nrow(par); Nt<-2	#Nt=no. trts

#Net Benefit based on current information
NB<-matrix(rep(0,Nsim*Nt),Nsim,Nt)
NB[,2]<- N*(1-par$a-par$b)*(par$M*par$e*(1-par$h) - T*(1-par$e*par$h))

#Expected Net Benefit
ENB<-apply(NB,2,mean)

#Optimal intervention based on current information
tstar<- which.max(ENB)

#Prob tstar is cost-effective
CE<-ifelse(NB[,2]>NB[,1],1,0)
probCE<-mean(CE)

### EVPI ###

#Find maximum NB for each simulation minus the NB using tstar
max.NBgain<-apply(NB,1,max) - NB[,tstar]

#Compute EVPI
EVPI<- 7.7217*mean(max.NBgain)
EVPI

#Running mean to assess convergence
EVPI.run<-c(rep(0,150))
for (i in 1:150){
  EVPI.run[i]<-7.7217*mean(max.NBgain[1:(i*1000)])
}

#Convergence check plot
plot(seq(1000,150000,1000), EVPI.run, type="l",lty=1, xlab="Simulation", ylab="EVPI")
plot(seq(101000,150000,1000), EVPI.run[101:150], type="l",lty=1, xlab="Simulation", ylab="EVPI")


### EVPI in BCEA ###

#Using BCEA to compute EVPI
effects<-matrix(rep(0,Nsim*Nt),Nsim,Nt)
costs<-matrix(rep(0,Nsim*Nt),Nsim,Nt)

effects[,2]<-N*(1-par$a-par$b)*(par$M*par$e*(1-par$h))
costs[,2]<-N*(1-par$a-par$b)*(T*(1-par$e*par$h))

bcea.out <- bcea (eff=effects, cost=costs, ref=2, wtp=1, 
                  interventions=c("targetted", "universal"),plot=F)
summary(bcea.out)

EVPI.bcea<- bcea.out$evi*7.7217

EVPI.bcea


