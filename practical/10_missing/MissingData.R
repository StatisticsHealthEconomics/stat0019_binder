# 0. Sets up environment
library(BCEA)
library(R2OpenBUGS)

# Loads the data list
data=readRDS("missing_data.rds")
# Inspect the relevant data --- for example by producing a histogram of the costs in the control arm (t=1)
hist(data$c[[1]])

# Run the Bivariate Normal model
# Defines: 
# 1. model file
filein="Normal_Normal.txt"

# 2. data list
datalist=list(N1=data$n[[1]],eff1=data$e[[1]],cost1=data$c[[1]],u1=data$u[[1]],
          N2=data$n[[2]],eff2=data$e[[2]],cost2=data$c[[2]],u2=data$u[[2]])

# 3. parameters to monitor
params<-c("mu.e","sd.e","alpha0","alpha1","beta0","beta1","Delta_e","Delta_c","mu.c","sd.c","eff1","eff2","cost1","cost2")

# 4. number of iterations
n.iter<-10000

# 5. sets up initial values for crucial parameters
inits=function(){
   list(alpha0=rnorm(2),alpha1=rnorm(2),beta0=rnorm(2),beta1=rnorm(2),mu.u=runif(2))
}

# 6. runs the model
NN<-bugs(data=datalist,inits=inits,parameters.to.save=params,model.file=filein,n.chains=2,n.iter=n.iter,n.thin=1,DIC=TRUE)

# Now runs the CEA using BCEA
e<-cbind(NN$sims.list$mu.e[,1],NN$sims.list$mu.e[,2])
c<-cbind(NN$sims.list$mu.c[,1],NN$sims.list$mu.c[,2])
CEA_NN<-bcea(e=e,c=c,ref = 2)
ceplane.plot(CEA_NN)
ceac.plot(CEA_NN)

# Checks the imputed values (can try other individuals too...)
hist(NN$sims.list$eff1[,1],xlab="QALYs for the first missing individual in t=1",main="")
abline(v=1,lwd=2)
sum(NN$sims.list$eff1[,1]>1)/NN$n.sims

# Run the Beta/Gamma model
# NB: The sampler in OpenBUGS cannot make this model/data work. The initial values are difficult to set and the model can't initialize.
#     On the other hand, JAGS can compile and run this very quickly. So to run the Beta/Gamma model, you'll have to install JAGS
#     https://sourceforge.net/projects/mcmc-jags/files/
#     and then install R2jags
#     install.packages("R2jags")
#     The syntax of the model and the R script is effectively identical

# 1. model file
filein="Beta_Gamma.txt"

# 2. data list
datalist=list(N1=data$n[[1]],eff1=data$e[[1]],cost1=data$c[[1]],u1=data$u[[1]],
              N2=data$n[[2]],eff2=data$e[[2]],cost2=data$c[[2]],u2=data$u[[2]])

# 3. parameters to monitor
params<-c("mu.e","sd.e","alpha0","alpha1","beta0","beta1","Delta_e","Delta_c","mu.c","sd.c","eff1","eff2","cost1","cost2")
# 4. number of iterations
n.iter<-10000

# 5. runs the model
library(R2jags)
BG <- jags(data=datalist,inits=inits,parameters.to.save=params,model.file=filein,n.chains=2,n.iter=n.iter,n.thin=1,DIC=TRUE)
# Can show the output --- the print function in R2jags has some nice extra features
print(BG,interval=c(.025,.975),digits=3)

# Can do the CEA again
# NB: The JAGS object has an extra "layer". You can access the simulation list in the object $BUGSoutput$sims.list
#     so you have to modify the previous code to add `$BUGSoutput`
e<-cbind(BG$BUGSoutput$sims.list$mu.e[,1],BG$BUGSoutput$sims.list$mu.e[,2])
c<-cbind(BG$BUGSoutput$sims.list$mu.c[,1],BG$BUGSoutput$sims.list$mu.c[,2])
CEA_BG<-bcea(e=e,c=c,ref = 2)
ceplane.plot(CEA_BG)
ceac.plot(CEA_BG)
eib.plot(CEA_BG)


### This is an attempt at setting initial values to run the model in OpenBUGS
### Unfortunately, they fail --- probably to do with the sampler selected by OpenBUGS to run this model
### and the specific dataset at hand.

# 5. sets up initial values for crucial parameters
inits=function(){
   cost1=datalist$cost1
   cost1[is.na(datalist$cost1)]=median(cost1,na.rm=T); cost1[!is.na(datalist$cost1)]=NA
   cost2=datalist$cost2
   cost2[is.na(datalist$cost2)]=median(cost2,na.rm=T); cost2[!is.na(datalist$cost2)]=NA
   eff1=datalist$eff1
   eff1[is.na(datalist$eff1)]=mean(eff1,na.rm=T); eff1[!is.na(datalist$eff1)]=NA
   eff2=datalist$eff2
   eff2[is.na(datalist$eff2)]=mean(eff2,na.rm=T); eff2[!is.na(datalist$eff2)]=NA
   u1=datalist$u1; u2=datalist$u2
   u1[is.na(datalist$u1)]=mean(u1,na.rm=T); u1[!is.na(datalist$u1)]=NA
   u2[is.na(datalist$u2)]=mean(u2,na.rm=T); u2[!is.na(datalist$u2)]=NA
   sd.e=sd.u=mu.e=limit.e=limit.u=numeric()
   alpha0=runif(2,1,2)
   alpha1=runif(2,1,4)
   beta0=runif(2,0,2)
   beta1=runif(2,0,3)
   mu.u=runif(2,.8,1)
   mu.e=exp(alpha0)/(1+exp(alpha0))
   for(t in 1:2) {
      limit.e[t]=sqrt(mu.e[t]*(1-mu.e[t]))
      limit.u[t]=sqrt(mu.u[t]*(1-mu.u[t]))
      sd.e[t]=mean(c(0,limit.e[t]))
      sd.u[t]=mean(c(0,limit.u[t]))
   }
   sd.c=runif(2,10,30)
   list(alpha0=alpha0,alpha1=alpha1,beta0=beta0,beta1=beta1,mu.u=mu.u,sd.e=sd.e,sd.u=sd.u,
        cost1=cost1,cost2=cost2,eff1=eff1,eff2=eff2,u1=u1,u2=u2,sd.c=sd.c
   )
}
# 6. runs the model
BG<-bugs(data=datalist,inits=inits,parameters.to.save=params,model.file=filein,n.chains=2,n.iter=n.iter,n.thin=1,DIC=TRUE)

