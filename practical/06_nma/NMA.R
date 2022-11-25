library(R2OpenBUGS)
library(BCEA)

### Smoking cessation network meta-analysis data in format obtained
### from Lu & Ades tutorial "Introduction to Mixed Treatment Comparisons"
load(here::here("practical","06_nma","smoke.Rdata"))

### Explores the data object just loaded
names(smoke.list)
lapply(smoke.list,class)

#############    FIT MODELS IN OPENBUGS via R2OpenBUGS   #############
### Initial values
inits <- list(list(mu=rep(0,24), d=c(NA,0,0,0)),
              list(mu=rep(-1,24), d=c(NA,1,1,1)))

### FIXED EFFECTS MODEL

### Pilot run with no burn-in, to illustrate convergence using traceplots
res <- bugs(
  model=here::here("practical","06_nma","smokefix_model.txt"), 
  data=smoke.list, inits=inits,
  parameters=c("d"),
  n.chains=2, n.burnin=0, n.iter=10000
)
# Traceplots
par(mfrow=c(2,2))
plot(res$sims.array[,1,1],t="l",xlab="Iterations",ylab="d[2]",col="blue")
points(res$sims.array[,2,1],t="l",col="red")
plot(res$sims.array[,1,2],t="l",xlab="Iterations",ylab="d[3]",col="blue")
points(res$sims.array[,2,2],t="l",col="red")
plot(res$sims.array[,1,3],t="l",xlab="Iterations",ylab="d[4]",col="blue")
points(res$sims.array[,2,3],t="l",col="red")

# Can also use 'bmhe::traceplot' to automate this, for instance
bmhe::traceplot(res,"d")

# In Rstudio do par(mfrow=c(1,1)) to "reset" the graphical interface
## It's certainly converged by 1000, so no need for any more fancy diagnostics
res <- bugs(
  model=here::here("practical","06_nma","smokefix_model.txt"), 
  data=smoke.list, inits=inits,
  parameters=c("d","L","pq"),
  n.chains=2, n.burnin=1000, n.iter=5000
)
res

## How many further iterations do we need after convergence?  Consider
## the "effective sample size".  Rule of thumb is you need n.eff>4000
## to ensure 95\% credible limits have 94.5-95.5% coverage (true in
## this case).  

### RANDOM EFFECTS MODEL.  Check convergence of random effects SD in particular
inits <- list(list(mu=rep(0,24), d=c(NA,0,0,0), sd=1),
              list(mu=rep(-1,24), d=c(NA,1,1,1), sd=2))

res2 <- bugs(
  model=here::here("practical","06_nma","smokere_model.txt"), 
  data=smoke.list, inits=inits,
  parameters=c("or", "d", "sd", "pq", "L"),
  n.chains=2, n.burnin=1000, n.iter=20000
)
print(res2,digits=3)

### Cost-effectiveness analysis
unit.cost <- c(0,200,6000,600)
ints <- c("No contact","Self help","Individual counselling","Group counselling")
e <- c <- matrix(NA,res2$n.sims,4)
L <- res2$sims.list$L # MCMC sample from distribution of life-years gained by quitting
pq <- res2$sims.list$pq # ...and from distributions of probability of quitting for each of 4 interventions

for (t in 1:4) {
    e[,t] <- L*pq[,t]
    c[,t] <- unit.cost[t]
}
colnames(e) <- colnames(c) <- ints
round(apply(e, 2, quantile, c(0.025, 0.5, 0.975)), 1) # results on slide
m <- bcea(e,c,interventions=ints,Kmax=1000,ref=4)
summary(m)
plot(m)


### Appendix: Original pre-processing of dataset to format the data in a suitable way
smoke <- read.table(
  here::here("practical","06_nma","smoke_data_orig.txt"), 
  header=TRUE, nrow=24
)
names(smoke) <- gsub("\\.", "", names(smoke))
ns <- nrow(smoke)
nt <- 4
r <- smoke[,c("r1","r2","r3")]
n <- smoke[,c("n1","n2","n3")]
t <- smoke[,c("t1","t2","t3")]
rc <- nc <- matrix(nrow=ns, ncol=nt)
for (i in 1:2) {   # ith col should be treatment i
  rc[cbind(1:ns, t[,i])] <- r[,i]
  nc[cbind(1:ns, t[,i])] <- n[,i]
}
rc[cbind(1:2, t[1:2,3])] <- r[1:2,3]
nc[cbind(1:2, t[1:2,3])] <- n[1:2,3]
dc <- matrix(paste(rc, nc, sep="/"), nrow=ns, ncol=nt)
dc[dc=="NA/NA"] <- ""
comp <- paste(t[,1], t[,2], ifelse(is.na(t[,3]),"",t[,3]), sep="")
tnames <- c("A: None","B: Self-help","C: Individual","D: Group")

## Order by comparison
rc <- rc[order(comp),]
nc <- nc[order(comp),]
t <- as.matrix(t)[order(comp),]
na <- smoke$na[order(comp)]

### Format data for WinBUGS or OpenBUGS
## r[s,t] should be data for tth treatment in study s
smoke.list <- list(r=rc, n=nc, t=t, na=na, NS=ns, NT=nt)