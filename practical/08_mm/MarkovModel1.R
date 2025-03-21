#' Asthma cost-effectiveness analysis based on a Markov model
#' source: Price & Briggs (2002); BMHE (5.4)
#' 
#' Load all the relevant data
S <- 5			# number of clinical states
			# s=1: STW
			# s=2: UTW
			# s=3: Pex
			# s=4: Hex
			# s=5: TF

J <- 12 		# number of weeks in the follow up

#' Observed data on transitions among the states for the two treatments
r.0 <- (matrix(c(
66,32,0,0,2, 
42,752,0,5,20, 
0,4,0,1,0, 
0,0,0,0,0, 
0,0,0,0,156),c(S,S),byrow=TRUE))

r.1 <- (matrix(c(
210,60,0,1,1,
88,641,0,4,13,
1,0,0,0,1,
0,0,0,0,0,
0,0,0,0,81),c(S,S),byrow=TRUE))

#' Includes only the first three rows of data (Hex has no data and 
#' TP is modelled deterministically)
r.0 <- r.0[1:3,]
r.1 <- r.1[1:3,]

n.0 <- apply(r.0,1,sum)		# total number of patients in each state for t=0
n.1 <- apply(r.1,1,sum)		# total number of patients in each state for t=1

#' Parameters of the Dirichlets prior distributions
#' scale factor to use different version of minimally informative prior
#' Reducing this makes the prior less and less informative
scale <- 10
alpha.0 <- rep(scale,S)	
alpha.1 <- rep(scale,S)

#' Run the MCMC model
library(R2jags)
#' Defines the data and points to the model code
dataJags <- list("n.0","n.1","r.0","r.1","alpha.0","alpha.1","S")
filein <- here::here("08_mm","MarkovModel1.txt")
#'
#' Defines the model parameters and then the inital values
params <- c("lambda.0","lambda.1")
inits <- function(){
	temp.0 <- matrix(rgamma(4*S,scale,1),4,S)
  sum.temp.0 <- apply(temp.0,1,sum)
  mat.0 <- temp.0/sum.temp.0
  temp.1 <- matrix(rgamma(4*S,scale,1),4,S)
  sum.temp.1 <- apply(temp.1,1,sum)
  mat.1 <- temp.1/sum.temp.1
  list(lambda.0=rbind(mat.0,rep(NA,S)),lambda.1=rbind(mat.1,rep(NA,S)))
}
#' NB: lambda.0 and lambda.1 need: 4 random rows (even Hex has 
#' random parameters) 1 row of NAs (because TP is deterministic)

#' Number of simulations, burnin and thinning
n.iter <- 10000
n.burnin <- 5000
n.thin <- floor((n.iter-n.burnin)/500)

#' Finally runs the model using JAGS & prints the results
mm1 <- jags(
  data=dataJags, inits=inits, parameters.to.save=params, model.file=filein,
  n.chains=2, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC=TRUE
)
print(mm1,digits=3,intervals=c(0.025, 0.975))

attach.jags(mm1)
 
#' ## Alternatively, could use BUGS --- code basically identical
# library(R2OpenBUGS)
# dataBugs <- list("n.0","n.1","r.0","r.1","alpha.0","alpha.1","S")
# filein <- here::here("08_mm","MarkovModel1.txt")
# params <- c("lambda.0","lambda.1")
# # lambda.0 and lambda.1 need: 4 random rows (even Hex has random parameters)
# #	                      1 row of NAs (because TP is deterministic)
# inits <- function(){
# 	temp.0 <- matrix(rgamma(4*S,scale,1),4,S)
# 	sum.temp.0 <- apply(temp.0,1,sum)
# 	mat.0 <- temp.0/sum.temp.0
# 	temp.1 <- matrix(rgamma(4*S,scale,1),4,S)
# 	sum.temp.1 <- apply(temp.1,1,sum)
# 	mat.1 <- temp.1/sum.temp.1
# 	list(lambda.0=rbind(mat.0,rep(NA,S)),lambda.1=rbind(mat.1,rep(NA,S)))
# }
# 
# n.iter <- 10000
# n.burnin <- 5000
# n.thin <- floor((n.iter-n.burnin)/500)
# mm1 <- bugs(data=dataBugs,inits=inits,parameters.to.save=params,model.file=filein,
# 	n.chains=2, n.iter, n.burnin, n.thin, DIC=TRUE)
# print(mm1,digits=3)
# attach.bugs(mm1)


#' Now run the Markov model from R
#' NB: this is an analysis for 1 single patient, but can create a "virtual" 
#'     population of N individuals and allocate them across the states at 
#'     time j=0
start <- c(1,0,0,0,0)

#' Markov transitions
m.0 <- m.1 <- array(NA,c(n.sims,S,(J+1)))
for (s in 1:S){
	m.0[,s,1] <- start[s]
	m.1[,s,1] <- start[s]
}

#' When using JAGS, your output lambda.0 and lambda.1 will include the last 
#' row of the transition matrices for each simulation (so all the transitions
#' off the absorbing state, whose probabilities are 0 to go away and 1 to stay, 
#' with no variance whatsoever). So we can create the Markov model rather simply
#' as with the following code
for (i in 1:n.sims) {
  for (j in 2:(J+1)){
    for (s in 1:S){
      # Use the new matrices lam0 and lam1 to do the matrix multiplication
      m.0[i,s,j] <- sum(m.0[i,,j-1]*lambda.0[i,,s])
      m.1[i,s,j] <- sum(m.1[i,,j-1]*lambda.1[i,,s])
    }
  }
}

#' NB: Conversely, BUGS output the matrices lambda.0 and lambda.1 with 
#'     simulations for the "random" part, ie only the first 4 rows, as the 
#'     last one is deterministically defined as c(0,0,0,0,1) because once 
#'     a patient is in TF, they can't move away. So need to reconstruct a 
#'     full matrix with S rows and S columns for each MCMC simulations. 
#'     We do this by defining new arrays lam0 and lam1 and then stacking up the 
#'     simulated values for the first (S-1) rows saved in lambda.0[i,,] 
#'     and lambda.1[i,,] for MCMC simulations i with a row vector containing 
#'     (S-1) 0s and then a 1, ie c(0,0,0,0,1)
# lam0=lam1=array(NA,c(n.sims,S,S))
# for (i in 1:n.sims) {
#   lam0[i,,]=rbind(lambda.0[i,,],c(0,0,0,0,1))
#   lam1[i,,]=rbind(lambda.1[i,,],c(0,0,0,0,1))
# 	for (j in 2:(J+1)){
# 		for (s in 1:S){
# 		   # Use the new matrices lam0 and lam1 to do the matrix multiplication
# 			m.0[i,s,j] <- sum(m.0[i,,j-1]*lam0[i,,s])
# 			m.1[i,s,j] <- sum(m.1[i,,j-1]*lam1[i,,s])
# 		}
# 	}
# }



#' Barplot of the number of people in each state at each time point
#' First gets the aggregated data from the simulations and then 
#' wrangles the resulting dataset to a long format. This plot is often 
#' termed a "Markov trace"
library(tidyverse)
apply(m.0,c(2,3),sum) |> t() |> as_tibble() |>  
  setNames(c("STW","UTW","Pex","Hex","TF")) |> 
  mutate(Intervention="FP",time=seq(0,12)) |> 
  bind_rows(
    apply(m.1,c(2,3),sum) |> t() |> as_tibble() |>  
      setNames(c("STW","UTW","Pex","Hex","TF")) |> 
      mutate(Intervention="SFC",time=seq(0,12)) 
  ) |> 
  pivot_longer(!c(time,Intervention),names_to="state",values_to="npeople") |>
  #' Turns the variable 'state' into a factor, which we order according to
  #' the various states in our Markov model
  mutate(
    state=factor(state),
    state=forcats::fct_relevel(state,c("STW","UTW","Pex","Hex","TF"))
  ) |>
  #' Try and run up to here (removing the last pipe '|>' in the above row) to
  #' check what the resulting data look like...
  #' 
  #' Then uses ggplot to make the barplot & customise it...
  ggplot(aes(x=time,y=npeople,fill=as_factor(state))) + 
    geom_bar(position="stack",stat="identity") + 
  xlab("Virtual follow up") + ylab("Proportion of patients in each state") +
  theme_bw() + facet_wrap(~Intervention) +
  theme(legend.position = "bottom") + labs(fill="State") + scale_fill_grey()


#' Now runs the economic analysis
#' Defines the costs --- see BMHE page 193
unit.cost.0 <- c(2.38,2.38,1815.58,95.21)
unit.cost.1 <- c(7.96,7.96,1821.17,100.79)

cost0 <- cost1 <- matrix(NA,n.sims,J)
for (i in 1:n.sims) {
	for (j in 2:(J+1)) {
		cost0[i,j-1] <- m.0[i,S,j]*
		  (unit.cost.0%*%m.0[i,1:(S-1),j])/sum(m.0[i,1:(S-1),j]) + 
				unit.cost.0%*%m.0[i,1:(S-1),j]
		cost1[i,j-1] <- m.1[i,S,j]*
		  (unit.cost.0%*%m.0[i,1:(S-1),j])/sum(m.0[i,1:(S-1),j]) + 
				unit.cost.1%*%m.1[i,1:(S-1),j]
	}
}

#' General formulation to apply discount
delta.b <- 0.035	# discount rate for benefits (3.5%)
delta.c <- 0.035	# discount rate for costs (3.5%)

#' Defines the discount factors
disc.b <- numeric(); disc.c <- numeric()
disc.b[1] <- 1; disc.c[1] <- 1
for (j in 2:J) {
	disc.b[j] <- (1+delta.b)^(j-1)
	disc.c[j] <- (1+delta.c)^(j-1)
}
disc.cost0 <- disc.eff0 <- disc.cost1 <- disc.eff1 <- matrix(NA,n.sims,J)
for (j in 1:J) {
	disc.cost0[,j] <- cost0[,j]/disc.c[j]
	disc.cost1[,j] <- cost1[,j]/disc.c[j]
	disc.eff0[,j] <- m.0[,1,j]/disc.b[j]
	disc.eff1[,j] <- m.1[,1,j]/disc.b[j]
}

#' Sums the values across all time points and creates matrix of costs
c <- matrix(NA,n.sims,2)
c[,1] <- apply(cost0,1,sum)
c[,2] <- apply(cost1,1,sum)

#' Effectiveness
e <- matrix(NA,n.sims,2)
e[,1] <- apply(m.0[,1,],1,sum)
e[,2] <- apply(m.1[,1,],1,sum)


#' Cost-effectiveness analysis
library(BCEA)
ints <- c("FP","SFC")
m <- bcea(e,c,ref=2,interventions=ints,Kmax=300)

contour2(m,wtp = 300,graph="gg")
plot(m,graph="gg")
summary(m)
