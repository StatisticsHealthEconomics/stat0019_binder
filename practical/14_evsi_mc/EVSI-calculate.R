##########################################################
### Model Code for model published in Brennan and     ####
### Kharroubi and modified by Menzies included in the ####
### paper "Estimating the Expected Value of Sample    #### 
### Information across Different Sample Sizes using   ####
### Moment Matching and Non-Linear Regression"        ####
###                                                   ####
### Modified for use within the EVSI package.         ####
### Copyright: Anna Heath 2018                        ####
##########################################################

####Required Pacakges
install.packages("devtools")
library(devtools)
install_github("annaheath/EVSI")
library(R2OpenBUGS)
library(BCEA)
library(EVSI)

####Model Code for BUGS
Model<-function(){
  #PSA distributions for the parameters
  #All parameters are normal distributions
  #Parameterised in terms of mean and precision
  theta.1~dnorm(10000,0.01)
  theta.11~dnorm(15000,0.01)
  theta.2~dnorm(0.1,2500)
  theta.12~dnorm(.08,2500)
  theta.3~dnorm(5.2,1)
  theta.13~dnorm(6.1,1)
  theta.4~dnorm(4000,2.5E-07)
  theta.8~dnorm(.25,100)
  theta.17~dnorm(.20,400)
  theta.9~dnorm(-.1,2500)
  theta.18~dnorm(-.1,2500)
  theta.10~dnorm(.5,25)
  theta.19~dnorm(.5,25)
  
  #Correlated PSA distributions requires multivariate normal distribution
  #Parameters of the multivariate normal distribution must be specified as data.
  #Parameterised in terms of a mean vector and a precision matrix
  cor.theta.1[1:4]~dmnorm(Mu.1[],Tau.1[,])
  theta.5<-cor.theta.1[1]
  theta.7<-cor.theta.1[3]
  theta.14<-cor.theta.1[2]
  theta.16<-cor.theta.1[4]
  
  cor.theta.2[1:2]~dmnorm(Mu.2[],Tau.2[,])
  theta.6<-cor.theta.2[1]
  theta.15<-cor.theta.2[2]
  
  ##Sampling distribution for the data
  #Sampling from a multivariate normal distribution
  #Within EVSI package - cannot have multivariate data collection
  #Decompose multivariate normal into two scalar data points with conditional mean and variance
  #Decomposition found in "Bayesian Methods in Health Economics" Baio (2012), Page 156
  tau.14<-tau.X.14/(1-rho.X*rho.X)
  for(i in 1:N){
    X.5[i]~dnorm(theta.5,tau.X.5)
    mu.14[i]<-theta.14+sqrt(tau.X.5/tau.X.14)*rho.X*(X.5[i]-theta.5)
    X.14[i]~dnorm(mu.14[i],tau.14)
  }
}

#This writes the model we just defined as a text file for use within BUGS 
filein.model <- "Model_File.txt"
write.model(Model,filein.model)

####Define the parameters of the multivariate normal priors
#Within BUGS, we define a mean vector and a precision matrix.
#Mean Vectors
Mu.1<-c(.7,.8,3,3)

Mu.2<-c(.3,.3)

#Precision Matrices
#The original model worked with standard deviations so first
#specify all the standard deviations
sd.5<-0.1
sd.6<-0.1
sd.7<-0.5
sd.14<-0.1
sd.15<-0.05
sd.16<-1
rho<-0.6

#Start by defining variance matrix using the matrix command
S.1<-matrix(
  #Vector containing all the values to enter into the matrix
  c(sd.5^2,rho*sd.14*sd.5,rho*sd.5*sd.7,rho*sd.5*sd.16,
    rho*sd.5*sd.14,sd.14^2,rho*sd.14*sd.7,rho*sd.14*sd.16,
    rho*sd.5*sd.7,sd.14*rho*sd.7,sd.7^2,rho*sd.16*sd.7,
    rho*sd.5*sd.16,sd.14*rho*sd.16,rho*sd.7*sd.16,sd.16^2),
  #Specify the size of the matrix
  nrow=4)


S.2<-matrix(
  #Vector containing all the values to enter into the matrix
  c(sd.6^2,rho*sd.6*sd.15,
    rho*sd.6*sd.15,sd.15^2),
  #Specify the size of the matrix
  nrow=2)

#Precision matrix is the inverse of variance matrix
#solve finds the matrix inverse.
Tau.1<-solve(S.1)
Tau.2<-solve(S.2)

###Define the parameters for the sampling distribution of data.
#This calculates the precision from the standard deviation
tau.X.5<-1/0.2^2
tau.X.14<-1/0.2^2
rho.X<-0.6

###Define Health Economic model
#Functions to calculate the health economic outcomes from the model parameters
#The R syntax to create a function uses the command function()
#Within the brackets you write the arguments of the new function, these are
#the values you give to the function
#Inside the curly {} brackets you write what you want the function to do.
#Here we calculate the effects measure for the treatments from the parameters
#The return() function tells R what it should output - in this case a vector of the
#effects across the two treatments.

#These functions must take ONE value for each of the model parameters and output
#either the cost or effects for each treatment as a vector
effects<-function(theta.5,theta.6,theta.7,theta.8,theta.9,theta.10,
                  theta.14,theta.15,theta.16,theta.17,theta.18,theta.19){
  #This is a decision tree model so the effects are in the sum/product form.
  e.treatment1<-(theta.5*theta.6*theta.7)+(theta.8*theta.9*theta.10)
  e.treatment2<-(theta.14*theta.15*theta.16)+(theta.17*theta.18*theta.19)
    e<-c(e.treatment1,e.treatment2)
  return(e)
}

costs<-function(theta.1,theta.2,theta.3,theta.4,
                theta.11,theta.12,theta.13){
  #This is a decision tree model so the costs are in the sum/product form.
  c.treatment1<-(theta.1+theta.2*theta.3*theta.4)
  c.treatment2<-(theta.11+theta.12*theta.13*theta.4)
  c<-c(c.treatment1,c.treatment2)
  return(c)
}

#EVSI
#Specify a list of data that needs to be used by the BUGS model
#In this example, we need the mean vectors and precision matrices 
#for the multivaraite normal priors and the precision for the 
#sampling distribution of the data.
data.evsi<-list(Mu.1=Mu.1,
               Mu.2=Mu.2,
               Tau.1=Tau.1,
               Tau.2=Tau.2,
               tau.X.5=tau.X.5,
               tau.X.14=tau.X.14,
               rho.X=rho.X)

#This function uses BCEA to perform the standard PSA analysis and to calculate the EVPPI
#Then it calculates the variance of the posterior distribution across different sample sizes
#Finally, it fits the non-linear model to calculate the EVSI across sample size.
g <- mm.post.var(
  #The place where we saved the BUGS model.
  filein.model,
  #The variable names for the data that we will collect in the future trial
  c("X.5","X.14"),
  #The variable name of the sample size (i.e. how many data points we will collect)
  "N",
  #The smallest and largest sample sizes for which we calculate the EVSI
  c(5,200),
  #The function that calculates the effects from the parameters
  effects,
  #The function that calculates the costs from the parameters
  costs,
  #The variable names of the parameters that are updated by the future trial
  #This is used to calculate the EVPPI for the focal parameters
  parameters=c("theta.5","theta.14"),
  #Any additional data that is needed for the BUGS model to run
  #In this case, the mean vectors and precision matrices
  #If you had data from a previous trial it was be included here.
  data.stats=data.evsi,
  #The Bayesian updating should be done with BUGS (rather than JAGS).
  update="bugs",
  #The number of iterations taken from each posterior distribution to
  #calculate the posterior variance.
  n.iter=1000,
  #The number of "future" posterior samples that should be taken
  Q=50
)

#Using the previous function, evsi.calc calculates the EVSI across different
#sample sizes and willingness-to-pay parameters by rescaling the fitted values
#from the EVPPI calculations.
evsi.BK<-evsi.calc(g)

#Launches the WebApp so you can see the results of the EVSI calculation
launch.App(evsi.BK)


####################################################  
#### Section 2                                  ####
####################################################


###Simulate from the BUGS model to perform PSA
#This means that you can have more PSA simulations than used for the 
#EVSI calculations
#You can also check the model fitting of the EVPPI or use non-standard
#calculation methods
#You can also control the PSA analysis - here we take the willingness-to-pay
#up to 200000 as this model is from an American persepective.

#Create a list of prior parameters for BUGS
prior.spec<-list(Mu.1=Mu.1,
                 Mu.2=Mu.2,
                 Tau.1=Tau.1,
                 Tau.2=Tau.2,
                 tau.X.5=tau.X.5,
                 tau.X.14=tau.X.14,
                 rho.X=rho.X,
                 #Here you need to specify N as a parameter of the model
                 #but we set it to 0 as we don't have any data.
                 N=1)

#Parameters to control the MCMC for the PSA
n.chains<-3        #Number of chains for MCMC
n.burnin <- 1000   #Number of burn in iterations
n.thin<-20         #Thinning
n.iter <- ceiling(10000/n.chains + n.burnin) #Total number of iterations for each chain
#Using these MCMC controls gives 10000 PSA simulations

#We need to monitor all the model parameters so we can perform the health economic analysis
parameters.to.save <- c("theta.1","theta.2","theta.3","theta.4","theta.5",
                        "theta.6","theta.7","theta.8","theta.9","theta.10",
                        "theta.11","theta.12","theta.13","theta.14","theta.15",
                        "theta.16","theta.17","theta.18","theta.19")

#Use BUGS to sample from the PSA distributions
model.output <- bugs(
  #Set the data as our prior parameters
  data =  prior.spec,
  inits=NULL,
  #monitor all model parameters
  parameters.to.save = parameters.to.save,
  #use the model file created earlier
  model.file = filein.model,
  #control the MCMC
  n.chains = n.chains, 
  n.iter = n.iter, 
  n.thin = n.thin, 
  n.burnin = n.burnin,
  #do not DIC
  DIC=F) 

###Perform PSA using the parameter simulations and the model functions
S<-length(model.output$sims.list$theta.1)
e<-array(dim=c(S,2))
#The function take 1 value for each of the model parameters and gives one
#value of the effects for ecah treatment.

#Here we loop through all the PSA simulations to calculate the effects for 
#each simulation.
for(i in 1:S){
  e[i,]<-effects(model.output$sims.list$theta.5[i],
                 model.output$sims.list$theta.6[i],
                 model.output$sims.list$theta.7[i],
                 model.output$sims.list$theta.8[i],
                 model.output$sims.list$theta.9[i],
                 model.output$sims.list$theta.10[i],
                 model.output$sims.list$theta.14[i],
                 model.output$sims.list$theta.15[i],
                 model.output$sims.list$theta.16[i],
                 model.output$sims.list$theta.17[i],
                 model.output$sims.list$theta.18[i],
                 model.output$sims.list$theta.19[i])
}


c<-array(dim=c(S,2))
for(i in 1:S){
  c[i,]<-costs(model.output$sims.list$theta.1[i],
               model.output$sims.list$theta.2[i],
               model.output$sims.list$theta.3[i],
               model.output$sims.list$theta.4[i],
               model.output$sims.list$theta.11[i],
               model.output$sims.list$theta.12[i],
               model.output$sims.list$theta.13[i])
}

###Present PSA results using the BCEA package
#Here we change the maximum willingness-to-pay
he<-bcea(e,c,Kmax=200000,ref=2,interventions=c("SoC","Treatment"))
#Here we plot the cost-effectiveness results for a threshold WTP of 100000
plot(he,wtp=100000)


###Obtain Value of Information results
#EVPI
evi.plot(he)

#EVPPI - using evppi
#To start we need the parameter simulations to be in the right form
#They are saved in model.output$sims.matrix but are not named
#For example:
model.output$sims.matrix$theta.1
#gives an error.

#This command formally names each column using the column names
input.matrix<-as.data.frame(model.output$sims.matrix)
names(input.matrix)<-colnames(model.output$sims.matrix)
#Try
input.matrix$theta.1
#Now we can call all parameters by their name and use these names in the EVPPI function
evi<-evppi(c("theta.5","theta.14"),input.matrix,he,method="gam")
plot(evi)
evi$evppi[which(he$k==100000)]

#Previously we have a warning that the data cannot be done within the EVSI package
#We now generate the data externally to the mm.post.var function.
#gen.quantiles selects the quantiles and sample sizes that can be used to generate the
#future data
theta.gen <- gen.quantiles(c("theta.5","theta.14"),input.matrix,50,c(5,200))
#Using the updated code for the BUGS models to ensure that the data are compatible
tau.14<-tau.X.14/(1-rho.X*rho.X)
data.list <- list()
for(i in 1:50){
  X.5 <- rnorm(theta.gen[i,"N"],theta.gen[i,"theta.5"],tau.X.5)
  mu.14 <- theta.gen[i,"theta.14"]+sqrt(tau.X.5/tau.X.14)*rho.X*
    (X.5 - theta.gen[i,"theta.5"])
  X.14 <- rnorm(theta.gen[i,"N"],mu.14,tau.14)
  #Data must be stored as a list of named lists
  #Each list element is a set of data for a BUGS model.
  data.list[[i]] <- list(X.5 = X.5,X.14 = X.14)
}

###Calculating the EVSI using our updated cost-effectiveness results.
#This function calculates the variance of the posterior distribution 
#across different sample sizes
g<-mm.post.var(
  #The place where we saved the BUGS model.
  filein.model,
  #The list containing the future data that we generated above
  data.list,
  #The variable name of the sample size (i.e. how many data points we will collect)
  "N",
  #The sample sizes we used to generate the future data
  theta.gen[,"N"],
  #The function that calculates the effects from the parameters
  effects,
  #The function that calculates the costs from the parameters
  costs,
  #The bcea object containing all the cost-effectiveness analysis
  he=he,
  #The evppi object containing the EVPPI results
  #Because we have included the evppi object, we don't need to specify
  #the parameters that will be informed by the future data collection.
  evi=evi,
  #Any additional data that is needed for the BUGS model to run
  #In this case, the mean vectors and precision matrices
  #If you had data from a previous trial it was be included here.
  data.stats=data.evsi,
  #The Bayesian updating should be done with BUGS (rather than JAGS).
  update="bugs",
  #The number of iterations taken from each posterior distribution to
  #calculate the posterior variance.
  n.iter=1000,
  #The number of "future" posterior samples that should be taken
  Q=50
)

#Using the previous function, evsi.calc calculates the EVSI across different
#sample sizes and willingness-to-pay parameters by rescaling the fitted values
#from the EVPPI calculations.
evsi.BK<-evsi.calc(g)

#Launches the WebApp so you can see the results of the EVSI calculation
launch.App(evsi.BK)
