# Regression analysis with Galton's data 
# See: http://www.statistica.it/gianluca/teaching/intro-stats/regression.html & following sections

# This script loads the data into the R workspace, makes some pre-processing (to plot the raw data), then runs OpenBUGS remotely
# using the model stored in the file 'model_galton.txt' and finally post-process the output of the MCMC analysis 

# Loads the data (they are in a tabular format as a .txt file)
galton=read.table("Galton.txt",header=TRUE)

##################
# Pre-processing #
##################
# Visualise the data
head(galton)

# NB Can use fancier and more modern tools (eg based on the 'tidyverse' approach and a suite of packages aimed at providing
#    a much better ability to explore data, eg dplyr, ggplot etc). You can do ***much*** more than the following command,
#    but this is just to get you going... See: https://www.tidyverse.org/. 
#    Of course, if you have not installed the packages yet, you need to do this first, e.g. 'install.packages("dplyr")'
library(dplyr)
as_tibble(galton)

# Uses ggplot (a package for advanced graphing in R) to summarise graphically the data --- recreates Fig 6.1
library(ggplot2)
ggplot(galton,aes(Father,Height))+
  labs(x="Father's height (inches)",y="Child's height (inches)") +
  geom_point(size=1,aes(colour=Family)) + 
  theme_classic() +
  theme(axis.line = element_line(size=.7, colour = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title=element_text(size = 20),
        text=element_text(size = 12),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12),
        axis.ticks.length=unit(.2, "cm")
  ) 

# Recreates Fig 6.3 showing the data on the centered covariate (father's age)
plot(scale(galton$Father,scale=F),galton$Height,axes=FALSE,xlab="Centred father's height (inches)",ylab="Child's height (inches)",main="",cex=.7,col="grey40")
axis(1,at=pretty(range(scale(galton$Father,scale=FALSE)))); axis(2,at=pretty(range(galton$Height)))
m2=lm(galton$Height~scale(galton$Father,scale=FALSE))
abline(m2,lwd=2)
points(rep(0,3),c(-10,0,m2$coefficients[1]),t="l",lwd=1.5,lty=2)
points(c(-100,-10,0),rep(m2$coefficients[1],3),t="l",lwd=1.5,lty=2)
legend("topleft",legend=c(paste("$\\hat\\beta_0=",format(m2$coefficients[1],digits=3,nsmall=2),"$"),
                          paste("$\\hat\\beta_1=",format(m2$coefficients[2],digits=3,nsmall=2),"$")
                          ),
       bty="n",lty=1,col=c("white","white",cex=1.2)
)

#################
# Model running #
#################
# Now runs the Bayesian linear regression model
# Defines the matrix of covariates
X=cbind(rep(1,nrow(galton)),             # vector of ones (to obtain the intercept)
        scale(galton$Father,scale=F),    # "scaled" (= centered) father's height
        scale(galton$Mother,scale=F)     # "scaled" (= centered) mother's height
)

# You can visualise the covariates matrix, eg by typing
head(X)
# Also you can convince yourself that the covariates have been rescaled to have 0 mean by doing the following.
# This computes the mean along each of the columns in the matrix (that's the second dimension, hence the argument '2').
# Essentially, this instructs R to 'apply' the function 'mean' to the object 'X', along its dimension '2' (columns)
apply(X,2,mean)
# As you can see, the mean for the second and third columns are so small that they are exactly 0...

library(R2OpenBUGS)
data=list(
  y=galton$Height,               # children's height
  X=X,                           # matrix of covariates (father's & mother's height)
  n=nrow(X),                     # number of observations (= number of rows in X)
  K=ncol(X),                     # number of covariates (= number of columns in X)
  a=0.1,b=0.1, 
  mu.beta=c(65,0,0),
  tau.beta=c(20,10,10)^-2
)

# Runs the BUGS model by calling OpenBUGS remotely
model=bugs(data=data,                                     # defines the data list
           inits=NULL,                                    # for simplicity does not specify initial values (BUGS will do it)
           parameters.to.save=c("beta","sigma"),          # defines the parameters to monitor
           model.file="model_galton.txt",                 # specifies the file containing the model code
           n.chains=2,                                    # specifies the number of chains to run in parallel
           n.iter=10000,                                  # total number of iterations
           n.burnin=5000                                  # number of iterations in the burn-in phase (to be discarded)
)

###################
# Post-processing #
###################
# Prints a summary table with the model results
print(model,digits=3)

# Plots traceplots for the model parameters (beta[1], in this case)
# Uses the elements in the object 'sims.array' stored in the model output
# First plots the results for 'beta[1]' (intercept) for the first chain (indicated by the fact that the second index in 'sims.array' is '1')
plot(model$sims.array[,1,"beta[1]"],col="blue",t="l",bty="l",ylab="beta[1]",xlab="Iterations")
# Then adds on the plot the results for 'beta[1]' (intercept) for the second chain (indicated by the fact that the second index is '2')
points(model$sims.array[,2,"beta[1]"],col="red",t="l")

# Autocorrelation plot for beta[3] (uses the R function 'acf')
# Uses the elements in the object 'sims.list' stored in the model output
acf(model$sims.list$beta[,3],main="beta[3]")

# Uses the pacakge 'coda' (which is automatically loaded with R2OpenBUGS) to produce Gelman-Rubin plots
coda::gelman.plot(model)
