# Loads the R2OpenBUGS package to connect R to OpenBUGS
library(R2OpenBUGS)

### Exercise 1.
# Loads the data on costs only into R from the txt file (list originally prepared for BUGS)
cost.data=source("cost-data.txt")$value
# Can inspect the resulting object
names(cost.data)
# NB: 'cost.data' is an object containing different variables so need to use the '$' notation to access them!
head(cost.data$cost1)		# This for example shows the first 6 values of the variable 'cost1' inside 'cost.data'
							
# Runs the BUGS model from R
model.file="normal-mod.txt"				# Specifies the path to the model file
params <- c("mu","ss","ls","delta.c","dev","total.dev")	# Defines the list of parameters
n.chains=2						# How many parallel chains?
n.burnin=1000						# Number of simulations to discard before convergence
n.iter=2000						# Total number of simulations
n.thin=1						# Saves one simulation every n.thin (may need 
							#   to increase to reduce autocorr)
debug=FALSE						# If set to FALSE then does not show OpenBUGS

# Finally runs the model by calling OpenBUGS in the background
m=bugs(data=cost.data,inits=NULL,model.file=model.file,parameters.to.save=params,
	n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T,debug=debug)

# Now can print & manipulate the output (stored in the object 'm')
print(m)
# NB: the object 'm' contains all the relevant simulations & information from the BUGS model
#     For example, the object 'sims.list' is a R list with the simulation values for all the monitored node
#     So you can access them by simply using the $ notation in R --- see below
#     Also, notice that BUGS increases all variables by one dimension (because it adds the simulation values): 
#     so scalars become vectors, vectors become matrices, matrices become arrays.
plot(
	m$sims.list$mu[,1],					# What to plot on the x-axis?
	m$sims.list$mu[,2],					# What to plot on the y-axis?
	xlab="Population average cost in arm 1 (x£1000)",	# Label for the x-axis
	ylab="Population average cost in arm 2 (x£1000)",	# Label for the y-axis
	pch=20,							# Type of symbol used in the plot 									# (see http://www.endmemo.com/program/R/pchsymbols.php)
	cex=.8,							# Amount by which plotting text and symbols should be scaled 
								# relative to the default (here 80%)
	main="Joint distribution of mean costs"			# Title for the plot
)
# Rescales the difference in cost by £1000
delta.c=m$sims.list$delta.c*1000
hist(
	delta.c,						# variable to use for the histogram 
	xlab="Cost differential (£)",				# Label for the x-axis
	main="Posterior distribution"				# Title of the graph
)
# Adds a vertical line at 0
abline(
	v=0,		# Plots a vertical line at 0
	lwd=2,		# line width (=2pt)
	lty=2		# line type (=2 = a line)
)
# Computes the posterior probability that delta.c is positive (increase in average costs) as the 
# number of times that delta.c is positive divided by the total number of simulations m$n.sims
sum(delta.c>0)/m$n.sims


## Exercise 3.
# Now loads the data for costs & utilities into R from the txt file (list originally prepared for BUGS)
cost.utility=source("cost-util-data.txt")$value		
model.file="cgeg-mod.txt"				# Specifies the new file with the cost-effectiveness model
inits1=source("cgeg-inits1.txt")$value			# Loads the initial values for the 1st chain
inits2=source("cgeg-inits2.txt")$value			# Loads the initial values for the 2nd chain
inits3=source("cgeg-inits3.txt")$value			# Loads the initial values for the 3rd chain
inits=list(inits1,inits2,inits3)			# Combines them into a single list --- can be lazy and only use 2...
params=c("mu.c","mu.e","delta.c","delta.e",		# Defines the list of parameters to be monitored
	"shape.c","shape.e","beta","INB","CEAC")
# Finally runs the model by calling OpenBUGS in the background
n.burnin=1000
n.iter=4000			# NB: this adds 3000 simulations to the 1000 of burnin
n.chains=3
m2=bugs(data=cost.utility,inits=inits,model.file=model.file,parameters.to.save=params,
	n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T,debug=debug)

# Now can print & manipulate the output (stored in the object 'm')
print(m2)

# Plots the cost-effectiveness plane
plot(m2$sims.list$delta.e,m2$sims.list$delta.c,pch=20,cex=.8,xlab="Effectiveness differential",
	ylab="Cost differential",main="Cost-effectiveness plane")

## Exercise 4.
# Figures out which of the indices is associated with willingness to pay = £500
K=numeric()			# defines an empty numeric vector, named 'K'
K.space=0.1			# in line with the BUGS model
for (j in 1:11) {		# Loop to compute the K vector for all j values between 1 and 11
	K[j]=(j-1)*K.space
}
idx=which(K==0.5)	# Now finds the index for which K is equal to 0.5 (=£500, as all costs are in £1000)
			# NB: When dealing with logical operations (eg if, while, absolute equality), R wants the "double ="
# Now can compute the statistics for INB at K=£500
mean(m2$sims.list$INB[,idx])	# NB: the object INB is a matrix with as many rows as you've saved simulations and 11 columns
				#     so you have to pick all the rows in column 'idx'
quantile(m2$sims.list$INB[,idx],.025)	# 2.5% quantile of the posterior distribution = lower end of the 95% interval
quantile(m2$sims.list$INB[,idx],.975)	# 97.5% quantile of the posterior distribution = upper end of the 95% interval

# Now can identify the probability of cost-effectiveness for the new treatment (arm 2)
# First a graphical representation --- "visual inspection"
hist(m2$sims.list$INB[,idx],xlab="Incremental Net Benefit for k=£500",main="Posterior distribution")
abline(v=0,lwd=2,lty=2)
# Then can actually compute the probability as the proportion of simulations that are positive
prob.ce=sum(m2$sims.list$INB[,idx]>0)/m2$n.sims

# Finally plots the CEAC
# First extracts the mean values of the CEAC for each willingness to pay
CEAC=numeric()
for (i in 1:length(K)) {		# length(K) is the length of the vector K = number of willingness to pay points
	CEAC[i]=mean(m2$sims.list$CEAC[,i])	# takes the mean of each column for the matrix CEAC
}
plot(K,CEAC,xlab="Willingness to pay (x£1000)",ylab="Cost-effectiveness acceptability curve",main="",
	t="l"	# This instructs R to plot using a solid curve (t="type", "l"=line)
)
