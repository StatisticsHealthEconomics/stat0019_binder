### 1. Laplace model
# Sets the data
y=241945
n=493527
data=list(y=y,n=n)

# And inspect the values and the whole object
names(data)
data$y

# Runs the model
# Use 'here::here' to make sure the correct path is specified to get the model
# code file
filein=here::here("practical","03_bcea","ModelLaplace.txt")
params="theta"
# Two versions of the initial values --- one is deterministic, the other is random
inits_det=list(list(theta=.1),list(theta=.9))
inits_ran=function(){list(theta=runif(1))}

library(R2OpenBUGS)
model= bugs(data=data,inits=inits_det,parameters.to.save=params,
            model.file=filein,n.chains=2,n.iter=10000,
            n.burnin=4500,n.thin=1,DIC=TRUE) 

# Inspect the outcome of the BUGS model
names(model)
model$n.iter

print(model,digits=3)

# Makes *all* of the BUGS model output available to the main R workspace
attach.bugs(model)

# Visual representation of the output
hist(theta)


### 2. BCEA
# Loads the data and shows what's inside the object - again using 
# 'here::here' to ensure we get the correct full path to the file, 
# irrespective of what the current working directory is
load(here::here("practical","03_bcea","Vaccine.RData"))
names(he)

# Uses BCEA to post-process the object
library(BCEA)
# Computes the ICER and shows its value
ICER=mean(he$delta.c)/mean(he$delta.e) 
ICER

# Computes the EIB for k=30000 and shows its value
k=30000
EIB=k*mean(he$delta.e)-mean(he$delta.c)
EIB

# Graphical representation of the economic results
ceplane.plot(he,wtp=10000)
eib.plot(he, plot.cri=FALSE)
contour(he)