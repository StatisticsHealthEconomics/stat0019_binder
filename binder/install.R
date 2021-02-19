# Removes the installation of survHE, rstan and INLA, which take forever. However, 
# all the relevant background packages are either installed or installable, so 
# these can always be brought in with a simple install.package...

#install.packages("survHE")
install.packages("tidyverse")
install.packages("rmarkdown")
install.packages('BCEA')
install.packages('Rjags')
install.packages("R2OpenBUGS")
install.packages("R2jags")
#install.packages("rstan")
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)


