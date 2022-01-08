install.packages("rmarkdown")
install.packages('BCEA')
install.packages('rjags')
install.packages("R2OpenBUGS")
install.packages("R2jags")
remotes::install_github("giabaio/bmhe_utils")  
install.packages("rstan",dep=TRUE)
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
# For now only manages to install the CRAN version of survHE
# That's not ideal as there are some changes in 'devel' that 'main' doesn't have
# but at least gets it going...
install.packages("survHE")
