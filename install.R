options(timeout=600)
install.packages("BCEA", repos = c("https://giabaio.r-universe.dev", "https://cloud.r-project.org"))
install.packages("tidyverse")
install.packages("here")
install.packages("remotes")
install.packages("rmarkdown")
install.packages("cli")
install.packages("rjags")
install.packages("R2OpenBUGS")
install.packages("R2jags")
install.packages("survHE", repos = c("https://giabaio.r-universe.dev", "https://cloud.r-project.org"))
install.packages('INLA', repos='https://inla.r-inla-download.org/R/stable')
install.packages("MatrixModels",repos="https://r-forge.r-project.org/")
install.packages("sn")
install.packages("voi")
install.packages("multinma")
remotes::install_github("giabaio/survHEinla")
remotes::install_github("StatisticsHealthEconomics/outstandR")
remotes::install_github("n8thangreen/simcovariates")
## For some reason, packages that are *only* on r-universe.dev break the mybinder
## installation... They work OK if installed through remotes. And in fact, other
## repos don't seem to interact badly with bspm (which looks for packages on r2u
## and installs it using apt install r-cran-<pkg>). But r-universe.dev has something
## that messes up with that system...
#install.packages("outstandR", repos = c("https://statisticshealtheconomics.r-universe.dev", "https://cloud.r-project.org"))
#install.packages("simcovariates", repos = c("https://n8thangreen.r-universe.dev", "https://cloud.r-project.org"))
#install.packages("survHEinla", repos = c("https://giabaio.r-universe.dev", "https://cloud.r-project.org"))
