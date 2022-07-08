# This files essentially does the same analysis as OpenBUGS

# NB: Remember to set the working directory to the one where the files are stored.
# You can do this by going on the bottom-right panel, selecting the tab "Files"
# (which is the default/first one), then navigate to "practicals", then
# "01_monte-carlo" and then "R-files". Then click on the right-most button
# "More" and select from the slider menu "Set As Working Directory"
#
# Alternatively, you can do this in the R terminal by typing the command
# setwd("PATH-TO-YOUR-FOLDER")
# The path changes depending on the operating system. In the case of the 
# Binder Virtual Machine, this would be
# setwd("~/practical/01_monte-carlo/R-files")
# or in a more verbose way
# setwd("/home/rstudio/practical/01_monte-carlo/R-files")
# (note that '~' is a shortcut for the home folder '/home/rstudio', in this case)

# Utility package
library(bmhe)

# Set up 
n = 10 # number of "trials"
theta = 0.5 # probability of "success" (heads) - the coin is "balanced"

# Simulates n tosses of a coin
nsim = 1000 # number of simulations
y = rbinom(nsim,n,theta)  # NB: the built-in R command 'rxxxx' (eg 'rbinom')
                          # is used to simulate random values ('r') from a
                          # given distribution ('binom'=Binomial). Type 
                          # '?Distributions' (without ''!) in your R terminal
                          # to get some help

# Computes the probability of 8 or more heads in each simulation
threshold = 7.5 # y is only integers, so if the threshold is 7.5, it means
                # there are 8 or more heads in a single simulation
P8 = y > threshold

# First creates a matrix by stacking together the two columns (one for P8
# and the other for y), using the 'cbind' R built-in command (type '?cbind'
# to open the help window)
sims = cbind(P8,y)
# Then uses the function 'stats' (which is made available by 'utils.R') to show
# summary statistics like BUGS would do
stats(sims)

# Now uses R plots to show graphical summaries of the distributions
hist(as.numeric(P8))  # note that R defines 'P8' as a "logical" variable
                      # ie it takes values TRUE/FALSE. In order to plot its 
                      # distribution, we need to turn it into a "numeric"
                      # variable, which we do by using the R built-in 
                      # 'as.numeric' function
hist(y)               # We don't need any transformation for y, because it is
                      # already a numeric variable



#### OPTIONAL
# NB: you can use more modern tools for visualisation and data-wrangling
# For this, you need the package 'tidyverse', which includes several other 
# packages. It takes a bit to get used to it, but once you do, it is *very*
# powerful
library(tidyverse)

# This creates a histogram of the distribution for y
# First, you need to turn your vector of simulations (y) into a different
# type of object (it needs to be either a 'data.frame' or a 'tibble'). You
# can also use 'piping' ('%>%'), which takes the argument to the left and 
# applies the function to the right. So 
# y %>% as_tibble
# is equivalent to 
# as_tibble(y)
# and what it does is to create a 'tibble' (a structured data frame) with the 
# variable (column) y.
y %>% as_tibble() %>% 
  # Then you start creating a 'ggplot' object. 'ggplot2' is a very powerful
  # package for advanced graphics. It's a lot more structured than the simpler
  # 'base' graphical engine and the syntax is a lot more complex. BUT: it can
  # be very, very versatile and powerful, so worth playing a bit with it...
  # Here, we define the "aesthetic", eg what are the variables involved in the 
  # plot. In this case, we want to do a histogram, or more precisely a barplot,
  # and so the only variable involved is 'y'.
  ggplot(aes(y)) + 
  # Then you start adding components (layers) to the graph. In this case, we 
  # specify the "geometry", that is the type of graph we want to make. In this
  # case, we will do a barplot, so the relevant command is 'geom_bar'
  geom_bar() + 
  # We can also specify a set "theme" for how the graph will look like. You 
  # can check '?ggtheme' to see which ones are available by default
  theme_bw() + 
  labs(title="Histogram of the distribution for y")

# Here's another example, with yet more structure to it...
y %>% as_tibble() %>% mutate(simulation=1:length(y)) %>% 
  ggplot(aes(simulation,y)) + geom_line() + theme_classic() + 
  labs(title="Traceplot for y") + 
  theme(
    plot.title = element_text(color="red", size=22, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  )
