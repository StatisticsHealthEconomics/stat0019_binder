##### Example To Test Plotting #####
### Function to translate the output from the evsi command (evsi.obj) into
### an object to plot
library(dplyr)
evsi.plot.adapt <- function (ouputs, inputs, pars, evsi.obj, method) 
{
  if (is.null(evsi.obj)) {
    stop("You have not provided the EVSI. Please include an evsi.obj object from the evsi function")
  }

  evsi <- evsi.obj %>% tidyr::pivot_wider(names_from = "k", values_from = "evsi")
  evpi <- voi::evpi(ouputs)
  evppi <- evppi(ouputs, inputs, pars, method)
  to.return <- list(evsi = evsi, attrib = list(k = as.numeric(colnames(evsi))[-1], 
                                                   N = evsi$n), 
                    evppi = evppi, evpi = evpi)
  class(to.return) <- "evsi.plot"
  return(to.return)
}

### Plots evsi across willingess to pay. Originally plot.evsi in EVSI package
evsi.wtp.plot <- function (evsi, pos = c(0, 0.8), N = NULL) 
{
  if (class(evsi) != "evsi.plot") {
    stop("plot.evsi must be used with an evsi.plot object created using the evsi.plot.adapt function.")
  }
  alt.legend <- pos
  if (is.numeric(alt.legend) & length(alt.legend) == 2) {
    temp <- ""
    if (alt.legend[2] == 0) 
      temp <- paste0(temp, "bottom")
    else if (alt.legend[2] != 0.5) 
      temp <- paste0(temp, "top")
    if (alt.legend[1] == 1) 
      temp <- paste0(temp, "right")
    else temp <- paste0(temp, "left")
    alt.legend <- temp
    if (length(grep("^((bottom|top)(left|right)|right)$", 
                    temp)) == 0) 
      alt.legend <- FALSE
  }
  if (is.logical(alt.legend)) {
    if (!alt.legend) 
      alt.legend = "topright"
    else alt.legend = "topleft"
  }
  N.length <- length(evsi$attrib$N)
  if (class(N) == "numeric") {
    select.length <- length(N)
    select <- rep(NA, select.length)
    for (i in 1:select.length) {
      select[i] <- which.min((evsi$attrib$N - N[i])^2)
    }
    warning("The EVSI is only calculated for sample sizes listed in evsi$attrib$N. If N is not in this list then EVSI calculated for closest possible sample size.")
    N <- evsi$attrib$N[select]
  }
  if (is.null(N)) {
    select <- 1:N.length
    N <- evsi$attrib$N
    select.length = length(select)
  }

  EVSI <- evsi$evsi[select, -1]
  x.range <- range(range(evsi$evpi$k), range(evsi$evppi$k), 
                   range(evsi$attrib$k))
  x.range[2] <- x.range[2] * 1.05
                   
  plot(evsi$evpi$k, evsi$evpi$evpi, t = "l", xlab = "Willingness to pay", 
       ylab = "", main = "Expected Value of Sample Information", 
       lwd = 2, ylim = range(range(evsi$evpi$evpi), range(evsi$evppi$evppi), 
                             range(EVSI)),
       xlim = x.range,
  col = "black")
    points(evsi$evppi$k, evsi$evppi$evppi, t = "l", col = "black", 
           lty = 1)
    
    text(evsi$evpi$k[length(evsi$evpi$k)], evsi$evpi$evpi[length(evsi$evpi$k)], "EVPI", pos = 4,cex = 0.6)
    text(evsi$evppi$k[length(evsi$evppi$k)], evsi$evppi$evppi[length(evsi$evppi$k)], "EVPPI", pos = 4, cex = 0.6)
    colours <- colorRampPalette(colors = c("grey90", "grey60", "grey20"))(select.length)
    if (length(evsi$attrib$k) < 30) {
      for (s in 1:select.length) {
        points(evsi$attrib$k, EVSI[s, ], pch = 19, col = colours[s])
      }
    }
    if (length(evsi$attrib$k) >= 30) {
      for (s in 1:select.length) {
        points(evsi$attrib$k, EVSI[s, ], type = "l", col = colours[s])
      }
    }
    if (select.length == 1) {
      legend(alt.legend, c("EVPI", "EVPPI for focal parameters", 
                           paste("EVSI for sample size of", N)), col = c("black", 
                                                                                "black", colours), cex = 0.7, box.lwd = 0, box.col = "white", 
                                                                                bg = "white", lty = c(1, 1, 1), lwd = c(2, 1))
    }
    if (select.length > 1) {
      legend(alt.legend, legend = c(min(N), rep(NA, max(0, 
                                                        select.length - 2)), max(N)), fill = colours, border = colours, 
             cex = 0.75, y.intersp = max(0.1, 1.2/select.length), 
             box.lwd = 0, box.col = "white", bg = "white")
    }
    box()
}

### Plots evsi across sample size for a fixed willingness to pay. 
### Originally plot.samplesize in EVSI package
evsi.ss.plot <- function (evsi, k = NULL, pos = c("bottomright")) 
{
  alt.legend <- pos
  if (is.numeric(alt.legend) & length(alt.legend) == 2) {
    temp <- ""
    if (alt.legend[2] == 0) 
      temp <- paste0(temp, "bottom")
    else if (alt.legend[2] != 0.5) 
      temp <- paste0(temp, "top")
    if (alt.legend[1] == 1) 
      temp <- paste0(temp, "right")
    else temp <- paste0(temp, "left")
    alt.legend <- temp
    if (length(grep("^((bottom|top)(left|right)|right)$", 
                    temp)) == 0) 
      alt.legend <- FALSE
  }
  if (is.logical(alt.legend)) {
    if (!alt.legend) 
      alt.legend = "topright"
    else alt.legend = "topleft"
  }
  if (!(class(k) %in% c("numeric" , "integer"))) {
    k.select <- which.max(evsi$evpi$evpi)[1]
    k <- evsi$attrib$k[k.select]
  }
  
  if (class(k) %in% c("numeric" , "integer")) {
    k.select <- which.min(abs(evsi$attrib$k - k))[1]
    k <- evsi$attrib$k[k.select]
  }

  if (length(evsi$attrib$N) == 1) {
    stop("This plot gives the EVSI across sample size. Do not use on a single design.")
  }
  EVSI <- evsi$evsi[k.select + 1]
  evppi <- evsi$evppi$evppi[k.select] 
  
  plot(1, 1, ylim = c(min(EVSI) * 0.95, evppi), 
       xlim = c(min(evsi$attrib$N), max(evsi$attrib$N)), col = "white", 
       xlab = expression("Sample Size"), ylab = "Per Person EVSI", 
       oma = c(0, 0, -1, 0), main = "Expected Value of Sample Information across Sample Size")

  if (length(evsi$attrib$N) < 15) {
      points(evsi$attrib$N, t(EVSI), pch = 19, lwd = 2, 
             lty = 1)
  }
  if (length(evsi$attrib$N) >= 15) {
      points(evsi$attrib$N, t(EVSI), type = "l", lwd = 2, 
             lty = 1)
  }
  
  abline(h = evppi, col = "darkgray", lwd = 3, lty = 2)
  legend(alt.legend, c("EVSI", "EVPPI"), col = c("black", "darkgray"), 
         lwd = c(2,3), lty = c(1,2), box.lwd = 0, 
         box.col = "white", bg = "white")
  box()
}

### Calculates the ENBS (internal function)
ENBS.fun <- function(evsi, Pop, Time, Dis, cost) {
  enbs <- evsi * Pop/Dis * (1 - exp(-Dis * Time)) - cost
  return(enbs)
}

### Calculates the standard devation for the ENBS to account for unknown costs
### (internal function)
ENBS.sd.calc <- function(evsi.sd, Pop, Time, Dis, cost.sd) {
  var <- (Pop/Dis * (1 - exp(-Dis * Time)))^2 * evsi.sd^2 + 
    cost.sd^2
  return(sqrt(var))
}

### Plots the probability of a cost-effective trial as a sentivity analysis
### to decision horizon and population size. Originally plot.prob.ce in EVSI package
evsi.prob.plot <- function (evsi, trial.cost = NULL, setup = NULL, pp = NULL, 
          Pop = c(0, 10000), Time = c(1, 20), Dis = 0.035, k = NULL, N = NULL, 
          pos = c("topright")) 
{
  alt.legend <- pos
  if (is.numeric(alt.legend) & length(alt.legend) == 2) {
    temp <- ""
    if (alt.legend[2] == 0) 
      temp <- paste0(temp, "bottom")
    else if (alt.legend[2] != 0.5) 
      temp <- paste0(temp, "top")
    if (alt.legend[1] == 1) 
      temp <- paste0(temp, "right")
    else temp <- paste0(temp, "left")
    alt.legend <- temp
    if (length(grep("^((bottom|top)(left|right)|right)$", 
                    temp)) == 0) 
      alt.legend <- FALSE
  }
  if (is.logical(alt.legend)) {
    if (!alt.legend) 
      alt.legend = "topright"
    else alt.legend = "topleft"
  }
  if (!(class(k) %in% c("numeric" , "integer"))) {
    k.select <- which.max(evsi$evpi$evpi)[1]
    k <- evsi$attrib$k[k.select]
  }
  
  if (class(k) %in% c("numeric" , "integer")) {
    k.select <- which.min(abs(evsi$attrib$k - k))[1]
    k <- evsi$attrib$k[k.select]
  }
  
  if (class(evsi$attrib$N) != "numeric") {
    N.select <- 1
    if (class(N) != "numeric") {
      N <- evsi$attrib$N
    }
  }
  if (class(evsi$attrib$N) == "numeric") {
    if (class(N) != "numeric") {
      N.select <- ceiling(length(evsi$attrib$N)/2)
      N <- evsi$attrib$N[N.select]
    }
    if (class(N) == "numeric") {
      N.select <- which.min(abs(evsi$attrib$N - N))
    }
  }
  
    type.evsi <- "det"
    evsi.focal <- as.numeric(c(evsi$evsi[N.select, k.select + 1]))
    evsi.params <- c(as.numeric(evsi$evsi[N.select, k.select + 1]), 
                     0)

  if (is.null(trial.cost)) {
    if (class(N) == "character") {
      stop("Please define the trial costs using trial.costs or the sample size of experiment using N=")
    }
    if (is.null(setup) || is.null(pp)) {
      stop("Please give the trial costs using either trial.costs for the full costs\n                                         or setup and pp to give the set up and per person costs ")
    }
    setup.params <- c(mean(setup), (range(setup)[2] - range(setup)[1])/4)
    pp.params <- c(mean(pp), (range(pp)[2] - range(pp)[1])/4)
    trial.cost <- c(setup.params[1] + pp.params[1] * N, sqrt(setup.params[2]^2 + 
                                                               N^2 * pp.params[2]^2))
  }
  colours <- colorRampPalette(c("black", "grey20", "grey40", "grey60", 
                                       "grey90", "white"))(100)

                                       if (class(trial.cost) == "numeric") {
                                         if (length(trial.cost) == 1) {
                                           trial.cost.params <- c(trial.cost, 0)
                                         }
                                         if (length(trial.cost) == 2) {
                                           trial.cost.params <- c(mean(trial.cost), (range(trial.cost)[2] - 
                                                                                       range(trial.cost)[1])/2)
                                         }
                                         Time.min <- min(Time)
                                         Time.max <- max(Time)
                                         Pop.min <- min(Pop)
                                         Pop.max <- max(Pop)
                                         dens.points <- 100
                                         Time.seq <- seq(Time.min, Time.max, length.out = dens.points)
                                         Pop.seq <- seq(Pop.min, Pop.max, length.out = dens.points)
                                         Prob.mat <- matrix(NA, nrow = length(Time.seq), ncol = length(Pop.seq))
                                         for (i in 1:length(Time.seq)) {
                                           for (j in 1:length(Pop.seq)) {
                                             ENBS.mean <- ENBS.fun(evsi.params[1], Pop.seq[j], 
                                                               Time.seq[i], Dis, trial.cost.params[1])
                                             ENBS.sd <- ENBS.sd.calc(evsi.params[2], Pop.seq[j], 
                                                                     Time.seq[i], Dis, trial.cost.params[2])
                                             Prob.mat[i, j] <- pnorm(0, ENBS.mean, ENBS.sd, 
                                                                     lower.tail = FALSE)
                                           }
                                         }
                                       }
                                       image(x = Time.seq, y = Pop.seq, z = Prob.mat, col = colours, 
                                             main = "Probability of Cost-Effective Trial", xlab = "Decision Horizon", 
                                             ylab = "Incidence Population", xlim = c(Time.min, Time.max), 
                                             ylim = c(Pop.min, Pop.max), breaks = seq(0, 1, length.out = 101))
                                       legend(alt.legend, c("Prob=0", rep(NA, 98/2), "Prob=.5", 
                                                            rep(NA, 96/2), "Prob=1"), fill = colours, border = colours, 
                                              cex = 0.75, y.intersp = 0.15, box.lwd = 0, box.col = "white", 
                                              bg = "white")
                                       box()
}

### Determines the optimal sample size of the study/
### Originally optim.samplesize in EVSI package
optim.ss <- function (evsi, setup, pp, Pop, Time, k = NULL, Dis = 0.035) 
{
  if (!(class(k) %in% c("numeric" , "integer"))) {
    k.select <- which.max(evsi$evpi$evpi)[1]
    k <- evsi$attrib$k[k.select]
  }
  
  if (class(k) %in% c("numeric" , "integer")) {
    k.select <- which.min(abs(evsi$attrib$k - k))[1]
    k <- evsi$attrib$k[k.select]
  }
  
  EVSI <- evsi$evsi[, k.select + 1]
  if ((length(setup) > 1) || (length(pp) > 1)) {
    setup <- mean(setup)
    pp <- mean(pp)
  }
  ENBS <- as.numeric(as.matrix(Pop * EVSI/Dis * (1 - exp(-Dis * Time)) - setup - 
                                 pp * evsi$attrib$N))
  max.select <- which.max(ENBS)
  max.less <- max.select - 1
  max.greater <- max.select + 1
  if ((max.less < 1) | (max.greater > length(ENBS))) {
    N.max <- evsi$attrib$N[max.select]
    ENBS.max <- ENBS[max.select]
    warning("Optimal sample size is at the limit of the considered values for N. An alternative sample size may be optimal,\n            please consider alternative values of N in the evsi.calc function.")
  }
  else {
    N.fit <- evsi$attrib$N[c(max.less, max.select, max.greater)]
    N2 <- N.fit^2
    ENBS.fit <- ENBS[c(max.less, max.select, max.greater)]
    model.optim <- lm(ENBS.fit ~ N.fit + N2)
    N.max <- round(-model.optim$coefficients[2]/(2 * model.optim$coefficients[3]))
    ENBS.max <- predict(model.optim, list(N.fit = N.max, 
                                          N2 = N.max^2))
  }
  tol <- ENBS.max - abs(ENBS.max * 0.05)
  limits <- which(ENBS > tol)
  if(length(limits) > 0){
    N.range <- range(evsi$attrib$N[limits])
  }
  if(length(limits) == 0){
    c.quad <- model.optim$coefficients[1] - tol
    b.quad <- model.optim$coefficients[2]
    a.quad <- model.optim$coefficients[3]
    
    N.range <- c((-b.quad - sqrt(b.quad^2 - 4 * a.quad * c.quad)) / (2 * a.quad),
                 (-b.quad + sqrt(b.quad^2 - 4 * a.quad * c.quad)) / (2 * a.quad))
  }
  return(list(SS.max = N.max, ENBS = ENBS.max, SS.I = N.range))
}

### Plots ENBS across sample size for fixed WTP/population size and decision horizon
### Originally plot.enbs in EVSI package
evsi.enbs.plot <- function (evsi, setup, pp, Pop = 10000, Time = 10, 
                            Dis = 0.035, k = NULL, N = NULL, pos = c("bottomright"),
                            ylim = NULL, legend = TRUE) 
{
  alt.legend <- pos
  if (is.numeric(alt.legend) & length(alt.legend) == 2) {
    temp <- ""
    if (alt.legend[2] == 0) 
      temp <- paste0(temp, "bottom")
    else if (alt.legend[2] != 0.5) 
      temp <- paste0(temp, "top")
    if (alt.legend[1] == 1) 
      temp <- paste0(temp, "right")
    else temp <- paste0(temp, "left")
    alt.legend <- temp
    if (length(grep("^((bottom|top)(left|right)|right)$", 
                    temp)) == 0) 
      alt.legend <- FALSE
  }
  if (is.logical(alt.legend)) {
    if (!alt.legend) 
      alt.legend = "topright"
    else alt.legend = "topleft"
  }
  if (!(class(k) %in% c("numeric" , "integer"))) {
    k.select <- which.max(evsi$evpi$evpi)[1]
    k <- evsi$attrib$k[k.select]
  }
  
  if (class(k) %in% c("numeric" , "integer")) {
    k.select <- which.min(abs(evsi$attrib$k - k))[1]
    k <- evsi$attrib$k[k.select]
  }
  if (class(N) != "numeric") {
    N <- evsi$attrib$N
  }
  if (class(N) == "character") {
    stop("Please define the sample size of your experiment using N=")
  }
  length.N <- length(N)
  N.select <- array(NA, dim = length.N)
  for (i in 1:length.N) {
    N.select[i] <- which.min((evsi$attrib$N - N[i])^2)
  }
  
  prob <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  length.prob <- length(prob)
  type.evsi <- "det"
  evsi.params <- cbind(evsi$evsi[N.select, k.select + 1], 0)
  if (length(setup) != 2 || length(pp) != 2) {
    stop("Please give minimum and maximum values for the setup and per person trial costs.")
  }
  setup.params <- c(mean(setup), (range(setup)[2] - range(setup)[1])/4)
  pp.params <- c(mean(pp), (range(pp)[2] - range(pp)[1])/4)
  trial.cost <- array(NA, dim = c(length.N, 2))
  for (i in 1:length.N) {
    trial.cost[i, ] <- c(setup.params[1] + pp.params[1] * 
                           N[i], sqrt(setup.params[2]^2 + N[i]^2 * pp.params[2]^2))
  }
  colours <- colorRampPalette(c("black", "grey20", "grey40", "grey60", "grey90", "white"))(100)
                                       ENBS.mat <- array(NA, dim = c(length.N, 5))
                                       for (j in 1:length.N) {
                                         ENBS.mean <- ENBS.fun(evsi.params[j, 1], Pop, Time, Dis, 
                                                               trial.cost[j, 1])
                                         ENBS.sd <- ENBS.sd.calc(evsi.params[j, 2], Pop, Time, 
                                                                 Dis, trial.cost[j, 2])
                                         ENBS.mat[j, ] <- qnorm(prob, ENBS.mean, ENBS.sd)
                                       }
                                       if (length.prob%%2 == 1) {
                                         lwd <- c(1:ceiling(length.prob/2), (ceiling(length.prob/2) - 
                                                                               1):1, 1)
                                         lty <- c(ceiling(length.prob/2):1, 2:ceiling(length.prob/2), 
                                                  1)
                                       }
                                       if (length.prob%%2 == 0) {
                                         lwd <- c(1:(length.prob/2), (length.prob/2):1, 1)
                                         lty <- c((length.prob/2):1, 2:(length.prob/2), 1)
                                       }
                                       plot.new()
                                       if(is.null(ylim)){
                                         plot.window(xlim = c(min(N), max(N)), ylim = c(min(ENBS.mat), 
                                                                                        max(ENBS.mat)))
                                       }
                                       if(!is.null(ylim)){
                                         plot.window(xlim = c(min(N), max(N)), ylim = ylim)
                                       }

                                       title(main = "Expected Net Benefit of Sampling by Sample Size", 
                                             xlab = "Sample Size", ylab = "ENBS")
                                       axis(side = 2)
                                       if (length.N < 15) {
                                         for (l in 1:length.prob) {
                                           points(N, ENBS.mat[, l], pch = 19, lwd = lwd[l])
                                           points(N, ENBS.mat[, l], type = "l", lwd = lwd[l])
                                         }
                                         legend(alt.legend, c(as.character(prob), "ENBS=0"), 
                                                col = c(rep("black", length.prob), "darkgray"), 
                                                lwd = lwd, box.lwd = 0, pch = 19,
                                                box.col = "white", bg = "white")
                                       }
                                       if (length.N >= 15) {
                                         for (l in 1:length.prob) {
                                           points(N, ENBS.mat[, l], type = "l", lwd = lwd[l], 
                                                  lty = lty[l])
                                         }
                                         if(legend){
                                         legend(alt.legend, c(as.character(prob), "ENBS=0"), col = c(rep("black", 
                                                                                                                length.prob), "darkgray"), lwd = lwd, lty = lty, box.lwd = 0, 
                                                box.col = "white", bg = "white")
                                         }
                                       }
                                       abline(h = 0, col = "darkgray", lwd = lwd[length.prob + 
                                                                                      1], lty = lty[length.prob + 1])
                                       
                                       box()
                                       optimal <- optim.ss(evsi, setup, pp, Pop, Time, k = k, 
                                                           Dis = Dis)
                                       axis(side = 1)
                                       a <- 0.04
                                       poi <- (1 + a) * min(ENBS.mat) - a * max(ENBS.mat)
                                       points(c(optimal$SS.I), c(poi, poi), type = "l", col = "red", 
                                              lwd = 3)
                                       points(c(optimal$SS.max), min(poi), pch = 9, col = "red", 
                                              lwd = 3)
}

### Plots the curve of optimal sample size.
### Not in the EVSI package but function structure and argument names match to
### the other functions
coss <-  function (evsi, setup, pp, Pop = 10000, Time = 10, 
                   Dis = 0.035, N = NULL, pos = c("bottomright")) 
{
  alt.legend <- pos
  if (is.numeric(alt.legend) & length(alt.legend) == 2) {
    temp <- ""
    if (alt.legend[2] == 0) 
      temp <- paste0(temp, "bottom")
    else if (alt.legend[2] != 0.5) 
      temp <- paste0(temp, "top")
    if (alt.legend[1] == 1) 
      temp <- paste0(temp, "right")
    else temp <- paste0(temp, "left")
    alt.legend <- temp
    if (length(grep("^((bottom|top)(left|right)|right)$", 
                    temp)) == 0) 
      alt.legend <- FALSE
  }
  if (is.logical(alt.legend)) {
    if (!alt.legend) 
      alt.legend = "topright"
    else alt.legend = "topleft"
  }
  
  if (class(N) != "numeric") {
    N <- evsi$attrib$N
  }
  length.N <- length(N)
  N.select <- array(NA, dim = length.N)
  for (i in 1:length.N) {
    N.select[i] <- which.min((evsi$attrib$N - N[i])^2)
  }
  
  length.k <- length(evsi$attrib$k)
  k <- evsi$attrib$k
  
  type.evsi <- "det"
  evsi.params <- evsi$evsi[N.select, ]
  
  oss.mat <- array(NA, dim = c(length.k, 3))
  for (j in 1:length.k) {
    optimal <- optim.ss(evsi, setup, pp, Pop, Time, k = k[j], 
                        Dis = Dis)
    oss.mat[j, 1] <- optimal$SS.max
    oss.mat[j, 2:3] <- optimal$SS.I
  }
  
  
  plot.new()
  plot.window(xlim = c(min(k), max(k)), ylim = c(min(oss.mat), 
                                                 max(oss.mat)))
  title(main = "Curve of Optimal Sample Size", 
        xlab = "Willingness to Pay", ylab = "Optimal Sample Size")
  axis(side = 2)
  
  points(k, oss.mat[, 1], type = "l", lwd = 2)
  points(k, oss.mat[, 2], type = "l", lwd = 1, lty = 2)
  points(k, oss.mat[, 3], type = "l", lwd = 1, lty = 2)
  box()
  
  axis(side = 1)
  
}

