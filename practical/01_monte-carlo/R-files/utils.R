# Utility functions to produce summary statistics like BUGS would do

bugs.stats=function (x, na.rm = FALSE) 
{
    if (na.rm == FALSE) {
        c(mean(x), sd(x), quantile(x, 0.025), median(x), quantile(x, 
            0.975))
    }
    else {
        c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE), quantile(x, 
            0.025, na.rm = TRUE), median(x, na.rm = TRUE), quantile(x, 
            0.975, na.rm = TRUE))
    }
}

stats=function (x, dim = 2, na.rm = FALSE) 
{
    if (is.null(dim(x)) == TRUE) {
        tab <- bugs.stats(x, na.rm)
        names(tab) <- c("mean", "sd", "2.5%", "median", "97.5%")
    }
    if (is.null(dim(x)) == FALSE) {
        tab <- t(apply(x, dim, function(x) bugs.stats(x, na.rm)))
        colnames(tab) <- c("mean", "sd", "2.5%", "median", "97.5%")
    }
    list(tab = tab)
}
