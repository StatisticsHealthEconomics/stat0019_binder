# This code replicates the plots in the lecture slides for lecture 6
library(bmhe)

# Loads data list
load(here::here("practical","06_nma","smoke.Rdata"))
# Treatment names
tnames <- c("A: None","B: Self-help","C: Individual","D: Group")

# Format data for table in the slides
rc <- smoke.list$r
nc <- smoke.list$n
dc <- matrix(paste(rc, nc, sep="/"), nrow=smoke.list$NS, ncol=smoke.list$NT)
dc[dc=="NA/NA"] <- ""

# FIXED EFFECTS MODEL
# Initial values
inits <- list(list(mu=rep(0,24), d=c(NA,0,0,0)),
              list(mu=rep(-1,24), d=c(NA,1,1,1)))
res <- bugs(
  model=here::here("practical","06_nma","smokefix_model.txt"), 
  data=smoke.list, inits=inits,
  parameters=c("d","or","L","pq"),
  n.chains=2, n.burnin=1000, n.iter=20000
)
st <- res$summary  # Gives posterior median or mean and 95% CIs for odds ratios, for plotting

#  Organise results for plotting
or <- st[grep("or", rownames(st)),c("2.5%", "50%", "97.5%")]
colnames(or) <- c("l95","est","u95")
or <- as.data.frame(or)
or$act <- as.numeric(sub("or\\[([0-9]+),([0-9]+)\\]", "\\1", rownames(or)))
or$com <- as.numeric(sub("or\\[([0-9]+),([0-9]+)\\]", "\\2", rownames(or)))
or <- or[or$act>or$com,]
or <- or[order(or$com,or$act),]
 

# RANDOM EFFECTS MODEL.  Check convergence of random effects SD in particular
inits <- list(list(mu=rep(0,24), d=c(NA,0,0,0), sd=1),
             list(mu=rep(-1,24), d=c(NA,1,1,1), sd=2))

res2 <- bugs(
  model=here::here("practical","06_nma","smokere_model.txt"), 
  data=smoke.list, inits=inits,
  parameters=c("or", "d", "sd", "pq", "L"),
  n.chains=2, n.burnin=1000, n.iter=20000
)
st <- res2$summary

#  Organise results for plotting
orre <- st[grep("or", rownames(st)),c("2.5%", "50%", "97.5%")]
colnames(orre) <- c("l95","est","u95")
orre <- as.data.frame(orre)
orre$act <- as.numeric(sub("or\\[([0-9]+),([0-9]+)\\]", "\\1", rownames(orre)))
orre$com <- as.numeric(sub("or\\[([0-9]+),([0-9]+)\\]", "\\2", rownames(orre)))
orre <- orre[orre$act>orre$com,]
orre <- orre[order(orre$com,orre$act),]
 
#  Calculate direct evidence - classical ORs and CIs.
#  These exist for all comparisons here
fe <- dati <- ntrials <- numeric()
datij <- array(dim=c(4,4,4)); dimnames(datij) <- list(c("r1","r2","n1","n2"),NULL,NULL)
for (i in 1:3) { # comparator
   for (j in (i+1):4){ # active
       pw <- !is.na(rc[,i]) & !is.na(rc[,j])
       ntrials <- c(ntrials, sum(pw))
       r1 <- sum(rc[pw,i]); n1 <- sum(nc[pw,i])
       r2 <- sum(rc[pw,j]); n2 <- sum(nc[pw,j])
       lori <- log(r2) - log(n2-r2) - log(r1) + log(n1-r1)
       sei <- sqrt(1/r2 - 1/(n2-r2) + 1/r1 - 1/(n1-r1))
       fe <- rbind(fe, c(j, i, exp(lori+qnorm(c(0.025, 0.5, 0.975))*sei)))
       colnames(fe) <- c("act","com","l95","est","u95")
       datij[,i,j] <- c(r1,r2,n1,n2)
       dati <- rbind(dati, c(r1,r2,n1,n2))
   }
}
#  PLOT:  direct against fixed effects mixed
par(lwd=2, mar=c(6.1, 0.1, 0.1, 2.1))
plot(or[,"est"], 1:6, type="n", xlim=c(-2, 5), ylim=c(0,6.5), ylab="", bty="n", xlab="Odds ratio", axes=FALSE)
points(or[,"est"], 1:6, pch=19, col="red")
axis(1, at=0:5)
abline(v=c(1), col="lightgray", lwd=1)
segments(or[,"l95"], 1:6, or[,"u95"], 1:6, col="red")
points(fe[,"est"], 1:6-0.2, pch=19)
segments(fe[,"l95"], 1:6-0.2, fe[,"u95"], 1:6-0.2)
text(-0.2,6.5,"Pooled N", pos=4)
text(-2, 1:6, paste(tnames[fe[,"act"]],tnames[fe[,"com"]],sep=" / "), pos=4, cex=0.8)
N <- dati[,3]+dati[,4]
text(-0.1, 1:6, N, pos=4)
legend("bottomright", col=c("red","black"), c("Fixed effects mixed","Fixed effects direct"), lwd=c(2,2,2), bty="n")

hres <- numeric()
for (i in 1:3) { # comparator
   for (j in (i+1):4){ # active
       pw <- !is.na(rc[,i]) & !is.na(rc[,j])
       r1 <- rc[pw,i]; n1 <- nc[pw,i]
       r2 <- rc[pw,j]; n2 <- nc[pw,j]
       lori <- log(r2) - log(n2-r2) - log(r1) + log(n1-r1)
       sei <- sqrt(1/r2 - 1/(n2-r2) + 1/r1 - 1/(n1-r1))
       hresi <- exp(cbind(lori, lori+qnorm(0.025)*sei, lori+qnorm(0.975)*sei))
       hres <- rbind(hres, cbind(j, i, hresi))
   }
}
colnames(hres)[3:5] <- c("or","l95","u95")
hres[c(7,20),"or"] <- -100
hres[c(7,20),"l95"] <- 1
hres[c(7,20),"u95"] <- 1
 
par(bg="white", xlog=TRUE)
comps <- paste(tnames[hres[,"j"]],tnames[hres[,"i"]],sep=" / ")
ypos <- 1:nrow(hres) + rep(1:6*2, table(comps))
plot(0, ylim=range(ypos), xlim=c(0.1, 10), axes=FALSE, bty="n", ylab="", xlab="Odds ratio (log scale axis)", type="n", log="x",cex.lab=1.3)
apos <- cumsum(table(comps)+2)+1
del <- 0.1
lrec <- 0.1; urec <- 10
rect(lrec, 0+del, urec, apos[1], col="lightgray", border=NA)
rect(lrec, apos[2]+del, urec, apos[3], col="lightgray", border=NA)
rect(lrec, apos[4]+del, urec, apos[5], col="lightgray", border=NA)
rect(lrec, apos[1]+del, urec, apos[2], col="lightblue", border=NA)
rect(lrec, apos[3]+del, urec, apos[4], col="lightblue", border=NA)
rect(lrec, apos[5]+del, urec, apos[6], col="lightblue", border=NA)
points(hres[,"or"], ypos, pch=19)
off.inds <- c(7,8,9,20)
hres[off.inds,"u95"] <- 10.5
segments(hres[,"l95"], ypos, hres[,"u95"], ypos)
rinf <- hres[off.inds,]
arrows(rinf[,"l95"], ypos[off.inds], rinf[,"u95"], ypos[off.inds], length=0.1)
axis(1, at=c(0.2, 0.5, 0:5, 10))
abline(v=1, col="purple")
amid <- c(0, apos[1:5]) + diff(c(0,apos))/2
text(0.1, amid, unique(comps), pos=4, cex=1.4)

par(lwd=2, mar=c(6.1, 0.1, 0.1, 2.1))
plot(or[,"est"], 1:6, type="n", xlim=c(-2, 5), ylim=c(0,6.5), ylab="", bty="n", xlab="Odds ratio", axes=FALSE, cex.lab=1.2)
points(or[,"est"], 1:6, pch=19, col="red")
axis(1, at=0:5)
abline(v=c(1), col="lightgray", lwd=1)
segments(or[,"l95"], 1:6, or[,"u95"], 1:6, col="red")
points(fe[,"est"], 1:6-0.2, pch=19)
segments(fe[,"l95"], 1:6-0.2, fe[,"u95"], 1:6-0.2)
points(orre[,"est"], 1:6+0.2, pch=19, col="blue")
segments(orre[,"l95"], 1:6+0.2, orre[,"u95"], 1:6+0.2, col="blue")
arrows(orre[,"l95"], 3.2, 5.3, 3.2, col="blue", length=0.1)
text(-0.2,6.5,"Pooled N", pos=4, cex=1.1)
text(-2, 1:6, paste(tnames[fe[,"act"]],tnames[fe[,"com"]],sep=" / "), pos=4, cex=1.1)
N <- dati[,3]+dati[,4]
text(-0.1, 1:6, N, pos=4,cex=1.1)
legend("bottomright", col=c("blue","red","black"), c("Random effects mixed","Fixed effects mixed","Fixed effects direct"), lwd=c(2,2,2), bty="n")
text(x=orre[,"est"], y=1:6+0.29,  col="blue",
    labels=paste(round(orre[,"est"],2), " (",round(orre[,"l95"],2), ", ", round(orre[,"u95"],2), ")", sep=""), cex=.8)
text(x=2.2, y=6.29, col="blue", labels="Posterior median (95\\% credible interval)", pos=4, cex=.8)
