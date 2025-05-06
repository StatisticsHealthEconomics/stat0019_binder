#' Smoking cessation network meta-analysis data in format obtained
#' from Lu & Ades tutorial "Introduction to Mixed Treatment Comparisons"
library(tidyverse)
smoking_data=readRDS(here::here("06_nma/smoke.rds"))

#' Loads R2jags and creates the fixed effect model as an R function
library(R2jags)

smoking_fe=function(){
  for (i in 1:N) {
    r[i] ~ dbin(pi[i],n[i])
    logit(pi[i]) <- alpha[s[i]] + delta[s[i],t[i]]
    delta[s[i],t[i]] <- d[t[i]] - d[b[i]]
  }
  
# Prior for baseline log-odds
  for (i in 1:NS){
    alpha[i] ~ dnorm(0,0.01)
  }
# Prior for the incremental effects (log-OR)  
  d[1] <- 0  # log odds ratio compared to treatment 1 (e.g. placebo)
  for (i in 2:NT) {
    d[i] ~ dnorm(0,1)
  }

# Maps odds ratios for *all* treatment comparisons
  for (c in 1:(NT-1)) {
    for (k in (c+1):NT)  {
      or[c,k] <- exp(d[c] - d[k])
      or[k,c] <- 1/or[c,k]
    }
  }
# Log odds of quitting successfully under no intervention (from published data)
  rho ~ dnorm(-2.6, 6.925208) # = SD 0.38
# Absolute probability of quitting successfully under each intervention
  for (i in 1:NT) {
    logit(p[i]) <- rho + d[i]
  }
}

#' Creates the data list
data.list=list(
  t=smoking_data$t,s=smoking_data$s,r=smoking_data$r,
  n=smoking_data$n,b=smoking_data$b,
  # Number of rows in the data
  N=smoking_data |> nrow(),
  # Total number of treatment arms
  NT=smoking_data$t |> max(),
  # Total number of studies
  NS=smoking_data$s |> max()
)

#' Runs model
m_fe=jags(
  data=data.list,parameters.to.save=c("d","or","p"),n.thin=10,
  inits=NULL,n.chains=2,n.burnin=5000,n.iter=25000,model.file=smoking_fe
)

#' Results & diagnostics
print(m_fe)
bmhe::traceplot(m_fe)
bmhe::acfplot(m_fe)

#' Analysis of heterogeneity
#' First creates a function to select the relevant data to do specialised
#' analyses. The function takes two possible inputs. When type='direct',
#' then does pooled ORs for *only direct pairwise comparisons*. When 
#' type='heterogeneity', then does *individual estimates* for ORs 
#' (study-by-study). This function is used to select the relevant data so that
#' the two types of models can be run
make_data_no_pooling=function(type="direct"){
  # Creates a tibble with the relevant pairwise comparisons
  comparisons=tibble(
    A=c(2,3,3,4,4,4),B=c(1,1,2,1,2,3)
  )
  ylabs=c(
    "B: Self-help / A: None","C: Individual / A: None",
    "C: Individual / B: Self-help","D: Group / A: None",
    "D: Group / B: Self-help","D: Group / C: Individual"
  )

  # Initialise the output tibble
  out=list()
  for (i in 1:nrow(comparisons)) {
    A=comparisons$A[i]
    B=comparisons$B[i]
    
    data_direct=smoking_data |> group_by(s) |> 
      mutate(comparison=paste(t,collapse=",")) |> ungroup() |> 
      dplyr::filter(
        (grepl(A,comparison)),(grepl(B,comparison)),
        t %in% c(A,B)
      ) |>
      select(r,n,t,b,s) |> group_by(s) |> mutate(ss=cur_group_id()) |> 
      ungroup() |> select(-s) |> rename(s=ss) |> as.list()
    data_direct$N=length(data_direct$r)
    data_direct$NS=max(data_direct$s)
    data_direct$Baseline=B
    data_direct$Active=A
    # Runs the model
    if (type=="direct") {
      m=jags(
        data=data_direct,parameters.to.save=c("or"),n.thin=4,
        inits=NULL,n.chains=2,n.burnin=1000,n.iter=10000,
        model.file=direct_evidence
      )
    }
    if (type=="heterogeneity") {
      m=jags(
        data=data_direct,parameters.to.save=c("or"),n.thin=4,
        inits=NULL,n.chains=2,n.burnin=1000,n.iter=10000,
        model.file=heterogeneity
      )
    }
    
    # Saves the output
    out[[i]]=m$BUGSoutput$summary[,c("mean","sd","2.5%","97.5%")] |>
      as_tibble() |> mutate(
        Parameter=m$BUGSoutput$summary[,c("mean","sd","2.5%","97.5%")] |> rownames()
      ) |> dplyr::filter(Parameter!="deviance")
    
    if (type=="heterogeneity"){
      out[[i]]=out[[i]] |> mutate(
        comparison=ylabs[i],
        N=data_direct |> as_tibble() |> select(n,s) |> group_by(s) |> 
          summarise(N=sum(n)) |> ungroup() |> pull(N),
        s=smoking_data |> group_by(s) |> 
          mutate(comparison=paste(t,collapse=",")) |> ungroup() |> 
          dplyr::filter(
            (grepl(A,comparison)),(grepl(B,comparison)),
            t %in% c(A,B)
            ) |> select(s) |> unique() |> pull(s)
      )
    }
  }
  if (type=="direct") {
    out=out |> bind_rows() |> mutate(
      comparison=ylabs[i],
      Parameter=paste0("or[",comparisons$A,",",comparisons$B,"]"),
      num=row_number()
    ) |> select(Parameter,mean,sd,"2.5%","97.5%",num)
  } 
  if (type=="heterogeneity") {
    out=out |> bind_rows()
  }
  return(out)
}

#' Now defines the model code for pooling over studies doing 
#' direct pairwise comparisons
direct_evidence=function() {
  for (i in 1:N) {
    r[i] ~ dbin(pi[i],n[i])
    logit(pi[i]) <- alpha[s[i]] + delta[s[i],t[i]]
    delta[s[i],t[i]] <- d[t[i]]-d[b[i]]
  }
  for (i in 1:NS) {
    alpha[i] ~ dnorm(0,0.01)
  }
  d[1] <- 0  # log odds ratio compared to treatment 1 (e.g. placebo)
  for (i in 2:4) {
    d[i] ~ dnorm(0,1)
  }
  or <- exp(d[Active] - d[Baseline])
}


#' Make graphs for the direct evidence analysis
#' 
#' First defines the comparisons (ORs) to select for plotting 
select=c("or[2,1]","or[3,1]","or[3,2]","or[4,1]","or[4,2]","or[4,3]")
#' Then runs 'bmhe::coefplot' to extract the relevant data from the 
#' fixed effect model ('m_fe')
toplot=bmhe::coefplot(m_fe,parameter=select)$data

#' Now selects the labels to add to the y-axis
ylabs=c(
  "B: Self-help / A: None (N=2867)",
  "C: Individual / A: None (N=12846)",
  "C: Individual / B: Self-help (N=255)",
  "D: Group / A: None (N=318)",
  "D: Group / B: Self-help (N=441)",
  "D: Group / C: Individual (N=764)"
)

#' Generate the data for the 'direct' comparison using the 
#' make_data_no_pooling' function defined above
out1=make_data_no_pooling()
#'
#' You should check the output 'out1' to see what it looks like... 
out1
#'
#' Now adds the relevant data to the 'toplot' object to combine all the results
#' This creates a 'coefplot' with the comparison between the fixed effect and 
#' the "direct evidence only" models
toplot |> mutate(model="Fixed effect (indirect)",Parameter=ylabs) |> 
  bind_rows(out1 |> mutate(model="No pooling effect (direct)",Parameter=ylabs)) |> 
  ggplot(aes(mean,Parameter)) + 
  geom_point(aes(col=model),position=position_dodge2(width=c(.25))) +
  geom_linerange(
    aes(xmin =`2.5%`, xmax =`97.5%`,col=model),
    position=position_dodge2(width=c(.25))
  ) + geom_vline(xintercept=1, col="gray50") + 
  labs(x="Odds Ratio") +  
  theme(
    legend.position="inside",legend.position.inside=c(.75,.1),
    legend.background=element_blank()
  ) + labs(col="") + ylab("") + xlim(0,6) +
  scale_color_manual(values=c("black","red"))


#' Now defines the model code for no-pooling individual studies
heterogeneity=function() {
  for (i in 1:N) {
    r[i] ~ dbin(pi[i],n[i])
    logit(pi[i]) <- alpha[s[i]] + delta[s[i],t[i]]
    delta[s[i],t[i]] <- d[s[i],t[i]] - d[s[i],b[i]]
  }
  for (i in 1:NS) {
    alpha[i] ~ dnorm(0,0.01)
    d[i,1] <- 0
    for (j in 2:4) {
      d[i,j] ~ dnorm(0,1)
    }
    or[i] <- exp(d[i,Active] - d[i,Baseline])
  }
}


#' Generate the data for the 'heterogeneity' comparison using the 
#' make_data_no_pooling' function defined above
out2=make_data_no_pooling(type="heterogeneity") 
#'
#' Check this out...
out2

#' Now creates a 'coefplot' with all the individual studies estimates, 
#' grouped by comparison using 'ggplot' "facets"
out2 |> mutate(Parameter=as.factor(paste0("s=",s," N=",N))) |>
  # mutate(id=row_number()) |> arrange(comparison,desc(id)) |> 
  # mutate(ttt=paste0(id,": N=",N)) |> mutate(Parameter=as.factor(ttt)) |>
  ggplot(aes(mean,fct_reorder(Parameter,s,.desc=TRUE))) +
  geom_linerange(
    aes(xmin =`2.5%`, xmax =`97.5%`)
  ) + geom_vline(xintercept=1, linetype = "dashed") + 
  geom_point() +
  labs(x="Odds Ratio (log scale axis)") + 
  facet_grid(comparison~.,scales="free",space="free") +
  ylab("") +
  theme(strip.text.y = element_text(angle = 0)) +
 scale_x_continuous(
    trans='log',labels=c(.2,.5,1,2,5,10,25),
    breaks=c(.2,.5,1,2,5,10,25)
  ) 


#' Now defines the model with random effects. In line with the lecture slides,
#' we select a Half-Cauchy prior for the standard deviation of the random 
#' effects. This is coded up using the property that 
#' HC(mu,lambda) = t(mu, lambda^(-2), 1) I(0,\infty)
#' ie a *Truncated* t with 1 degree of freedom is equivalent to a Half-Cauchy.
#' NB: We code the truncation in JAGS using the '%_%T(0,)' notation
smoking_re=function(){
  for (i in 1:N) {
    r[i] ~ dbin(pi[i],n[i])
    logit(pi[i]) <- alpha[s[i]] + delta[s[i],t[i]]
    delta[s[i],t[i]] ~ dnorm(mu[s[i],t[i]],tau)
    mu[s[i],t[i]] <- d[t[i]] - d[b[i]]
  }
  
  # Prior for baseline log-odds
  for (i in 1:NS){
    alpha[i] ~ dnorm(0,0.01)
  }
  # Prior for the incremental effects (log-OR)  
  d[1] <- 0  # log odds ratio compared to treatment 1 (e.g. placebo)
  for (i in 2:NT) {
    d[i] ~ dnorm(0,1)
  }
  # Half-Cauchy prior for common sd for the random effects. NB: in JAGS, the 
  # t distribution is indexed in terms of the *precision* (=1/variance)!
  sigma ~ dt(0,0.25,1)%_%T(0,)
  tau <- pow(sigma,-2)
  
  # Maps odds ratios for all treatment comparisons
  for (c in 1:(NT-1)) {
    for (k in (c+1):NT)  {
      or[c,k] <- exp(d[c] - d[k])
      or[k,c] <- 1/or[c,k]
    }
  }
  # Log odds of quitting successfully under no intervention (from published data)
  rho ~ dnorm(-2.6,6.925208) # = SD 0.38
  # Absolute probability of quitting successfully under each intervention
  for (i in 1:NT) {
    logit(p[i]) <- rho + d[i]
  }
}

#' And runs the model
m_re=jags(
  data=data.list,parameters.to.save=c("d","or","p","sigma"),n.thin=10,
  inits=NULL,n.chains=2,n.burnin=5000,n.iter=25000,model.file=smoking_re
)
print(m_re)

#' Diagnostics
bmhe::diagplot(m_re)
bmhe::diagplot(m_re,what="n.eff",label = TRUE) + ylim(0,4250)

#' Compares all the models (fixed effects, no pooling and random effects)
#' using data extracted by the output of the 'bmhe::coefplot' function
#' 
#' Selects the relevant comparisons in terms of ORs
select=c("or[2,1]","or[3,1]","or[3,2]","or[4,1]","or[4,2]","or[4,3]")
#'
#' Extracts the data to plot from 'bmhe::coefplot' using the facilities of
#' 'ggplot' (the resulting output would be a 'ggplot' object...)
toplot_fe=bmhe::coefplot(m_fe,parameter=select)$data
toplot_re=bmhe::coefplot(m_re,parameter=select)$data

#' Then brings everything together and plots the three models...
toplot_fe |> mutate(model="Fixed effect (indirect)",Parameter=ylabs) |> 
  bind_rows(out1 |> mutate(model="No pooling (direct)",Parameter=ylabs)) |> 
  bind_rows(toplot_re |> mutate(model="Random effects",Parameter=ylabs)) |>
  ggplot(aes(mean,Parameter)) + 
  geom_point(aes(col=model),position=position_dodge2(width=c(.35))) +
  geom_linerange(
    aes(xmin =`2.5%`, xmax =`97.5%`,col=model),
    position=position_dodge2(width=c(.35))
  ) + geom_vline(xintercept=1, col="gray50") + 
  labs(x="Odds Ratio") +  
  theme(
    legend.position=c(.75,.12),legend.background=element_blank()
  ) + labs(col="") + ylab("") + xlim(0,6) +
  scale_color_manual(values=c("black","red","blue"))


#' Now compares the estimates in terms of the absolute probability of 
#' quitting smoking in the fixed and random effects models
#' 
#' First extracts the data, this time only for the parameter 'p'
toplot_fe=bmhe::coefplot(m_fe,parameter = "p")$data
toplot_re=bmhe::coefplot(m_re,parameter = "p")$data

#' And then combines them all together to make the plots
toplot=toplot_fe |> mutate(model="Fixed effects") |> 
  bind_rows(toplot_re |> mutate(model="Random effects")) |> 
  mutate(
    Parameter=stringr::str_replace_all(Parameter,"p\\[1\\]","No intervention"),
    Parameter=stringr::str_replace_all(Parameter,"p\\[2\\]","Self-help"),
    Parameter=stringr::str_replace_all(Parameter,"p\\[3\\]","Individual counselling"),
    Parameter=stringr::str_replace_all(Parameter,"p\\[4\\]","Group counselling"),
    Parameter=factor(Parameter,levels=c("No intervention","Self-help","Individual counselling","Group counselling"))
  ) 
#' Makes the plot
toplot |> ggplot(aes(mean,Parameter)) + 
  geom_point(aes(col=model),position=position_dodge2(width=c(.25))) +
  geom_linerange(
    aes(xmin =`2.5%`, xmax =`97.5%`,col=model),
    position=position_dodge2(width=c(.25))
  ) + labs(x="Probablity of quitting smoking") +  
  theme(
    legend.position=c(.75,.15),legend.background=element_blank()
  ) + labs(col="") + ylab("") +
  annotate(
    "text",toplot|> dplyr::filter(model=="Fixed effects") |> pull(`2.5%`),y=1:4,
    label=round(toplot|> dplyr::filter(model=="Fixed effects") |> pull(`2.5%`),3),
    vjust=2.25
  ) +
  annotate(
    "text",toplot|> dplyr::filter(model=="Fixed effects") |> pull(`97.5%`),y=1:4,
    label=round(toplot|> dplyr::filter(model=="Fixed effects") |> pull(`97.5%`),3),
    vjust=2.25
  ) +
  annotate(
    "text",toplot|> dplyr::filter(model=="Fixed effects") |> pull(mean),y=1:4,
    label=round(toplot|> dplyr::filter(model=="Fixed effects") |> pull(mean),3),
    vjust=2.25
  ) +
  annotate(
    "text",toplot|> dplyr::filter(model=="Random effects") |> pull(`2.5%`),y=1:4,
    label=round(toplot|> dplyr::filter(model=="Random effects") |> pull(`2.5%`),3),
    vjust=-1.15
  ) +
  annotate(
    "text",toplot|> dplyr::filter(model=="Random effects") |> pull(`97.5%`),y=1:4,
    label=round(toplot|> dplyr::filter(model=="Random effects") |> pull(`97.5%`),3),
    vjust=-1.15
  ) +
  annotate(
    "text",toplot|> dplyr::filter(model=="Random effects") |> pull(mean),y=1:4,
    label=round(toplot|> dplyr::filter(model=="Random effects") |> pull(mean),3),
    vjust=-1.15
  ) 

#' Cost-effectiveness analysis
#' 
#' Defines fixed unit costs for each intervention
unit.cost=c(0,200,6000,600)
#'
#' Labels for the four interventions
ints=c("No contact","Self help","Individual counselling","Group counselling")
#' Extracts effects and costs from the random effect model
e=c=matrix(NA,m_re$BUGSoutput$n.sims,4)

#' Monte Carlo samples from distribution of life-years gained by quitting
#' Here uses the same number of simulations as for the random effects models
#' and sets mean=15 and sd=4 (see lecture slides)
L=rnorm(m_re$BUGSoutput$n.sims,15,4)

#' Extracts the absolute probability of quitting from the random effect model.
#' This is a matrix with 'm_re$BUGSoutput$n.sims' simulations and 4 columns,
#' one for each intervention
p=m_re$BUGSoutput$sims.list$p 

#' Now computes benefits (L*p) and costs (unit.cost)
for (t in 1:4) {
  e[,t] <- L*p[,t]
  c[,t] <- unit.cost[t]
}
#' And adds column names (interventions)
colnames(e) <- colnames(c) <- ints

#' Runs BCEA to perform the economic evaluation
library(BCEA)
m <- bcea(e,c,interventions=ints,Kmax=1000,ref=4)
summary(m)
plot(m,graph="gg")
