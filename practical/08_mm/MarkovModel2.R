# Loads 'dplyr', which is part of the 'tidyverse' --- this is very useful for data "wrangling" and preparation
# If you don't have 'dplyr' installed on your computer, you need to run 'install.packages("tidyverse")' to install
# the full set of packages included in the tidyverse environment. Or you can simply run 'install.packages("dplyr")'
# and 'install.packages("ggplot2")' (which you need to produce the graphs, later in the practical)
library(dplyr)

# Set the working directory, for example something like 
setwd(here::here("practical","08_mm"))

# Now loads the data. Notice the 'pipe' operator '%>%' (https://cran.r-project.org/web/packages/magrittr/vignettes/magrittr.html) 
data=read.table("data.txt",header=T,sep="\t") %>% 
  # Structure the data as a "tibble" for easier visualisation and manipulation
  as_tibble() %>%
  # Now create new variables in the dataset
  mutate(
    # Create a patient ID for convenience 
    patid=row_number(),
    # Convert times from months to years --- you may do the analysis on this scale
    prog_ty=prog_t/12,
    death_ty=death_t/12,
    progdeath_ty=progdeath_t/12,
    # Create competing risk indicator
    # 0=did not die or have a progression, 1=had a progression, 2=died without progression
    crstatus=case_when(
      prog==1~1,
      (prog==0 & death)==1~2,
      TRUE~0
    )
  )

# Visualise the data
data

# Now modify the dataset to create the relevant subsets for the analysis
# 'tidyverse' code to re-arrange the data in "multistate" format
msmdata=
  # Transition Pre to Post
  data %>% mutate(                                                        # use the original data and start making changes to them
    id=patid,                                                               # patient ID
    from=1,                                                                 # starting state (1="Pre-progression")
    to=2,                                                                   # arriving state (2="Progressed")
    trans=1,                                                                # transition ID (1="Pre-progression -> Progressed")
    Tstart=0,                                                               # entry time
    Tstop=prog_t,                                                           # exit time (time at which the event happens)
    time=Tstop-Tstart,                                                      # observed time
    status=case_when(                                                       # event indicator
      prog==1~1,                                                            #   if progressed then 1
      TRUE~0                                                                #   otherwise 0 (censored for progression)
    ),
    treat=treat                                                             # treatment arm
  ) %>% select(id,from,to,trans,Tstart,Tstop,time,status,treat) %>%         # selects only the relevant columns (for simplicity)
  bind_rows(                                                              # stack these new rows below those selected above
    # Transition Pre to Death
    data %>% mutate(                                                      # use the original data and start making changes to them
      id=patid,                                                           # patient ID
      from=1,                                                             # starting state (1="Pre-progression")
      to=3,                                                               # arriving state (3="Death")
      trans=2,                                                            # transition ID (2="Pre-progression -> Death")
      Tstart=0,                                                           # entry time
      Tstop=prog_t,                                                       # exit time (time at which the event happens)
      time=Tstop-Tstart,                                                  # observed time
      status=case_when(                                                   # event indicator
        (death==1 & prog_t==death_t)~1,                                   #   if death then 1
        TRUE~0                                                            #   otherwise 0 (censored for death)
      ),
      treat=treat                                                         # treatment arm
    ) %>% select(id,from,to,trans,Tstart,Tstop,time,status,treat)         # selects only the relevant columns (for simplicity)
  ) %>% 
  bind_rows(                                                              # stack these new rows below those selected above
    # Transition Post to Death
    data %>% filter(prog==1) %>% mutate(                                  # use the original data, but **filter only those who have progressed**
      id=patid,                                                           # patient ID
      from=2,                                                             # starting state (2="Progressed")
      to=3,                                                               # arriving state (3="Death")
      trans=3,                                                            # transition ID (3="Progressed -> Death")
      Tstart=prog_t,                                                      # entry time. NB: this time is the time of progression!
      Tstop=death_t,                                                      # exit time (time at which the event happens). NB: this time it's death!
      time=Tstop-Tstart,                                                  # observed time
      status=case_when(                                                   # event indicator
        death==1~1,                                                       #   if death then 1
        TRUE~0                                                            #   otherwise 0 (censored for death **after progression**)
      ),
      treat=treat                                                         # treatment arm
    ) %>% select(id,from,to,trans,Tstart,Tstop,time,status,treat)         # selects only the relevant columns (for simplicity)
  ) %>% arrange(id,trans)

# Visualise the data
msmdata

# Now runs the survival analysis for the three subsets. You need 'survHE' to do this. 
# Loads survHE
library(survHE)

# Runs survival models on the specific subsets to obtain estimate of the various transition probabilities
m_12=fit.models(Surv(time,status)~as.factor(treat),              # model 'formula': defines the time and censoring indicator and the covariates
                data=msmdata %>% filter(trans==1),               # subsets the msmdata by filtering transition number 1 (pre-progression->progressed)
                distr="gom"#,                                    # selects the Gompertz model
# You can run the Bayesian version of this model by uncommenting these two lines + installing 'survHEhmc'
#                method="hmc",                                    # instructs R to use HMC/Bayesian modelling
#                priors=list(gom=list(a_alpha=1.5,b_alpha=1.5))   # specifies the informative prior
)

m_13=fit.models(Surv(time,status)~as.factor(treat),              # model 'formula': defines the time and censoring indicator and the covariates
                data=msmdata %>% filter(trans==2),               # subsets the msmdata by filtering transition number 2 (pre-progression->death)
                distr="gom"#,                                    # selects the Gompertz model
# You can run the Bayesian version of this model by uncommenting these two lines + installing 'survHEhmc'
#                method="hmc",                                    # instructs R to use HMC/Bayesian modelling
#                priors=list(gom=list(a_alpha=1.5,b_alpha=1.5))   # specifies the informative prior
)

m_23=fit.models(Surv(time,status)~as.factor(treat),              # model 'formula': defines the time and censoring indicator and the covariates
                data=msmdata %>% filter(trans==3),               # subsets the msmdata by filtering transition number 3 (progressed->death)
                distr="gom"#,                                    # selects the Gompertz model
# You can run the Bayesian version of this model by uncommenting these two lines + installing 'survHEhmc'
#                method="hmc",                                    # instructs R to use HMC/Bayesian modelling
#                priors=list(gom=list(a_alpha=1.5,b_alpha=1.5))   # specifies the informative prior
)

# Now run the function 'three_state_mm', which in turns calls the function 'make.transition.probs' (which computes
# the transition probabilities from the survival curves) and then 'make_state_occupancy' (to compute the number of 
# people in each state at each time point in the follow up). These functions are "hidden" inside 'survHE', so you
# need to use the ":::" notation to access them!
mm=survHE:::three_state_mm(
  m_12,m_13,m_23,         # these are the three objects containing the parameters estimates
  t=seq(0,130),           # specifies that the Markov model needs to be run for discrete times from 0 to 130 (months)
  nsim=1,                 # only uses 1 simulation from the distribution of the model parameters (the mean)
  start=c(1000,0,0)       # initial population: 1000 in pre-progression, 0 in progressed and 0 in death
)

# Visualises the results of the Markov model in terms of state occupancy & transition probabilities
mm

# Of course, in reality you would run the analysis for, say, 1000 simulations to fully characterise uncertainty in the 
# underlying parameters and then the Markov model.

# Finally, can visualise the Markov trace, using the specialised function 'markov_trace'
markov_trace(mm)
