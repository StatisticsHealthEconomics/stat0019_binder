#' make.transition.probs
#' 
#' Computes the transition probabilities (to be passed to a Markov model) from
#' the survival curves obtained using \code{fit.models}, using the formula 
#' p(t)=1-S(t+k)/S(t), where k is the Markov model cycle length and t is a 
#' generic time
#' 
#' @aliases make.transition.probs
#' @param fit an object obtained as output of the call to \code{fit.models}
#' @param ...  additional arguments. Includes the standard inputs to the 
#' call to \code{make.surv}, so \code{mod} (the index of the possibly many
#' models stored in the 'survHE' object), \code{t} (the vector of times 
#' over which to compute the survival curves), \code{newdata} (a list that
#' defines the profile of covariates) and \code{nsim} (the number of 
#' simulations to use - default is \code{nsim}=1)
#' @return A tibble 'lambda' with an indicator for the treatment arm,
#' the times at which the probabilities have been computed and \code{nsim}
#' columns each with a simulation of the transition probabilities for 
#' all the times specified by the user
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Something will go here
#' @keywords Transition probabilities Markov models
#' @examples
#' \dontrun{
#' # Something will go here
#' }
#' 
#' @export make.transition.probs
make.transition.probs <- function(fit,...) {
  exArgs <- list(...)
  
  # Makes default for parameters to the call to 'make.surv' (which are overwritten if the user has specified them
  # separately)
  if(exists("mod",exArgs)) {mod=exArgs$mod} else {mod=1}
  if(exists("t",exArgs)) {t=exArgs$t} else {t=NULL}
  if(is.null(t)) {
    t <- sort(unique(fit$misc$km$time))
    # Add an extra time=0 at the beginning. This ensures the computation can be done for all the actual times
    # specified by the user (because the first one has no lag and so the ratio gives NA...)
    if(t[1]>0) {t=c(0,t)}
  }
  if(exists("newdata",exArgs)) {newdata=exArgs$newdata} else {newdata=NULL}
  if(exists("nsim",exArgs)) {nsim=exArgs$nsim} else {nsim=1}
  # Now computes the simulations using 'make.surv'
  s=make.surv(fit,mod=mod,t=t,newdata=newdata,nsim=nsim)
  
  # Now retrieves the transition probabilities 'lambda' applying the formula
  # lambda(t)=1-S(t+k)/S(t) where k is the MM cycle length and t is a generic time
  lambda=s$mat %>% 
    # First stacks together all the matrices with the simulation(s) for the survival curves
    bind_rows() %>% 
    # Then creates a treatment indicator (1,2,...) & place it as the first column
    mutate(treat=rep(1:length(s$mat),each=nrow(s$mat[[1]]))) %>% select(treat,everything()) %>% 
    # Then group by treatment and computes the probabilities using '1-S/lag(S)'
    group_by(treat) %>% mutate(across(starts_with("S"),~1-./lag(.))) %>% 
    # Then removes the first row (in each treatment) for which the transition probability can't be computed
    # because lag(S) doesn't exist
    slice(-1) %>% ungroup()
  ###########################################################################################################
  # This replaces all the NaNs with a 1 - when you have a NaN it's because the survival curve has gone to 
  # absolutely 0 and so there's no people left to make that transition. This is relevant because with this fix
  # it can extrapolate indefinitely without problems (or removing NaN values). What it means is that 
  # essentially there's a sudden drop in the number of people in a given state, though, which may be a bit weird...
  # NB: If use
  #replace(is.na(.),0) %>% 
  # then both NAs (which happen for the first time in each treatment arm) *and* NaNs (which happen for cases
  # when the survival curve is exactly 0, so that for two consecutive points you get 0/0=NaN) are set to 0
  # which means that the first time (which should be discarded) isn't removed and is artificially set to 0
  ###########################################################################################################
  # Then removes the NAs (which occur if the time vector includes 0) and ungroups
  # If this is removed then can do extrapolation to whatever time --- NaNs will exist (and possibly negative 
  # transition probabilities when one of them hits 1; but this won't matter as when there are NaNs the extrapolation
  # will effectively stop and nothing will be shown. That's beyond the point when everybody has died...)
  #filter(across(starts_with("S"),~!is.na(.))) %>% ungroup() 
  # And now renames the columns from S(_1,S_2,...,S_nsim) to lambda(_1,lambda_2,...,lambda_nsim)
  if (nsim==1) {
    lambda=lambda %>% rename_with(starts_with("S"),.fn=~"lambda")
  } else {
    lambda=lambda %>% rename_with(starts_with("S"),.fn=~paste0("lambda_",1:nsim))
  }
  
  # Returns the object containing the probabilities
  return(lambda)
}

#' three_state_mm
#' 
#' General purpose function to run a standard three-state Markov model
#' (typically used in cancer modelling). The states are typically
#' 'Pre-progression', 'Progressed' and 'Death'. No backward transition
#' from 'Progressed' to 'Pre-progression' is allowed and 'Death' is 
#' obviously an absorbing state. All other transitions are possible.
#' The crucial assumption is that *individual-level data* are available
#' recording an indicator and the time of progression and death for each
#' individual. The function returns the full transition matrix
#' 
#' @aliases three_state_mm
#' @param m_12 A 'survHE' object (output to a call to \code{fit.models})
#' estimating the parameters of a model for the transition from 
#' 'Pre-progression' (state 1) to 'Progressed' (state 2). Given the 
#' individual level data with the complete event history (in the object
#' 'data'), can be done with a call like 'x=make_data_multi_state(data)'
#' and then 'fit.models(Surv(time,status)~...,data=x %>% filter(trans==1),...)'
#' @param m_13 A 'survHE' object (output to a call to \code{fit.models})
#' estimating the parameters of a model for the transition from 
#' 'Pre-progression' (state 1) to 'Death' (state 3).  Given the 
#' individual level data with the complete event history (in the object
#' 'data'), can be done with a call like 'x=make_data_multi_state(data)'
#' and then 'fit.models(Surv(time,status)~...,data=x %>% filter(trans==2),...)'
#' @param m_23 A 'survHE' object (output to a call to \code{fit.models})
#' estimating the parameters of a model for the transition from 
#' 'Progressed' (state 2) to 'Death' (state 3).  Given the 
#' individual level data with the complete event history (in the object
#' 'data'), can be done with a call like 'x=make_data_multi_state(data)'
#' and then 'fit.models(Surv(time,status)~...,data=x %>% filter(trans==3),...)'
#' @param nsim The number of simulations for the model parameters that are 
#' used to compute the survival curves. Defaults to \code{nsim}=1,
#' which simply creates one survival curve for each treatment arm.
#' @param start A vector of initial state occupancy. By default assumes 1000 
#' individuals, all initially allocated to 'Pre-progression'
#' @param basecase Should the base case be computed as well, based on the 
#' point estimate of the underlying model parameters? (Default=FALSE)
#' @param ...  additional arguments. 
#' @return A list including the state occupancy simulations in an object 'm'.
#' This is a tibble with the number of individuals in each of the 3 states 
#' at each of the times specified by the user. If \code{nsim}>1, then the tibble
#' also contains a simulation index to keep track of that. The list also 
#' includes the computation time to obtain the state occupancy tibble (in the
#' object 'running_time'). If \code{basecase==TRUE}, then the function also
#' computes the "base case scenario" (based on 1 simulation from of the 
#' underlying survival curves, i.e. the point estimate of the model parameters)
#' and stores it in the object 'base_case'
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso make.transition.probs make_data_multi_state
#' @references Something will go here
#' @keywords Transition probabilities Markov models Three-state cancer model
#' @examples
#' \dontrun{
#' # Something will go here
#' }
#' 
#' @export three_state_mm
three_state_mm = function(m_12,m_13,m_23,nsim=1,start=c(1000,0,0),basecase=FALSE,...){
  
  exArgs <- list(...)
  
  # Initialises the base_case object
  base_case=NULL
  
  # Makes default for parameters to the call to 'make.surv' (which are overwritten if the user has specified them
  # separately)
  if(exists("mod",exArgs)) {mod=exArgs$mod} else {mod=1}
  if(exists("t",exArgs)) {t=exArgs$t} else {t=NULL}
  if(is.null(t)) {
    t <- sort(unique(m_12$misc$km$time))
  }
  if(exists("newdata",exArgs)) {newdata=exArgs$newdata} else {newdata=NULL}
  
  # Computes the transition probabilities for the transitions that are directly
  # identifiable from the observed data
  lambda_12=make.transition.probs(m_12,mod=mod,t=t,newdata=newdata,nsim=nsim)
  lambda_13=make.transition.probs(m_13,mod=mod,t=t,newdata=newdata,nsim=nsim)
  lambda_23=make.transition.probs(m_23,mod=mod,t=t,newdata=newdata,nsim=nsim)
  
  # Derives lambda_11 by subtraction (as all transition probs out of state 1 must sum to 1)
  lambda_11=(lambda_12 %>% select(starts_with("lambda")) + lambda_13 %>% select(starts_with("lambda"))) %>% 
    as_tibble() %>% bind_cols(lambda_12 %>% select(treat,t)) %>% select(treat,t,everything()) %>% 
    mutate(across(starts_with("lambda"),~1-.))
  # Derives lambda_22 by subtraction (as all transition probs out of state 2 must sum to 1)
  lambda_22=(1-lambda_23 %>% select(starts_with("lambda"))) %>% as_tibble() %>% 
    bind_cols(lambda_23 %>% select(treat,t)) %>% select(treat,t,everything())
  
  # Computes the state occupancy for all the simulations
  tic=Sys.time()
  m=make_state_occupancy(nsim=nsim,lambda_11,lambda_12,lambda_13,lambda_22,lambda_23,start)
  toc=Sys.time()
  running_time=toc-tic
  
  if(basecase) {
    # Makes base-case Markov model (by considering the point estimate of the model parameters)
    # Computes the transition probabilities for the transitions that are directly
    # identifiable from the observed data
    lambda_12=make.transition.probs(m_12,mod=mod,t=t,newdata=newdata,nsim=1)
    lambda_13=make.transition.probs(m_13,mod=mod,t=t,newdata=newdata,nsim=1)
    lambda_23=make.transition.probs(m_23,mod=mod,t=t,newdata=newdata,nsim=1)
    
    # Derives lambda_11 by subtraction (as all transition probs out of state 1 must sum to 1)
    lambda_11=(lambda_12 %>% select(starts_with("lambda"))+lambda_13 %>% select(starts_with("lambda"))) %>% 
      as_tibble() %>% bind_cols(lambda_12 %>% select(treat,t)) %>% select(treat,t,everything()) %>% 
      mutate(across(starts_with("lambda"),~1-.))
    # Derives lambda_22 by subtraction (as all transition probs out of state 2 must sum to 1)
    lambda_22=(1-lambda_23 %>% select(starts_with("lambda"))) %>% as_tibble() %>% 
      bind_cols(lambda_23 %>% select(treat,t)) %>% select(treat,t,everything())
    
    base_case=make_state_occupancy(nsim=1,lambda_11,lambda_12,lambda_13,lambda_22,lambda_23,start)
  }
  
  # Outputs of the function
  list(m=m,running_time=running_time,base_case=base_case)
}


#' make_state_occupancy
#' 
#' Utility function to compute the state occupancy in a three state MM
#' 
#' @param nsim The number of simulations for the model parameters that are 
#' used to compute the survival curves. 
#' @param lambda_11 the tibble containing the transition probabilities for
#' the transition "Pre-progression -> Pre-progression" 
#' @param lambda_12 the tibble containing the transition probabilities for
#' the transition "Pre-progression -> Progression" 
#' @param lambda_13 the tibble containing the transition probabilities for
#' the transition "Pre-progression -> Death" 
#' @param lambda_22 the tibble containing the transition probabilities for
#' the transition "Progression -> Progression" 
#' @param lambda_23 the tibble containing the transition probabilities for
#' the transition "Progression -> Death" 
#' @return A list including the state occupancy simulations in an object 'm'.
#' This is a tibble with the number of individuals in each of the 3 states 
#' at each of the times specified by the user. If \code{nsim}>1, then the tibble
#' also contains a simulation index to keep track of that.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso make.transition.probs make_data_multi_state
#' @references Something will go here
#' @keywords Transition probabilities Markov models Three-state cancer model
#' @examples
#' \dontrun{
#' # Something will go here
#' }
#' 
#' @noRd 
make_state_occupancy=function(nsim,lambda_11,lambda_12,lambda_13,lambda_22,lambda_23,start) {
  
  m=list()
  # Initialises the lists
  m=lapply(1:nsim,function(i) {
    # Creates the tibbles (one per each of the nsim simulations)
    m[[i]]=tibble(treat=lambda_11$treat,t=lambda_11$t,`Pre-progressed`=NA,Progressed=NA,Death=NA)
    # Now adds in the relevant transition probabilities in the correct rows
    m[[i]]=m[[i]] %>% left_join(lambda_11 %>% select(treat,t,starts_with("lambda")[[i]]) %>% rename("lambda_11"=starts_with("lambda")), by=c("treat","t")) %>% 
      left_join(lambda_12 %>% select(treat,t,starts_with("lambda")[[i]]) %>% rename("lambda_12"=starts_with("lambda")), by=c("treat","t")) %>% 
      left_join(lambda_13 %>% select(treat,t,starts_with("lambda")[[i]]) %>% rename("lambda_13"=starts_with("lambda")), by=c("treat","t")) %>% 
      left_join(lambda_22 %>% select(treat,t,starts_with("lambda")[[i]]) %>% rename("lambda_22"=starts_with("lambda")), by=c("treat","t")) %>% 
      left_join(lambda_23 %>% select(treat,t,starts_with("lambda")[[i]]) %>% rename("lambda_23"=starts_with("lambda")), by=c("treat","t"))
    # Initialise the tibbles with the start values
    m[[i]]=m[[i]] %>% group_by(treat) %>% mutate(
      `Pre-progressed`=replace(`Pre-progressed`,row_number()==1,start[1]),
      Progressed=replace(Progressed,row_number()==1,start[2]),
      Death=replace(Death,row_number()==1,start[3]),
    ) %>% ungroup()
  })
  
  # Convert the list into a massive tibble with a simulation index
  m=m %>% bind_rows() %>% mutate(sim_idx=rep(1:nsim,each=nrow(m[[1]])))
  
  # Loops over times to compute the state occupancy
  for (j in 2:nrow(m)) {
    # This is a simple trick to restart the computation when the treatment indicator changes
    if (m$treat[j]==m$treat[j-1]) {
      m$`Pre-progressed`[j]=m$`Pre-progressed`[j-1]*m$lambda_11[j-1]
      m$Progressed[j]=m$`Pre-progressed`[j-1]*m$lambda_12[j-1] + m$Progressed[j-1]*m$lambda_22[j-1]
      m$Death[j]=m$`Pre-progressed`[j-1]*m$lambda_13[j-1] + m$Progressed[j-1]*m$lambda_23[j-1] + m$Death[j-1]
    }
  }
  m=m %>% select(treat,t,`Pre-progressed`,Progressed,Death,sim_idx,everything())
  # Need to set to 0 all negative values!
  m=m %>% mutate(
    `Pre-progressed`=if_else(`Pre-progressed`<0,0,`Pre-progressed`),
    Progressed=if_else(Progressed<0,0,Progressed),
    Death=if_else(Death<0,0,Death)
  )
  
  return(m)
}


#' markov_trace
#' 
#' Computes the transition probabilities (to be passed to a Markov model) from
#' the survival curves obtained using \code{fit.models}, using the formula 
#' p(t)=1-S(t+k)/S(t), where k is the Markov model cycle length and t is a 
#' generic time
#' 
#' @aliases markov_trace
#' @param mm 
#' @interventions 
#' @param ...  additional arguments. 
#' @return 
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Something will go here
#' @keywords Transition probabilities Markov models Markov trace
#' @examples
#' \dontrun{
#' # Something will go here
#' }
#' 
#' @export markov_trace
markov_trace=function(mm,interventions=NULL,...) {
  # First reshape the data
  if(!is.null(interventions)) {
    # Figures out how many observations there are in each treatment & replaces the values passed 
    # as arguments in 'interventions'
    mm$m=mm$m %>% mutate(treat=rep(interventions,each=mm$m %>% count(treat) %>% slice(1) %>% pull(n)))
  }
  
  pl=mm$m %>% select(treat,t,`Pre-progressed`) %>% rename(npeople=`Pre-progressed`) %>% mutate(group="Pre-progressed") %>%
    bind_rows(mm$m %>% select(treat,t,`Progressed`) %>% rename(npeople=Progressed) %>% mutate(group="Progressed")) %>%
    bind_rows(mm$m %>% select(treat,t,Death) %>% rename(npeople=Death) %>% mutate(group="Death")) %>%
    # Create a numeric/factor group label to help manage the apparance of the graph
    mutate(
      grp_lab=as.factor(case_when(
        group=="Pre-progressed"~3,
        group=="Progressed"~2,
        TRUE~1
      ))
    ) %>% 
    ggplot(aes(x=t,y=npeople,fill=grp_lab))+geom_bar(position="stack",stat="identity") +
    labs(x="Cycle",y="Number of people",title="Markov trace",fill="State") +
    facet_wrap(~treat) + theme_bw() +
    # Add control to the legend label
    scale_fill_discrete(breaks=c(3,2,1),labels=c("Pre-progressed","Progressed","Death"))
  
  return(pl)
}

