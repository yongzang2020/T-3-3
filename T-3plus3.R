#######################################################################################################################
## ## This file contains two functions for implementation of the proposed time-to-event 3+3 (T-3+3) design to find the maximum
## tolerated dose (MTD) with pending dose-limiting toxicity (DLT).
##
## t3plus3.oc() is the function used to generate operating characteristics for the proposed design
## t3plus3.df() is the function used for dose-finding of the actual trial.
## t3plus3.select() is a function used for determine the MTD of the actual trial.
##
## At the end of this file, an example is provided to illustrate how to use the proposed designs to conduct a clinical trial.
########################################################################################################################


##########################################################################################################
## Function to generate operating characteristics of the T-3+3 design for single agent trials with 
## pending DLT by simulating trials. 
##
## Arguments:
## target: the target toxicity rate
## p.true: a vector containing the true toxicity probabilities of the investigational dose levels.
## cohortsize: the cohort size
## maxt: the maximum follow-up time 
## prior.p: a vector of length 3, which specifies the prior probability that the time to toxicity 
##          lies inside the time interval (0,\code{maxt}/3), (\code{maxt}/3,\code{2maxt}/3), 
##          (\code{2maxt}/3,1). The default value is \code{prior.p=c(1/3,1/3,1/3)}. 
## accrual: the accrual rate, i.e., the number of patients accrued in 1 unit of time
## phiE: a number from (0,1) that controls the minimum probability required for escalation. 
##       The default is \cod{phiE=0.5}.
## phiR: a number from (0,1) that controls the minimum probability required for retention. 
##       The default is \cod{phiR=0.5}.
## phiD: a number from (0,1) that controls the minimum probability required for de-escalation. 
##       The default is \cod{phD=0.75}.
## dist1: the underlying distribution of the time to toxicity outcomes; \code{dist1=1} is the
##        uniform distribution, \code{dist1=2} corresponds to the Weibull distribution,
##        \code{dist1=3} is the log-logistic distribution, \code{dist1=4} is the piecewise 
##        exponential distribution.
## dist2: the underlying distribution of patient arrival time; \code{dist2=1} is the 
##        uniform distribution, \code{dist2=2} is the exponential distribution
## alpha: a number from (0,1) that controls alpha*100% events in (0, 1/2T). 
##        The default is \code{alpha=0.5}.             
## startdose: the starting dose level for the trial
## ntrial:  the total number of trials to be simulated.
#########################################################################################################

t3plus3.oc <- function(target,
                       p.true, 
                       cohortsize=3, 
                       maxt=1, 
                       prior.p=NA, 
                       accrual=1, 
                       phiE=0.5, 
                       phiR=0.5, 
                       phiD=0.75,  
                       dist1=1, 
                       dist2=1, 
                       alpha=0.5, 
                       startdose = 1, 
                       ntrial=1000, 
                       seed=123){
  
  select.mtd <- function(target, npts, ntox, verbose = TRUE) {
    ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
    pava <- function(x, wt = rep(1, length(x))) {
      n <- length(x)
      if (n <= 1)
        return(x)
      if (any(is.na(x)) || any(is.na(wt))) {
        stop("Missing values in 'x' or 'wt' not allowed")
      }
      lvlsets <- (1:n)
      repeat {
        viol <- (as.vector(diff(x)) < 0)
        if (!(any(viol)))
          break
        i <- min((1:(n - 1))[viol])
        lvl1 <- lvlsets[i]
        lvl2 <- lvlsets[i + 1]
        ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
        x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
        lvlsets[ilvl] <- lvl1
      }
      x
    }
    
    ## no dose should be selected (i.e., selectdose=0) if the first dose is already very toxic
    y = ntox
    n = npts
    if (y[1]>=2) {
      selectdose = 0
    } else {
      adm.set = (npts != 0)
      adm.index = which(adm.set == T)
      y.adm = y[adm.set]
      n.adm = n[adm.set]
      
      ## poster mean and variance of toxicity probabilities using beta(0.05, 0.05) as the prior
      phat = (y.adm + 0.05)/(n.adm + 0.1)
      phat.var = (y.adm + 0.05) * (n.adm - y.adm + 0.05)/((n.adm + 0.1)^2 * (n.adm + 0.1 + 1))
      
      ## perform the isotonic transformation using PAVA
      phat = pava(phat, wt = 1/phat.var)
      phat = phat + (1:length(phat)) * 1e-10  ## break ties by adding an increasingly small number
      selectd = sort(abs(phat - target), index.return = T)$ix[1]  ## select dose closest to the target as the MTD
      selectdose = adm.index[selectd]
    }
    
    if (verbose == TRUE) {
      if (selectdose == 0) {
        out = list(target = target, MTD = selectdose,
                   p_est = data.frame(cbind('dose'=1:length(npts), 'phat'=rep("----",length(npts)),
                                            'CI'=paste("(", rep("----",length(npts)),",",rep("----",length(npts)),")",sep="")))
        )
      } else {
        trtd = (n != 0)
        poverdose = pava(1 - pbeta(target, y[trtd] + 0.05, n[trtd] - y[trtd] + 0.05))
        phat.all = pava((y[trtd] + 0.05)/(n[trtd] + 0.1), wt = 1/((y[trtd] + 0.05) * (n[trtd] - y[trtd] + 0.05)/((n[trtd] + 0.1)^2 * (n[trtd] + 0.1 + 1))))
        
        A1 = A2 = NA
        A3 = NA
        A4 = A5 = NA
        ## output summary statistics
        for (i in 1:ndose) {
          if (n[i] > 0) {
            A1 = append(A1, formatC(phat.all[i], digits = 2, format = "f"))
            A2 = append(A2, formatC(qbeta(0.025, y[i] + 0.05, n[i] - y[i] + 0.05), digits = 2, format = "f"))
            A3 = append(A3, formatC(qbeta(0.975, y[i] + 0.05, n[i] - y[i] + 0.05), digits = 2, format = "f"))
            A4 = append(A4, formatC(poverdose[i], digits = 2, format = "f"))
          } else {
            # no estimate output for doses never used to treat patients
            A1 = append(A1, "----")
            A2 = append(A2, "----")
            A3 = append(A3, "----")
            A4 = append(A4, "----")
          }
        }
        p_est = data.frame(cbind('dose'=1:length(npts), 'phat'=A1[-1], 'CI'=paste("(", A2[-1],",",A3[-1],")",sep="")))
        
        out = list(target = target, MTD = selectdose, p_est=p_est, p_overdose = A4[-1])
      }
    } else {
      out = list(target = target, MTD = selectdose)
    }
    return(out)
  }
  
  
  gen.tite<-function(dist=1, n, pi, alpha=0.5, Tobs=1)
  {
    ############ subroutines ############
    weib<-function(n, pi, pihalft)
    {
      ## solve parameters for Weibull given pi=1-S(T) and phalft=1-S(T/2)
      alpha = log(log(1-pi)/log(1-pihalft))/log(2);
      lambda = -log(1-pi)/(Tobs^alpha);
      t = (-log(runif(n))/lambda)^(1/alpha);
      return(t);
    }
    
    llogit<-function(n, pi, pihalft)
    {
      ## solve parameters for log-logistic given pi=1-S(T) and phalft=1-S(T/2)
      alpha = log((1/(1-pi)-1)/(1/(1-pihalft)-1))/log(2);
      lambda = (1/(1-pi)-1)/(Tobs^alpha);
      t = ((1/runif(n)-1)/lambda)^(1/alpha);
      return(t);
    }
    
    pwexp<-function(n, pi, pihalft)
    {
      ## solve parameters for piecewise exponential given pi=1-S(T) and phalft=1-S(T/2)
      lambda1 = -log(1-pihalft)*2/Tobs;
      lambda2 = -log((1-pi)/(1-pihalft))*2/Tobs;
      t = NULL;
      U = runif(n);
      for (i in seq_along(U)){
        if (U[i]<=pihalft) {x = -log(1-U[i])/lambda1}
        else {x = Tobs/2-log((1-U[i])/(1-pihalft))/lambda2}
        t = c(t,x)
      }
      return(t);
    }
    
    ############ end of subroutines ############
    
    
    tox = rep(0, n);
    t.tox = rep(0, n);
    
    #### uniform
    if(dist==1) {  # 50% event in (0, 1/2T)
      tox = rbinom(n, 1, pi);
      ntox.st = sum(tox);
      t.tox[tox==0]=Tobs;
      t.tox[tox==1]=runif(ntox.st, 0, Tobs);
    }
    #### Weibull
    if(dist==2)
    {
      pihalft = alpha*pi;  # alpha*100% event in (0, 1/2T)
      t.tox = weib(n, pi, pihalft);
      tox[t.tox<=Tobs]=1;
      ntox.st = sum(tox);
      t.tox[tox==0]=Tobs;
    }
    #### log-logistic
    if(dist==3)
    {
      pihalft = alpha*pi;  # alpha*100% event in (0, 1/2T)
      t.tox = llogit(n, pi, pihalft);
      tox[t.tox<=Tobs]=1;
      ntox.st = sum(tox);
      t.tox[tox==0]=Tobs;
    }
    #### piecewise exponential
    if(dist==4)
    {
      pihalft = alpha*pi;  # alpha*100% event in (0, 1/2T)
      t.tox = pwexp(n, pi, pihalft);
      tox[t.tox<=Tobs]=1;
      ntox.st = sum(tox);
      t.tox[tox==0]=Tobs;
    }
    return(list(tox=tox, t.tox=t.tox, ntox.st=ntox.st));
  }
  
  
  ### simple error checking
  if(target<0.2) {cat("Error: the target is too low! \n"); return();}
  if(target>0.4)  {cat("Error: the target is too high! \n"); return();}
  if(!is.na(prior.p[1])){if(length(prior.p)!=3){cat("Error: The length of the prior probabilities should be 3! \n"); return();}}
  
  set.seed(seed);
  if(is.na(prior.p[1])){prior.p = rep(1/3,3)}
  prior.p = prior.p/sum(prior.p)
  ndose = length(p.true);
  Y = matrix(rep(0, ndose * ntrial), ncol = ndose);
  N = matrix(rep(0, ndose * ntrial), ncol = ndose);
  dselect = rep(0, ntrial);
  durationV = rep(0, ntrial);
  npendV = rep(0, ntrial);
  
  
  for(trial in 1:ntrial)
  {
    y=NULL;  #toxicity indicator for each subject
    dv=NULL;  #dose for each subject
    n.d = rep(0, ndose);  # number of toxicity at each dose
    y.d = rep(0, ndose);  # number of patient at each dose
    t.enter=NULL; # time enter the study
    t.event=NULL; # time to event
    t.decision = 0; # decision making time
    d = startdose;  # current dose level
    earlystop = 0; #indicate if trial stops early
    stop = 0;
    npend = 0;
    mtd = rep(0, ndose)
    
    while(stop==0)
    {
      # generate data for the new patient
      for(j in 1:cohortsize)
      {
        if(j==1) { t.enter = c(t.enter, t.decision); 
        }else { 
          if(dist2==1){ t.enter = c(t.enter, t.enter[length(t.enter)] + runif(1, 0, 2/accrual))}
          if(dist2==2){ t.enter = c(t.enter, t.enter[length(t.enter)] + rexp(1, rate=accrual))}
        }
      }
      obscohort = gen.tite(dist1, cohortsize, p.true[d], alpha=alpha, Tobs=maxt);
      t.event = c(t.event, obscohort$t.tox);
      y = c(y, obscohort$tox);
      dv = c(dv, rep(d, cohortsize));
      t.decision = t.enter[length(t.enter)];
      nobs=-1; pending=1;
      d.curr=d;
      npend = npend-1;
      while(pending==1)
      {
        npend = npend+1;
        pending = 0;
        if(dist2==1){t.decision = t.decision + runif(1, 0, 2/accrual)}
        if(dist2==2){t.decision = t.decision + rexp(1, rate=accrual)}
        
        # determine which observation are observed
        delta = ((t.enter+t.event)<=t.decision);
        t = pmin(t.event, t.decision-t.enter, maxt);  ## used for recording potential censoring time
        cset = (dv==d);
        # pending when no patient completes the follow-up
        if (sum(delta[cset])==0) {pending=1;next;}
        delta.curr = delta[cset];
        t.curr = t[cset];
        ntox.curr = sum((y[cset])[delta.curr==1]);
        #totalt = sum(t.curr[delta.curr==0])/maxt;
        totalt = t.curr[delta.curr==0]
        totalt = 3*prior.p[1]*totalt*(totalt<=maxt/3)+
          ((prior.p[1]-prior.p[2])*maxt+3*prior.p[2]*totalt)*(maxt/3<totalt & totalt<=2*maxt/3)+
          ((prior.p[1]+prior.p[2]-2*prior.p[3])*maxt+3*prior.p[3]*totalt)*(2*maxt/3<totalt & totalt<=maxt)
        totalt = sum(totalt)/maxt
        n.curr = sum(cset);
        n.pend = sum(delta[cset]==0);
        nobs = sum(delta[cset]);
        n.d[d] = n.curr
        n1 = nobs
        r1 = ntox.curr
        n2 = n.pend
        n2t = totalt 
        
        #set l1, l2 based on number of patients treated at current dose
        if (n.curr==3) {l1=0; l2=2}
        if (n.curr==6) {l1=1; l2=2}
        
        #calculate probability of escalation, de-escalation and retention
        if (n2==0) {if (r1<=l1) {pe=1;pd=0} else if (r1>=l2) {pd=1;pe=0} else {pd=0;pe=0}} 
        if (n2>0) {if (l1>=r1) {pe=pbb(l1-r1,n2,r1+1,n1-r1+n2t+1)} else {pe=0}}
        if (n2>0) {if (l2-r1>=1) {pd=1-pbb(l2-r1-1,n2,r1+1,n1-r1+n2t+1)} else {pd=1}}
        pr=1-pe-pd
        
        #make dose assignment decisions
        maxp=max(pe,pd,pr)
        if ((maxp==pe & maxp<phiE) || (maxp==pr & maxp<phiR) || (maxp==pd & maxp<phiD)) {pending=1
        } else if (maxp==pe) {if (mtd[d+1]==2||d==ndose){if (n.curr==6){stop=1;break} else {d=d}} else {d=d+1}
        } else if (maxp==pr) {if (n.curr==3) {d=d} else {stop=1;break}
        } else if (maxp==pd) {if (d==1||n.d[d-1]==6) {stop=1;break} else {mtd[d]=2;d=d-1}}
      }
    }
    
    for(k in 1:ndose){
      y.d[k] = sum(y[dv==k]);
      n.d[k] = sum(dv==k);
    }
    
    npendV[trial]= npend;
    Y[trial, ] = y.d
    N[trial, ] = n.d
    durationV[trial] = t.decision
    dselect[trial] = select.mtd(target, n.d, y.d, verbose=FALSE)$MTD
  }
  
  selpercent = rep(0, ndose)
  nptsdose = apply(N,2,mean);
  pptsdose = round(nptsdose/sum(nptsdose)*100,2)
  ntoxdose = apply(Y,2,mean);
  npts = rowSums(N)
  
  for(i in 0:ndose) {selpercent[i+1]=sum(dselect==i)/ntrial*100; }
  
  out=list(selpercent=selpercent, ptpercent=pptsdose, ntox=ntoxdose, 
           totaltox=sum(Y)/ntrial, totaln=sum(N)/ntrial, duration=mean(durationV),
           simu.setup=data.frame(target=target, p.true=c('-', p.true), ntrial = ntrial, dose=0:ndose, 
                                 selected=selpercent, ppatients=c('-',pptsdose), ntox=c('-',ntoxdose)));
  
  return(out);
}


############################################################################
## Function for real trial dose-finding procedure of T-3+3
##
## ndose: number of dose levels 
## toxv: toxicity event for each individual by interim time
## enter_time: time for individuals entering the trial
## event_time: time from entering the trial to toxicity event, 
##           if do not experience the toxicity by interim time, simply plug in the time assessment window for toxicity
## decision_time: interim decison time
## dv: dose level for each individual
## maxt: time assess window for toxicity
## prior.p: a vector of length 3, which specifies the prior probability that the time to toxicity 
##          lies inside the time interval (0,\code{maxt}/3), (\code{maxt}/3,\code{2maxt}/3), 
##          (\code{2maxt}/3,1). The default value is \code{prior.p=c(1/3,1/3,1/3)}. 
## phiE: a number from (0,1) that controls the minimum probability required for escalation. 
##       The default is \cod{phiE=0.5}.
## phiR: a number from (0,1) that controls the minimum probability required for retention. 
##       The default is \cod{phiR=0.5}.
## phiD: a number from (0,1) that controls the minimum probability required for de-escalation. 
##       The default is \cod{phiD=0.75}.
############################################################################

t3plus3.df = function(ndose,
                      toxv,
                      enter_time,
                      event_time,
                      decision_time,
                      dv,                       
                      maxt,
                      prior.p=NA,
                      phiE=0.5, 
                      phiR=0.5, 
                      phiD=0.75){

  pending=0
  stop=0
  if(is.na(prior.p[1])){prior.p = rep(1/3,3)}
  prior.p = prior.p/sum(prior.p)
  n.d = rep(0, ndose) #number of patients at each dose by decision time
  for (dd in 1:length(dv)) {
    n.d[dv[dd]]=n.d[dv[dd]]+1
  }
  mtd = rep(0, ndose) #set the dose of toxicity to 2
  for (dd in 1:(length(dv)-1)) {
    if (dv[dd+1]<dv[dd]) {mtd[dv[dd]]=2}
  }
  d = dv[(length(dv))] #current dose
  delta = ((enter_time+event_time)<=decision_time);
  t = pmin(event_time, decision_time-enter_time, maxt);  ## used for recording potential censoring time
  cset = (dv==d);
  # pending when no patient completes the follow-up
  if (sum(delta[cset])==0) {pending=1
  }else {
    delta.curr = delta[cset];
    t.curr = t[cset];
    ntox.curr = sum((toxv[cset])[delta.curr==1]);
    totalt = t.curr[delta.curr==0]
    totalt = 3*prior.p[1]*totalt*(totalt<=maxt/3)+
      ((prior.p[1]-prior.p[2])*maxt+3*prior.p[2]*totalt)*(maxt/3<totalt & totalt<=2*maxt/3)+
      ((prior.p[1]+prior.p[2]-2*prior.p[3])*maxt+3*prior.p[3]*totalt)*(2*maxt/3<totalt & totalt<=maxt)
    totalt = sum(totalt)/maxt
    n.curr = sum(cset);
    n.pend = sum(delta[cset]==0);
    nobs = sum(delta[cset]);
    n.d[d] = n.curr
    n1 = nobs
    r1 = ntox.curr
    n2 = n.pend
    n2t = totalt 
    
    #set l1, l2 based on number of patients treated at current dose
    if (n.curr==3) {l1=0; l2=2}
    if (n.curr==6) {l1=1; l2=2}
    
    #calculate probability of escalation, de-escalation and retention
    if (n2==0) {if (r1<=l1) {pe=1;pd=0} else if (r1>=l2) {pd=1;pe=0} else {pd=0;pe=0}} 
    if (n2>0) {if (l1>=r1) {pe=pbb(l1-r1,n2,r1+1,n1-r1+n2t+1)} else {pe=0}}
    if (n2>0) {if (l2-r1>=1) {pd=1-pbb(l2-r1-1,n2,r1+1,n1-r1+n2t+1)} else {pd=1}}
    pr=1-pe-pd
    
    #make dose assignment decisions
    maxp=max(pe,pd,pr)
    if ((maxp==pe & maxp<phiE) || (maxp==pr & maxp<phiR) || (maxp==pd & maxp<phiD)) {pending=1
    } else if (maxp==pe) {if (mtd[d+1]==2||d==ndose){if (n.curr==6){stop=1} else {d=d}} else {d=d+1}
    } else if (maxp==pr) {if (n.curr==3) {d=d} else {stop=1}
    } else if (maxp==pd) {if (d==1||n.d[d-1]==6) {stop=1} else {mtd[d]=2;d=d-1}}
    
    if (pending==1) {return(list("next dose"="suspend the trial"))
    }else if (stop==1) {return(list("next dose"="stop the trial"))
    }else {return(list("next dose"=d))}
  }
}


############################################################################
## Function for real trial MTD determination by isotonic regression
##
## Arguments:
## target: the target toxicity rate
## npts: number of patients treated at each dose level
## ntox: number of patients reporting toxicity at each dose level
############################################################################
t3plus3.select = function(target, npts, ntox, verbose = TRUE) {
  ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
  pava <- function(x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1)
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol)))
        break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }
  
  ## no dose should be selected (i.e., selectdose=0) if the first dose is already very toxic
  y = ntox
  n = npts
  if (y[1]>=2) {
    selectdose = 0
  } else {
    adm.set = (npts != 0)
    adm.index = which(adm.set == T)
    y.adm = y[adm.set]
    n.adm = n[adm.set]
    
    ## poster mean and variance of toxicity probabilities using beta(0.05, 0.05) as the prior
    phat = (y.adm + 0.05)/(n.adm + 0.1)
    phat.var = (y.adm + 0.05) * (n.adm - y.adm + 0.05)/((n.adm + 0.1)^2 * (n.adm + 0.1 + 1))
    
    ## perform the isotonic transformation using PAVA
    phat = pava(phat, wt = 1/phat.var)
    phat = phat + (1:length(phat)) * 1e-10  ## break ties by adding an increasingly small number
    selectd = sort(abs(phat - target), index.return = T)$ix[1]  ## select dose closest to the target as the MTD
    selectdose = adm.index[selectd]
  }
  
  if (verbose == TRUE) {
    if (selectdose == 0) {
      out = list(target = target, MTD = selectdose,
                 p_est = data.frame(cbind('dose'=1:length(npts), 'phat'=rep("----",length(npts)),
                                          'CI'=paste("(", rep("----",length(npts)),",",rep("----",length(npts)),")",sep="")))
      )
    } else {
      trtd = (n != 0)
      poverdose = pava(1 - pbeta(target, y[trtd] + 0.05, n[trtd] - y[trtd] + 0.05))
      phat.all = pava((y[trtd] + 0.05)/(n[trtd] + 0.1), wt = 1/((y[trtd] + 0.05) * (n[trtd] - y[trtd] + 0.05)/((n[trtd] + 0.1)^2 * (n[trtd] + 0.1 + 1))))
      
      A1 = A2 = NA
      A3 = NA
      A4 = A5 = NA
      ## output summary statistics
      for (i in 1:ndose) {
        if (n[i] > 0) {
          A1 = append(A1, formatC(phat.all[i], digits = 2, format = "f"))
          A2 = append(A2, formatC(qbeta(0.025, y[i] + 0.05, n[i] - y[i] + 0.05), digits = 2, format = "f"))
          A3 = append(A3, formatC(qbeta(0.975, y[i] + 0.05, n[i] - y[i] + 0.05), digits = 2, format = "f"))
          A4 = append(A4, formatC(poverdose[i], digits = 2, format = "f"))
        } else {
          # no estimate output for doses never used to treat patients
          A1 = append(A1, "----")
          A2 = append(A2, "----")
          A3 = append(A3, "----")
          A4 = append(A4, "----")
        }
      }
      p_est = data.frame(cbind('dose'=1:length(npts), 'phat'=A1[-1], 'CI'=paste("(", A2[-1],",",A3[-1],")",sep="")))
      
      out = list(target = target, MTD = selectdose, p_est=p_est)
    }
  } else {
    out = list(target = target, MTD = selectdose)
  }
  return(out)
}


########################################## an example #######################################

####### generate operating characteristics of the designs ##########
target=0.3
accrual<-2
maxt<-3
ntrial<-10000
P<-c(0.05,0.10,0.20,0.31,0.50,0.70)
t3plus3.oc(target,P,maxt=maxt,accrual=accrual,dist1=2,dist2=2,ntrial=ntrial)
#$selpercent
#[1]  2.05  5.21 18.39 33.34 31.30  9.24  0.47

#$ptpercent
#[1] 22.32 25.17 24.86 18.17  8.13  1.35

#$ntox
#[1] 0.1785 0.4103 0.8021 0.9049 0.6529 0.1541

#$totaltox
#[1] 3.1028

#$totaln
#[1] 16.2051

#$duration
#[1] 16.05361

#$simu.setup
#target p.true ntrial dose selected ppatients   ntox
#1    0.3      -  10000    0     2.05         -      -
#2    0.3   0.05  10000    1     5.21     22.32 0.1785
#3    0.3    0.1  10000    2    18.39     25.17 0.4103
#4    0.3    0.2  10000    3    33.34     24.86 0.8021
#5    0.3   0.31  10000    4    31.30     18.17 0.9049
#6    0.3    0.5  10000    5     9.24      8.13 0.6529
#7    0.3    0.7  10000    6     0.47      1.35 0.1541

## using the example in trial implementation t3plus3.df(6,toxv,enter_time,event_time,decision_time,dv,maxt)
## call dose-finding function to obtain the recommended dose assignment for the next cohort at day 106
toxv=c(0,0,0)
enter_time=c(1,16,31)
event_time=c(90,90,90)
decision_time=106
dv=c(1,1,1)
t3plus3.df(6,toxv,enter_time,event_time,decision_time,dv,90)
#$`next dose`
#[1] 2

## call dose-finding function to obtain the recommended dose assignment for the next cohort at day 211
toxv=c(0,0,0,1,0,0)
enter_time=c(1,16,31,106,121,136)
event_time=c(90,90,90,60,90,90)
decision_time=211
dv=c(1,1,1,2,2,2)
t3plus3.df(6,toxv,enter_time,event_time,decision_time,dv,90)
#$`next dose`
#[1] 2

## suppose we want to obtain the recommended dose assignment for the next cohort at day 300,
## we need to suspend the trail for more information to make the decision
toxv=c(0,0,0,1,0,0,0,0,0)
enter_time=c(1,16,31,106,121,136,211,226,241)
event_time=c(90,90,90,60,90,90,90,90,90)
decision_time=300
dv=c(1,1,1,2,2,2,2,2,2)
t3plus3.df(6,toxv,enter_time,event_time,decision_time,dv,90)
#$`next dose`
#[1] "suspend the trial"

## ....... repeat the above procedure for incoming cohorts................

## Suppose after the last cohort of patients have been treated, the updated data are
npts=c(3,6,6,3,0,0)
ntox=c(0,1,1,2,0,0)
t3plus3.select(target,npts,ntox)

#$target
#[1] 0.3

#$MTD
#[1] 3

#$p_est
#dose phat          CI
#1    1 0.02 (0.00,0.20)
#2    2 0.17 (0.01,0.53)
#3    3 0.17 (0.01,0.53)
#4    4 0.66 (0.16,0.99)
#5    5 ---- (----,----)
#6    6 ---- (----,----)
