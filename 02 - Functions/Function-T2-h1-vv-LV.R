# ----------------------------------------------------------------------------------------------------------------------------------
# Functional Response Surface with Trophic interaction modifications 
# ----------------------------------------------------------------------------------------------------------------------------------
# Fitting Function for Type II Functional Response with TIM on space clearance rate (VV) and handling time (linear) with Lotka-Voltera competition
# ----------------------------------------------------------------------------------------------------------------------------------
#
# The following R-code of the source functions is based on the method of:
#
# Rosenbaum, B. & Rall, B.C. (2018). Fitting functional responses: Direct parameter estimation by
# simulating differential equations. Methods in Ecology and Evolution.
#
# A third equation has been added for the modifier species - accounts for changes in modifier density and the 
# competitive interactions between the prey species and the modifier
#
# Single estimate for predator conversion efficiency (i.e. c does not change with paramecium density)
#
# Interference parameter w represents time taken away from the predator's time budget due to interactions with the modifier
# Original parameter used to describe interference of predators with each other Beddington (1975), DeAngelis et al. (1975) adapted by  
# van Veen et al. (2005) to describe time lost to non-susceptible host
#
# This is a type II functional response thus the space clearance rate a is constant across prey density 
#
# Linear constant for a allows for the density of the modifier to impact the rate at which the predator clears the environment of prey - linear scaling
# previously seen in Hauzy et al. 2010, and Mocq et al. 2021 (no code provided) - a' = (a + a12*R2)
#
# Linear constant for h allows for the density of the modifier to impact the handling time of the predator - h' = (h + h12*R2)
#
# Lotka Volterra competition coefficients used to capture the intra- and inter-specific interaction between the the prey and modifier species 
#
# ----------------------------------------------------------------------------------------------------------------------------------


library(odeintr)

#### Set ODE Equation to Solve ####

# x[0] = predator; x[1] = prey, x[2] = modifier

FRS.t2.h1.VV.lv = '

dxdt[0] = x[0] * c * ((a * x[1]) / (1 + a * (h + h12 * x[2]) * x[1] + w * x[2])) - x[0] * m;

dxdt[1] = - ((a * x[1]) / (1 + a * (h + h12 * x[2]) * x[1] + w * x[2])) * x[0] + r1 * x[1] * (1- alpha11*x[1] - alpha12*x[2]);

dxdt[2] = r2 * x[2] * (1- alpha22*x[2] - alpha21*x[1]);
'


#### Compile ####
compile_sys("FRS_t2_h1_VV_lv", FRS.t2.h1.VV.lv, pars = c("a", "h", "h12", "w", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")


#### Solve model eq.odeint ####

eq.odeint.t2.h1.VV.lv = function(start.v, a, h, h12, w, r1, alpha11, alpha12, r2, alpha22, alpha21, c, m, Tt, timesteplength){
  FRS_t2_h1_VV_lv_set_params(a=a, h=h, h12=h12, w=w, r1=r1, alpha11=alpha11, alpha12=alpha12, r2=r2, alpha22=alpha22, alpha21=alpha21, c=c, m=m)
  FRS_t2_h1_VV_lv(start.v,Tt,timesteplength)
}


#### Solve model eaten.odeint ####
eaten.odeint.t2.h1.VV.lv = function(N0, a, h, h12, w, r1, alpha11, alpha12, r2, alpha22, alpha21, c, m, Tt, M0, P, steps=100){
  
  predrep.est = vector()
  Neaten.est = vector()
  modgrowth.est = vector()
  
  for(i.eaten in 1:length(N0)){
    
    temp = eq.odeint.t2.h1.VV.lv(start.v = c(P[i.eaten], N0[i.eaten], M0[i.eaten]),
                                 a = a,
                                 h = h,
                                 h12 = h12,
                                 w = w, 
                                 r1 = r1, 
                                 alpha11 = alpha11, 
                                 alpha12 = alpha12, 
                                 r2 = r2, 
                                 alpha22 = alpha22, 
                                 alpha21 = alpha21, 
                                 c = c,
                                 m = m,
                                 Tt= Tt[i.eaten],
                                 timesteplength=Tt[i.eaten]/steps)
    
    predrep.est[i.eaten] = temp[steps+1,2] - P[i.eaten] # how many at the end - how many at the beginning = how many new
    Neaten.est[i.eaten] = N0[i.eaten] - temp[steps+1,3] # how many at the beginning - how many at the end = how many eaten
    modgrowth.est[i.eaten] = temp[steps+1,4] - M0[i.eaten] #M0 growth = how many at the end - how many at the beginning
  }
  
  ret = list(predrep.est, Neaten.est, modgrowth.est, temp)
  
  #print(tail(temp, n=1))
  return(ret)
  
}


#### Determine nll ####
nll.odeint.t2.h1.VV.lv = function(a_log, h_log, h12, log_w, r1_log, alpha11, alpha12, r2_log, alpha22, alpha21, c_log, m, N0, Ndead, P, P.end, M0, M0.end, Tt, steps=100, sigma, sigma2, verbose = F){ #log_w
  if(sigma <= 0) return(Inf)
  temp2 = eaten.odeint.t2.h1.VV.lv(N0 = N0,
                                   a = exp(a_log),
                                   h = exp(h_log),
                                   h12 = h12,
                                   w = exp(log_w),
                                   
                                   r1 = exp(r1_log), 
                                   alpha11 = alpha11, 
                                   alpha12 = alpha12, 
                                   
                                   r2 = exp(r2_log), 
                                   alpha22 = alpha22, 
                                   alpha21 = alpha21, 
                                   
                                   c = exp(c_log),
                                   m = m,
                                   Tt=Tt,
                                   P=P,
                                   M0 = M0,
                                   steps=steps)
  
  y = temp2[[2]] # prey - number of prey eaten
  z = temp2[[1]] # pred - number of new predators
  w = temp2[[3]] # modifier - number of new modifiers
  
  #print(c(Ndead, y))
  
  nll1 = -1*sum(dnorm(x = log(N0-Ndead), #OBSERVED number of prey at the end of the experiment
                      mean = log(N0-y), #PREDICTED number of prey at the end of the experiment
                      sd = sigma,
                      log=T))
  
  nll2 = -1*sum(dpois(x = P.end, #OBSERVED number of predators at the end of the experiment
                      lambda = P+z, #PREDICTED number of predators at the end of the experiment
                      log=T))
  
  nll3 = -1*sum(dnorm(x = log(M0.end), #OBSERVED number of modifiers at the end of the experiment
                      mean = log(M0+w), #PREDICTED number of modifiers at the end of the experiment
                      sd = sigma2,
                      log=T))
  
  nll = nll1 + nll2 + nll3
  
  return(nll)
}


