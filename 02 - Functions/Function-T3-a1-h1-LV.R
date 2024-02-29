# ----------------------------------------------------------------------------------------------------------------------------------
# Functional Response Surface with Trophic interaction modifications 
# ----------------------------------------------------------------------------------------------------------------------------------
# Fitting Function for Type III Functional Response with TIM on space clearance rate (linear) and handling time (constant) with Lotka-Voltera competition
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
# This is a type III functional response thus the space clearance rate a is dependent on prey density. To prevent biological unrealistic 
# estimates of a=bN^q that can occur with Real (1977), we use a form of a sigmoidal functional response in which a is asymptotic 
# (i.e. predator–prey interactions are suppressed at low prey abundances but converge on what might happen in a type II scenario as prey 
# levels increase) outlined in DeLong (2021) and originally proposed by Hassel et al. (1977). kr is a half-saturation constant (i.e. 
# the density of prey where the realized space clearance rate reaches half of the maximum value). This “asymptotic a” model recognizes 
# that space clearance rate has limits and cannot increase indefinitely with prey levels.   
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

FRS.t3.a1.h1.lv = '

dxdt[0] = x[0] * c * (((((a + a12 * x[2]) * x[1]) / (kr + x[1]))  * x[1]) / (1 + (((a + a12 * x[2]) * x[1]) / (kr + x[1])) * (h + h12 * x[2]) * x[1])) - x[0] * m;

dxdt[1] = - (((((a + a12 * x[2]) * x[1]) / (kr + x[1]))  * x[1]) / (1 + (((a + a12 * x[2]) * x[1]) / (kr + x[1])) * (h + h12 * x[2]) * x[1])) * x[0] + r1 * x[1] * (1- alpha11*x[1] - alpha12*x[2]);

dxdt[2] = r2 * x[2] * (1- alpha22*x[2] - alpha21*x[1]);
'


#### Compile ####
compile_sys("FRS_t3_a1_h1_lv", FRS.t3.a1.h1.lv, pars = c("a", "a12", "kr", "h", "h12", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")


#### Solve model eq.odeint ####

eq.odeint.t3.a1.h1.lv = function(start.v, a, a12, kr, h, h12, r1, alpha11, alpha12, r2, alpha22, alpha21, c, m, Tt, timesteplength){
  FRS_t3_a1_h1_lv_set_params(a=a, a12=a12, kr=kr, h=h, h12=h12, r1=r1, alpha11=alpha11, alpha12=alpha12, r2=r2, alpha22=alpha22, alpha21=alpha21, c=c, m=m)
  FRS_t3_a1_h1_lv(start.v,Tt,timesteplength)
}


#### Solve model eaten.odeint ####
eaten.odeint.t3.a1.h1.lv = function(N0, a, a12, kr, h, h12, r1, alpha11, alpha12, r2, alpha22, alpha21, c, m, Tt, M0, P, steps=100){
  
  predrep.est = vector()
  Neaten.est = vector()
  modgrowth.est = vector()
  
  for(i.eaten in 1:length(N0)){
    
    temp = eq.odeint.t3.a1.h1.lv(start.v = c(P[i.eaten], N0[i.eaten], M0[i.eaten]),
                                 a = a,
                                 a12 = a12,
                                 kr = kr,
                                 h = h,
                                 h12 = h12,
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
nll.odeint.t3.a1.h1.lv = function(a_log, a12, h_log, h12, kr_log, r1_log, alpha11, alpha12, r2_log, alpha22, alpha21, c_log, m, N0, Ndead, P, P.end, M0, M0.end, Tt, steps=100, sigma, sigma2, verbose = F){ #log_w
  if(sigma <= 0) return(Inf)
  temp2 = eaten.odeint.t3.a1.h1.lv(N0 = N0,
                                   a = exp(a_log),
                                   a12 = a12,
                                   h = exp(h_log),
                                   h12 = h12,
                                   kr = exp(kr_log),
                                   
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


