# ----------------------------------------------------------------------------------------------------------------------------------
# Response Surface for Competitive Interactions (Intra and Inter specific)
# ----------------------------------------------------------------------------------------------------------------------------------
# Fitting Function for Lotka-Voltera competition
# ----------------------------------------------------------------------------------------------------------------------------------
#
# The following R-code of the source functions is based on the method of:
#
# Rosenbaum, B. & Rall, B.C. (2018). Fitting functional responses: Direct parameter estimation by simulating differential equations. 
# Methods in Ecology and Evolution.
#
# A second equation has been added for the modifier species - accounts for changes in modifier density and the 
# competitive interactions between the prey species and the modifier
#
# Lotka Volterra competition coefficients used to capture the intra- and inter-specific interaction between the the prey and modifier species 
#
# ----------------------------------------------------------------------------------------------------------------------------------

library(odeintr)

#### Set ODE Equation to Solve ####

# x[0] = prey, x[1] = modifier

FRS.lv = '

dxdt[0] = r1 * x[0] * (1- alpha11*x[0] - alpha12*x[1]);

dxdt[1] = r2 * x[1] * (1- alpha22*x[1] - alpha21*x[0]);
'

#### Compile ####
compile_sys("FRS_lv", FRS.lv, pars = c("r1", "r2","alpha11", "alpha22", "alpha12", "alpha21"), method = "rk54")

#### Solve model eq.odeint ####

eq.odeint.lv = function(start.v, r1, alpha11, alpha12, r2, alpha22, alpha21, Tt, timesteplength){
  FRS_lv_set_params(r1=r1, alpha11=alpha11, alpha12=alpha12, r2=r2, alpha22=alpha22, alpha21=alpha21)
  FRS_lv(start.v,Tt,timesteplength)
}


#### Solve model eaten.odeint ####
eaten.odeint.lv = function(N0, r1, alpha11, alpha12, r2, alpha22, alpha21, Tt, M0, steps=100){
  
  Neaten.est = vector()
  modgrowth.est = vector()
  
  for(i.eaten in 1:length(N0)){
    
    temp = eq.odeint.lv(start.v = c(N0[i.eaten], M0[i.eaten]),
                        r1 = r1,
                        alpha11 = alpha11, 
                        alpha12 = alpha12, 
                        r2 = r2, 
                        alpha22 = alpha22, 
                        alpha21 = alpha21, 
                        Tt= Tt[i.eaten],
                        timesteplength=Tt[i.eaten]/steps)
                                    
    Neaten.est[i.eaten] = N0[i.eaten] - temp[steps+1,2] # how many at the beginning - how many at the end = how many eaten
    modgrowth.est[i.eaten] = temp[steps+1,3] - M0[i.eaten] #M0 growth = how many at the end - how many at the beginning
  }
  
  ret = list(Neaten.est, modgrowth.est, temp)
  
  #print(tail(temp, n=1))
  return(ret)
  
}


#### Determine nll ####
nll.odeint.lv = function(r1_log, alpha11, alpha12, r2_log, alpha22, alpha21, N0, Ndead, M0, M0.end, Tt, steps=100, sigma, sigma2, verbose = F){
  if(sigma <= 0) return(Inf)
  temp2 = eaten.odeint.lv(N0 = N0,
                                      
                                      r1 = exp(r1_log), 
                                      alpha11 = alpha11, 
                                      alpha12 = alpha12, 
                                      
                                      r2 = exp(r2_log), 
                                      alpha22 = alpha22, 
                                      alpha21 = alpha21, 
                                      
                                      Tt=Tt,
                                      M0 = M0,
                                      steps=steps)
  
  y = temp2[[1]] # prey - number of prey eaten
  w = temp2[[2]] # modifier - number of new modifiers
  
  #print(c(Ndead, y))
  
  nll1 = -1*sum(dnorm(x = log(N0-Ndead), #OBSERVED number of prey at the end of the experiment
                      mean = log(N0-y), #PREDICTED number of prey at the end of the experiment
                      sd = sigma,
                      log=T))
  
  nll2 = -1*sum(dnorm(x = log(M0.end), #OBSERVED number of modifiers at the end of the experiment
                      mean = log(M0+w), #PREDICTED number of modifiers at the end of the experiment
                      sd = sigma2,
                      log=T))
  
  nll = nll1 + nll2
  
  return(nll)
}
