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

# x[0] = prey

FRS.lvgrowth = '

dxdt[0] = r1 * x[0] * (1- alpha11*x[0]);
'

#### Compile ####
compile_sys("FRS_lvgrowth", FRS.lvgrowth, pars = c("r1", "alpha11"), method = "rk54")

#### Solve model eq.odeint ####

eq.odeint.lvgrowth = function(start.v, r1, alpha11, Tt, timesteplength){
  FRS_lvgrowth_set_params(r1=r1, alpha11=alpha11)
  FRS_lvgrowth(start.v,Tt,timesteplength)
}


#### Solve model eaten.odeint ####
eaten.odeint.lvgrowth = function(N0, r1, alpha11, Tt, steps=100){
  
  Neaten.est = vector()
  
  for(i.eaten in 1:length(N0)){
    
    temp = eq.odeint.lvgrowth(start.v = c(N0[i.eaten]),
                        r1 = r1,
                        alpha11 = alpha11, 
                        Tt= Tt[i.eaten],
                        timesteplength=Tt[i.eaten]/steps)
    
    Neaten.est[i.eaten] = N0[i.eaten] - temp[steps+1,2] # how many at the beginning - how many at the end = how many eaten
  }
  
  ret = list(Neaten.est, temp)
  
  #print(tail(temp, n=1))
  return(ret)
  
}


#### Determine nll ####
nll.odeint.lvgrowth = function(r1_log, alpha11, N0, Ndead, steps=100, sigma, verbose = F){
  if(sigma <= 0) return(Inf)
  temp2 = eaten.odeint.lvgrowth(N0 = N0,
                                
                                r1 = exp(r1_log), 
                                alpha11 = alpha11, 
                                
                                Tt=Tt,
                                steps=steps)
  
  y = temp2[[1]] # prey - number of prey eaten
 
  #print(c(Ndead, y))
  
  nll1 = -1*sum(dnorm(x = log(N0-Ndead), #OBSERVED number of prey at the end of the experiment
                      mean = log(N0-y), #PREDICTED number of prey at the end of the experiment
                      sd = sigma,
                      log=T))
  
  nll = nll1
  
  return(nll)
}
