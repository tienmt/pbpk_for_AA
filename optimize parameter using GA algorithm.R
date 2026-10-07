
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# function to simulate the model for a given parameterset theta and time points x
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# x is time
# rm( k_AAuptake,k_onAA_T,V_max_p450,V_max_EH,k_onAA_B,k_onAA_Li,k_onAA_Ki,k_onGA_T,k_onGA_B,k_exc_AAMA,k_exc_GAMA )
# rm(params)
simulator <- function(times, theta) {
  params <- c( Vmax_AA_GSH  = theta[1],  
              pAA_TB       = theta[2],  
              V_max_p450   = theta[3], 
              Vmax_GA_GSH  = theta[4],  
              pGA_TB      = theta[5], 
              KM_p450     = theta[6], 
              k_exc_AAMA  =  theta[7]  , 
              k_exc_GAMA  =  theta[8] ,
              K_FORM_AA_VAL = theta[9]
  )
  
  diet <- data.frame(var = "m_AA_dose", method = "add",
                     time = c(0.0),  value = 0.5 *BW/1000 )  # dose of 0.05 microg/kg bw
  out <- ode(y = yini, times = times, func = PBPKmodelAA, parms = params, 
             events = list( data = diet ) )
  return( list( 'm_AAMA_urinary' = out[  ,'m_AAMA'],
                'm_GAMA_urinary' = out[  ,'m_GAMA']
  )  )
}
###########################################################
## measured data 
### manual readout by Sophie from Kopp and Dekant 2009 Fig 3A
##amount exc in urine nmol
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
times <- seq(from = 0, to = 80, by = 0.1)

pars_init <- c(  Vmax_AA_GSH ,  
                 pAA_TB   , 
                 V_max_p450 , 
                 Vmax_GA_GSH  ,  
                 pGA_TB   , 
                 KM_p450  , 
                 k_exc_AAMA    ,  
                 k_exc_GAMA,
                 K_FORM_AA_VAL
)
lower <- 0.9 * pars_init
upper <- 1.1 * pars_init

simulator( yobs_urine$time , pars_init)
yobs_urine$AAMA

## fitness function
fitnessL2 <- function(theta, x, y){
  - max( 1000 * ( y$AAMA - simulator(x, theta)$m_AAMA_urinary )^2 ) -
    max( 1000 * ( y$GAMA - simulator(x, theta)$m_GAMA_urinary )^2 )
}  # ng scale





### Optimisation
library(GA) ;library(parallel)
my_n_cores <- detectCores()/2  # get number of cores for parallel computing

GA <- ga(type = "real-valued", 
         fitness = fitnessL2,
         x = yobs_urine$time,  y =  yobs_urine , 
         lower = lower,    # lower bounds for parameters
         upper = upper, # upper bounds for parameters
         popSize = 200 ,   maxiter = 200,  run = 20,   parallel = my_n_cores
)
(estimated_parameters <- apply( GA@solution,2,max ) )


