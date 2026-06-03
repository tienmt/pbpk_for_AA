
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# function to simulate the model for a given parameterset theta and time points x
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# x is time
# rm( k_AAuptake,k_onAA_T,V_max_p450,V_max_EH,k_onAA_B,k_onAA_Li,k_onAA_Ki,k_onGA_T,k_onGA_B,k_exc_AAMA,k_exc_GAMA )
# rm(params)
simulator <- function(times, theta) {
  params <-c( k_AAuptake  = theta[1],  # (0,1)
              k_onAA_T    = theta[2],  # (0,1)
              V_max_p450  = theta[3],  # (0,4)
              V_max_EH    = theta[4],  # (0,4)
              k_onAA_B    = theta[5],  # (0,1)
              k_onAA_Li   = theta[6],  # (0,1)
              k_onAA_Ki   = theta[7],  # (0,1)
              k_onGA_T    = theta[8],  # (0,1)
              k_onGA_B    = theta[9],   # (0,1)
              k_exc_AAMA  =  theta[10]  ,    # 0.049 #1/h or 0.13 # this will change the shape
              k_exc_GAMA  =  theta[11] 
  )
  
  diet <- data.frame(var = "m_AA_dose", method = "add",
                     time = c(0.0),  value = 0.5 *BW/1000 )  # dose of 0.05 microg/kg bw
  out <- ode(y = yini, times = times, func = PBPKmodelAA, parms = params, 
             events = list( data = diet ) )
  return( list( 'm_AAMA_urinary' = out[  ,'m_AAMA'],
                'm_GAMA_urinary' = out[  ,'m_GAMA']
  )  )
}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## measured data 
### manual readout by Sophie from Kopp and Dekant 2009 Fig 3A
##amount exc in urine nmol
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
times <- seq(from = 0, to = 80, by = 0.1)

pars_init <- c( k_AAuptake  ,  
                k_onAA_T   ,  
                V_max_p450  ,  
                V_max_EH   ,  
                k_onAA_B   ,  
                k_onAA_Li ,  
                k_onAA_Ki  ,  
                k_onGA_T   ,  
                k_onGA_B   ,  
                k_exc_AAMA  ,  
                k_exc_GAMA 
                )
lower <- 0.9 * pars_init
upper <- 1.1 * pars_init

simulator( yobs_urine$time , pars_init)
yobs_urine$AAMA
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## fitness function
fitnessL2 <- function(theta, x, y){
  - max( 1000 * ( y$AAMA - simulator(x, theta)$m_AAMA_urinary )^2 ) -
    max( 1000 * ( y$GAMA - simulator(x, theta)$m_GAMA_urinary )^2 )
}  #ng scale
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## optimisation
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++v
library(GA) ;library(parallel)
my_n_cores <- detectCores()/2  # get number of cores for parallel computing
GA <- ga(type = "real-valued", 
         fitness = fitnessL2,
         x = yobs_urine$time,  y =  yobs_urine , 
         lower = lower,    # lower bounds for parameters
         upper = upper, # upper bounds for parameters
         popSize = 200 ,   maxiter = 200,  run = 20,   parallel = my_n_cores
)
(estimated_parameters <- apply( GA@solution,2,max ) ); names(estimated_parameters) <- NULL

optimized_fit <- cumsum( simulator(yobs_urine$time,  estimated_parameters )$m_AAMA_urinary )
# plot for AAMA in urinary
par(mfrow=c(1,2))
plot(yobs_urine$time, optimized_fit , ylab = 'AAMA',
     type = 'l',col="lightblue", ylim = c(0,.1) , lwd = 2 ); grid()
points(yobs_urine$time, cumsum(yobs_urine$AAMA) , col="blue", lwd = 4)
lines(yobs_urine$time, cumsum( out[,'m_AAMA'][time_points_measure_unrine ] ),type = 'l', col = 'pink', lwd = 2, lty=2)
legend("topleft", legend = c( "pbpk_w_pre_specified",
                              'para_optimizations',
                              'observed data'),
       col = c('pink','lightblue','blue'), lty = c(2,1,1), lwd =3 )

optimized_fit_GAMA <- cumsum( simulator(yobs_urine$time,  estimated_parameters )$m_GAMA_urinary )

plot(yobs_urine$time, optimized_fit_GAMA ,type = 'l',col="lightblue", ylim = c(0,.02) , lwd = 2 ); grid()
points(yobs_urine$time, cumsum(yobs_urine$GAMA) , col="blue", lwd = 4)
lines(yobs_urine$time, cumsum( out[,'m_GAMA'][time_points_measure_unrine ] ),type = 'l', col = 'pink', lwd = 2, lty=2)
legend("topleft", legend = c( "pbpk_w_pre_specified",
                              'para_optimizations',
                              'observed data'),
       col = c('pink','lightblue','blue'), lty = c(2,1,1), lwd =3 )

names( estimated_parameters ) <- c( 'k_AAuptake','k_onAA_T','V_max_p450','V_max_EH','k_onAA_B','k_onAA_Li','k_onAA_Ki','k_onGA_T','k_onGA_B','k_exc_AAMA','k_exc_GAMA' )
estimated_parameters


