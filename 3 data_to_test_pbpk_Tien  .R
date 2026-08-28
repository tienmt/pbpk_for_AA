

################################################################
# manual readout from Kopp and Dekant 2009
################################################################
times <- seq(from = 0, to = 60, by = 0.1)
diet <- data.frame(var = "m_AA_dose", method = "add",
                   time = c(0.0),  value = 0.05 *BW )  # dose of 0.05 microg/kg bw
out <- ode(y = yini, times = times, func = PBPKmodelAA, parms = params, events = list(data = diet))
yobs_urine <- data.frame(
  time = c(0.0, 3.9, 8.3, 14, 19.5, 28, 37, 45.9),
  AAMA = c(0.0/1000000, 32.8/1000000, 69.5/1000000, 42.8/1000000, 45.2/1000000, 24.5/1000000, 15/1000000, 7/1000000)*234.28, #https://pubchem.ncbi.nlm.nih.gov/compound/N-Acetyl-S-_2-carbamoylethyl_-L-cysteine
  GAMA = c(0.0, 2.15*1e-6, 3.66*1e-6, 4.38*1e-6, 6.57*1e-6, 5.26*1e-6, 3.50*1e-6, 3.0*1e-6)*250.27
)
# plot for AAMA in urinary
par(mfrow=c(2,2))
plot(out[,'time'], out[,'m_AAMA_urinary'],type = 'l', ylim = c(0,.03) )
points(yobs_urine$time, yobs_urine$AAMA , col="blue", lwd = 4)
time_points_measure_unrine = c(1, 40, 84, 141, 196, 280 , 371, 460)
tamtam = out[,'m_AAMA_urinary'][time_points_measure_unrine ]
plot(yobs_urine$time, cumsum( tamtam ),type = 'l', ylim = c(0,.1) ); grid()
points(yobs_urine$time, cumsum(yobs_urine$AAMA) , col="blue", lwd = 4)

# plot for GAMA 
plot(out[,'time'], out[,'m_GAMA_urinary'],type = 'l', ylim = c(0,.002) )
points(yobs_urine$time, yobs_urine$GAMA , col="blue", lwd = 4); 
time_points_measure_unrine = c(1, 40, 84, 141, 196, 280 , 371, 460)
tamtam = out[,'m_GAMA_urinary'][time_points_measure_unrine ]
plot(yobs_urine$time, cumsum( tamtam ),type = 'l', ylim = c(0,.01) ); grid()
points(yobs_urine$time, cumsum(yobs_urine$GAMA) , col="blue", lwd = 4)







################################################################
# test new data Wang et al 2017 https://link.springer.com/article/10.1007/s00204-016-1869-6
# data 1, in women  ##############################################
yobs_urine <- data.frame(
  time = c(2, 3.5, 5.6, 7.9, 9.7, 14, 24.1, 36, 48.2),
  AAMA = c(220.4, 574.5, 1291.9, 1795.0, 1925.4, 2046.5, 1161.4, 630.4, 472.0)*234.28/1e6 ,
  GAMA = c( 29 , 33.3, 50.3, 79.2, 105.1, 118.5, 119.2, 70.3, 66.6)*250.27/1e6
) 
diet <- data.frame(var = "m_AA_dose",
                   time = c(1),  value = 12.6 *BW , method = "add")# dose of 12.6 microg/kg bw
times <- seq(from = 0, to = 50, by = 0.1)
out <- ode(y = yini, times = times, func = PBPKmodelAA, parms = params, events = list(data = diet), atol = 1e-6, rtol= 1e-8)
# plot for AAMA in urinary
par(mfrow=c(2,2))
plot(out[,'time'], out[,'m_AAMA_urinary'],type = 'l', ylim = c(0, 1 ) )
points(yobs_urine$time, yobs_urine$AAMA , col="blue", lwd = 4)
time_points_measure_unrine = c(2, 3.5, 5.6, 7.9, 9.7, 14, 24.1, 36, 48.2)*10
tamtam = out[,'m_AAMA_urinary'][time_points_measure_unrine ]
plot(yobs_urine$time, cumsum( tamtam ) , type = 'l', ylim = c(0, 4) ); grid()
points(yobs_urine$time, cumsum(yobs_urine$AAMA) , col="blue", lwd = 4)

# plot for GAMA 
plot(out[,'time'], out[,'m_GAMA_urinary'],type = 'l', ylim = c(0,.05) )
points(yobs_urine$time, yobs_urine$GAMA , col="blue", lwd = 4)
tamtam = out[,'m_GAMA_urinary'][time_points_measure_unrine ]
plot(yobs_urine$time, cumsum( tamtam ),type = 'l', ylim = c(0,.3) ); grid()
points(yobs_urine$time, cumsum(yobs_urine$GAMA) , col="blue", lwd = 4)






################################################################
# manual readout from Kopp and Dekant 2009
################################################################
yobs_urine <- data.frame(
  time = c(0.0, 4.2, 9.5, 14.1, 22.1, 30.5, 38.2),
  AAMA = c(0.0/1e6, 770.7/1e6, 2556.3/1e6, 1837.9/1e6, 1710.3/1e6, 892.2/1e6, 448.0/1e6)*234.28, 
  GAMA = c(0.0/1e6, 38.1/1e6, 129.7/1e6, 190.8/1e6, 297.7/1e6, 175.5/1e6, 129.7/1e6)*250.27
)
diet <- data.frame(var = "m_AA_dose", method = "add",
                   time = c(0),  value = 20 *BW )# dose of 20 microg/kg bw
times <- seq(from = 0, to = 50, by = 0.1)
out <- ode(y = yini, times = times, func = PBPKmodelAA, parms = params, events = list(data = diet), atol = 1e-6, rtol= 1e-8)
# plot for AAMA in urinary
par(mfrow=c(2,2))
plot(out[,'time'], out[,'m_AAMA_urinary'],type = 'l', ylim = c(0, 2 ) )
points(yobs_urine$time, yobs_urine$AAMA , col="blue", lwd = 4)
time_points_measure_unrine = c(1, 42, 84, 141, 196, 280 , 371)
tamtam = out[,'m_AAMA_urinary'][time_points_measure_unrine ]
plot(yobs_urine$time, cumsum( tamtam ) , type = 'l', ylim = c(0, 3) ); grid()
points(yobs_urine$time, cumsum(yobs_urine$AAMA) , col="blue", lwd = 4)

# plot for GAMA 
plot(out[,'time'], out[,'m_GAMA_urinary'],type = 'l', ylim = c(0,.1) )
points(yobs_urine$time, yobs_urine$GAMA , col="blue", lwd = 4)
time_points_measure_unrine = c(1, 40, 84, 141, 196, 280 , 371 )
tamtam = out[,'m_GAMA_urinary'][time_points_measure_unrine ]
plot(yobs_urine$time, cumsum( tamtam ),type = 'l', ylim = c(0,.3) ); grid()
points(yobs_urine$time, cumsum(yobs_urine$GAMA) , col="blue", lwd = 4)


