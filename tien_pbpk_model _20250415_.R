library(deSolve)
library(minpack.lm)

graphics.off()

#body weight in kg
BW = 70 

# fractions of BW
##https://journals.sagepub.com/doi/pdf/10.1177/ANIB_32_3-4
F_B = 0.079 
F_B_AB = 0.35 #fraction: arterial blood of blood volume
F_B_VB = 0.65 #fraction: venous blood of blood volume
F_Li = 0.026
F_Ki = 0.0044
F_T = 1-(F_B+F_Li+F_Ki)

#organ volumes in L
V_AB <- BW*F_B*F_B_AB 
V_VB <- BW*F_B*F_B_VB 
V_T <- BW*F_T
V_Ki <- BW*F_Ki
V_Li <- BW*F_Li

#cardiac output in L/(kg*h)
##https://journals.sagepub.com/doi/pdf/10.1177/ANIB_32_3-4
QCC = 16 #12.5 # maybe 14 or 6.5 check?
Q_C = QCC*BW^0.75
FQ_P = 2.5/100 
Q_P <- Q_C*FQ_P # 6.5*BW 

FQ_Li = 0.19
FQ_Ki = 0.255 #total - there is also an aterial flow rate 
FQ_T = 1 -( FQ_Li + FQ_Ki +FQ_P)
Q_Li <- Q_C*FQ_Li
Q_Ki <- Q_C*FQ_Ki
Q_T <- Q_C*FQ_T

# partition coefficient
pAA_TB = 0.4
pAA_KiB = 0.8
pAA_LiB = 0.3
pAA_LuAir = 31000000.0

pGA_TB = 0.97
pGA_LiB = 0.5
pGA_KiB = 1.0
pGA_LuAir = 98000000.0

#reaction rate constants 
k_AAuptake = 1 #1/h

k_onAA_T = 0.028 #1/h
k_onAA_B = 0.0036 #1/h
k_onAA_Li = 0.071 #1/h
k_onAA_Ki = 0.13 #1/h

k_syn_GSH = 0.025 #mmol/h 
k_cl_GSH = 0.35 #1/h
k_onAA_GSH = 0.55 #L/(mmol h)
k_onGA_GSH = 0.8 #L/(mmol h)

k_onGA_T = 0.089 #1/h
k_onGA_B = 0.0108 #1/h
k_onGA_Ki = k_onAA_T/2 #1/h
k_onGA_Li = 0.215 #1/h

k_exc_AAMA = 0.029      # 0.049 #1/h or 0.13 # this will change the shape
k_exc_GAMA = 0.017     # 0.077 #1/h  # or 0.027
k_exc_GAOH = 0.077     # 1/h
k_exc_GA <- 2.48       # 1/h Q_Ki*0.025
k_exc_AA <- 2


MW_aa = 71   # mg/mmol
MW_ga = 87   # mg/mmol

MW_GSH = 307.32 # mg/mmol

# Maximum velocity for enzymatic reaction 
V_max_p450 = 9 /MW_aa*BW^0.7  # 0.235 mg/h (Tien: AA to GA mmol/hr)
V_max_EH = 20 /MW_ga*BW^0.7  
KM_p450 = 10.0
KM_EH = 100.0

c_AA_inh = 0
c_GA_inh = 0

k_0_GSH = 7
### CHECK  WHY     m_GSH_Li = 7.0*V_Li, Tien: 7.0 as in Sweeney is the Initial GSH (mmol/l)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# create list of parameter
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
params <- unlist(c(data.frame(pAA_TB, pAA_LiB, pAA_LuAir, pAA_KiB,
                              pGA_TB, pGA_LiB, pGA_LuAir, pGA_KiB,
                              k_onAA_T, k_onAA_B, k_onAA_Li, k_onAA_GSH, k_onAA_Ki,
                              k_onGA_T, k_onGA_B, k_onGA_Li, k_onGA_Ki, k_exc_GA,
                              k_onGA_GSH, k_exc_AAMA, k_exc_GAMA, k_exc_GAOH,
                              k_cl_GSH, V_max_p450, KM_p450, V_max_EH, KM_EH,
                              c_AA_inh, c_GA_inh, k_AAuptake)
                          ))
yini <- c(m_AA_AB = 0.0, m_GA_AB = 0.0, m_AA_VB = 0.0, m_GA_VB = 0.0,
          m_AA_Ki = 0.0, m_GA_Ki = 0.0, m_AA_dose = 0 ,
          m_AA_Li = 0.0, m_AAMA = 0.0,  #m_AAMA_urinary = 0.0,
          m_GA_Li = 0.0, m_GAMA = 0.0,
          #m_GAMA_urinary = 0.0, 
          m_GSH_Li = k_0_GSH * V_Li * MW_GSH, 
          m_AA_T = 0.0, m_GA_T = 0.0,
          m_GAOH = 0.0, 
          m_AA_out = 0.0,
          m_AA_accum = 0.0,
          m_AA_free = 0.0,
          m_AA_in = 0.0,
          m_GA_out = 0.0,
          m_GA_accum = 0.0,
          m_GA_free = 0.0,
          m_GA_in = 0.0
)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PBPK model for Acrylamide 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PBPKmodelAA <- function(t,state,parameter){
  with(as.list(c(t,state,parameter)), {
    #concentrations blood in the organs in mg/L
    c_AA_AB <- m_AA_AB/V_AB
    c_AA_VB <- m_AA_VB/V_VB
    c_AA_T <- m_AA_T/V_T
    c_AA_Li <- m_AA_Li/V_Li
    c_AA_Ki <- m_AA_Ki/V_Ki
    
    c_GSH_Li <- m_GSH_Li/V_Li
    
    c_GA_AB <- m_GA_AB/V_AB
    c_GA_VB <- m_GA_VB/V_VB
    c_GA_T <- m_GA_T/V_T
    c_GA_Li <- m_GA_Li/V_Li
    c_GA_Ki <- m_GA_Ki/V_Ki
    
    #ODEs for amount of AA in compartment in mg
    c_AA_Lu <- (Q_C*c_AA_VB +Q_P*c_AA_inh)/(Q_C +(Q_P/pAA_LuAir))
    c_GA_Lu <- (Q_C*c_GA_VB +Q_P*c_GA_inh)/(Q_C +(Q_P/pGA_LuAir))
    # units checked -> mg/h   
    dm_AA_AB <- Q_C*(c_AA_Lu - c_AA_AB) -k_onAA_B*m_AA_AB
    dm_GA_AB <- Q_C*(c_GA_Lu - c_GA_AB) -k_onGA_B*m_GA_AB
    # units checked -> mg/h   
    dm_AA_VB <- Q_T*(c_AA_T/pAA_TB) +Q_Li*(c_AA_Li/pAA_LiB) +Q_Ki*(c_AA_Ki/pAA_KiB) - Q_C*c_AA_VB -k_onAA_B*m_AA_VB
    dm_GA_VB <- Q_T*(c_GA_T/pGA_TB) +Q_Li*(c_GA_Li/pGA_LiB) +Q_Ki*(c_GA_Ki/pGA_KiB) - Q_C*c_GA_VB -k_onGA_B*m_GA_VB
    # units checked -> mg/h   
    dm_AA_Ki <- Q_Ki*c_AA_AB -Q_Ki*(c_AA_Ki/pAA_KiB) -k_onAA_Ki*m_AA_Ki 
    dm_GA_Ki <- Q_Ki*c_GA_AB -Q_Ki*(c_GA_Ki/pGA_KiB) -k_onGA_Ki*m_GA_Ki 
    # units checked -> mg/h   
    dm_AA_dose <- -k_AAuptake*m_AA_dose
    
    metAA_GSH <- k_onAA_GSH*c_GSH_Li*m_AA_Li / MW_GSH
    metAA_P450 <- V_max_p450 * MW_aa*c_AA_Li/ (KM_p450+c_AA_Li)
    # units checked -> mg/h  
    #   RAL1 = QL*CA1 + KA*STOM - QL*CVL1 - RP450 - RGST1 - rpbl1
    dm_AA_Li <- Q_Li*(c_AA_AB - c_AA_Li/pAA_LiB) + k_AAuptake*m_AA_dose - 
                  k_onAA_Li*m_AA_Li  - metAA_P450 - metAA_GSH
    # Trine: I think the two last parts of the equation should be subtracted. 
    # Metabolism of AA by P450 and EH will remove AA
    
    # units checked -> mg/h   
    dm_AAMA <- metAA_P450  - m_AAMA*k_exc_AAMA
    # units checked -> mg/h
    metGA_GSH <- k_onGA_GSH*c_GSH_Li*m_GA_Li / MW_GSH
    metGA_EH <- (V_max_EH *MW_ga *c_GA_Li)/(KM_EH+c_GA_Li)
    # units checked -> mg/h 
    #  RAL2 = QL*CA2 + KA2*STOM2 - QL*CVL2 + RP450 - RGST2 - REH - rpbl2
    # Trine: I think the two last parts of the equation should be subtracted. 
    # Metabolism of AA by P450 and EH will remove AA
    dm_GA_Li <- Q_Li*(c_GA_AB - c_GA_Li/pGA_LiB) - k_onGA_Li*m_GA_Li + 
                        metAA_P450 - metGA_GSH -metGA_EH
    # units checked -> mg/h, 
    # old: dm_GAMA <-  metGA_GSH - m_GAMA*k_exc_GAMA 
    #Trine: I do not think metGA_EH should be here. That's the metabolism to glyceramide, and will not contribute to GAMA
    # Trine: However I think, maybe, that the metAA_EH shold be here, as the metabolism of AA by EH is producing GA. 
    # Trine: Maybe production of glyceramide should be a separate equation? Need to discuss!!
    dm_GAMA <-  metGA_GSH  -m_GAMA*k_exc_GAMA
    # units checked -> mg  
    # units checked -> mg/h
    # Trine: I think the metAA_P450 and metAA_GSH should not be in this equation. 
    # The metGA_EH and metGA_GSH should be subtracted, as they will remove GA from the liver.
    dm_GSH_Li <-  -k_cl_GSH*m_GSH_Li + metAA_GSH - metGA_GSH
    # units checked -> mg/h   
    dm_AA_T <- Q_T*(c_AA_AB - c_AA_T/pAA_TB) -k_onAA_T*m_AA_T
    dm_GA_T <- Q_T*(c_GA_AB - c_GA_T/pGA_TB) -k_onGA_T*m_GA_T

    dm_GAOH <- metGA_EH -m_GAOH*k_exc_GAMA

    dm_AA_out <- metAA_P450 + metAA_GSH
    dm_AA_accum <- k_onAA_B*m_AA_AB +k_onAA_B*m_AA_VB +k_onAA_T*m_AA_T +k_onAA_Li*m_AA_Li +k_onAA_Ki*m_AA_Ki
    dm_AA_free <- dm_AA_AB + dm_AA_VB +dm_AA_T +dm_AA_Li +dm_AA_Ki
    dm_AA_in <- Q_P*c_AA_inh + k_AAuptake*m_AA_dose
    
    dm_GA_out <- metGA_GSH + metGA_EH 
    dm_GA_accum <- k_onGA_B*m_GA_AB +k_onGA_B*m_GA_VB +k_onGA_T*m_GA_T +k_onGA_Li*m_GA_Li +k_onGA_Ki*m_GA_Ki
    dm_GA_free <- dm_GA_AB + dm_GA_VB +dm_GA_T +dm_GA_Li +dm_GA_Ki 
    dm_GA_in <- Q_P*c_GA_inh + metGA_GSH +metGA_EH
    
    return(list(c(dm_AA_AB,  dm_GA_AB,  dm_AA_VB,    dm_GA_VB,
                  dm_AA_Ki,  dm_GA_Ki,  dm_AA_dose,  dm_AA_Li,
                  dm_AAMA,  # dm_AAMA_urinary,  
                  dm_GA_Li,  dm_GAMA,
                  #dm_GAMA_urinary,    
                  dm_GSH_Li, dm_AA_T,   dm_GA_T,
                  dm_GAOH,
                  dm_AA_out,
                  dm_AA_accum,
                  dm_AA_free,
                  dm_AA_in,
                  dm_GA_out,
                  dm_GA_accum,
                  dm_GA_free,
                  dm_GA_in
    )))
  })
}
PBPKmodelAA <- compiler::cmpfun(PBPKmodelAA)



################################################################
# manual readout from Kopp and Dekant 2009
################################################################
times <- seq(from = 0, to = 60, by = 0.1)
diet <- data.frame(var = "m_AA_dose", method = "add",
                   time = c(0.0),  value = 0.5 *BW /1000 )  # dose of 0.05 microg/kg bw
out <- ode(y = yini, times = times, func = PBPKmodelAA, parms = params, events = list(data = diet))
yobs_urine <- data.frame(
  time = c(0.0, 3.9, 8.3, 14, 19.5, 28, 37, 45.9),
  AAMA = c(0.0/1000000, 32.8/1000000, 69.5/1000000, 42.8/1000000, 45.2/1000000, 24.5/1000000, 15/1000000, 7/1000000)*234.28, #https://pubchem.ncbi.nlm.nih.gov/compound/N-Acetyl-S-_2-carbamoylethyl_-L-cysteine
  GAMA = c(0.0, 2.15*1e-6, 3.66*1e-6, 4.38*1e-6, 6.57*1e-6, 5.26*1e-6, 3.50*1e-6, 3.0*1e-6)*250.27
)
# plot for AAMA in urinary
par(mfrow=c(2,2),mar=c(2,4,.5,.5))
plot(out[,'time'], out[,'m_AAMA']   ,type = 'l',xlab = 'time', ylab = 'aama', ylim = c(0,.03) )
points(yobs_urine$time, yobs_urine$AAMA , col="blue", lwd = 4) ; grid()
time_points_measure_unrine = c(1, 40, 84, 141, 196, 280 , 371, 460)
tamtam = out[,'m_AAMA'][time_points_measure_unrine ]
plot(yobs_urine$time, cumsum( tamtam ),type = 'l', ylab = 'aama', ylim = c(0,.1), xlab = 'time' ); grid()
points(yobs_urine$time, cumsum(yobs_urine$AAMA) , col="blue", lwd = 4)

# plot for GAMA 
plot(out[,'time'], out[,'m_GAMA']  ,type = 'l',xlab = '', ylab = 'GAMA', ylim = c(0,.002) )
points(yobs_urine$time, yobs_urine$GAMA , col="blue", cex = 1.5, pch = 17); grid()
time_points_measure_unrine = c(1, 40, 84, 141, 196, 280 , 371, 460)
tamtam = out[,'m_GAMA'][time_points_measure_unrine ]
plot(yobs_urine$time, cumsum( tamtam ),type = 'l',xlab = '', ylab = 'GAMA', ylim = c(0,.01) ); grid()
points(yobs_urine$time, cumsum(yobs_urine$GAMA) , col="blue",cex = 1.5, pch = 17)







