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
QCC = 16        #12.5 # maybe 14 or 6.5 check?
Q_C = QCC*BW^0.75
FQ_P = 2.5/100 
FQ_Li = 0.19
FQ_Ki = 0.255 #total - there is also an aterial flow rate 
FQ_T = 1 -( FQ_Li + FQ_Ki +FQ_P)
Q_Li <- Q_C*FQ_Li
Q_Ki <- Q_C*FQ_Ki
Q_T <- Q_C*FQ_T

# partition coefficient
pAA_TB = 0.2
pAA_KiB = 0.2
pAA_LiB = 1.5

pGA_TB = 0.97
pGA_LiB = 0.9
pGA_KiB = .3

#reaction rate constants 
k_AAuptake = 1 #1/h

k_onAA_T = 0.28 #1/h
k_onAA_B = 0.36 #1/h
k_onAA_Li = 0.71 #1/h
k_onAA_Ki = 0.83 #1/h

k_syn_GSH = 0.25 #mmol/h 
k_cl_GSH = 0.35 #1/h
k_onAA_GSH = 0.55 #L/(mmol h)
k_onGA_GSH = 0.8 #L/(mmol h)

k_onGA_T = 0.089 #1/h
k_onGA_B = 0.108 #1/h
k_onGA_Ki = k_onAA_T/2 #1/h
k_onGA_Li = 0.215 #1/h

k_exc_AAMA = 0.049      # 0.049 #1/h or 0.13 # this will change the shape
k_exc_GAMA = 0.027     # 0.077 #1/h  # or 0.027
k_exc_GAOH = 0.077     # 1/h
k_exc_GA <- 1.48       # 1/h Q_Ki*0.025
k_exc_AA <- 2

PL1 = 1.7     # !$'liver/blood partition AA'
PL2 = 0.9    # !$'liver/blood partition GA'

MW_aa = 71   # mg/mmol
MW_ga = 87   # mg/mmol

MW_GSH = 307.32 # mg/mmol

VMAXGC1 = 1   #!Vmax for GSH-AA conjugation mg/hr-kg^0.7
VMAXG1 = VMAXGC1/MW_aa*BW**0.7    # !$'Liver AA-GSH rate'

# Maximum velocity for enzymatic reaction 
V_max_p450 = 9 /MW_aa*BW^0.7  # 0.235 mg/h (Tien: AA to GA mmol/hr)
V_max_EH = 20 /MW_ga*BW^0.7  
KM_p450 = 10.0
KM_EH = 100.0

k_0_GSH = 7

AGSH0 = k_0_GSH*V_Li

KMG1 = 100/MW_aa  #!Km with respect to AA for GSH conjugation mM
KMGG = 0.1        # !KM with respect to GSH for AA or GA conjugation with GSH mM
KMG2 = 100/MW_ga  #!Km with respect to GA for GSH conjugation mM


KPT_Li = 0.015     # !'protein turnover rate in liver'
KPT_Ki = 0.013     #  !'protein turnover rate in kidney'
KPTR = 0.013     # !'protein turnover rate in rpt'
KPTS = 0.0039    # !'protein turnover rate in spt'

#Trine: Why don't we use the protein turnover rate for blood?
# Trine: Since KPTRB and KPTPL are the same, could we combine ethis as a turnover rate in blood?
KPTRB = 0.0039   # !'protein turnover rate in rbc'
KPTPL = 0.0039   # !'protein turnover rate in plasma'
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# create list of parameter
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
params <- unlist(c(data.frame(pAA_TB, pAA_LiB, pAA_KiB,
                              pGA_TB, pGA_LiB, pGA_KiB,
                              k_onAA_T, k_onAA_B, k_onAA_Li, k_onAA_GSH, k_onAA_Ki,
                              k_onGA_T, k_onGA_B, k_onGA_Li, k_onGA_Ki, k_exc_GA,
                              k_onGA_GSH, k_exc_AAMA, k_exc_GAMA, k_exc_GAOH,
                              k_cl_GSH, V_max_p450, KM_p450, V_max_EH, KM_EH,
                              k_AAuptake)
))
yini <- c(m_AA_AB = 0.0, m_GA_AB = 0.0, m_AA_VB = 0.0, m_GA_VB = 0.0,
          m_AA_Ki = 0.0, m_GA_Ki = 0.0, m_AA_dose = 0 ,
          m_AA_Li = 0.0, m_AAMA = 0.0,  #m_AAMA_urinary = 0.0,
          m_GA_Li = 0.0,
          a_pb_GA_Li   = 0.0,
          m_GAMA  = 0.0,
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
          m_GA_in = 0.0,
          a_pb_AA_Ki = 0,
          a_pb_AA_Li = 0,
          a_pb_GA_Ki = 0
)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PBPK model for Acrylamide 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PBPKmodelAA <- function(t,state,parameter){
  with(as.list(c(t,state,parameter)), {
    ########################################################################################  
    #----#---#  Model for aa  #---#---#
    #concentrations blood in the organs in mg/L
    c_AA_AB <- m_AA_AB/V_AB ;    c_AA_VB <- m_AA_VB/V_VB
    c_AA_T <- m_AA_T/V_T
    c_AA_Li <- m_AA_Li/(V_Li*PL1)
    c_AA_Ki <- m_AA_Ki/V_Ki
    
    c_GSH_Li <- m_GSH_Li/V_Li
    
    # Blood
    #-----------------------------
    # units checked -> mg/h
    #Trine: Why is the equation for the AB built up differently than VB? I suggest it to be simmilar, and according to the 
    # overall figure of our model.
    dm_AA_AB <- Q_C*c_AA_VB - Q_T*(c_AA_T/pAA_TB) - Q_Li*(c_AA_Li/pAA_LiB) - Q_Ki*(c_AA_Ki/pAA_KiB) -k_onAA_B*m_AA_AB
    # units checked -> mg/h   
    dm_AA_VB <- Q_T*(c_AA_T/pAA_TB) +Q_Li*(c_AA_Li/pAA_LiB) +Q_Ki*(c_AA_Ki/pAA_KiB) - Q_C*c_AA_VB -k_onAA_B*m_AA_VB
    
    # Trine: I added this. Can we make an overall turnover rate in blood? Is this correct? 
    da_pb_AA_B <- k_onAA_B*m_AA_AB + k_onAA_B*m_AA_VB - da_pb_AA_B*KPTRB
    
    # Kidney
    #---------------------------------
    # units checked -> mg/h   
    dm_AA_Ki <- Q_Ki*c_AA_AB -Q_Ki*(c_AA_Ki/pAA_KiB) -k_onAA_Ki*m_AA_Ki 
    
    # protein tunr over AA in Kidney
    da_pb_AA_Ki <- k_onAA_Ki*m_AA_Ki - a_pb_AA_Ki*KPT_Ki
    
    # Liver
    #--------------------------------
    # units checked -> mg/h   
    dm_AA_dose <- -k_AAuptake*m_AA_dose
    # units checked -> mg/h
    dm_GSH_Li <- AGSH0 - k_cl_GSH*m_GSH_Li 
    
    metAA_GSH <- VMAXG1 *c_AA_Li * c_GSH_Li /(c_AA_Li + KMG1)/(c_GSH_Li + KMGG)
    metAA_P450 <- V_max_p450 * MW_aa*c_AA_Li/ (KM_p450+c_AA_Li)
    
    # units checked -> mg/h  
    dm_AA_Li <- Q_Li*(c_AA_AB - c_AA_Li/pAA_LiB) + k_AAuptake*m_AA_dose - k_onAA_Li*m_AA_Li  - metAA_P450 - metAA_GSH 
    
    
    # protein tunr over AA in Liver
    da_pb_AA_Li <- k_onAA_Li*m_AA_Li - a_pb_AA_Li * KPT_Li
    
    # Tissue
    #--------------------------------------------------------------
    # units checked -> mg/h   
    dm_AA_T <- Q_T*(c_AA_AB - c_AA_T/pAA_TB) -k_onAA_T*m_AA_T
    
    # Urine
    #-------------------------------------------------------------
    # units checked -> mg/h   
    dm_AAMA <- metAA_P450 + metAA_GSH  - m_AAMA*k_exc_AAMA
    
    #Trine: We need to include the protein turnover also in urine. See Sweeney
    
    dm_AA_out <- metAA_P450 + metAA_GSH
    dm_AA_in <- k_AAuptake*m_AA_dose
    dm_AA_accum <- k_onAA_B*m_AA_AB +k_onAA_B*m_AA_VB +k_onAA_T*m_AA_T +k_onAA_Li*m_AA_Li +k_onAA_Ki*m_AA_Ki
    dm_AA_free <- dm_AA_AB + dm_AA_VB +dm_AA_T +dm_AA_Li +dm_AA_Ki
########################################################################################  
########################################################################################    
    #----#---#  Model for Ga  #---#---#    
    #concentrations blood in the organs in mg/L
    c_GA_AB <- m_GA_AB/V_AB ;   c_GA_VB <- m_GA_VB/V_VB
    c_GA_T <- m_GA_T/V_T
    c_GA_Li <- m_GA_Li/(V_Li*PL2)
    c_GA_Ki <- m_GA_Ki/V_Ki
    # units checked -> mg/h   
    dm_GA_AB <- Q_C*(c_GA_VB - c_GA_AB) -k_onGA_B*m_GA_AB
    # units checked -> mg/h   
    dm_GA_VB <- Q_T*(c_GA_T/pGA_TB) +Q_Li*(c_GA_Li/pGA_LiB) +Q_Ki*(c_GA_Ki/pGA_KiB) - Q_C*c_GA_VB -k_onGA_B*m_GA_VB
    # units checked -> mg/h   
    dm_GA_Ki <- Q_Ki*c_GA_AB -Q_Ki*(c_GA_Ki/pGA_KiB) -k_onGA_Ki*m_GA_Ki 
    
    
    # protein tunr over AA in Kidney
    da_pb_GA_Ki <- k_onGA_Ki*m_GA_Ki - a_pb_GA_Ki*KPT_Ki
    
    
    # units checked -> mg/h
    metGA_GSH <- k_onGA_GSH * c_GSH_Li *m_GA_Li / (c_GA_Li + KMG2)/(c_GSH_Li + KMGG)
    metGA_EH <- V_max_EH *MW_ga *c_GA_Li / (KM_EH + c_GA_Li)
    # units checked -> mg/h 
    dm_GA_Li <- Q_Li*(c_GA_AB - c_GA_Li/pGA_LiB) - k_onGA_Li*m_GA_Li + metAA_P450 - metGA_GSH -metGA_EH 
    
    
    # protein tunr over GA in Liver
    da_pb_GA_Li = k_onGA_Li*m_GA_Li - a_pb_GA_Li*KPT_Li 
    
    
    # units checked -> mg/h, 
    dm_GAMA <-  metGA_GSH  -m_GAMA*k_exc_GAMA
    
    dm_GA_T <- Q_T*(c_GA_AB - c_GA_T/pGA_TB) -k_onGA_T*m_GA_T
    
    dm_GAOH <- metGA_EH -m_GAOH*k_exc_GAMA
    
    dm_GA_out <- metGA_GSH + metGA_EH 
    dm_GA_in <- metGA_GSH +metGA_EH
    dm_GA_accum <- k_onGA_B*m_GA_AB +k_onGA_B*m_GA_VB +k_onGA_T*m_GA_T +k_onGA_Li*m_GA_Li +k_onGA_Ki*m_GA_Ki
    dm_GA_free <- dm_GA_AB + dm_GA_VB +dm_GA_T +dm_GA_Li +dm_GA_Ki 

    
    return(list(c(dm_AA_AB,  dm_GA_AB,  dm_AA_VB,    dm_GA_VB,
                  dm_AA_Ki,  dm_GA_Ki,  dm_AA_dose,
                  dm_AA_Li,
                  dm_AAMA,   
                  dm_GA_Li, da_pb_GA_Li, dm_GAMA,
                  dm_GSH_Li, dm_AA_T,   dm_GA_T,
                  dm_GAOH,
                  dm_AA_in,
                  dm_AA_out,
                  dm_AA_accum,
                  dm_AA_free,
                  dm_GA_out,
                  dm_GA_accum,
                  dm_GA_free,
                  dm_GA_in,
                  da_pb_AA_Ki,
                  da_pb_AA_Li,
                  da_pb_GA_Ki
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


plot(out)




