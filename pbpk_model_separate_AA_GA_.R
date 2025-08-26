library(deSolve)
library(minpack.lm)
graphics.off()

# body weight in kg
BW = 70 
# fractions of BW
##https://journals.sagepub.com/doi/pdf/10.1177/ANIB_32_3-4
F_B = 0.079 
F_B_AB = 0.35 #fraction: arterial blood of blood volume
F_B_VB = 0.65 #fraction: venous blood of blood volume
F_Li = 0.026
F_Ki = 0.0044
F_T = 1-(F_B+F_Li+F_Ki)
# organ volumes in L
V_AB <- BW*F_B*F_B_AB 
V_VB <- BW*F_B*F_B_VB 
V_T <- BW*F_T
V_Ki <- BW*F_Ki
V_Li <- BW*F_Li
# cardiac output in L/(kg*h)
##https://journals.sagepub.com/doi/pdf/10.1177/ANIB_32_3-4
QCC = 16        #12.5 # maybe 14 or 6.5 check?
Q_C = QCC*BW^0.75
#FQ_P = 2.5/100 Trine: We don't have a lung compartment in this model
FQ_Li = 0.255
FQ_Ki = 0.19 #total - there is also an aterial flow rate 
FQ_T = 1 -( FQ_Li + FQ_Ki) # +FQ_P) Trine: We don't have a lung compartment in this model
Q_Li <- Q_C*FQ_Li
Q_Ki <- Q_C*FQ_Ki
Q_T <- Q_C*FQ_T

MW_aa = 71   # mg/mmol
MW_ga = 87   # mg/mmol
MW_GSH = 307.32 # mg/mmol
k_0_GSH = 7   # from Sweeney paper.
MW_GSH = 307.32     # mg/mmol
AGSH0 = k_0_GSH * V_Li * MW_GSH

# Maria: I came across this paper and has a lot of interesting toxicokinetc data. Humans are
# exposed to acrylamide and measure urine metabolites. They measured 60% or possibly more.
# of ingested acrylamide is absorbed. 
# We can use half-lives, metabolite fractions, pathway ratios etc. to validate the PBPK in humans.
# 'Toxicokinetics of acrylamide in humans after ingestion of a defined dose in a test meal to improve risk assessment for acrylamide carcinogenicity'
# DOI: 10.1158/1055-9965.EPI-05-0647

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# partition coefficient 
pAA_TB = .2        # Tissue::blood partition AA, Sweeney et al (2010) paper
pAA_KiB = 0.8      # Kidney::blood partition AA, Sweeney et al (2010) paper
pAA_LiB = 0.4      # liver::blood partition AA, Sweeney et al (2010) paper
# Maybe we said it earlier and these parameters are not sensitive to the model, but from which paper did you get the values of these partition coefs., 
# in the paper of Sweeney and of Li are different pAA_LiB = 0.4 and pAA_KiB = 0.8.

# reaction rate constants 
k_AAuptake = 1 #1/h

k_onAA_T = 0.028  #1/h     # we took it from Sweenney
k_onAA_B = 0.036  #1/h
k_onAA_Li = 0.71 #1/h
k_onAA_Ki = 0.13 #1/h

# from 0.23 to 0.35
k_cl_GSH = 0.035   # 1/h   #  reference from Maria

# k_onAA_GSH = 0.55 # L/(mmol h)
# check, we do not use it any more.
k_exc_AAMA = 0.049   # 0.049 #1/h or 0.13 # this will change the shape

# Maximum velocity for enzymatic reaction 
V_max_p450 = 9 /MW_aa*BW^0.7  # mg/h  from Sweenney paper, checked!
# maximum rate of GA and AA conjugation with GSH (mg/(h BW^0.7))
Vmax_AA_GSH <- 24*BW^0.7 # from the Sweeny  fitted   mg/h




#????????????????????????????????
#????????????????????????????????
#????????????????????????????????
#????????????????????????????????
#????????????????????????????????
#????????????????????????????????
### PARAMETERS FOR GA
# partition coefficient
pGA_TB = 0.97
pGA_LiB = 0.9
pGA_KiB = .3

#k_onGA_GSH = 0.8 #L/(mmol h)
k_onGA_T = 0.69        #1/h
k_onGA_B = 0.108       #1/h
k_onGA_Ki = k_onAA_T/2 #1/h
k_onGA_Li = 0.115      #1/h

k_exc_GAMA = 0.027     # 0.077 #1/h  # or 0.027 # from sweeney
k_exc_GAOH = 0.077     # 1/h  # from sweeney
# (removed as no longer used) k_exc_GA <- .4         # 1/h Q_Ki*0.025


VMAXGC1 = 22   # maximum rate of AA conjugation with GSH mg/hr-kg^0.7, Sweeney et al (2010) paper
VMAXG1 = VMAXGC1/MW_aa*BW**0.7    # !$'Liver AA-GSH rate'

# Maximum velocity for enzymatic reaction 
V_max_EH = 20 /MW_ga*BW^0.7  
KM_p450 = 10.0
KM_EH = 100.0
# maximum rate of GA and AA conjugation with GSH (mg/(h BW^0.7))
Vmax_GA_GSH <- 20*BW^0.7 # from the Sweeny  fitted   mg/h

KMG1 = 20  #!Km with respect to AA for GSH conjugation mg/L (from Sweeny code) 
# Trine: This should be mg. Delete the term MW_aa
# Sweeney has this parameter for rats ans 100 mg/L why did you use this number and why divived by MW??

KMGG = .1/MW_GSH        # !KM with respect to GSH for AA or GA conjugation with GSH mmol/L. Trine: should be mg/L
####and that this is 307 mg/L #This comment is wrong ignore
KMG2 = 90  #!Km with respect to GA for GSH conjugation mM. Trine: This should be mg. Delete the term MW_aa
# The same comment as for KMG1

# from the code, not from paper
KPT_Li = 0.015    # !'protein turnover rate in liver'# this is confimred in the exp. values of Sweeney
KPT_Ki = 0.013    # !'protein turnover rate in kidney'# 1/h this is 0.012 in Sweeney, I guess its not much of difference
KPTR = 0.013      # !'protein turnover rate in rpt'# this is 0.012 in Li et.al
KPTS = 0.0039     # !'protein turnover rate in spt'# and this they have it 0.0051 Sweeney and Li 0.0012 refering to Sweeney
KPTRB = 0.0039    # !'protein turnover rate in rbc'# how did we get that?
KPTPL = 0.012    # !'protein turnover rate in plasma'# Sweeney has 0.012 for this value
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# create list of parameter
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
params <- unlist(c(data.frame(pAA_TB, pAA_LiB, pAA_KiB,
                              pGA_TB, pGA_LiB, pGA_KiB,
                              k_onAA_T, k_onAA_B, k_onAA_Li,  k_onAA_Ki,
                              k_onGA_T, k_onGA_B, k_onGA_Li, k_onGA_Ki, 
                              k_exc_AAMA, k_exc_GAMA,
                              k_cl_GSH, V_max_p450, KM_p450, V_max_EH, KM_EH,
                              k_AAuptake)
                              ))
yini <- c(m_AA_AB = 0.0, m_GA_AB = 0.0,      a_pb_AA_B = 0 ,
          m_AA_VB = 0.0, m_GA_VB = 0.0,
          m_AA_Ki = 0.0, m_GA_Ki = 0.0, m_AA_dose = 0 ,
          m_AA_Li = 0.0, m_AAMA = 0.0, 
          m_GA_Li = 0.0, a_pb_GA_Li = 0.0,
          m_GAMA  = 0.0,       
          m_GSH_Li = k_0_GSH * V_Li * MW_GSH, 
          m_AA_T = 0.0, m_GA_T = 0.0,
          m_GAOH = 0.0, 
          m_AA_out = 0.0, m_AA_in = 0.0, m_AA_accum = 0.0, m_AA_free = 0.0, 
          m_GA_out = 0.0, m_GA_in = 0.0, m_GA_accum = 0.0, m_GA_free = 0.0, 
          a_pb_AA_Ki = 0,  a_pb_AA_Li = 0,
          a_pb_GA_Ki = 0, a_pb_AA_T = 0, a_pb_GA_B = 0, a_pb_GA_T = 0
)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PBPK model for Acrylamide 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PBPKmodelAA <- function(t,state,parameter){
  with(as.list(c(t,state,parameter)), {
    #----#---#  Model for AA  #---#---#
    # concentrations blood in the organs in mg/L
    c_AA_AB <- m_AA_AB/V_AB   
    c_AA_VB <- m_AA_VB/V_VB
    c_AA_T <- m_AA_T/V_T
    c_AA_Li <- m_AA_Li/V_Li #Trine: Should we multiply with PL1 as we do for GA
    c_AA_Ki <- m_AA_Ki/V_Ki
    
    c_GSH_Li <- m_GSH_Li/V_Li
    # units checked -> mg/h   
    dm_AA_dose <- - k_AAuptake*m_AA_dose
    
    # Blood
    #-----------------------------
    # units checked -> mg/h   
    dm_AA_AB <- Q_C*(c_AA_VB - c_AA_AB) - k_onAA_B*m_AA_AB
    dm_AA_VB <- Q_T*(c_AA_T/pAA_TB) + Q_Li*(c_AA_Li/pAA_LiB) + Q_Ki*(c_AA_Ki/pAA_KiB) - Q_C*c_AA_VB -k_onAA_B*m_AA_VB
    # Trine: I added this. Can we make an overall turnover rate in blood? Is this correct? 
    da_pb_GA_B <- k_onAA_B*m_AA_AB + k_onAA_B*m_AA_VB - a_pb_GA_B*KPTRB
    
    # Kidney
    #---------------------------------
    # units checked -> mg/h   
    dm_AA_Ki <- Q_Ki*c_AA_AB -Q_Ki*(c_AA_Ki/pAA_KiB) - k_onAA_Ki*m_AA_Ki 
    # protein tunr over AA in Kidney
    da_pb_AA_Ki <- k_onAA_Ki*m_AA_Ki - a_pb_AA_Ki*KPT_Ki
   
    # unit chekced mg/h 
    metAA_GSH <- Vmax_AA_GSH * c_AA_Li * c_GSH_Li / ( (KMG1 + c_AA_Li) * (KMGG + c_GSH_Li) )

    metAA_P450 <- V_max_p450 * MW_aa * c_AA_Li/ (KM_p450 + c_AA_Li)
    
    # units checked -> mg/h  
    dm_AA_Li <- Q_Li*(c_AA_AB - c_AA_Li/pAA_LiB) + k_AAuptake*m_AA_dose - k_onAA_Li*m_AA_Li  - metAA_P450 - metAA_GSH 
    # protein turn over AA in Liver
    da_pb_AA_Li <- k_onAA_Li*m_AA_Li - a_pb_AA_Li * KPT_Li
    
    # Tissue
    #--------------------------------------------------------------
    # units checked -> mg/h   
    dm_AA_T <- Q_T*(c_AA_AB - c_AA_T/pAA_TB) - k_onAA_T*m_AA_T
    # protein turn over AA in Tussue
    da_pb_AA_T <- k_onAA_T*m_AA_T - a_pb_AA_T * KPTS
    
    # Urine
    #-------------------------------------------------------------
    # units checked -> mg/h   
    dm_AAMA <- metAA_GSH - m_AAMA*k_exc_AAMA 

    dm_AA_out <- metAA_P450 + metAA_GSH + k_onAA_Li*m_AA_Li
    dm_AA_in <- k_AAuptake * m_AA_dose
    dm_AA_accum <- k_onAA_B*m_AA_AB +k_onAA_B*m_AA_VB +k_onAA_T*m_AA_T +k_onAA_Ki*m_AA_Ki
    dm_AA_free <- dm_AA_AB + dm_AA_VB +dm_AA_T +dm_AA_Li +dm_AA_Ki

      #--------------------------------------------------------------
    ########################################################################################    
    #-GA-#-GA-#  Model for GA  #-GA-#-GA-#    
    #concentrations blood in the organs in mg/L
    c_GA_AB <- m_GA_AB/V_AB ;   
    c_GA_VB <- m_GA_VB/V_VB
    c_GA_T <- m_GA_T/V_T
    c_GA_Li <- m_GA_Li/V_Li
    c_GA_Ki <- m_GA_Ki/V_Ki
    
    # Blood
    #-----------------------------
    # units checked -> mg/h   
    dm_GA_AB <- Q_C*(c_GA_VB - c_GA_AB) -k_onGA_B*m_GA_AB
    # units checked -> mg/h   
    dm_GA_VB <- Q_T*(c_GA_T/pGA_TB) + Q_Li*(c_GA_Li/pGA_LiB) + Q_Ki*(c_GA_Ki/pGA_KiB) - Q_C*c_GA_VB - k_onGA_B*m_GA_VB
    # units checked -> mg/h   
    da_pb_AA_B <- k_onGA_B*m_AA_AB + k_onGA_B*m_AA_VB - a_pb_GA_B*KPTRB
    
    # Kidney
    #---------------------------------
    dm_GA_Ki <- Q_Ki*c_GA_AB -Q_Ki*(c_GA_Ki/pGA_KiB) -k_onGA_Ki*m_GA_Ki 
    # protein tunr over AA in Kidney
    da_pb_GA_Ki <- k_onGA_Ki*m_GA_Ki - a_pb_GA_Ki*KPT_Ki
    
    # Liver
    #--------------------------------
    # units checked -> mg/h
    metGA_GSH <- Vmax_GA_GSH * c_GSH_Li *c_GA_Li / (c_GA_Li + KMG2)/(c_GSH_Li + KMGG)
    metGA_EH <- V_max_EH *MW_ga *c_GA_Li / (KM_EH + c_GA_Li)
    # units checked -> mg/h 
    dm_GA_Li <- Q_Li*(c_GA_AB - c_GA_Li/pGA_LiB) - k_onGA_Li*m_GA_Li + metAA_P450 - metGA_GSH -metGA_EH 
    # protein tunr over GA in Liver
    da_pb_GA_Li = k_onGA_Li*m_GA_Li - a_pb_GA_Li*KPT_Li 
    
    # units checked -> mg/h, 
    dm_GAMA <-  metGA_GSH  -m_GAMA*k_exc_GAMA
    
    #????????????????????????????????????????????
    # we need a new k_exc_GAOH
    dm_GAOH <- metGA_EH - m_GAOH*k_exc_GAOH
    
    # Tissue    #-----------------------------------------------------
    dm_GA_T <- Q_T*(c_GA_AB - c_GA_T/pGA_TB) - k_onGA_T*m_GA_T
    da_pb_GA_T <- k_onGA_T*m_GA_T - a_pb_GA_T * KPTS
    
    # Liver
    #--------------------------------
    # units checked -> mg/h
    dm_GSH_Li <- k_0_GSH * V_Li * MW_GSH - k_cl_GSH*m_GSH_Li - metAA_GSH - metGA_GSH
    
    dm_GA_out <- metGA_GSH + metGA_EH 
    dm_GA_in <- metAA_P450
    dm_GA_accum <- k_onGA_B*m_GA_AB +k_onGA_B*m_GA_VB +k_onGA_T*m_GA_T +k_onGA_Li*m_GA_Li +k_onGA_Ki*m_GA_Ki
    dm_GA_free <- dm_GA_AB + dm_GA_VB +dm_GA_T +dm_GA_Li +dm_GA_Ki 
    
    return(list(c(dm_AA_AB,    dm_GA_AB,   da_pb_AA_B,
                  dm_AA_VB,    dm_GA_VB,
                  dm_AA_Ki,    dm_GA_Ki,    dm_AA_dose,
                  dm_AA_Li,    dm_AAMA,   
                  dm_GA_Li,   da_pb_GA_Li, dm_GAMA,
                  dm_GSH_Li,  dm_AA_T,     dm_GA_T,
                  dm_GAOH,
                  dm_AA_out,  dm_AA_in,   dm_AA_accum,  dm_AA_free,
                  dm_GA_out,  dm_GA_in,   dm_GA_accum,  dm_GA_free,
                  da_pb_AA_Ki,    da_pb_AA_Li,
                  da_pb_GA_Ki,    da_pb_AA_T, da_pb_GA_B , da_pb_GA_T   )))
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
plot(out[,'time'], out[,'m_GAMA']  ,type = 'l',xlab = '', ylab = 'GAMA', ylim = c(0,.005) )
points(yobs_urine$time, yobs_urine$GAMA , col="blue", cex = 1.5, pch = 17); grid()
time_points_measure_unrine = c(1, 40, 84, 141, 196, 280 , 371, 460)
tamtam = out[,'m_GAMA'][time_points_measure_unrine ]
plot(yobs_urine$time, cumsum( tamtam ),type = 'l',xlab = '', ylab = 'GAMA', ylim = c(0,.04) ); grid()
points(yobs_urine$time, cumsum(yobs_urine$GAMA) , col="blue",cex = 1.5, pch = 17)



