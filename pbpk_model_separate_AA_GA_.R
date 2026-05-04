library(deSolve) ; library(minpack.lm)
graphics.off()

BW = 70    # body weight in kg
# fractions of BW # https://journals.sagepub.com/doi/pdf/10.1177/ANIB_32_3-4
F_B = 0.079 
F_B_AB = 0.35 # fraction: arterial blood of blood volume

F_B_VB = 0.65 #fraction: venous blood of blood volume
F_Li = 0.026
F_Ki = 0.0044
F_T = 1 - ( F_B + F_Li + F_Ki )
# organ volumes in L
V_AB <- BW*F_B*F_B_AB 
V_VB <- BW*F_B*F_B_VB 

V_T <- BW*F_T
V_Ki <- BW*F_Ki
V_Li <- BW*F_Li

#mixed venous delay compartment volume
V_mix <- V_VB

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

# partition coefficient 
pAA_TB = .2        # Tissue::blood partition AA, Sweeney et al (2010) paper
pAA_KiB = 0.8      # Kidney::blood partition AA, Sweeney et al (2010) paper
pAA_LiB = 0.4      # liver::blood partition AA, Sweeney et al (2010) paper
# Maybe we said it earlier and these parameters are not sensitive to the model,
# but from which paper did you get the values of these partition coefs., 
# in the paper of Sweeney and of Li are different pAA_LiB = 0.4 and pAA_KiB = 0.8.

# reaction rate constants 
k_AAuptake = 1 #1/h

k_onAA_T = 0.028  #1/h     # we took it from Sweenney
k_onAA_B = 0.003 #0.036  #1/h # 0.036 h⁻¹ corresponds to tissue binding (KPBS1)in Sweeney, not blood
k_onAA_Li = 0.71 #1/h
k_onAA_Ki = 0.13 #1/h

# from 0.23 to 0.35
k_cl_GSH = 0.035   # 1/h   #  Reed et al. 2008 'A mathematical model of glutathione metabolism'
# doi:10.1186/1742-4682-5-8

k_exc_AAMA = 0.049  ## ##MK Furh et al. report on Elimination half-life of AAMA as ~17.4 h, then k_exc_AAMA=0.0398 h-1#
# 'Toxicokinetics of acrylamide in humans after ingestion of a defined dose in a test meal to improve risk assessment for acrylamide carcinogenicity'
# DOI: 10.1158/1055-9965.EPI-05-0647
# 0.049 #1/h or 0.13 # this will change the shape

# Maximum velocity for enzymatic reaction 
V_max_p450 = 0.5 * BW^0.7 # 9 /MW_aa*BW^0.7  # mg/h  from Sweenney paper, checked!
# maximum rate of GA and AA conjugation with GSH (mg/(h BW^0.7))
Vmax_AA_GSH <- 22 * BW^0.7 # 24*BW^0.7 # from the Sweeny  fitted   mg/h


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PARAMETERS FOR GA
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# partition coefficient
pGA_TB = 1.35 # Have been calculated from the Doerge et al 2005 paper on TK in mice
pGA_LiB = 0.63 # Have been calculated from the Doerge et al 2005 paper on TK in mice
pGA_KiB = .3

#k_onGA_GSH = 0.8 #L/(mmol h)
k_onGA_T = 0.69        #1/h
k_onGA_B = 0.108       #1/h
k_onGA_Ki = k_onAA_T/2 #1/h
k_onGA_Li = 0.115      #1/h

k_exc_GAMA = 0.027   ##MK Furh et al. report on Elimination half-life of GAMA: ~25.1 hours, then k_exc_GAMA = 0.0276 h-1##  
# 0.077 #1/h  # or 0.027 # from sweeney # 
k_exc_GAOH = 0.077     # 1/h  # from sweeney
# (removed as no longer used) k_exc_GA <- .4         # 1/h Q_Ki*0.025

# Maximum velocity for enzymatic reaction 
V_max_EH = 20 *BW^0.7  
KM_p450 = 10.0
KM_EH = 100.0
# maximum rate of GA and AA conjugation with GSH (mg/(h BW^0.7))
Vmax_GA_GSH <- 20*BW^0.7 # from the Sweeny ( Vmax_GC2 )  fitted   mg/h

KMG1 = 100  # In Sweeney: KMG1 = 100 mg/L (human)  #50  #!Km with respect to AA for GSH conjugation mg/L (from Sweeny code) 
# Trine: This should be mg. Delete the term MW_aa
# Sweeney has this parameter for rats ans 100 mg/L why did you use this number and why divived by MW??

KMGG <- 0.1 * MW_GSH # In Sweeney: KMGG = 0.1 mM # .1/MW_GSH        # !KM with respect to GSH for AA or GA conjugation with GSH mmol/L. Trine: should be mg/L
####and that this is 307 mg/L #This comment is wrong ignore
KMG2 = 100  #!Km with respect to GA for GSH conjugation mM. Trine: This should be mg. Delete the term MW_aa
# The same comment as for KMG1

# from the code, not from paper
KPT_Li = 0.015    # !'protein turnover rate in liver'# this is confimred in the exp. values of Sweeney
KPT_Ki = 0.013    # !'protein turnover rate in kidney'# 1/h this is 0.012 in Sweeney, I guess its not much of difference
#KPTR   = 1        # !'protein turnover rate in rpt'# this is 0.012 in Li et.al
KPTS   = 0.0039   # !'protein turnover rate in spt'# and this they have it 0.0051 Sweeney and Li 0.0012 refering to Sweeney
KPTRB  = 0.00035  # (correct by Maria)    # !'protein turnover rate in rbc'# how did we get that?
KPTPL  = 0.012    # !'protein turnover rate in plasma'# Sweeney has 0.012 for this value

kpb_Li = 0.1    #!protein binding in liver
reactratio = 2.0  #!reactivity ratio with protein, AA/GA
kpb_Li_2 = kpb_Li/reactratio    #!protein binding in liver

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CORRECTED HB ADDUCT PARAMETERS
# 2 Formation Constants (L/h)
K_FORM_AA_VAL =  6400  # !fmol AA-val/mg globin per mM AA-hr

K_FORM_GA_VAL <- 59000 # fmol GA-val/mg globin per mM GA-hr

# 3. Removal Constants (1/h)
# Adducts are removed with RBC turnover (Lifespan ~120 days)
# Should match KPTRB (Protein turnover in RBC)
K_REM_AA_VAL <- 0.00035 #AA-val adducts from RBC per hr
K_REM_GA_VAL <- 0.00035 #AA-val adducts from RBC per hr



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# create list of parameter
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
params <- unlist(c(data.frame(pAA_TB, pAA_LiB, pAA_KiB,
                              pGA_TB, pGA_LiB, pGA_KiB,
                              k_onAA_T, k_onAA_B, k_onAA_Li,  k_onAA_Ki,
                              k_onGA_T, k_onGA_B, k_onGA_Li, k_onGA_Ki, 
                              k_exc_AAMA, k_exc_GAMA,
                              k_cl_GSH, V_max_p450, KM_p450, V_max_EH, KM_EH,
                              k_AAuptake)  ))
yini <- c(m_AA_AB = 0.0, m_GA_AB = 0.0,      
          m_AA_VB = 0.0, m_GA_VB = 0.0,
          m_AA_mix = 0.0,   # AA amount in mixed venous delay compartment (mg)
          m_GA_mix = 0.0,   # GA amount in mixed venous delay compartment (mg)
          m_AA_Ki = 0.0, m_GA_Ki = 0.0, m_AA_dose = 0 ,
          m_AA_Li = 0.0, m_AAMA = 0.0, 
          m_GA_Li = 0.0, a_pb_GA_Li = 0.0,
          m_GAMA  = 0.0,       
          m_GSH_Li = k_0_GSH * V_Li * MW_GSH, 
          m_AA_T = 0.0, m_GA_T = 0.0,     m_GAOH = 0.0, 
          m_AA_Hb = 0.0,   # <-- AA-hemoglobin adduct (fmol adducts per mg globin)
          m_GA_Hb = 0.0    # <-- GA-hemoglobin adduct (fmol adducts per mg globin)
)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PBPK model for Acrylamide (AA) with GA and mass-balance diagnostics
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PBPKmodelAA <- function(t, state, parameter) {
  with(as.list(c(t, state, parameter)), {
    
    #----#---#  Model for AA  #---#---#
    # concentrations blood in the organs in mg/L
    c_AA_AB <- m_AA_AB / V_AB
    c_AA_VB <- m_AA_VB / V_VB
    c_AA_mix <- m_AA_mix / V_mix   # mixed compartment
    c_AA_T  <- m_AA_T  / V_T
    c_AA_Li <- m_AA_Li / V_Li   
    c_AA_Ki <- m_AA_Ki / V_Ki
    
    c_GSH_Li <- m_GSH_Li / V_Li
    
    # units checked -> mg/h
    dm_AA_dose <- - k_AAuptake * m_AA_dose
    
    # units checked -> mg/h
    # Arterial blood
    dm_AA_AB <- Q_C * (c_AA_VB - c_AA_AB) - k_onAA_B * c_AA_AB * V_AB
    
    # Mixed venous delay compartment: receives tissue outflows
    dm_AA_mix <- Q_T  * (c_AA_T  / pAA_TB) +
      Q_Li * (c_AA_Li / pAA_LiB) +
      Q_Ki * (c_AA_Ki / pAA_KiB) -
      Q_C  * c_AA_mix
    
    # Venous blood: receives flow from mixed compartment
    dm_AA_VB <- Q_C * (c_AA_mix - c_AA_VB) - k_onAA_B * c_AA_VB * V_VB
    
    # !terminal valine hemoglobin adducts
    # we define m_AA_Hb in unit as "fmol adducts per mg globin"
    dm_AA_Hb <- K_FORM_AA_VAL * (c_AA_VB / MW_aa) - K_REM_AA_VAL * m_AA_Hb
    
    
    
    
    # Kidney
    dm_AA_Ki <- Q_Ki * c_AA_AB - Q_Ki * (c_AA_Ki / pAA_KiB) - k_onAA_Ki * m_AA_Ki
   
    # liver metabolism (AA -> GA via P450 and AA -> AAMA via GSH conjugation)
    metAA_GSH  <- Vmax_AA_GSH * c_AA_Li * c_GSH_Li / ((KMG1 + c_AA_Li) * (KMGG + c_GSH_Li))
    metAA_P450 <- V_max_p450 * c_AA_Li / (KM_p450 + c_AA_Li)
    
    # Liver mass balance for AA
    dm_AA_Li <- Q_Li * (c_AA_AB - c_AA_Li / pAA_LiB) + k_AAuptake * m_AA_dose -
      k_onAA_Li * m_AA_Li - metAA_P450 - metAA_GSH

    
    # Tissue
    dm_AA_T <- Q_T * (c_AA_AB - c_AA_T / pAA_TB) - k_onAA_T * m_AA_T
    
    # Urine / metabolite AAMA (from AA via GSH conjugation)
    dm_AAMA <- metAA_GSH - m_AAMA * k_exc_AAMA
    
    
    ########################################################################################
    #-GA-#-GA-#  Model for GA  #-GA-#-GA-#    
    # concentrations blood in the organs in mg/L
    c_GA_AB <- m_GA_AB / V_AB
    c_GA_VB <- m_GA_VB / V_VB
    c_GA_T  <- m_GA_T  / V_T
    c_GA_Li <- m_GA_Li / V_Li
    c_GA_Ki <- m_GA_Ki / V_Ki
    
    
    c_GA_mix <- m_GA_mix / V_mix   # mixed compartment
    
    # Blood
    dm_GA_AB <- Q_C * (c_GA_VB - c_GA_AB) - k_onGA_B * c_GA_AB * V_AB
    
    # Mixed venous delay compartment: receives tissue outflows
    dm_GA_mix <- Q_T * (c_GA_T  / pGA_TB) +
      Q_Li * (c_GA_Li / pGA_LiB) +
      Q_Ki * (c_GA_Ki / pGA_KiB) -
      Q_C  * c_GA_mix
    
    # Venous blood: receives flow from mixed compartment
    dm_GA_VB <- Q_C * (c_GA_mix - c_GA_VB) - k_onGA_B * c_GA_VB * V_VB
    
    # hemoglobin
    dm_GA_Hb <- K_FORM_GA_VAL * ( c_GA_VB / MW_ga ) - K_REM_GA_VAL * m_GA_Hb
    
    
    # Kidney
    dm_GA_Ki <- Q_Ki * c_GA_AB - Q_Ki * (c_GA_Ki / pGA_KiB) - k_onGA_Ki * m_GA_Ki

    # Liver metabolism for GA (and formation from AA via p450)
    metGA_GSH <- Vmax_GA_GSH * c_GSH_Li * c_GA_Li / ((c_GA_Li + KMG2) * (c_GSH_Li + KMGG))
    metGA_EH  <- V_max_EH  * c_GA_Li / (KM_EH + c_GA_Li)
    
    
    metAA_P450_GA <- metAA_P450 * (MW_ga / MW_aa)
    dm_GA_Li <- Q_Li * (c_GA_AB - c_GA_Li / pGA_LiB) - k_onGA_Li * m_GA_Li +
      metAA_P450_GA - metGA_GSH - metGA_EH
    # protein turnover GA in Liver
    da_pb_GA_Li <- kpb_Li_2 * m_GA_Li - a_pb_GA_Li * KPT_Li
    
    # GA metabolites (excretable)
    dm_GAMA  <- metGA_GSH - m_GAMA * k_exc_GAMA
    dm_GAOH  <- metGA_EH  - m_GAOH * k_exc_GAOH
    
    # Tissue
    dm_GA_T <- Q_T * (c_GA_AB - c_GA_T / pGA_TB) - k_onGA_T * m_GA_T

    # GSH pool in liver (consumed by both AA and GA conjugation) # correction of the eqs.
    # earlier eq. 'production − clearance − consumption' 
    #This eq. 'GSH self-corrects to baseline, matches physiology, from Sweeney
    dm_GSH_Li <- k_cl_GSH * (k_0_GSH * V_Li * MW_GSH - m_GSH_Li) - metAA_GSH - metGA_GSH
    
    
    
    
    # -------------------------------------------------------------------------
    # Mass-balance diagnostics (do not change dynamics; returned as named outputs)
    # Total mass pools (AA, GA, and combined AA+GA + metabolites)
    total_AA_mass <- m_AA_AB + m_AA_VB + m_AA_T + m_AA_Li + m_AA_Ki + m_AAMA
    total_GA_mass <- m_GA_AB + m_GA_VB + m_GA_T + m_GA_Li + m_GA_Ki + m_GAMA + m_GAOH
    
    # Rate of change of total AA (sum of AA derivatives)
    dm_total_AA <- dm_AA_dose + dm_AA_AB + dm_AA_VB + dm_AA_T + dm_AA_Li + dm_AA_Ki + dm_AAMA
    # Rate of change of total GA (sum of GA derivatives)
    dm_total_GA <- dm_GA_AB + dm_GA_VB + dm_GA_T + dm_GA_Li + dm_GA_Ki + dm_GAMA + dm_GAOH
    
    # Combined pool derivative (AA + GA + excretable metabolites)
    combined_deriv <- dm_total_AA + dm_total_GA
    
    # Total excretion (rate at which mass leaves the modeled system)
    total_excretion_rate <- m_AAMA * k_exc_AAMA + m_GAMA * k_exc_GAMA + m_GAOH * k_exc_GAOH
    
    # For perfect mass conservation: combined_deriv + total_excretion_rate should be ~ 0
    # (because conversions AA->GA are internal transfers; only excretion removes mass)
    mass_balance_residual <- combined_deriv + total_excretion_rate
    
    # Flow/volume diagnostics
    flow_residual <- Q_C - (Q_T + Q_Li + Q_Ki)   # if non-zero, flows don't sum to cardiac output
    volume_sum    <- V_AB + V_VB + V_T + V_Li + V_Ki  # check if sum matches expected total blood/tissue volume
    
    # -------------------------------------------------------------------------
    # Return: first element must be vector of derivatives (in the same order as `state`)
    derivs <- c(
      dm_AA_AB,    dm_GA_AB,   
      dm_AA_VB,    dm_GA_VB,
      dm_AA_mix,         dm_GA_mix,
      dm_AA_Ki,    dm_GA_Ki,   dm_AA_dose,
      dm_AA_Li,    dm_AAMA,
      dm_GA_Li,    da_pb_GA_Li, dm_GAMA,
      dm_GSH_Li,   dm_AA_T,     dm_GA_T,
      dm_GAOH,
      dm_AA_Hb ,   # <-- AA-hemoglobin adduct
      dm_GA_Hb     # <-- GA-hemoglobin adduct
    )
    
    # Named diagnostics returned as additional items (deSolve will include them in output)
    return(list(    derivs,
      total_AA_mass = total_AA_mass,
      total_GA_mass = total_GA_mass,
      combined_total_mass = total_AA_mass + total_GA_mass,
      dm_total_AA = dm_total_AA,
      dm_total_GA = dm_total_GA,
      combined_deriv = combined_deriv,
      total_excretion_rate = total_excretion_rate,
      mass_balance_residual = mass_balance_residual,
      flow_residual = flow_residual,
      volume_sum = volume_sum
    ))
  })
}
PBPKmodelAA <- compiler::cmpfun(PBPKmodelAA)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# manual readout from Vikstrom 2011 InVitroAcrylamide
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

n_days = 1 # simulate for 5 days
times <- seq(from = 0, to = n_days * 24   , by = 0.1)
diet <- data.frame(var = "m_AA_dose", method = "add",   #  24* c(0:(n_days/2))
                   time = c(0 ),  value = 1000 *BW /1000 ) #  0.5 *BW /1000   # dose of 0.05 microg/kg bw
out <- ode(y = yini, times = times, func = PBPKmodelAA, parms = params, events = list(data = diet))


par( mfrow = c(1,2) , mar = c(4,4,1,1) )
# plot for hemoglobin adducts AA 
plot(out[,'time'], 
     out[,'m_AA_Hb'],
     type = 'l', lwd=2,  xlab = 'Time (hours)',    ylab = 'Hb Adducts (pmol/g globin)', 
     main = 'AA Hemoglobin Adducts',   
     ylim = c(0, 3000)); grid() # Adjusted ylim to typical range (0-150)
# Add reference points (ensure these values are also in pmol/g globin)
#points(c(24, 120), c(51, 110), col="blue", cex = 1.5, pch = 17) ; grid()

plot(out[,'time'],  out[,'m_GA_Hb'], 
     type = 'l', lwd=2,  xlab = 'Time (hours)',    ylab = 'Hb Adducts (pmol/g globin)', 
     main = 'GA Hemoglobin Adducts',   
     ylim = c(0, 300)); grid() # Adjusted ylim to typical range (0-150)
# Add reference points (ensure these values are also in pmol/g globin)
#points(c(24, 120), c(51, 100), col="blue", cex = 1.5, pch = 17) ; grid()





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#     Sensitivity analysis using the pksensi package in R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 1. First, we must ensure all tested parameters are actually in the numeric 'params' vector.
# Let's rebuild the numeric params vector to include Q_C and k_exc_GAOH
params <- unlist(c(data.frame(Q_C, pAA_TB, pAA_LiB, pAA_KiB,
                              pGA_TB, pGA_LiB, pGA_KiB,
                              k_onAA_T, k_onAA_B, k_onAA_Li,  k_onAA_Ki,
                              k_onGA_T, k_onGA_B, k_onGA_Li, k_onGA_Ki, 
                              k_exc_AAMA, k_exc_GAMA, k_exc_GAOH,
                              k_cl_GSH, V_max_p450, KM_p450, V_max_EH, KM_EH,
                              k_AAuptake)))

## Define the distribution of the parameters that you will analyse in the sensitivity test 
LL <- 0.9 # 10% lower limit
UL <- 1.1 # 10% upper limit

# 2. Define a character vector of ONLY the parameters that actually exist in the model
sens_params <- c('BW' ,'Q_C', 'pAA_TB', 'pAA_LiB',  'pAA_KiB'
                 # ,
                 # 'pGA_TB', 'pGA_LiB',  'pGA_KiB',
                 # 'k_onAA_T', 'k_onAA_B', 'k_onAA_Li', 'k_onAA_Ki',
                 # 'k_onGA_T', 'k_onGA_B', 'k_onGA_Li', 'k_onGA_Ki', 
                 # 'k_exc_AAMA', 'k_exc_GAMA', 'k_exc_GAOH',
                 # 'k_cl_GSH', 'V_max_p450', 'KM_p450', 'V_max_EH'
                 )

# 3. Dynamically build the q.arg list so you don't have to type them all out manually
q.arg <- lapply(sens_params, function(p) {
  list(min = params[p] * LL, max = params[p] * UL)
})

q <- rep("qunif", length(sens_params))

## Create parameter matrix ##  
set.seed(124)

# Use the sens_params character vector here
x <- rfast99(params = sens_params, n = 200, q = q, q.arg = q.arg, rep = 1)

library(pksensi)

## Conduct simulation pksensi ##
# 4. The output must perfectly match a state variable name (m_GAMA, not m_GAMA_urinary)
outputs <- c("m_AA_Hb")

out_sensitivity <- solve_fun(x, time = times, func = PBPKmodelAA, 
                             events = list(data = diet), initState = yini, outnames = outputs)

pksim(out_sensitivity)

## Output from the sensitivity analysis ##
plot(out_sensitivity)






