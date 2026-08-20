# Recalculate physiology from BW inside the model
# Keep original parameter architecture: V_max_p450, Vmax_AA_GSH, Vmax_GA_GSH, V_max_EH
# Add BW to params and sens_params
# Recalculate actual Vmax values from sampled baseline Vmax values using (BW/BW_ref)^0.7

library(deSolve)
library(minpack.lm)
library(pksensi)
library(dplyr)
library(tidyr)
library(ggplot2)

graphics.off()


# Baseline constants
BW_ref <- 70
BW <- 70

# Fractions of BW
F_B <- 0.079
F_B_AB <- 0.35
F_B_VB <- 0.65
F_Li <- 0.026
F_Ki <- 0.0044
F_T <- 1 - (F_B + F_Li + F_Ki)

# Cardiac output and flow fractions
QCC <- 16
FQ_Li <- 0.255
FQ_Ki <- 0.19
FQ_T <- 1 - (FQ_Li + FQ_Ki)

# Molecular weights
MW_aa <- 71
MW_ga <- 87
MW_GSH <- 307.32

# Baseline GSH
k_0_GSH <- 7

# Partition coefficients
pAA_TB <- 0.55
pAA_KiB <- 0.8
pAA_LiB <- 0.4

pGA_TB <- 1.35
pGA_LiB <- 0.5
pGA_KiB <- 0.3

# Reaction and excretion constants
k_AAuptake <- 1

k_onAA_T <- 0.028
k_onAA_B <- 0.003
k_onAA_Li <- 0.71
k_onAA_Ki <- 0.13

k_cl_GSH <- 0.035

k_exc_AAMA <- 0.049
k_exc_GAMA <- 0.027
k_exc_GAOH <- 0.077

# Original actual Vmax parameters at BW_ref = 70 kg
V_max_p450 <- 0.5 * BW_ref^0.7
Vmax_AA_GSH <- 22 * BW_ref^0.7
Vmax_GA_GSH <- 20 * BW_ref^0.7
V_max_EH <- 20 * BW_ref^0.7

KM_p450 <- 10.0
KM_EH <- 100.0
KMG1 <- 100
KMGG <- 0.1 * MW_GSH
KMG2 <- 100

# GA parameters
k_onGA_T <- 0.69
k_onGA_B <- 0.108
k_onGA_Ki <- k_onAA_T / 2
k_onGA_Li <- 0.115

# Protein turnover and binding
KPT_Li <- 0.015
KPT_Ki <- 0.013
KPTS <- 0.0039
KPTRB <- 0.00035
KPTPL <- 0.012

kpb_Li <- 0.1
kpbk_ki <- 0.08
kpbr1 <- 0.08
kpbs1 <- 0.08
kpbrb1 <- 0.01
kpbpl1 <- 0.01
reactratio <- 2.0

kpb_Li_2 <- kpb_Li / reactratio
kpbk2 <- kpbk_ki / reactratio
kpbr2 <- kpbr1 / reactratio
kpbs2 <- kpbs1 / reactratio
kpbrb2 <- kpbrb1 / reactratio
kpbpl2 <- kpbpl1 / reactratio

# Hb adduct parameters
K_FORM_AA_VAL <- 6400
K_FORM_GA_VAL <- 59000
K_REM_AA_VAL <- 0.00035
K_REM_GA_VAL <- 0.00035


# Parameters vector
params <- unlist(c(data.frame(
  BW,
  pAA_TB, pAA_LiB, pAA_KiB,
  pGA_TB, pGA_LiB, pGA_KiB,
  k_onAA_T, k_onAA_B, k_onAA_Li, k_onAA_Ki,
  k_onGA_T, k_onGA_B, k_onGA_Li, k_onGA_Ki,
  k_exc_AAMA, k_exc_GAMA, k_exc_GAOH,
  k_cl_GSH,
  V_max_p450, KM_p450, Vmax_AA_GSH,
  V_max_EH, KM_EH, Vmax_GA_GSH,
  KMG1, KMG2, KMGG,
  k_AAuptake,
  K_FORM_AA_VAL, K_FORM_GA_VAL,
  K_REM_AA_VAL, K_REM_GA_VAL
)))


# Initial conditions
# Initial GSH uses baseline BW. During simulation, the target GSH level is recalculated
# from the sampled BW inside PBPKmodelAA.

V_Li_baseline <- BW_ref * F_Li

yini <- c(
  m_AA_AB = 0.0, m_GA_AB = 0.0,
  m_AA_VB = 0.0, m_GA_VB = 0.0,
  m_AA_mix = 0.0,
  m_GA_mix = 0.0,
  m_AA_Ki = 0.0, m_GA_Ki = 0.0, m_AA_dose = 0.0,
  m_AA_Li = 0.0, m_AAMA = 0.0,
  m_GA_Li = 0.0, a_pb_GA_Li = 0.0,
  m_GAMA = 0.0,
  m_GSH_Li = as.numeric(k_0_GSH * V_Li_baseline * MW_GSH),
  m_AA_T = 0.0, m_GA_T = 0.0,
  m_GAOH = 0.0,
  a_pb_AA_Ki = 0.0, a_pb_AA_Li = 0.0,
  a_pb_GA_Ki = 0.0, a_pb_AA_T = 0.0,
  a_pb_GA_T = 0.0,
  m_AA_Hb = 0.0,
  m_GA_Hb = 0.0
)

# PBPK model
PBPKmodelAA <- function(t, state, parameter) {
  with(as.list(c(t, state, parameter)), {
    # BW-scaled physiology
    V_AB <- BW * F_B * F_B_AB
    V_VB <- BW * F_B * F_B_VB
    V_T <- BW * F_T
    V_Ki <- BW * F_Ki
    V_Li <- BW * F_Li
    V_mix <- V_VB

    Q_C <- QCC * BW^0.75
    Q_Li <- Q_C * FQ_Li
    Q_Ki <- Q_C * FQ_Ki
    Q_T <- Q_C * FQ_T
    # BW-scaling of actual Vmax parameters
    # For a sampled BW, they are scaled at BW_ref = 70 kg as:
    # Vmax(BW) = Vmax(70 kg) * (BW / 70)^0.7
    V_max_p450_BW <- V_max_p450 * (BW / BW_ref)^0.7
    Vmax_AA_GSH_BW <- Vmax_AA_GSH * (BW / BW_ref)^0.7
    Vmax_GA_GSH_BW <- Vmax_GA_GSH * (BW / BW_ref)^0.7
    V_max_EH_BW <- V_max_EH * (BW / BW_ref)^0.7

    #----#---#  Model for AA  #---#---#
    c_AA_AB <- m_AA_AB / V_AB
    c_AA_VB <- m_AA_VB / V_VB
    c_AA_mix <- m_AA_mix / V_mix
    c_AA_T <- m_AA_T / V_T
    c_AA_Li <- m_AA_Li / V_Li
    c_AA_Li_free <- c_AA_Li / pAA_LiB
    c_AA_Ki <- m_AA_Ki / V_Ki
    
    m_GSH_Li_BW <- m_GSH_Li * (BW / BW_ref)
    c_GSH_Li <- m_GSH_Li_BW / V_Li

   
    dm_AA_dose <- -k_AAuptake * m_AA_dose

    dm_AA_AB <- Q_C * (c_AA_VB - c_AA_AB) - k_onAA_B * c_AA_AB * V_AB

    dm_AA_mix <- Q_T * (c_AA_T / pAA_TB) +
      Q_Li * (c_AA_Li / pAA_LiB) +
      Q_Ki * (c_AA_Ki / pAA_KiB) -
      Q_C * c_AA_mix

    dm_AA_VB <- Q_C * (c_AA_mix - c_AA_VB) - k_onAA_B * c_AA_VB * V_VB

    dm_AA_Hb <- K_FORM_AA_VAL * (c_AA_VB / MW_aa) - K_REM_AA_VAL * m_AA_Hb

    dm_AA_Ki <- Q_Ki * c_AA_AB - Q_Ki * (c_AA_Ki / pAA_KiB) - k_onAA_Ki * m_AA_Ki

    da_pb_AA_Ki <- kpbk_ki * m_AA_Ki - a_pb_AA_Ki * KPT_Ki

    metAA_GSH <- Vmax_AA_GSH_BW * c_AA_Li_free * c_GSH_Li /
      ((KMG1 + c_AA_Li_free) * (KMGG + c_GSH_Li))

    metAA_P450 <- V_max_p450_BW * c_AA_Li_free /
      (KM_p450 + c_AA_Li_free)

    dm_AA_Li <- Q_Li * (c_AA_AB - c_AA_Li / pAA_LiB) +k_AAuptake * m_AA_dose -
      k_onAA_Li * m_AA_Li -metAA_P450 -metAA_GSH

    da_pb_AA_Li <- kpb_Li * m_AA_Li - a_pb_AA_Li * KPT_Li

    dm_AA_T <- Q_T * (c_AA_AB - c_AA_T / pAA_TB) - k_onAA_T * m_AA_T

    da_pb_AA_T <- kpbs1 * m_AA_T - a_pb_AA_T * KPTS

    dm_AAMA <- metAA_GSH - m_AAMA * k_exc_AAMA

    # GA model
    c_GA_AB <- m_GA_AB / V_AB
    c_GA_VB <- m_GA_VB / V_VB
    c_GA_T <- m_GA_T / V_T
    c_GA_Li <- m_GA_Li / V_Li
    c_GA_Li_free <- c_GA_Li / pGA_LiB
    c_GA_Ki <- m_GA_Ki / V_Ki
    c_GA_mix <- m_GA_mix / V_mix

    dm_GA_AB <- Q_C * (c_GA_VB - c_GA_AB) - k_onGA_B * c_GA_AB * V_AB

    dm_GA_mix <- Q_T * (c_GA_T / pGA_TB) +Q_Li * (c_GA_Li / pGA_LiB) +
      Q_Ki * (c_GA_Ki / pGA_KiB) - Q_C * c_GA_mix

    dm_GA_VB <- Q_C * (c_GA_mix - c_GA_VB) - k_onGA_B * c_GA_VB * V_VB

    dm_GA_Hb <- K_FORM_GA_VAL * (c_GA_VB / MW_ga) - K_REM_GA_VAL * m_GA_Hb

    dm_GA_Ki <- Q_Ki * c_GA_AB - Q_Ki * (c_GA_Ki / pGA_KiB) - k_onGA_Ki * m_GA_Ki

    da_pb_GA_Ki <- kpbk_ki * m_GA_Ki - a_pb_GA_Ki * KPT_Ki

    metGA_GSH <- Vmax_GA_GSH_BW * c_GSH_Li * c_GA_Li_free /
      ((KMG2 + c_GA_Li_free) * (KMGG + c_GSH_Li))

    metGA_EH <- V_max_EH_BW * c_GA_Li_free / (KM_EH + c_GA_Li_free)

    metAA_P450_GA <- metAA_P450 * (MW_ga / MW_aa)

    dm_GA_Li <- Q_Li * (c_GA_AB - c_GA_Li / pGA_LiB) -k_onGA_Li * m_GA_Li +
      metAA_P450_GA -metGA_GSH -metGA_EH

    da_pb_GA_Li <- kpb_Li_2 * m_GA_Li - a_pb_GA_Li * KPT_Li

    dm_GAMA <- metGA_GSH - m_GAMA * k_exc_GAMA

    dm_GAOH <- metGA_EH - m_GAOH * k_exc_GAOH

    dm_GA_T <- Q_T * (c_GA_AB - c_GA_T / pGA_TB) - k_onGA_T * m_GA_T

    da_pb_GA_T <- kpbs2 * m_GA_T - a_pb_GA_T * KPTS

    dm_GSH_Li <- (k_cl_GSH * (k_0_GSH * V_Li * MW_GSH - m_GSH_Li_BW) -
        metAA_GSH - metGA_GSH) * (BW_ref / BW)
    
    
    # dm_GSH_Li <- k_cl_GSH * (k_0_GSH * V_Li * MW_GSH - m_GSH_Li) -
    #  metAA_GSH -  metGA_GSH

    # Diagnostics
    total_AA_mass <- m_AA_AB + m_AA_VB + m_AA_mix + m_AA_T + m_AA_Li + m_AA_Ki + m_AAMA + m_AA_dose
    total_GA_mass <- m_GA_AB + m_GA_VB + m_GA_mix + m_GA_T + m_GA_Li + m_GA_Ki + m_GAMA + m_GAOH

    dm_total_AA <- dm_AA_dose + dm_AA_AB + dm_AA_VB + dm_AA_mix + dm_AA_T + dm_AA_Li + dm_AA_Ki + dm_AAMA
    dm_total_GA <- dm_GA_AB + dm_GA_VB + dm_GA_mix + dm_GA_T + dm_GA_Li + dm_GA_Ki + dm_GAMA + dm_GAOH

    combined_deriv <- dm_total_AA + dm_total_GA

    total_excretion_rate <- m_AAMA * k_exc_AAMA + m_GAMA * k_exc_GAMA + m_GAOH * k_exc_GAOH

    mass_balance_residual <- combined_deriv + total_excretion_rate

    flow_residual <- Q_C - (Q_T + Q_Li + Q_Ki)
    volume_sum <- V_AB + V_VB + V_mix + V_T + V_Li + V_Ki

    derivs <- c(
      dm_AA_AB, dm_GA_AB,
      dm_AA_VB, dm_GA_VB,
      dm_AA_mix, dm_GA_mix,
      dm_AA_Ki, dm_GA_Ki, dm_AA_dose,
      dm_AA_Li, dm_AAMA,
      dm_GA_Li, da_pb_GA_Li, dm_GAMA,
      dm_GSH_Li, dm_AA_T, dm_GA_T,
      dm_GAOH,
      da_pb_AA_Ki, da_pb_AA_Li,
      da_pb_GA_Ki, da_pb_AA_T,
      da_pb_GA_T,
      dm_AA_Hb,
      dm_GA_Hb
    )

    return(list(
      derivs,
      total_AA_mass = total_AA_mass,
      total_GA_mass = total_GA_mass,
      combined_total_mass = total_AA_mass + total_GA_mass,
      dm_total_AA = dm_total_AA,
      dm_total_GA = dm_total_GA,
      combined_deriv = combined_deriv,
      total_excretion_rate = total_excretion_rate,
      mass_balance_residual = mass_balance_residual,
      flow_residual = flow_residual,
      volume_sum = volume_sum,
      V_max_p450_BW = V_max_p450_BW,
      Vmax_AA_GSH_BW = Vmax_AA_GSH_BW,
      Vmax_GA_GSH_BW = Vmax_GA_GSH_BW,
      V_max_EH_BW = V_max_EH_BW,
      metAA_GSH = metAA_GSH,
      metAA_P450 = metAA_P450,
      metGA_GSH = metGA_GSH,
      metGA_EH = metGA_EH))})
}

PBPKmodelAA <- compiler::cmpfun(PBPKmodelAA)



# Dose is defined as mg/kg and scaled by the sampled BW
# within the deSolve dosing event.
dose_mgkg <- 0.5 / 1000  # 0.5 µg/kg = 0.0005 mg/kg
dose_times <- c(0, 24, 48)


# Dosing scenario
n_days <- 5
times <- seq(from = 0, to = n_days * 24, by = 1)

# Function for instantaneous Oral-dose event.
# For every FAST run, parms["BW"] is that run's sampled BW.
dose_event <- function(t, y, parms) {
  BW_i <- as.numeric(parms["BW"])
  y["m_AA_dose"] <- y["m_AA_dose"] + dose_mgkg * BW_i
  return(y)}

#  Sensitivity analysis
run_model <- function(parms) {
  out <- ode(
    y = yini,
    times = times,
    func = PBPKmodelAA,
    parms = parms,
    events = list(
      func = dose_event,
      time = dose_times))
  return(as.data.frame(out))
}

LL <- 0.9
UL <- 1.1

sens_params <- c(
  'BW',
  'pAA_TB', 'pAA_LiB', 'pAA_KiB',
  'k_onAA_B', 'k_onAA_T', 'k_onAA_Li',
  'k_onAA_Ki', 'k_exc_AAMA', 'V_max_EH', 'KM_EH',
  'K_FORM_AA_VAL', 'K_REM_AA_VAL',
  'k_cl_GSH', 'V_max_p450', 'KM_p450',
  'pGA_TB', 'pGA_LiB', 'pGA_KiB',
  'k_onGA_T', 'k_onGA_B', 'k_onGA_Li', 'k_onGA_Ki',
  'Vmax_AA_GSH', 'Vmax_GA_GSH', 'k_exc_GAMA', 'k_exc_GAOH',
  'K_FORM_GA_VAL', 'K_REM_GA_VAL'
)

sens_params <- unique(sens_params)

print(table(sens_params))
print(setdiff(sens_params, names(params)))
print(setdiff(names(params), sens_params))

if (length(setdiff(sens_params, names(params))) > 0) {
  stop("Some sens_params are missing from params.")
}

q.arg <- lapply(sens_params, function(p) {
  list(
    min = params[p] * LL,
    max = params[p] * UL
  )
})

q <- rep("qunif", length(sens_params))

set.seed(124)

x <- rfast99(params = sens_params,n = 200,q = q,q.arg = q.arg,rep = 1)

if (any(!is.finite(x$a))) {
  print(head(x$a))
  stop("rfast99 produced non-finite values.")
}

par(mfrow = c(1, 1), mar = c(4, 4, 1, 1))

outputs <- c('m_AA_Hb')  # m_AA_Hb ; m_AAMA ; m_GAMA ; m_GA_Hb

# The event definition is fixed, but dose_event() reads the sampled BW
# from parms separately for every FAST model run.
out_sensitivity <- solve_fun(
  x,
  time = times,
  func = PBPKmodelAA,
  initState = yini,
  outnames = outputs,
  events = list(
    func = dose_event,
    time = dose_times))

tSI <- out_sensitivity$tSI[, , 1]
tSI_df <- as.data.frame(tSI)
tSI_df$time <- as.numeric(rownames(tSI_df))

tSI_long <- tSI_df %>%
  pivot_longer(cols = -time, names_to = "Parameter",values_to = "Sensitivity")

ggplot(
  tSI_long,
  aes(x = time, y = Parameter, fill = Sensitivity)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "black",
    midpoint = median(tSI_long$Sensitivity, na.rm = TRUE)) +
  theme_bw(base_size = 14) +
  labs(
    title = paste("Total Sensitivity:", outputs[1]),
    x = "Time",
    y = "Parameter")

