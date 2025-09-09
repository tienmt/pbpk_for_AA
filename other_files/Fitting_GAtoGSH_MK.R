### Simple one compartment model trying to fit the Vmax and Km
### of GA to GSH conjugation using urinary data from Kopp paper###

###The values are far from what Sweeney report as fitted values

###Limitations are that in the Kopp they do not have GA plasma data
###and that we use one compartment model not the whole PBPK as they 
##claim that they did in Sweeney

library(deSolve)
library(minpack.lm)

#  Model definition 
mm_model_urine <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    rate_GA_GSH <- Vmax * GA / (Km + GA)
    rate_excr <- k_exc * GAMA
    
    dGA <- -rate_GA_GSH
    dGAMA <- rate_GA_GSH - rate_excr
    dGAMA_urine <- rate_excr
    
    list(c(dGA, dGAMA, dGAMA_urine))
  })
}

# Data from Kopp et al. urinary GAMA data (mg) 
time_gama <- c(0, 3.9, 8.3, 14, 19.5, 28, 37, 45.9)
gama_umol <- c(0.0, 2.15e-6, 3.66e-6, 4.38e-6, 6.57e-6, 5.26e-6, 3.50e-6, 3.0e-6)
MW_GAMA <- 250.27
gama_mg_obs <- gama_umol * MW_GAMA

#Initial conditions 
GA0 <- 0.418  # mg (30% of a 1.4 mg acrylamide dose)
state0 <- c(GA = GA0, GAMA = 0, GAMA_urine = 0)

# --- Fixed parameter ---
k_exc <- 0.693 / 26.3  # h⁻¹, based on GAMA half-life

# Simulation function
simulate_model <- function(par) {
  names(par) <- c("Vmax", "Km")
  times <- sort(unique(c(time_gama, seq(0, 50, by = 0.1))))
  out <- ode(y = state0, times = times, func = mm_model_urine, parms = c(par, k_exc = k_exc))
  pred <- approx(out[, "time"], out[, "GAMA_urine"], xout = time_gama)$y
  return(pred)
}

#  Residuals for fit
residuals_fn <- function(par) {
  pred <- simulate_model(par)
  return(pred - gama_mg_obs)
}

# Starting values and fit
start_par <- c(Vmax = 0.01, Km = 10)
fit <- nls.lm(par = start_par, fn = residuals_fn, lower = c(0.0001, 0.1), upper = c(1, 500))
fitted_pars <- fit$par
print(fitted_pars)

# Caclculating allometric params
Vmax_fit <- fitted_pars["Vmax"] # mg/h
Km_fit <- fitted_pars["Km"] # mg/L
BW <- 70
Vmax_allometric <- Vmax_fit / (BW^0.7)


