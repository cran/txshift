## ----setup--------------------------------------------------------------------
library(data.table)
library(haldensify)
library(txshift)
set.seed(11249)

## ----make_data----------------------------------------------------------------
# parameters for simulation example
n_obs <- 200  # number of observations

# baseline covariate -- simple, binary
W <- rbinom(n_obs, 1, 0.5)

# create treatment based on baseline W
A <- rnorm(n_obs, mean = 2 * W, sd = 0.5)

# create outcome as a linear function of A, W + white noise
Y <- A + W + rnorm(n_obs, mean = 0, sd = 1)

# two-phase sampling based on covariates
Delta_samp <- rbinom(n_obs, 1, plogis(W))

# treatment shift parameter
delta <- 0.5

## ----est_simple---------------------------------------------------------------
est_shift <- txshift(
  W = W, A = A, Y = Y, delta = delta,
  g_exp_fit_args = list(
    fit_type = "hal", n_bins = 5,
    grid_type = "equal_mass",
    lambda_seq = exp(seq(-1, -10, length = 100))
  ),
  Q_fit_args = list(
    fit_type = "glm",
    glm_formula = "Y ~ .^2"
  )
)
est_shift

## ----setup_sl, eval = FALSE---------------------------------------------------
#  library(sl3)
#  
#  # SL learners to be used for most fits (e.g., IPCW, outcome regression)
#  mean_learner <- Lrnr_mean$new()
#  glm_learner <- Lrnr_glm$new()
#  rf_learner <- Lrnr_ranger$new()
#  Q_lib <- Stack$new(mean_learner, glm_learner, rf_learner)
#  sl_learner <- Lrnr_sl$new(learners = Q_lib, metalearner = Lrnr_nnls$new())
#  
#  # SL learners for fitting the generalized propensity score fit
#  hose_learner <- make_learner(Lrnr_density_semiparametric,
#    mean_learner = glm_learner
#  )
#  hese_learner <- make_learner(Lrnr_density_semiparametric,
#    mean_learner = rf_learner,
#    var_learner = glm_learner
#  )
#  g_lib <- Stack$new(hose_learner, hese_learner)
#  sl_learner_density <- Lrnr_sl$new(
#    learners = g_lib,
#    metalearner = Lrnr_solnp_density$new()
#  )

## ----est_with_sl, eval = FALSE------------------------------------------------
#  est_shift_sl <- txshift(
#    W = W, A = A, Y = Y, delta = delta,
#    g_exp_fit_args = list(
#      fit_type = "sl",
#      sl_learners_density = sl_learner_density
#    ),
#    Q_fit_args = list(
#      fit_type = "sl",
#      sl_learners = sl_learner
#    )
#  )
#  est_shift_sl

## ----ipcw_est_shift-----------------------------------------------------------
est_shift_ipcw <- txshift(
  W = W, A = A, Y = Y, delta = delta,
  C_samp = Delta_samp, V = c("W", "Y"),
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(
    fit_type = "hal", n_bins = 5,
    grid_type = "equal_mass",
    lambda_seq = exp(seq(-1, -10, length = 100))
  ),
  Q_fit_args = list(
    fit_type = "glm",
    glm_formula = "Y ~ .^2"
  ),
  eif_reg_type = "glm"
)
est_shift_ipcw

## ----confint------------------------------------------------------------------
(ci_est_shift <- confint(est_shift))

## ----fit_external, eval=FALSE-------------------------------------------------
#  # compute censoring mechanism and produce IPC weights externally
#  pi_mech <- plogis(W)
#  ipcw_out <- pi_mech
#  
#  # compute treatment mechanism (propensity score) externally
#  ## first, produce the down-shifted treatment data
#  gn_downshift <- dnorm(A - delta, mean = tx_mult * W, sd = 1)
#  ## next, initialize and produce the up-shifted treatment data
#  gn_upshift <- dnorm(A + delta, mean = tx_mult * W, sd = 1)
#  ## now, initialize and produce the up-up-shifted (2 * delta) treatment data
#  gn_upupshift <- dnorm(A + 2 * delta, mean = tx_mult * W, sd = 1)
#  ## then, initialize and produce the un-shifted treatment data
#  gn_noshift <- dnorm(A, mean = tx_mult * W, sd = 1)
#  ## finally, put it all together into an object like what's produced internally
#  gn_out <- as.data.table(cbind(gn_downshift, gn_noshift, gn_upshift,
#                                gn_upupshift))[C == 1, ]
#  setnames(gn_out, c("downshift", "noshift", "upshift", "upupshift"))
#  
#  # compute outcome regression externally
#  # NOTE: transform Y to lie in the unit interval and bound predictions such that
#  #       no values fall near the bounds of the interval
#  Qn_noshift <- (W + A - min(Y)) / diff(range(Y))
#  Qn_upshift <- (W + A + delta - min(Y)) / diff(range(Y))
#  Qn_noshift[Qn_noshift < 0] <- 0.025
#  Qn_noshift[Qn_noshift > 1] <- 0.975
#  Qn_upshift[Qn_upshift < 0] <- 0.025
#  Qn_upshift[Qn_upshift > 1] <- 0.975
#  Qn_out <- as.data.table(cbind(Qn_noshift, Qn_upshift))[C == 1, ]
#  setnames(Qn_out, c("noshift", "upshift"))
#  
#  # construct efficient estimator by applying wrapper function
#  est_shift_spec <- txshift(
#    W = W, A = A, Y = Y, delta = delta,
#   samp_fit_args = NULL,
#   samp_fit_ext = ipcw_out,
#   g_exp_fit_args = list(fit_type = "external"),
#   Q_fit_args = list(fit_type = "external"),
#   gn_exp_fit_ext = gn_out,
#   Qn_fit_ext = Qn_out
#  )

