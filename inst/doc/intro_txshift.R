## ----setup--------------------------------------------------------------------
library(data.table)
library(haldensify)
library(txshift)
set.seed(429153)

## -----------------------------------------------------------------------------
# simulate simple data for tmle-shift sketch
n_obs <- 1000  # number of observations
n_w <- 1  # number of baseline covariates
tx_mult <- 2  # multiplier for the effect of W = 1 on the treatment

## baseline covariate -- simple, binary
W <- as.numeric(replicate(n_w, rbinom(n_obs, 1, 0.5)))

## create treatment based on baseline W
A <- as.numeric(rnorm(n_obs, mean = tx_mult * W, sd = 1))

# create outcome as a linear function of A, W + white noise
Y <- A + W + rnorm(n_obs, mean = 0, sd = 1)

# shift parameter
delta <- 0.5

## -----------------------------------------------------------------------------
tmle_hal_shift_1 <- txshift(
  W = W, A = A, Y = Y, delta = delta,
  fluctuation = "standard",
  g_exp_fit_args = list(fit_type = "hal", n_bins = 5,
                        grid_type = "equal_mass",
                        lambda_seq = exp(seq(-1, -12, length = 500))),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .")
)
tmle_hal_shift_1

## -----------------------------------------------------------------------------
tmle_hal_shift_2 <- txshift(
  W = W, A = A, Y = Y, delta = delta,
  fluctuation = "weighted",
  g_exp_fit_args = list(fit_type = "hal", n_bins = 5,
                        grid_type = "equal_mass",
                        lambda_seq = exp(seq(-1, -12, length = 500))),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .")
)
tmle_hal_shift_2

## ---- eval = FALSE------------------------------------------------------------
#  # SL learners to be used for most fits (e.g., IPCW, outcome regression)
#  mean_learner <- Lrnr_mean$new()
#  glm_learner <- Lrnr_glm$new()
#  rf_learner <- Lrnr_ranger$new()
#  Q_lib <- Stack$new(mean_learner, glm_learner, rf_learner)
#  sl_learner <- Lrnr_sl$new(learners = Q_lib, metalearner = Lrnr_nnls$new())
#  
#  # SL learners for fitting the generalized propensity score fit
#  hse_learner <- make_learner(Lrnr_density_semiparametric,
#    mean_learner = glm_learner
#  )
#  mvd_learner <- make_learner(Lrnr_density_semiparametric,
#    mean_learner = rf_learner,
#    var_learner = glm_learner
#  )
#  g_lib <- Stack$new(hse_learner, mvd_learner)
#  sl_learner_density <- Lrnr_sl$new(learners = g_lib,
#                                    metalearner = Lrnr_solnp_density$new())

## ---- eval = FALSE------------------------------------------------------------
#  tmle_sl_shift_1 <- txshift(
#    W = W, A = A, Y = Y, delta = delta,
#    fluctuation = "standard",
#    g_exp_fit_args = list(fit_type = "sl",
#                          sl_learners_density = sl_learner_density),
#    Q_fit_args = list(fit_type = "sl", sl_learners = sl_learner)
#  )
#  tmle_sl_shift_1

## ---- eval = FALSE------------------------------------------------------------
#  tmle_sl_shift_2 <- txshift(
#    W = W, A = A, Y = Y, delta = delta,
#    fluctuation = "weighted",
#    g_exp_fit_args = list(fit_type = "sl",
#                          sl_learners_density = sl_learner_density),
#    Q_fit_args = list(fit_type = "sl", sl_learners = sl_learner)
#  )
#  tmle_sl_shift_2

## -----------------------------------------------------------------------------
(ci_param <- confint(tmle_hal_shift_1))

## -----------------------------------------------------------------------------
# compute treatment mechanism (propensity score) externally
## first, produce the down-shifted treatment data
gn_downshift <- dnorm(A - delta, mean = tx_mult * W, sd = 1)
## next, initialize and produce the up-shifted treatment data
gn_upshift <- dnorm(A + delta, mean = tx_mult * W, sd = 1)
## now, initialize and produce the up-up-shifted (2 * delta) treatment data
gn_upupshift <- dnorm(A + 2 * delta, mean = tx_mult * W, sd = 1)
## then, initialize and produce the un-shifted treatment data
gn_noshift <- dnorm(A, mean = tx_mult * W, sd = 1)
## finally, put it all together into an object like what's produced internally
gn_out <- as.data.table(cbind(gn_downshift, gn_noshift, gn_upshift,
                              gn_upupshift))
setnames(gn_out, c("downshift", "noshift", "upshift", "upupshift"))

# compute outcome regression externally
Qn_noshift <- (W + A - min(Y)) / diff(range(Y))
Qn_upshift <- (W + A + delta - min(Y)) / diff(range(Y))
Qn_noshift[Qn_noshift < 0] <- .Machine$double.neg.eps
Qn_noshift[Qn_noshift > 1] <- 1 - .Machine$double.neg.eps
Qn_upshift[Qn_upshift < 0] <- .Machine$double.neg.eps
Qn_upshift[Qn_upshift > 1] <- 1 - .Machine$double.neg.eps
Qn_out <- as.data.table(cbind(Qn_noshift, Qn_upshift))
setnames(Qn_out, c("noshift", "upshift"))

# invoke the wrapper function only to apply the targeting step
tmle_shift_spec <- txshift(W = W, A = A, Y = Y, delta = delta,
                           fluctuation = "standard",
                           samp_fit_args = NULL,
                           g_exp_fit_args = list(fit_type = "external"),
                           Q_fit_args = list(fit_type = "external"),
                           gn_exp_fit_ext = gn_out,
                           Qn_fit_ext = Qn_out)
tmle_shift_spec

