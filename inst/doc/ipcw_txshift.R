## ----setup--------------------------------------------------------------------
library(data.table)
library(haldensify)
library(txshift)
set.seed(429153)

## -----------------------------------------------------------------------------
# simulate simple data for tmle-shift sketch
n_obs <- 500  # sample size
n_w <- 1  # just one baseline covariate for this example
tx_mult <- 2  # multiplier for the effect of W = 1 on the treatment

# baseline covariate -- simple, binary
W <- as.numeric(replicate(n_w, rbinom(n_obs, 1, 0.5)))

# set and organize treatment based on baseline W
A <- as.numeric(rnorm(n_obs, mean = tx_mult * W, sd = 1))

# create outcome as linear function of A, W + white noise
Y <- A + W + rnorm(n_obs, mean = 0, sd = 1)

# censoring based on covariates
C <- rbinom(n_obs, 1, plogis(W))

# treatment shift parameter
delta <- 0.5

## ----ipcw_tmle_ghal_Qglm_1----------------------------------------------------
tmle_hal_shift_1 <- txshift(W = W, A = A, Y = Y,
                            C = C, V = c("W", "Y"),
                            delta = delta,
                            fluctuation = "standard",
                            max_iter = 2,
                            ipcw_fit_args = list(fit_type = "glm"),
                            g_fit_args = list(fit_type = "hal",
                                              n_bins = 5,
                                              grid_type = "equal_mass",
                                              lambda_seq =
                                                exp(seq(-1, -9,
                                                        length = 300))),
                            Q_fit_args = list(fit_type = "glm",
                                              glm_formula = "Y ~ ."),
                            eif_reg_type = "glm"
                           )
summary(tmle_hal_shift_1)

## ----ipcw_tmle_ghal_Qglm_2----------------------------------------------------
tmle_hal_shift_2 <- txshift(W = W, A = A, Y = Y,
                            C = C, V = c("W", "Y"),
                            delta = delta,
                            fluctuation = "weighted",
                            max_iter = 2,
                            ipcw_fit_args = list(fit_type = "glm"),
                            g_fit_args = list(fit_type = "hal",
                                              n_bins = 5,
                                              grid_type = "equal_mass",
                                              lambda_seq =
                                                exp(seq(-1, -9,
                                                        length = 300))),
                            Q_fit_args = list(fit_type = "glm",
                                              glm_formula = "Y ~ ."),
                            eif_reg_type = "glm"
                           )
summary(tmle_hal_shift_2)

## ----make_sl, eval = FALSE----------------------------------------------------
#  # SL learners to be used for most fits (e.g., IPCW, outcome regression)
#  lrnr_lib <- make_learner_stack("Lrnr_mean", "Lrnr_glm", "Lrnr_ranger")
#  sl_learner <- Lrnr_sl$new(learners = lrnr_lib, metalearner = Lrnr_nnls$new())
#  
#  # SL learners for conditional densities to be used for the propensity score fit
#  lrnr_dens_lib <- make_learner_stack(list("Lrnr_haldensify", n_bins = 3,
#                                            grid_type = "equal_range",
#                                            lambda_seq = exp(seq(-1, -7,
#                                                                 length = 100))),
#                                      list("Lrnr_haldensify", n_bins = 3,
#                                            bin_method = "equal_mass",
#                                            lambda_seq = exp(seq(-1, -7,
#                                                                 length = 100)))
#                                     )
#  sl_learner_density <- Lrnr_sl$new(learners = lrnr_dens_lib,
#                                    metalearner = Lrnr_solnp_density$new())

## ----ipcw_tmle_sl_1_not_run, eval=FALSE---------------------------------------
#  tmle_sl_shift_1 <- txshift(W = W, A = A, Y = Y,
#                             C = C, V = c("W", "Y"),
#                             delta = delta,
#                             fluctuation = "standard",
#                             max_iter = 2,
#                             ipcw_fit_args = list(fit_type = "sl",
#                                                  sl_learnrs = sl_learner),
#                             g_fit_args = list(fit_type = "sl",
#                                               sl_learners_density =
#                                                 sl_learner_density),
#                             Q_fit_args = list(fit_type = "sl",
#                                               sl_learners = sl_learner),
#                             eif_reg_type = "glm"
#                            )
#  summary(tmle_sl_shift_1)

## ----ipcw_tmle_sl_2_not_run, eval=FALSE---------------------------------------
#  tmle_sl_shift_2 <- txshift(W = W, A = A, Y = Y,
#                             C = C, V = c("W", "Y"),
#                             delta = delta,
#                             fluctuation = "weighted",
#                             max_iter = 2,
#                             ipcw_fit_args = list(fit_type = "sl",
#                                                  sl_learners = sl_learner),
#                             g_fit_args = list(fit_type = "hal",
#                                               n_bins = 5,
#                                               grid_type = "equal_mass",
#                                               lambda_seq =
#                                                 exp(seq(-1, -9,
#                                                         length = 300))),
#                             Q_fit_args = list(fit_type = "sl",
#                                               sl_learners = sl_learner),
#                             eif_reg_type = "glm"
#                            )
#  summary(tmle_sl_shift_2)

## ----ipcw_tmle_ghal_sl_inefficient, eval = FALSE------------------------------
#  tmle_sl_ineffic <- txshift(W = W, A = A, Y = Y,
#                             C = C, V = c("W", "Y"),
#                             delta = delta,
#                             fluctuation = "standard",
#                             ipcw_fit_args = list(fit_type = "sl",
#                                                  sl_learners = sl_learner),
#                             g_fit_args = list(fit_type = "hal",
#                                               n_bins = 5,
#                                               grid_type = "equal_mass",
#                                               lambda_seq =
#                                                 exp(seq(-1, -9,
#                                                         length = 300))),
#                             Q_fit_args = list(fit_type = "sl",
#                                               sl_learners = sl_learner),
#                             ipcw_efficiency = FALSE
#                            )
#  summary(tmle_sl_ineffic)

## -----------------------------------------------------------------------------
(ci_shift <- confint(tmle_hal_shift_1))

## ----manual_ipcw_tmle---------------------------------------------------------
# compute the censoring mechanism and produce IPC weights externally
pi_mech <- plogis(W)
ipc_weights_out <- (as.numeric(C == 1) / pi_mech)[C == 1]
ipcw_out <- list(pi_mech = pi_mech, ipc_weights = ipc_weights_out)

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
                              gn_upupshift))[C == 1, ]
setnames(gn_out, c("downshift", "noshift", "upshift", "upupshift"))

# compute outcome regression externally
Qn_noshift <- (W + A - min(Y)) / diff(range(Y))
Qn_upshift <- (W + A + delta - min(Y)) / diff(range(Y))
Qn_noshift[Qn_noshift < 0] <- .Machine$double.neg.eps
Qn_noshift[Qn_noshift > 1] <- 1 - .Machine$double.neg.eps
Qn_upshift[Qn_upshift < 0] <- .Machine$double.neg.eps
Qn_upshift[Qn_upshift > 1] <- 1 - .Machine$double.neg.eps
Qn_out <- as.data.table(cbind(Qn_noshift, Qn_upshift))[C == 1, ]
setnames(Qn_out, c("noshift", "upshift"))

# invoke the wrapper function only to apply the targeting step
tmle_shift_spec <- txshift(W = W, A = A, Y = Y, delta = delta,
                           C = C, V = c("W", "Y"),
                           fluctuation = "standard",
                           max_iter = 2,
                           ipcw_fit_args = list(fit_type = "external"),
                           g_fit_args = list(fit_type = "external"),
                           Q_fit_args = list(fit_type = "external"),
                           ipcw_fit_ext = ipcw_out,
                           gn_fit_ext = gn_out,
                           Qn_fit_ext = Qn_out,
                           eif_reg_type = "glm")
summary(tmle_shift_spec)

