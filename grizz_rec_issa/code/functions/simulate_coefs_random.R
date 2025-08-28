
########################################################################
#
# Title: Function for simulating iSSF coefficients from a glmmTMB 
# mixed-effects model and including random slope variance

# Author: Eric Palm
#
#######################################################################

require(glmmTMB)
require(stringr)
require(dplyr)

# Need the model and the number of sets of coefficients you want as arguments
sim_coef_ran <-
  function(model, n) {

  # Get the fixed effects point estimates  
  mu_fixed <- glmmTMB::fixef(model)$cond
  
  # Get the full variance-covariance matrix from the model
  sig_temp <- vcov(model, full = T)
  
  # Extract the names of all the random effects
  names_rand <- attr(sig_temp, "dimnames")[[1]][-c(1:length(mu_fixed))] %>% 
    as.character()
  
  # Create a vector of mu values where mu is set to 0 for all random effects
  # Random slopes should be offsets from fixed effect coefficients
  mu_all <- c(mu_fixed, rep(0, length(names_rand)) %>% `names<-`(names_rand))
  
  # Find the index of the random intercept, which we will drop
  index_int <- which(names(mu_all) == "theta_1|stratum.1")
  
  # Drop the random intercept from the vector of mus
  mu_final <- mu_all[-index_int]
  
  # Now drop the random intercept from the vcov matrix
  sig <- sig_temp[-index_int, -index_int]

  # Simulate n sets of coefficients from a multivariate normal distribution
  coefs_raw <-
    MASS::mvrnorm(n = n, mu = mu_final, Sigma = sig) %>% 
    dplyr::as_tibble()

  # Rename the random slopes to match the fixed effects
  names_clean <- 
    coefs_raw %>% 
    dplyr::select(contains("theta")) %>% 
    colnames() %>% 
    stringr::str_sub(., 9, -6)  
  
  # Sum the coefficients for fixed effects and random slopes (offsets) for each variable
  # within each set (simulated animal) of coefficients 
  coefs_raw %>% 
    dplyr::select(!contains("theta")) %>% 
    dplyr::mutate(index = dplyr::row_number()) %>% 
    dplyr::bind_rows(., coefs_raw %>%
                dplyr::select(tidyselect::contains("theta")) %>%
                dplyr::rename_with(~ names_clean) %>% 
                dplyr::mutate(index = dplyr::row_number())) %>% 
    dplyr::mutate(dplyr::across(tidyselect::everything(), ~tidyr::replace_na(., 0))) %>% 
    tidyr::pivot_longer(cols = -index, names_to = "covariate") %>% 
    dplyr::group_by(index, covariate) %>% 
    dplyr::summarize(coef = sum(value)) %>% 
    tidyr::pivot_wider(names_from = "covariate", values_from = "coef") %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-index)
}
