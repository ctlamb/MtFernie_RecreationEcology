#####################################################################
#
# Title: Fit a mixed-effects iSSF in 'glmmTMB'
# Author: Eric Palm
#
#####################################################################

library(glmmTMB)
library(dplyr)
library(purrr)
library(tibble)
library(furrr)
library(stringr)
library(broom.mixed)
library(DHARMa)
library(ggplot2)
library(tidyverse)
plan(multicore)

# Import annotated ssf data frame
df_unscaled <- readRDS(file.path("data", "derived", "ssf_2_hrs_unscaled.rds"))


##summary stats on crossings
df_unscaled%>%
  group_by(case_==1)%>%
  summarise_at(vars(barrier_hwy:barrier_river),  sum)%>%
  pivot_longer(barrier_hwy:barrier_river)%>%
  pivot_wider(names_from = `case_ == 1`, values_from="value")%>%
  mutate(s=`TRUE`/(`FALSE`)*10)

# Center and scale continuous covariates within seasons, but not movement parameters
df_scaled <-
  df_unscaled %>%
  dplyr::mutate(
    dplyr::across(
      -c(id:step_id_, year, stratum, tidyselect::contains(c(
        "_01", "barrier", "sl_", "cos_ta",
        "Broadleaf", "Coniferous", "Exposed_Barren_Land", "Herb", "Mixedwood", "Rock_Rubble", "Shrubland", "Wetland", "bldg.builtup_end", "bldg.builtup_start"
      ))),
      ~ as.numeric(scale(.)[, 1])
    ),
    month = lubridate::month(t1_, abbr = TRUE, label = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  data.table::setDT()


# Plot correlation matrix
df_unscaled %>%
  # dplyr::filter(season == "winter") %>%
  dplyr::select(bldg.builtup_end:tri_end, delta_elev:log_tri_end, -contains("log")) %>%
  tidyr::drop_na() %>%
  cor() %>%
  ggcorrplot::ggcorrplot(.,
    hc.order = TRUE,
    type = "lower",
    lab = TRUE
  )

sum_mean_dat <- df_unscaled %>%
  group_by(pop) %>%
  dplyr::select(bldg.builtup_end:tri_end, delta_elev:log_tri_end) %>%
  summarise_all(mean, na.rm = T)

sum_min_dat <- df_unscaled %>%
  group_by(pop) %>%
  dplyr::select(bldg.builtup_end:tri_end, delta_elev:log_tri_end) %>%
  summarise_all(min, na.rm = T)

colnames(sum_min_dat)
gc()


### Fit models

# Shared parts of your model
base_formula_parts <- c(
  "-1",
  
##habitat and topography
  "evi_summer",
"log_sl_km:evi_summer",
  "log_tri_end",
"log_sl_km:log_tri_start",
  "elevation_end",
"log_sl_km:elevation_start",
  "Coniferous", 
"log_sl_km:Coniferous",
  "Broadleaf",
"log_sl_km:Broadleaf",
  "Shrubland",
"log_sl_km:Shrubland",
  "Mixedwood",
"log_sl_km:Mixedwood",

##distance to (town, hwy, trails added below in custom)
  "log_dist_to_rd_no_hwy_end",
  "log_sl_km:log_dist_to_rd_no_hwy_start",

##movement
  "delta_elev",
"log_sl_km:delta_elev",
  "barrier_rd",
"log_sl_km:barrier_rd",
  "barrier_river",
"log_sl_km:barrier_river",
  "cos_ta", "log_sl_km",
  "cos_ta:log_sl_km",


#random intercept (always put this first of random effects)
  "(1 | stratum)",

#random slopes
  "(0 + log_tri_end | id)",
  "(0 + evi_summer | id)",
  "(0 + Coniferous | id)", "(0 + Broadleaf | id)", "(0 + Shrubland | id)", "(0 + Mixedwood | id)",
  "(0 + log_dist_to_rd_no_hwy_end | id)",
  "(0 + delta_elev | id)",
  "(0 + cos_ta | id)", "(0 + log_sl_km | id)",
  "(0 + cos_ta:log_sl_km | id)",
  "(0 + barrier_rd | id)",
  "(0 + barrier_river | id)",
  "(0 + elevation_end | month)"
)


model_variants <- tribble(
  ~model_name,         ~custom_terms,
  "hwy",            c("log_dist_to_hwy_end", "log_sl_km:log_dist_to_hwy_start", "(0 + log_dist_to_hwy_end | id)",
                      "barrier_hwy", "log_sl_km:barrier_hwy", "(0 + barrier_hwy | id)"),
  
  "town",           c("log_dist_to_bldg_end",  "(0 + log_dist_to_bldg_end | id)",
                      "barrier_builtup", "log_sl_km:barrier_builtup", "(0 + barrier_builtup | id)"),
  
  "trail+hwy",      c("log_dist_to_trail_end", "log_sl_km:log_dist_to_trail_start", "(0 + log_dist_to_trail_end | id)",
                         "barrier_trail", "log_sl_km:barrier_trail", "(0 + barrier_trail | id)",
                      
                         "log_dist_to_hwy_end", "log_sl_km:log_dist_to_hwy_start", "(0 + log_dist_to_hwy_end | id)",
                         "barrier_hwy", "log_sl_km:barrier_hwy", "(0 + barrier_hwy | id)"),
  
  "trail+town",     c("log_dist_to_trail_end", "log_sl_km:log_dist_to_trail_start", "(0 + log_dist_to_trail_end | id)",
                         "barrier_trail", "log_sl_km:barrier_trail", "(0 + barrier_trail | id)",
                      
                         "log_dist_to_bldg_end", "log_sl_km:log_dist_to_bldg_start", "(0 + log_dist_to_bldg_end | id)",
                         "barrier_builtup", "log_sl_km:barrier_builtup", "(0 + barrier_builtup | id)"),
  
  "trail+towndist+hwy",     c("log_dist_to_trail_end", "log_sl_km:log_dist_to_trail_start", "(0 + log_dist_to_trail_end | id)",
                      "barrier_trail", "log_sl_km:barrier_trail", "(0 + barrier_trail | id)",
                      
                      "log_dist_to_bldg_end", "log_sl_km:log_dist_to_bldg_start", "(0 + log_dist_to_bldg_end | id)",
                      "barrier_builtup", "log_sl_km:barrier_builtup", "(0 + barrier_builtup | id)",
                      
                      "barrier_hwy", "log_sl_km:barrier_hwy", "(0 + barrier_hwy | id)"),

  
  "trail+town+hwydist",     c("log_dist_to_trail_end", "log_sl_km:log_dist_to_trail_start", "(0 + log_dist_to_trail_end | id)",
                      "barrier_trail", "log_sl_km:barrier_trail", "(0 + barrier_trail | id)",
                      
                      "barrier_builtup", "log_sl_km:barrier_builtup", "(0 + barrier_builtup | id)",
                      
                      "log_dist_to_hwy_end", "log_sl_km:log_dist_to_hwy_start", "(0 + log_dist_to_hwy_end | id)",
                      "barrier_hwy", "log_sl_km:barrier_hwy", "(0 + barrier_hwy | id)")

)




hist(log(df_unscaled$dist_to_bldg_start+1))
hist(df_unscaled$dist_to_builtup_end)
hist(sqrt(df_unscaled$dist_to_builtup_end))


# Function to count number of random slope terms
count_random_slopes <- function(base_parts, custom_parts) {
  n <- sum(grepl("\\|", c(base_parts, custom_parts))) - 1  # minus 1 fixed intercept
  return(n)
}

# Function to build and fit a model
build_model <- function(base_parts, extras, data, map_ids) {
  all_terms <- paste(c(base_parts, extras), collapse = " + ")
  full_formula <- as.formula(paste("case_ ~", all_terms))
  
  glmmTMB(
    formula = full_formula,
    family = poisson,
    data = data[, `:=`(id = as.factor(id), stratum = as.factor(stratum))],
    map = list(theta = factor(c(NA, map_ids))),
    start = list(theta = c(log(1e3), rep(0, length(map_ids)))),
    verbose = FALSE
  )
}

# Attach map_ids (based on # of random slopes) to model variants
model_variants <- model_variants %>%
  rowwise() %>%
  mutate(
    n_random_slopes = count_random_slopes(base_formula_parts, custom_terms),
    map_ids = list(1:n_random_slopes)
  ) %>%
  ungroup()

# Fit all models in parallel
model_results <- model_variants %>%
  mutate(model = future_pmap(list(custom_terms, map_ids),
                             ~ build_model(base_formula_parts, ..1, df_scaled, ..2)))

# Assess models
# Load necessary packages

# STEP 1: Extract top model by AIC --------------------------------------------
model_aics <- model_results %>%
  mutate(aic = map_dbl(model, AIC)) %>%
  arrange(aic)

# View sorted AIC table (lower is better)
print(model_aics %>% select(model_name, aic))

# Extract top model and name
top_model_name <- model_aics$model_name[1]
top_model <- model_aics$model[[1]]
cat("Top model selected:", top_model_name, "\n")

# STEP 2: View summary of model -----------------------------------------------
# Includes coefficient estimates, standard errors, and random effect structure
summary(top_model)

# STEP 3: Extract and inspect fixed effects with 95% confidence intervals -----
fixed_effects <- broom.mixed::tidy(top_model, effects = "fixed", conf.int = TRUE)
print(fixed_effects)

# Plot fixed effects with confidence intervals
# Watch for large error bars crossing zero → weak or uncertain predictors
fixed_effects %>%
  filter(!grepl("(Intercept)", term)) %>%
  ggplot(aes(x = estimate, y = reorder(term, estimate))) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high)) +
  theme_minimal() +
  labs(title = paste("Fixed Effects:", top_model_name), x = "Estimate", y = "Term")

# STEP 4: Residual diagnostics with DHARMa ------------------------------------
hist(predict(top_model, type = "link"), breaks = 100)
range(predict(top_model, type = "response"))


# STEP 6: Save models  -------------------------------
 saveRDS( model_results, file = "results/model.results.rds")

 



#######################################################
# Brings in the function to simulate coefficients that account for random effects
# (for those variables with random slopes)
source("code/functions/simulate_coefs_random.R")

# Create a quick plot to see how the distribution of coefficients changes when drawing
# just from fixed effects versus when including uncertainty from random effects
MASS::mvrnorm(n = 10000, Sigma = vcov(top_model)$cond, mu = glmmTMB::fixef(top_model)$cond) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(type = "fixed only") %>%
  dplyr::bind_rows(., sim_coef_ran(top_model, n = 10000) %>%
    dplyr::mutate(type = "random_slopes")) %>%
  tidyr::pivot_longer(cols = -type, names_to = "covariate") %>%
  ggplot(., aes(x = value, color = type)) +
  geom_line(stat = "density") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~covariate, scales = "free") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  theme(strip.background = element_blank())


#### Effect plot


library(dplyr)
library(broom.mixed)
library(ggplot2)
library(patchwork)

# 1. Tidy model output and compute p-values
tidy_mod <- broom.mixed::tidy(top_model, effects = "fixed") %>%
  mutate(
    p_value = 2 * pnorm(-abs(estimate / std.error)),
    significant = p_value < 0.05
  )

# 2. Panel A: Barrier effects (transform to % RSS decrease)
barriers <- tidy_mod %>%
  filter(term %in% c("barrier_trail", "barrier_trail_multi", "barrier_rd", "barrier_hwy", "barrier_builtup", "barrier_river")) %>%
  transmute(
    feature = recode(term,
                     "barrier_trail"    = "Trail crossing",
                     "barrier_trail_multi"    = "Trail crossing",
                     "barrier_rd"       = "Road crossing",
                     "barrier_hwy"      = "Highway crossing",
                     "barrier_builtup"  = "Town crossing",
                     "barrier_river" = "River crossing"
     ),
    avoid       = 100 * (1 - exp(estimate)),
    avoid_lower = 100 * (1 - exp(estimate + 1.96 * std.error)),
    avoid_upper = 100 * (1 - exp(estimate - 1.96 * std.error)),
    significant,
    p_value
  )

# 3. Panel B: Distance effects (raw β)
distance <- tidy_mod %>%
  filter(term %in% c(
    "log_dist_to_trail_end",
    "log_dist_to_trail_multi_end",
    "log_dist_to_rd_no_hwy_end",
    "log_dist_to_hwy_end",
    "log_dist_to_bldg_end"
  )) %>%
  transmute(
    feature = recode(term,
                     "log_dist_to_trail_multi_end"= "Farther from trail",
                     "log_dist_to_trail_end"      = "Farther from trail",
                     "log_dist_to_rd_no_hwy_end"  = "Farther from road",
                     "log_dist_to_hwy_end"        = "Farther from highway",
                     "log_dist_to_bldg_end"    = "Farther from town"
    ),
    estimate,
    lower       = estimate - 1.96 * std.error,
    upper       = estimate + 1.96 * std.error,
    significant,
    p_value
  )

# Step length interaction effects, interpreted as relative movement speed change
step_speed <- tidy_mod %>%
  mutate(
    feature = case_when(
      term == "barrier_river:log_sl_km" ~ "River crossing",
      term == "log_sl_km:log_dist_to_rd_no_hwy_start" ~ "Farther from road",
      term == "barrier_rd:log_sl_km" ~ "Road crossing",
      term == "log_sl_km:barrier_river" ~ "River crossing",
      term == "log_sl_km:log_dist_to_trail_start" ~ "Farther from trail",
      term == "log_sl_km:barrier_trail" ~ "Trail crossing",
      term == "log_sl_km:log_dist_to_bldg_start" ~ "Farther from town",
      term == "log_sl_km:barrier_builtup" ~ "Town crossing",
      term == "log_sl_km:barrier_hwy" ~ "Highway crossing",
      TRUE ~ NA_character_  # unmatched terms get dropped
    )
  ) %>%
  filter(!is.na(feature)) %>%
  mutate(
    percent_change = 100 * (exp(estimate) - 1),
    percent_change_low = 100 * (exp(estimate - 1.96 * std.error) - 1),
    percent_change_up = 100 * (exp(estimate + 1.96 * std.error) - 1)
  )

# 4. Plot Panel A, coloring bars by significance and adding ‘*’
p1 <- ggplot(barriers, aes(x = feature, y = avoid, fill = significant)) +
  geom_col() +
  geom_errorbar(aes(ymin = avoid_lower, ymax = avoid_upper), width = 0.2) +
  geom_text(
    aes(label = ifelse(significant, "*", "")),
    vjust = -0.5, hjust=-0.2, size = 6
  ) +
  scale_fill_manual(
    values = c("FALSE" = "grey70", "TRUE" = "firebrick")
  ) +
  labs(
    y = "Decrease in Selection (%)",
    x = NULL,
    title = "A) Barrier Crossing Effects",
    fill = "p < 0.05"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1))

# 5. Plot Panel B, coloring points by significance and adding ‘*’
p2 <- ggplot(distance, aes(x = feature, y = estimate, color = significant)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_text(
    aes(label = ifelse(significant, "*", "")),
    vjust = -0.5, hjust=-0.2, size = 6
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "firebrick")
  ) +
  labs(
    y = "Effect of Distance (β)\n(+ = avoidance, – = attraction)",
    x = NULL,
    title = "B) Distance to Feature Effects",
    color = "p < 0.05"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1))


p3 <- ggplot(step_speed, aes(x = reorder(feature, percent_change), y = percent_change, fill = significant)) +
  geom_col() +
  geom_errorbar(aes(ymin = percent_change_low, ymax = percent_change_up), width = 0.2) +
  geom_text(aes(label = ifelse(significant, "*", "")), vjust = -0.5, hjust=-0.2, size = 6) +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
  labs(
    y = "Change in Movement Speed (%)",
    x = NULL,
    title = "C) Speed Change with Environment",
    fill = "p < 0.05"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1))

# 6. Combine panels
p1 + p2 +p3

##save plot
ggsave("results/issa_rss.png", dpi=300, width=11, height=6, unit="in")

## summary stats for paper
## bears
df_scaled$id %>%
  unique() %>%
  length()

## days per bear average
df_scaled %>%
  filter(case_ == TRUE) %>%
  group_by(id) %>%
  summarise(
    n_locs = n(),
    n_days = n_distinct(t1_ %>% date())
  ) %>%
  ungroup() %>%
  summarise(
    n_locs_mean = mean(n_locs),
    locs_min=min(n_locs),
    locs_max=max(n_locs),
    n_days_mean = mean(n_days)
  )


#save outputs
# mod.output <- tibble()
# for(i in 1:nrow(model_aics)){
# # 1. Tidy model output and compute p-values
# tidy_mod <- broom.mixed::tidy(model_aics$model[[i]], effects = "fixed") %>%
#   mutate(
#     p_value = 2 * pnorm(-abs(estimate / std.error)),
#     significant = p_value < 0.05,
#     mod=model_aics$model_name[i],
#     buff=buf_dist/1000,
#     aic=AIC(model_aics$model[[i]]),
#     convergence=model_aics$model[[i]]$fit$convergence,
#     message=model_aics$model[[i]]$fit$message
#   )
# 
# # append to the rolling tibble
# mod.output <- bind_rows(mod.output, tidy_mod)
# 
# }

# 3) At this point, mod.output contains one row per fixed‐effect coefficient,
#    for every model in model_aics, all tagged with buffer_km = buf_dist/1000.

# # 4) Now write it to a CSV whose name embeds the buffer distance:
# out_dir <- file.path("results", paste0("issa_buf_", buf_dist/1000, "km"))
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 
# outfile <- file.path(out_dir, paste0("fixed_effects_summary_buf", buf_dist/1000, "km.csv"))
# readr::write_csv(mod.output, outfile)
# 
# message("Wrote summaries for buffer=", buf_dist/1000, " km to:\n  ", outfile)
# 
# ##save data used here for sanity too
# readr::write_csv(df_unscaled, file.path(out_dir, paste0("data_used", buf_dist/1000, "km.csv")))



# library(ggplot2)
# library(dplyr)
# 
# # Extract coefficients and variance-covariance components
# b_trail <- fixef(m.full.town2)$cond["log_dist_to_trail_end"]
# b_inter <- fixef(m.full.town2)$cond["log_dist_to_trail_end:dist_to_builtup_end"]
# 
# vc <- vcov(m.full.town2)$cond
# var_trail <- vc["log_dist_to_trail_end", "log_dist_to_trail_end"]
# var_inter <- vc["log_dist_to_trail_end:dist_to_builtup_end", "log_dist_to_trail_end:dist_to_builtup_end"]
# covar <- vc["log_dist_to_trail_end", "log_dist_to_trail_end:dist_to_builtup_end"]
# 
# # Sequence of distances from town (not log-transformed)
# town_seq <- seq(min(df_scaled$dist_to_builtup_end, na.rm = TRUE),
#   max(df_scaled$dist_to_builtup_end, na.rm = TRUE),
#   length.out = 100
# )
# 
# # Compute linear predictor, variance, RSS and CI
# rss_df <- tibble(
#   dist_to_town = town_seq,
#   eta = b_trail + b_inter * dist_to_town,
#   se_eta = sqrt(
#     var_trail +
#       dist_to_town^2 * var_inter +
#       2 * dist_to_town * covar
#   ),
#   RSS = exp(eta),
#   RSS_lower = exp(eta - 1.96 * se_eta),
#   RSS_upper = exp(eta + 1.96 * se_eta)
# )
# 
# # Plot with clear labels
# ggplot(rss_df, aes(x = dist_to_town, y = RSS)) +
#   geom_ribbon(aes(ymin = RSS_lower, ymax = RSS_upper), alpha = 0.3, fill = "steelblue") +
#   geom_line(color = "steelblue", size = 1.2) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
#   annotate("text", x = max(town_seq) * 0.95, y = 1.05, label = "Neutral", hjust = 1, size = 4, color = "grey30") +
#   annotate("text",
#     x = max(town_seq) * 0.95, y = max(rss_df$RSS_upper) * 0.9,
#     label = "Selection for greater\ndistance from trail", hjust = 1, size = 4, color = "forestgreen"
#   ) +
#   annotate("text",
#     x = max(town_seq) * 0.95, y = min(rss_df$RSS_lower) * 1.1,
#     label = "Selection for proximity\nto trail", hjust = 1, size = 4, color = "firebrick"
#   ) +
#   labs(
#     x = "Distance to Town (km)",
#     y = "Relative Selection Strength (RSS)",
#     title = "Trail Selection Across Distance to Town",
#     subtitle = "Interaction: log(trail dist) × town dist"
#   ) +
#   theme_minimal(base_size = 14)
# 
# 
# ## again for barrier model
# 
# # Pull the real observed range from your data
# own_seq <- seq(
#   min(df_scaled$dist_to_builtup_end, na.rm = TRUE),
#   max(df_scaled$dist_to_builtup_end, na.rm = TRUE),
#   length.out = 100
# )
# 
# # Coefficients and standard errors from m.full.town
# b_barrier_trail <- -0.327593
# se_barrier_trail <- 0.058799
# b_interaction <- 0.020753
# se_interaction <- 0.048816
# covar <- 0 # Can be extracted via vcov if needed
# 
# # Linear predictor and SE for RSS
# eta <- b_barrier_trail + b_interaction * own_seq
# se_eta <- sqrt(se_barrier_trail^2 +
#   (own_seq^2) * se_interaction^2 +
#   2 * own_seq * covar) # If covar ≠ 0
# 
# # Calculate RSS and CI
# rss <- exp(eta)
# rss_lo <- exp(eta - 1.96 * se_eta)
# rss_hi <- exp(eta + 1.96 * se_eta)
# 
# # Data frame for plotting
# df_rss <- tibble(
#   dist_to_builtup_end = own_seq,
#   RSS = rss,
#   RSS_lo = rss_lo,
#   RSS_hi = rss_hi
# )
# 
# # Plot
# library(ggplot2)
# 
# ggplot(df_rss, aes(x = dist_to_builtup_end, y = RSS)) +
#   geom_ribbon(aes(ymin = RSS_lo, ymax = RSS_hi), fill = "lightblue", alpha = 0.4) +
#   geom_line(color = "steelblue", size = 1.2) +
#   geom_hline(yintercept = 1, linetype = "dashed") +
#   labs(
#     x = "Distance to Town (km)",
#     y = "Relative Selection Strength (RSS)",
#     title = "Grizzly Bear Avoidance of Trails as a Function of Distance to Town",
#     subtitle = "RSS < 1 = Avoidance | RSS > 1 = Selection"
#   ) +
#   theme_minimal(base_size = 14)
