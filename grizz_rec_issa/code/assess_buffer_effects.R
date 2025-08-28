library(tidyverse)

# 1) Find all “fixed_effect…” CSVs under results/
fixed_files <- list.files(
  path       = "/Users/claytonlamb/Dropbox/Documents/University/Work/RecreationMonitoring/MtFernie/grizz_rec_issa/results",
  pattern    = "^fixed.*\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

# If your fixed‐effect files always live directly in “results/issa_buf_<X>km/” and
# start with “fixed_…buf<X>km.csv”, you could also do:
# fixed_files <- list.files("results", pattern = "fixed.*buf.*\\.csv$", recursive = TRUE, full.names = TRUE)

# 2) Read + bind them, adding a “buffer_km” column by parsing the folder name:
combined_fixed <- fixed_files %>%
  # For each path, read it and attach a buffer‐size column
  map_dfr(function(fp) {
    # Extract “5km”, “10km”, etc., from something like ".../issa_buf_5km/fixed_effects_...buf5km.csv"
    buf_label <- str_extract(fp, "(?<=issa_buf_)[0-9]+(?=km)")
    
    read_csv(fp, show_col_types = FALSE) %>%
      mutate(buffer_km = as.integer(buf_label))
  })


a <- combined_fixed%>%
  group_by(buff, mod)%>%
  summarise(aic=mean(aic))%>%
  group_by(buff)%>%
  arrange(-aic)

# 3) Inspect the combined data
glimpse(combined_fixed)

##plot

library(dplyr)
library(stringr)
library(ggplot2)

# 1) Subset to only those rows whose term contains "trail" (case‐insensitive)
trail_df <- combined_fixed %>%
  filter(str_detect(term, regex("barrier", ignore_case = TRUE)))

# 2) Compute 95% confidence intervals around each estimate
trail_df <- trail_df %>%
  mutate(
    conf_low  = estimate - 1.96 * std.error,
    conf_high = estimate + 1.96 * std.error
  )

# 3) Now plot:
ggplot(trail_df, aes(x = buffer_km, y = estimate, color = mod)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high),
                width = 0.2, alpha = 0.7) +
  facet_wrap(~ term) +
  labs(
    x = "Buffer (km)",
    y = "Coefficient estimate",
    color = "Model",
    title = "Trail‐related fixed‐effects across buffers"
  ) +
  theme_minimal() +
  theme(
    strip.text        = element_text(face = "bold"),
    panel.spacing     = unit(0.5, "lines"),
    axis.title        = element_text(size = 12),
    plot.title        = element_text(size = 14, face = "bold"),
    legend.position   = "bottom"
  )+
  geom_hline(yintercept = 0, linetype="dashed")



library(tidyverse)
library(patchwork)

# 0) If you haven't already, make sure `combined_fixed` is in memory:
#    combined_fixed <- read_csv("results/fixed_effects_all_buffers.csv")

# 1) Choose the single model you want to visualize
selected_model <- "trail+towndist+hwy"  # <-- replace with the exact mod you want

# 2) Filter to just that model
one_model_df <- combined_fixed %>%
  filter(mod == selected_model)

# 3) Define which terms are barriers vs. distances
barrier_terms <- c(
  "barrier_trail",
  "barrier_trail_multi",
  "barrier_trail_hu",
  "barrier_rd",
  "barrier_hwy",
  "barrier_builtup",
  "barrier_river"
)

distance_terms <- c(
  "log_dist_to_trail_end",
  "log_dist_to_trail_multi_end",
  "log_dist_to_trail_hu_end",
  "log_dist_to_rd_no_hwy_end",
  "log_dist_to_hwy_end",
  "log_dist_to_builtup_end"
)

# 4) Build a small tibble for barrier effects (percent avoidance):
barrier_1model <- one_model_df %>%
  filter(term %in% barrier_terms) %>%
  transmute(
    buffer_km = factor(buffer_km, levels = sort(unique(buffer_km))),
    feature = recode(
      term,
      "barrier_trail"       = "Trail crossing",
      "barrier_trail_multi" = "Trail crossing",
      "barrier_trail_hu" = "Trail crossing",
      "barrier_rd"          = "Road crossing",
      "barrier_hwy"         = "Highway crossing",
      "barrier_builtup"     = "Town crossing",
      "barrier_river"       = "River crossing"
    ),
    avoid     = 100 * (1 - exp(estimate)),
    avoid_low = 100 * (1 - exp(estimate + 1.96 * std.error)),
    avoid_high= 100 * (1 - exp(estimate - 1.96 * std.error))
  ) %>%
  mutate(feature = factor(feature, levels = c(
    "Trail crossing",
    "Road crossing",
    "Highway crossing",
    "Town crossing",
    "River crossing"
  )))

# 5) Build a small tibble for distance effects (raw β):
distance_1model <- one_model_df %>%
  filter(term %in% distance_terms) %>%
  transmute(
    buffer_km = factor(buffer_km, levels = sort(unique(buffer_km))),
    feature = recode(
      term,
      "log_dist_to_trail_end"       = "To trail",
      "log_dist_to_trail_multi_end" = "To trail",
      "log_dist_to_trail_hu_end" = "To trail",
      "log_dist_to_rd_no_hwy_end"   = "To road",
      "log_dist_to_hwy_end"         = "To highway",
      "log_dist_to_builtup_end"     = "To town"
    ),
    estimate = estimate,
    lower    = estimate - 1.96 * std.error,
    upper    = estimate + 1.96 * std.error
  ) %>%
  mutate(feature = factor(feature, levels = c(
    "To trail",
    "To road",
    "To highway",
    "To town"
  )))

# 6) Plot Panel A: Barrier‐crossing bars, faceted by feature
p_barriers_1model <- ggplot(barrier_1model, aes(x = buffer_km, y = avoid)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_errorbar(aes(ymin = avoid_low, ymax = avoid_high),
                width = 0.2, color = "black") +
  facet_wrap(~ feature, nrow = 1, scales = "free_y") +
  labs(
    x     = "Buffer radius (km)",
    y     = "Decrease in Selection (%)",
    title = paste0("Barrier‐crossing effects for ", selected_model)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    axis.text.x     = element_text(size = 10),
    axis.title.y    = element_text(size = 11),
    strip.text      = element_text(face = "bold")
  )

# 7) Plot Panel B: Distance bars, faceted by feature
p_distance_1model <- ggplot(distance_1model, aes(x = buffer_km, y = estimate)) +
  geom_col(fill = "coral", width = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ feature, nrow = 1, scales = "free_y") +
  labs(
    x     = "Buffer radius (km)",
    y     = expression("Effect of distance ("~beta~")"),
    title = paste0("Distance effects for ", selected_model)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    axis.text.x     = element_text(size = 10),
    axis.title.y    = element_text(size = 11),
    strip.text      = element_text(face = "bold")
  )

# 8) Combine the two panels vertically
combined_1model <- p_barriers_1model / p_distance_1model + plot_layout(heights = c(1, 1))
print(combined_1model)





library(tidyverse)
library(patchwork)
library(viridis)

# 1) Filter to buffer_km == 10, dropping NAs
fixed_10km <- combined_fixed %>%
  filter(buffer_km == 10, !is.na(estimate))

# 2) Define barrier and distance term sets
barrier_terms <- c(
  "barrier_trail",
  "barrier_trail_multi",
  "barrier_trail_hu",
  "barrier_rd",
  "barrier_hwy",
  "barrier_builtup",
  "barrier_river"
)

distance_terms <- c(
  "log_dist_to_trail_end",
  "log_dist_to_trail_multi_end",
  "log_dist_to_trail_hu_end",
  "log_dist_to_rd_no_hwy_end",
  "log_dist_to_hwy_end",
  "log_dist_to_builtup_end"
)

# 3) Build Panel A data (Barriers)
barrier_10km <- fixed_10km %>%
  filter(term %in% barrier_terms) %>%
  transmute(
    mod,
    feature = recode(
      term,
      "barrier_trail"       = "Trail crossing",
      "barrier_trail_multi" = "Trail crossing",
      "barrier_trail_hu" = "Trail crossing",
      "barrier_rd"          = "Road crossing",
      "barrier_hwy"         = "Highway crossing",
      "barrier_builtup"     = "Town crossing",
      "barrier_river"       = "River crossing"
    ),
    avoid       = 100 * (1 - exp(estimate)),
    avoid_low   = 100 * (1 - exp(estimate + 1.96 * std.error)),
    avoid_high  = 100 * (1 - exp(estimate - 1.96 * std.error)),
    significant,
    p_value
  ) %>%
  mutate(feature = factor(feature, levels = c(
    "Trail crossing",
    "Road crossing",
    "Highway crossing",
    "Town crossing",
    "River crossing"
  )))

# 4) Build Panel B data (Distance)
distance_10km <- fixed_10km %>%
  filter(term %in% distance_terms) %>%
  transmute(
    mod,
    feature = recode(
      term,
      "log_dist_to_trail_end"       = "To trail",
      "log_dist_to_trail_multi_end" = "To trail",
      "log_dist_to_trail_hu_end"    = "To trail",
      "log_dist_to_rd_no_hwy_end"   = "To road",
      "log_dist_to_hwy_end"         = "To highway",
      "log_dist_to_builtup_end"     = "To town"
    ),
    estimate,
    lower       = estimate - 1.96 * std.error,
    upper       = estimate + 1.96 * std.error,
    significant,
    p_value
  ) %>%
  mutate(feature = factor(feature, levels = c(
    "To trail",
    "To road",
    "To highway",
    "To town"
  )))

# 5) Plot Panel A: Wider dodge + viridis palette
p_barriers_10km <- ggplot(barrier_10km, 
                          aes(x = feature, y = avoid, fill = mod)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5) +
  geom_errorbar(
    aes(ymin = avoid_low, ymax = avoid_high),
    position = position_dodge(width = 0.9),
    width = 0.2, alpha = 0.8
  ) +
  geom_text(
    aes(label = ifelse(significant, "*", "")),
    position = position_dodge(width = 0.9),
    vjust = -0.5, size = 5
  ) +
  scale_fill_brewer(palette = "Set3")+
  labs(
    x     = NULL,
    y     = "Decrease in Selection (%)",
    fill  = "Model",
    title = "A) Barrier‐crossing effects at 10 km buffer"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 13, face = "bold"),
    axis.text.x     = element_text(angle = 25, vjust = 1, hjust = 1),
    axis.title.y    = element_text(size = 11),
    legend.position = "bottom"
  )

# 6) Plot Panel B: Wider dodge + viridis palette
p_distance_10km <- ggplot(distance_10km, 
                          aes(x = feature, y = estimate, color = mod)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(position = position_dodge(width = 0.9), size = 3) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.9),
    width = 0.2, alpha = 0.8
  ) +
  geom_text(
    aes(label = ifelse(significant, "*", "")),
    position = position_dodge(width = 0.9),
    vjust = -0.5, size = 5
  ) +
  scale_color_brewer(palette = "Set3")+
  labs(
    x     = NULL,
    y     = expression("Effect of distance ("~beta~")"),
    color = "Model",
    title = "B) Distance‐to‐feature effects at 10 km buffer"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 13, face = "bold"),
    axis.text.x     = element_text(angle = 25, vjust = 1, hjust = 1),
    axis.title.y    = element_text(size = 11),
    legend.position = "bottom"
  )

# 7) Combine panels
combined_10km_plot <- p_barriers_10km / p_distance_10km + plot_layout(heights = c(1, 1))
print(combined_10km_plot)

