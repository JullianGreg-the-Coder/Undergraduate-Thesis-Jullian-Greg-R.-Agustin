library(rlang)
library(cli)
library(writexl)
library(ggplot2)
library(dplyr)
library(lme4)
library(emmeans)
library(tidyverse)
library(performance)
library(DHARMa)
library(pROC)
library(glmmTMB)
library(arm)
library(readxl)
library(bbmle)
library(MuMIn)
library(performance)
library(nlme)
library(patchwork)
library(partR2)
library(logistf)
library(logistf)
library(ggrepel)

## A. Nucleotide Diversity
#1. Load Nucleotide Diversity with Basin dataset. Factorize, change names, add columns.
nd_withBasin <- read_excel("justND.xlsx")
View(nd_withBasin)
str(nd_withBasin)

colnames(nd_withBasin)[12] <- "Mitochondrial_Region"
nd_withBasin$IUCN_Status <- factor(nd_withBasin$IUCN_Status,
                                   levels = c("LC","NT","VU","EN","CR"))
nd_withBasin$Mitochondrial_Region <- factor(nd_withBasin$Mitochondrial_Region)
nd_withBasin$Region <- factor(nd_withBasin$Region)
nd_withBasin$Species_Name <- factor(nd_withBasin$Species_Name)
nd_withBasin$Basin <- factor(nd_withBasin$Basin)
nd_withBasin$IUCN_grouped <- factor(case_when(
  nd_withBasin$IUCN_Status %in% c("EN", "CR") ~ "More Threatened",
  nd_withBasin$IUCN_Status %in% c("NT", "VU") ~ "Less Threatened",
  nd_withBasin$IUCN_Status == "LC" ~ "Least Concern"))


#2. Filter the ND dataset
nd_filtered <- nd_withBasin %>% filter(Number >= 5)
View(nd_filtered)
str(nd_filtered)

#3. Data Distribution: 
#Overall Distribution: with zeroes; everything else is positive, right-skewed
nd_distribution_overall <- ggplot(nd_filtered, aes(x = Nucleotide_Diversity)) +
  geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
  labs(x = "", y = "Count")
nd_distribution_overall

#By Mitochondrial Region: with zeroes; everything else is positive, right-skewed
nd_distribution_mitoRegion <- ggplot(nd_filtered, aes(x = Nucleotide_Diversity)) +
  geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
  facet_wrap(~ Mitochondrial_Region) +
  labs(x = "Nucleotide Diversity", y = "Count")
nd_distribution_mitoRegion

combined_nd_distribution <- nd_distribution_overall / nd_distribution_mitoRegion
combined_nd_distribution
ggsave("ND Distribution.png", dpi = 600, plot = combined_nd_distribution, 
       width = 6, height = 7)
#Decision = Tweedie GLM (Handles Zeroes, Gamma GLMM)

#4.Tweedie GLMM of ND Data
model_tw_filtered <- glmmTMB(Nucleotide_Diversity ~ IUCN_grouped + Mitochondrial_Region + 
                               (1|Species_Name) + (1|Basin),
                             family = tweedie(link = "log"),
                             data = nd_filtered)
summary(model_tw_filtered)

check_collinearity(model_tw_filtered)
sim_tw_filtered <- simulateResiduals(model_tw_filtered, n = 1000)
plot(sim_tw_filtered)
testDispersion(sim_tw_filtered)
testZeroInflation(sim_tw_filtered)





## B. Haplotype Diversity
#1. Load Haplotype Diversity with Basin dataset. Factorize, change names, add columns.
hd_withBasin <- read_excel("justHD.xlsx")
View(hd_withBasin)
str(hd_withBasin)

colnames(hd_withBasin)[10] <- "Number"
colnames(hd_withBasin)[11] <- "Mitochondrial_Region"
hd_withBasin$IUCN_Status <- factor(hd_withBasin$IUCN_Status,
                                   levels = c("LC","NT","VU","EN","CR"))
hd_withBasin$Mitochondrial_Region <- factor(hd_withBasin$Mitochondrial_Region)
hd_withBasin$Species_Name <- factor(hd_withBasin$Species_Name)
hd_withBasin$Basin <- factor(hd_withBasin$Basin)
hd_withBasin$IUCN_grouped <- factor(case_when(
  hd_withBasin$IUCN_Status %in% c("EN", "CR") ~ "More Threatened",
  hd_withBasin$IUCN_Status %in% c("NT", "VU") ~ "Less Threatened",
  hd_withBasin$IUCN_Status == "LC" ~ "Least Concern"))

#2. Filter the HD Dataset
hd_filtered <- hd_withBasin %>% filter(Number >= 5)
str(hd_filtered)
View(hd_filtered)

#3. Data Distribution: 
#Overall Distribution: left-skewed; everything else is positive, with data at Zero and One
hd_distribution_overall <- ggplot(hd_filtered, aes(x = Haplotype_Diversity)) +
  geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
  labs(x = "", y = "Count")
hd_distribution_overall

#By Mitochondrial Region: left-skewed, positive, with data at Zero and at One
hd_distribution_mitoRegion <- ggplot(hd_filtered, aes(x = Haplotype_Diversity)) +
  geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
  facet_wrap(~ Mitochondrial_Region) +
  labs(x = "Haplotype Diversity", y = "Count")
hd_distribution_mitoRegion

hd_distribution_combined <- hd_distribution_overall / hd_distribution_mitoRegion
hd_distribution_combined
ggsave("HD Distribution.png", dpi = 600, plot = hd_distribution_combined, 
       width = 6, height = 7)
#Decision: pursue an ordered beta GLMM because it can handle zeroes and ones in dataset.

#4. Ordered Beta GLMM of HD Data
model_HD <- glmmTMB(
  Haplotype_Diversity ~ IUCN_grouped + Mitochondrial_Region + (1|Species_Name) + (1|Basin),
  family = ordbeta(),
  data   = hd_filtered)
summary(model_HD)

#Finding: significant predictive value of "more threatened" for higher HD
#Decision: proceed to collapsing to "threatened" and "non-threatened" to see if
#the significant relationship can withstand the extra noise

#5. Ordered Beta GLMM of HD Data with a collapse to Threatened vs. Non-threatened
hd_filtered$Threatened <- num(case_when(
  hd_filtered$IUCN_Status %in% c("EN", "CR", "VU") ~ 1,
  hd_filtered$IUCN_Status %in% c("NT", "LC") ~ 0))
View(hd_filtered)

model_HD2 <- glmmTMB(
  Haplotype_Diversity ~ Threatened + Mitochondrial_Region + (1|Species_Name) + (1|Basin),
  family = ordbeta(link = "logit"),
  data   = hd_filtered)
summary(model_HD2)
#Finding: significant predictive value of "more threatened" for higher HD
#Decision: proceed to logistic regression of Threatened~Haplotype_Diversity; that is,
#can haplotype diversity predict threatened status as a diagnostic tool?

#6. What percent of the variance is explained by fixed and random effects?
r2(model_HD2)
# Conditional R2: 0.881; Marginal R2: 0.152

#Fixed EffectS: Use Nested Models to get proportion of the variance
#captured by threat status and mitochondrial region
m_region_only <- glmmTMB(Haplotype_Diversity ~ Mitochondrial_Region +
                           (1|Species_Name) + (1|Basin),
                         family = ordbeta(), data = hd_filtered)

m_threatened_only <- glmmTMB(Haplotype_Diversity ~ Threatened +
                               (1|Species_Name) + (1|Basin),
                             family = ordbeta(), data = hd_filtered)

r2_full       <- r2(model_HD2)$R2_marginal #0.1519899
r2_region     <- r2(m_region_only)$R2_marginal #0.06923112 
r2_threatened <- r2(m_threatened_only)$R2_marginal #0.09582857

unique_threatened <- r2_full - r2_region
unique_region     <- r2_full - r2_threatened

cat("Unique R2 - Threatened:          ", round(unique_threatened, 4), "\n") #0.0828
cat("Unique R2 - Mitochondrial_Region:", round(unique_region, 4),     "\n") #0.0562

#Random Effects: get the variance of species and basin in the conditional model.
# get the fixed effects variance. divide FE variance by marginal R2 to get total R2.
# divide variance of species and basin by the total variance to get the
# proportion of the variance that they capture.
vc             <- VarCorr(model_HD2)$cond
sigma2_species <- as.numeric(vc$Species_Name)
sigma2_basin   <- as.numeric(vc$Basin)
sigma2_f       <- var(predict(model_HD2, re.form = NA, type = "link")) #fixed effects var
sigma2_total   <- sigma2_f / r2(model_HD2)$R2_marginal

R2_species  <- sigma2_species / sigma2_total
R2_basin    <- sigma2_basin   / sigma2_total
R2_residual <- 1 - r2(model_HD2)$R2_conditional

cat("Species:  ", round(R2_species, 4), "\n")
cat("Basin:    ", round(R2_basin, 4),   "\n")
cat("Residual: ", round(R2_residual, 4),"\n")

#7. Logistic Regression of Haplotype Diversity
## a. HD Logistic GLMM 1: Both Random Effects Added
hd_threshold1 <- glmmTMB(Threatened ~ Haplotype_Diversity + Mitochondrial_Region + 
                        (1|Species_Name) + (1|Basin),
                    family = binomial(link = "logit"),
                    data = hd_filtered)
summary(hd_threshold1)
# Decision: species = complete separation --> Hauck-Donner Effect --> Remove species
# Complete separation is the cause for the failed convergence.

## b.Did some threshold-setting with GLMM 1. Just exploratory.
# Threshold averaged based on the weights of the mitoregions in the dataset
# higher # of entries = higher proportion.
# This code for this section was generated using Claude Opus 4.7
fe1 <- fixef(hd_threshold1)$cond
w  <- prop.table(table(hd_filtered$Mitochondrial_Region))
intercept_avg <- fe1["(Intercept)"] +
  w["CR"]    * fe1["Mitochondrial_RegionCR"] +
  w["cyt b"] * fe1["Mitochondrial_Regioncyt b"]
threshold_all1 <- -intercept_avg / fe1["Haplotype_Diversity"]
threshold_all1 #6.580579

# Also tracked probabilities of threat status based on HD from 0 --> 1
# Some weightedness applied based on the number of entries by mitoregion
fe1 <- fixef(hd_threshold1)$cond
hd_vals <- c(0, 0.5, 1)

w <- prop.table(table(hd_filtered$Mitochondrial_Region))

p_coi1  <- plogis(fe1["(Intercept)"]                                   + fe1["Haplotype_Diversity"] * hd_vals)
p_cr1   <- plogis(fe1["(Intercept)"] + fe1["Mitochondrial_RegionCR"]    + fe1["Haplotype_Diversity"] * hd_vals)
p_cytb1 <- plogis(fe1["(Intercept)"] + fe1["Mitochondrial_Regioncyt b"] + fe1["Haplotype_Diversity"] * hd_vals)

p_avg1  <- w["COI"] * p_coi1 + w["CR"] * p_cr1 + w["cyt b"] * p_cytb1

data.frame(
  HD          = hd_vals,
  COI     = round(as.numeric(p_coi1), 10),
  CR      = round(as.numeric(p_cr1), 10),
  cyt_b   = round(as.numeric(p_cytb1), 10),
  Avg_pct     = round(as.numeric(p_avg1), 10))

## b. HD Logistic GLMM 2: Removed species, retain basin
hd_threshold2 <- glmmTMB(Threatened ~ Haplotype_Diversity + Mitochondrial_Region + 
                           + (1|Basin),
                         family = binomial(link = "logit"),
                         data = hd_filtered)
summary(hd_threshold2)

#How much of the basins have only 1 threat status?
hd_filtered %>%
  group_by(Basin) %>%
  summarise(n_statuses = n_distinct(Threatened)) %>%
  pull(n_statuses) %>%
  table()
# 90 = 1 threat status
# 4 = 2 threat statuses

#Decision: No convergence, but the Hauck-Donner Effect persists --> high SE
# --> really-bad for threshold-setting when using the Intercept (Threshold = Bo / B_HD)
# --> remove basin

# c. HD Logistic  GLM 3: Remove both species, basin. Only retain fixed effects
hd_threshold3 <- glm(Threatened ~ Haplotype_Diversity + Mitochondrial_Region,
                         family = binomial(link = "logit"),
                         data = hd_filtered)
summary(hd_threshold3)
#Results: SE is still really high --> why? --> investigate MitoRegion Fixed Effect

# Is there complete separation of threat status in the MitoRegions? Yes!
hd_filtered %>%
  distinct(Species_Name, Mitochondrial_Region, Threatened) %>%
  with(table(Mitochondrial_Region, Threatened))

#Results: complete separation in COI: all entries are non-threatened
# --> Hauck-Donner Effect --> inflated SE, intercept --> bad for threshold-setting
#Dilemma: a) Remove mitoregion or b) use Firth's regression
#Decision: Try Firth's regression to accomodate the complete separation happening
#at the COI.

# c. HD Logistic GLM 4: Modified logistic regression: Firth's Regression
hd_threshold4 <- logistf(Threatened ~ Haplotype_Diversity + Mitochondrial_Region, 
                         data = hd_filtered)
summary(hd_threshold4)
#Results: HD significantly predicts threat status. SE not overinflated anymore.
#Decision: Can proceed with Threshold-setting

#Effect of HD on Threatened Probability
fe4 <- coef(hd_threshold4)
w4  <- setNames(as.vector(prop.table(table(hd_filtered$Mitochondrial_Region))),
               names(table(hd_filtered$Mitochondrial_Region)))

hd_vals <- c(0, 0.5, 1)

p_coi4  <- plogis(fe4["(Intercept)"]                                   + fe4["Haplotype_Diversity"] * hd_vals)
p_cr4   <- plogis(fe4["(Intercept)"] + fe4["Mitochondrial_RegionCR"]    + fe4["Haplotype_Diversity"] * hd_vals)
p_cytb4 <- plogis(fe4["(Intercept)"] + fe4["Mitochondrial_Regioncyt b"] + fe4["Haplotype_Diversity"] * hd_vals)
p_avg4 <- w4["COI"] * p_coi4 + w4["CR"] * p_cr4 + w4["cyt b"] * p_cytb4

data.frame(
  HD        = hd_vals,
  COI_pct   = round(as.numeric(p_coi4)  * 100, 3),
  CR_pct    = round(as.numeric(p_cr4)   * 100, 3),
  cyt_b_pct = round(as.numeric(p_cytb4) * 100, 3),
  Avg_pct   = round(as.numeric(p_avg4)  * 100, 3))

#Threshold-setting
intercept_avg4 <- fe4["(Intercept)"] +
  w4["CR"]    * fe4["Mitochondrial_RegionCR"] +
  w4["cyt b"] * fe4["Mitochondrial_Regioncyt b"]
threshold_all4 <- -intercept_avg4 / fe4["Haplotype_Diversity"]
threshold_all4

#7. Subsampling Simulations:
## a. Tried subsampling simulations with GLMM1: Threat + Mito + 1|Spec + 1|Basin
# Mainly did this as a side exploratory step in deciding whether I should keep Species
# in the threshold-setting model.
# This section was generated using Claude Opus 4.7
# i. Set Parameters
iterations <- 1000
ratios     <- c(0.25, 0.5, 0.75, 1, 1.5, 2)
results_list <- list()

# ii. Species-level subsampling
threatened_species    <- hd_filtered %>% filter(Threatened == 1) %>% pull(Species_Name) %>% unique()
nonthreatened_species <- hd_filtered %>% filter(Threatened == 0) %>% pull(Species_Name) %>% unique()
n_threatened          <- length(threatened_species)  # 18 species

set.seed(1000)
counter <- 1

for (ratio in ratios) {
  n_sample <- round(n_threatened * ratio)
  cat(sprintf("Running ratio: %sx (n_species = %d)\n", ratio, n_sample))
  for (i in 1:iterations) {
    sampled_sp <- sample(nonthreatened_species, size = n_sample, replace = FALSE)
    sim_data   <- hd_filtered %>%
      filter(Species_Name %in% c(threatened_species, sampled_sp))
    
    sim_model <- tryCatch(suppressWarnings(
      glmmTMB(Threatened ~ Haplotype_Diversity + Mitochondrial_Region +
                (1|Species_Name) + (1|Basin),
              family = binomial(link = "logit"),
              data = sim_data)),
      error = function(e) NULL)
    
    if (is.null(sim_model)) {
      thresh <- NA_real_
    } else {
      fe <- fixef(sim_model)$cond
      if (!is.finite(fe["(Intercept)"]) ||
          !is.finite(fe["Haplotype_Diversity"]) ||
          fe["Haplotype_Diversity"] == 0) {
        thresh <- NA_real_
      } else {
        thresh <- -fe["(Intercept)"] / fe["Haplotype_Diversity"]
      }
    }
    
    results_list[[counter]] <- data.frame(
      Ratio     = paste0(ratio, "x"),
      Iteration = i,
      Threshold = as.numeric(thresh))
    counter <- counter + 1
  }
}
bootstrap_results <- bind_rows(results_list)

# Factor levels
ratio_levels <- c("0.25x", "0.5x", "0.75x", "1x", "1.5x", "2x")
bootstrap_results$Ratio <- factor(bootstrap_results$Ratio, levels = ratio_levels)

# iii. Summary statistics
summary_stats <- bootstrap_results %>%
  group_by(Ratio) %>%
  summarise(
    Mean_Threshold   = mean(Threshold, na.rm = TRUE),
    Median_Threshold = median(Threshold, na.rm = TRUE),
    Mode_Threshold   = {d <- density(Threshold[is.finite(Threshold)]); d$x[which.max(d$y)]},
    CI_low           = quantile(Threshold, 0.025, na.rm = TRUE),
    CI_high          = quantile(Threshold, 0.975, na.rm = TRUE),
    n_failed         = sum(is.na(Threshold)),
    .groups = "drop"
  )
summary_stats$Ratio <- factor(summary_stats$Ratio, levels = ratio_levels)
print(summary_stats)

ratio_to_n <- c("0.25x" = "4", "0.5x" = "9", "0.75x" = "14",
                "1x"    = "18", "1.5x" = "27", "2x" = "36")

plot_data <- bootstrap_results %>%
  group_by(Ratio) %>%
  filter(is.finite(Threshold),
         Threshold >= quantile(Threshold, 0.05, na.rm = TRUE),
         Threshold <= quantile(Threshold, 0.95, na.rm = TRUE)) %>%
  ungroup()

# iv. Visualize
threshold_combined <- ggplot(plot_data, aes(x = Threshold, fill = Ratio, colour = Ratio)) +
  geom_density(alpha = 0.5, linewidth = 0.8) +
  geom_vline(data = summary_stats,
             aes(xintercept = Median_Threshold, colour = Ratio, linetype = "Median"),
             linewidth = 1) +
  geom_text_repel(data = summary_stats,
                  aes(x = Median_Threshold, y = Inf,
                      label = round(Median_Threshold, 3), colour = Ratio),
                  vjust = 1.5, 
                  direction = "y", 
                  nudge_y = -0.02, 
                  size = 4, fontface = "bold",
                  show.legend = FALSE) +
  geom_vline(xintercept = 1, colour = "firebrick", linewidth = 0.8) +
  coord_cartesian(xlim = c(-5, 30)) + 
  scale_fill_viridis_d(option = "mako", begin = 0.2, end = 0.8,
                       name = "Non-threatened\nspecies sampled", labels = ratio_to_n) +
  scale_colour_viridis_d(option = "mako", begin = 0.2, end = 0.8,
                         name = "Non-threatened\nspecies sampled", labels = ratio_to_n) +
  scale_linetype_manual(name = "Statistic",
                        values = c("Median" = "dashed")) +
  labs(title = "Distribution of 50% Probability Thresholds at Different\nSampling Ratios of Non-Threatened to Threatened",
       subtitle = "1,000 iterations per sampling ratio",
       x = "Calculated HD Threshold",
       y = "Density") +
  theme_minimal(base_size = 14) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        legend.position = "right",
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 10))
threshold_combined

ggsave("Combined Thresholds.png", plot = threshold_combined,
       width = 12, height = 6, dpi = 600, type = "cairo")

## b. Subsampling simulations with Firth's Logistic Reg: Threat + Mito
# i. Set Parameters
iterations <- 1000
ratios     <- c(0.25, 0.5, 0.75, 1, 1.5, 2)
results_list_firth <- list()

threatened_species    <- hd_filtered %>% filter(Threatened == 1) %>% pull(Species_Name) %>% unique()
nonthreatened_species <- hd_filtered %>% filter(Threatened == 0) %>% pull(Species_Name) %>% unique()
n_threatened          <- length(threatened_species)

set.seed(1000)
counter <- 1

# ii. Subsampling
for (ratio in ratios) {
  n_sample <- round(n_threatened * ratio)
  cat(sprintf("Running ratio: %sx (n_species = %d)\n", ratio, n_sample))
  for (i in 1:iterations) {
    sampled_sp <- sample(nonthreatened_species, size = n_sample, replace = FALSE)
    sim_data   <- hd_filtered %>%
      filter(Species_Name %in% c(threatened_species, sampled_sp))
    
    sim_model <- tryCatch(
      logistf(Threatened ~ Haplotype_Diversity + Mitochondrial_Region,
              data = sim_data,
              control = logistf.control(maxit = 100)),
      error = function(e) NULL)
    
    if (is.null(sim_model)) {
      thresh <- NA_real_
    } else {
      fe      <- coef(sim_model)
      w_sim   <- prop.table(table(sim_data$Mitochondrial_Region))
      int_avg <- fe["(Intercept)"] +
        as.numeric(w_sim["CR"]    %||% 0) * (fe["Mitochondrial_RegionCR"]    %||% 0) +
        as.numeric(w_sim["cyt b"] %||% 0) * (fe["Mitochondrial_Regioncyt b"] %||% 0)
      if (!is.finite(int_avg) || !is.finite(fe["Haplotype_Diversity"]) ||
          fe["Haplotype_Diversity"] == 0) {
        thresh <- NA_real_
      } else {
        thresh <- -int_avg / fe["Haplotype_Diversity"]
      }
    }
    results_list_firth[[counter]] <- data.frame(
      Ratio     = paste0(ratio, "x"),
      Iteration = i,
      Threshold = as.numeric(thresh))
    counter <- counter + 1
  }
}
bootstrap_results_firth <- bind_rows(results_list_firth)

ratio_levels <- c("0.25x", "0.5x", "0.75x", "1x", "1.5x", "2x")
bootstrap_results_firth$Ratio <- factor(bootstrap_results_firth$Ratio, levels = ratio_levels)

# iii. Summary Stats
summary_stats_firth <- bootstrap_results_firth %>%
  group_by(Ratio) %>%
  summarise(
    Mean_Threshold   = mean(Threshold, na.rm = TRUE),
    Median_Threshold = median(Threshold, na.rm = TRUE),
    Mode_Threshold   = {d <- density(Threshold[is.finite(Threshold)]); d$x[which.max(d$y)]},
    CI_low           = quantile(Threshold, 0.025, na.rm = TRUE),
    CI_high          = quantile(Threshold, 0.975, na.rm = TRUE),
    n_failed         = sum(is.na(Threshold)),
    .groups = "drop"
  )
summary_stats_firth$Ratio <- factor(summary_stats_firth$Ratio, levels = ratio_levels)
print(summary_stats_firth)

# iv. Visualization
ratio_to_n <- c("0.25x" = "4", "0.5x" = "9", "0.75x" = "14",
                "1x"    = "18", "1.5x" = "27", "2x"   = "36")

plot_data_firth <- bootstrap_results_firth %>%
  group_by(Ratio) %>%
  filter(is.finite(Threshold),
         Threshold >= quantile(Threshold, 0.05, na.rm = TRUE),
         Threshold <= quantile(Threshold, 0.95, na.rm = TRUE)) %>%
  ungroup()

threshold_firth <- ggplot(plot_data_firth, aes(x = Threshold, fill = Ratio, colour = Ratio)) +
  geom_density(alpha = 0.5, linewidth = 0.8) +
  geom_vline(data = summary_stats_firth,
             aes(xintercept = Median_Threshold, colour = Ratio, linetype = "Median"),
             linewidth = 1) +
  geom_text(data = summary_stats_firth,
            aes(x = Median_Threshold, y = Inf,
                label = round(Median_Threshold, 3), colour = Ratio),
            vjust = 1.5, hjust = -0.2, size = 4, fontface = "bold",
            show.legend = FALSE) +
  geom_vline(xintercept = 1, colour = "firebrick", linewidth = 0.8) +
  coord_cartesian(xlim = c(
    quantile(plot_data_firth$Threshold, 0.02, na.rm = TRUE),
    quantile(plot_data_firth$Threshold, 0.98, na.rm = TRUE)
  )) +
  scale_fill_viridis_d(option = "mako", begin = 0.2, end = 0.8,
                       name = "Non-threatened\nspecies sampled", labels = ratio_to_n) +
  scale_colour_viridis_d(option = "mako", begin = 0.2, end = 0.8,
                         name = "Non-threatened\nspecies sampled", labels = ratio_to_n) +
  scale_linetype_manual(name = "Statistic", values = c("Median" = "dashed")) +
  labs(title = "Distribution of 50% Probability Thresholds at Different\nSampling Ratios of Non-Threatened to Threatened",
       subtitle = "1,000 iterations per sampling ratio",
       x = "Calculated HD Threshold",
       y = "Density") +
  theme_minimal(base_size = 14) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        legend.position = "right",
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 10))
threshold_firth

ggsave("Firth Thresholds.png", plot = threshold_firth, width = 12, height = 6.5, dpi = 600)

# Finding: Threshold increases with increasing class imbalance of Non-Threatened:Threatened
# Decision: There is a significant HD-Threat Status relationship, but it is not large enough
# to entail a strong structuring between high and low HD.





# C. dS GLM
#1. Data Distribution
# Distribution: no zero values, positive, right-skewed
ds_distribution <- ggplot(dNdS_cytB2, aes(x = dS)) +
  geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
  labs(x = "dS", y = "Count")
ds_distribution
ggsave("dS distribution.png", dpi = 600, plot = ds_distribution, height = 6,
       width = 7)
# Decision: Beta GLMM. Positive + right skew + no 1s/0s

#2. Gamma GLMM
model_dS <- glmmTMB(dS ~ IUCN_grouped + (1|Species) + (1|Basin),
                    family = Gamma(link = "log"),
                    data = dNdS_cytB2)
summary(model_dS)

sim_model_dS <- simulateResiduals(model_dS, n = 1000)
plot(sim_model_dS)
testDispersion(sim_model_dS)
testZeroInflation(sim_model_dS)

#variance Partitioning
r2(model_dS)
# Conditional R2: 0.649
# Marginal R2: 0.015 --> Variance in dS explained by IUCN Status

#Random Effects
vc_dS             <- VarCorr(model_dS)$cond
sigma2_species_dS <- as.numeric(vc_dS$Species)
sigma2_basin_dS   <- as.numeric(vc_dS$Basin)
sigma2_f_dS       <- var(predict(model_dS, re.form = NA, type = "link"))
sigma2_total_dS   <- sigma2_f_dS / r2(model_dS)$R2_marginal

R2_species_dS  <- sigma2_species_dS / sigma2_total_dS
R2_basin_dS    <- sigma2_basin_dS   / sigma2_total_dS
R2_residual_dS <- 1 - r2(model_dS)$R2_conditional

cat("Species:  ", round(R2_species_dS, 5), "\n")
cat("Basin:    ", round(R2_basin_dS, 5),   "\n")
cat("Residual: ", round(R2_residual_dS, 5),"\n")

# Finding: despite the downward trend in dS with a more threatened status,
# no significant predictive relationship can be found between threat status and dS.
# Could be due to the fact that we're only using a small dataset for modelling
# Fixed effects and random effects tend to fight over variance in this case.




View(nd_filtered)
View(hd_filtered)

write_xlsx(nd_filtered, "justND_filtered.xlsx")
write_xlsx(hd_filtered, "justHD_filtered.xlsx")

table(nd_filtered$'Body of Water')
table(nd_filtered$Mitochondrial_Region)
table(nd_filtered$IUCN_Status)
View(table(nd_filtered$Basin))

table(hd_filtered$'Body of Water')
table(hd_filtered$Mitochondrial_Region)
table(hd_filtered$IUCN_Status)
View(table(hd_filtered$Basin))