options("install.lock" = FALSE)
library(readxl)
library(rlang)
library(cli)
library(writexl)
library(ggplot2)
library(dplyr)

#1. Disprove Grant & Bowen's arbitrary thresholds by plotting ND and HD against each
# other. The result should be no clear clustering into high/low HD/ND by threat status.
View(both_div)
str(both_div)

iucn_colors <- c(
  "LC" = "#60C659",  
  "NT" = "#CCE226",  
  "VU" = "#F9E814",  
  "EN" = "#FC7F3F",  
  "CR" = "#D81E05")

ndhd <- both_div %>%
  ggplot(aes(
    x = ln_ND,
    y = Haplotype_Diversity,
    color = IUCN_Status
  )) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_vline(xintercept = log(0.005), linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.7, color = "black") +
  #geom_vline(xintercept = log(0.000001), linetype = "solid", linewidth = 0.7, color = "black") +
  labs(x = "ln(Nucleotide diversity)", y = "Haplotype diversity") +
  scale_color_manual(values = iucn_colors, name = "IUCN Status") +
  theme_classic()
ndhd
ggsave("ND vs HD.png", dpi = 600, plot = ndhd, width = 6, height = 5)


#2. Separate ND vs HD plot by Mitochondrial Region
#COI
ndhd_COI <- both_div %>%
  filter(Mitochondrial_Region == "COI") %>%
  ggplot(aes(
    x = ln_ND,
    y = Haplotype_Diversity,
    color = IUCN_Status
  )) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_vline(xintercept = log(0.005), linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.7, color = "black") +
  #geom_vline(xintercept = log(0.000001), linetype = "solid", linewidth = 0.7, color = "black") +
  labs(x = "ln(Nucleotide diversity)", y = "Haplotype diversity") +
  scale_color_manual(values = iucn_colors, name = "IUCN Status") +
  theme_classic()
ndhd_COI
ggsave("ND vs HD_COI.png", dpi = 600, plot = ndhd_COI, width = 6, height = 5)

#cyt b
ndhd_cytb <- both_div %>%
  filter(Mitochondrial_Region == "cyt b") %>%
  ggplot(aes(
    x = ln_ND,
    y = Haplotype_Diversity,
    color = IUCN_Status
  )) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_vline(xintercept = log(0.005), linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.7, color = "black") +
 # geom_vline(xintercept = log(0.000001), linetype = "solid", linewidth = 0.7, color = "black") +
  labs(x = "ln(Nucleotide diversity)", y = "Haplotype diversity") +
  scale_color_manual(values = iucn_colors, name = "IUCN Status") +
  theme_classic()
ndhd_cytb
ggsave("ND vs HD_cyt b.png", dpi = 600, plot = ndhd_cytb, width = 6, height = 5)

#CR
ndhd_CR <- both_div %>%
  filter(Mitochondrial_Region == "CR") %>%
  ggplot(aes(
    x = ln_ND,
    y = Haplotype_Diversity,
    color = IUCN_Status
  )) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_vline(xintercept = log(0.005), linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.7, color = "black") +
#  geom_vline(xintercept = log(0.000001), linetype = "solid", linewidth = 0.7, color = "black") +
  labs(x = "ln(Nucleotide diversity)", y = "Haplotype diversity") +
  scale_color_manual(values = iucn_colors, name = "IUCN Status") +
  theme_classic()
ndhd_CR
ggsave("ND vs HD_CR.png", dpi = 600, plot = ndhd_CR, width = 6, height = 5)


#3 Plot ln_ND vs. IUCN Status. Use the collapsed IUCN statuses because they represent
#a more "honest version" of the graphs by accounting for extreme class imbalance:
#ND: LC->NT->VU->EN->CR: 392, 13, 9, 48, 5
#HD: 377, 13, 14, 47, 5
iucn_grouped_colors <- c(
  "More Threatened" = "#FC7F3F",
  "Less Threatened" = "#F9E814",
  "Least Concern" = "#60C659")

ND <- ggplot(Nucleotide_Diversity, aes(x = IUCN_grouped, y = ln_ND, fill = IUCN_grouped)) +
  geom_boxplot(alpha = 1,                 
               outlier.alpha = 0.6,         
               outlier.size = 2) +          
  scale_fill_manual(values = iucn_grouped_colors) +
  labs(x = "IUCN Status", 
       y = "ln(Nucleotide Diversity)", 
       title = "Trend of Nucleotide Diversity across IUCN Status") +
  theme_classic(base_size = 12) +         
  theme(legend.position = "none",         
        plot.title = element_text(face = "bold", hjust = 0.5))
ND
ggsave("ND vs IUCN Status.png", dpi = 600, plot = ND, width = 6, height = 6)

#Split the plots by Mitochondrial Region
ND_MitoRegion <- ggplot(Nucleotide_Diversity, aes(x = IUCN_grouped, y = ln_ND, fill = IUCN_grouped)) +
  geom_boxplot(alpha = 1,
               outlier.alpha = 0.6,
               outlier.size = 2) +
  scale_fill_manual(values = iucn_grouped_colors) +
  facet_wrap(~ Mitochondrial_Region, scales = "free_y") +
  labs(x = "IUCN Status",
       y = "ln(Nucleotide Diversity)",
       title = "Nucleotide Diversity across IUCN Status\nby Mitochondrial Region") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 20, hjust = 1))
ND_MitoRegion
ggsave("ND vs IUCN Status_MitoRegion.png", dpi = 600, plot = ND_MitoRegion, width = 6, height = 6)

#4 Plot HD vs IUCN Status
HD <- ggplot(Haplotype_Diversity, aes(x = IUCN_grouped, y = Haplotype_Diversity, fill = IUCN_grouped)) +
  geom_boxplot(alpha = 1,                 
               outlier.alpha = 0.6,         
               outlier.size = 2) +          
  scale_fill_manual(values = iucn_grouped_colors) +
  labs(x = "IUCN Status", 
       y = "Haplotype Diversity", 
       title = "Trend of Haplotype Diversity across IUCN Status") +
  theme_classic(base_size = 12) +         
  theme(legend.position = "none",         
        plot.title = element_text(face = "bold", hjust = 0.5))
HD
ggsave("HD vs IUCN Status.png", dpi = 600, plot = HD, width = 6, height = 6)

#Split the plots by Mitochondrial Region
HD_MitoRegion <- ggplot(Haplotype_Diversity, aes(x = IUCN_grouped, y = Haplotype_Diversity, fill = IUCN_grouped)) +
  geom_boxplot(alpha = 1,
               outlier.alpha = 0.6,
               outlier.size = 2) +
  scale_fill_manual(values = iucn_grouped_colors) +
  facet_wrap(~ Mitochondrial_Region, scales = "free_y") +
  labs(x = "IUCN Status",
       y = "Haplotype Diversity",
       title = "Haplotype Diversity across IUCN Status\nby Mitochondrial Region") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 20, hjust = 1))
HD_MitoRegion
ggsave("HD vs IUCN Status_MitoRegion.png", dpi = 600, plot = HD_MitoRegion, width = 6, height = 4)

#5 Plot dS vs. IUCN Status
dS <- ggplot(dNdS_cytB2, aes(x = IUCN_grouped, y = ln_dS, fill = IUCN_grouped)) +
  geom_boxplot() +
  scale_fill_manual(values = iucn_grouped_colors) +
  labs(x = "IUCN Status", y = "ln(dS)", title = "Trends of dS across IUCN Status") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))
dS
ggsave("dS vs IUCN Status.png", plot = dS, dpi = 600, height = 5, width = 5)
