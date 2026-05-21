options("install.lock" = FALSE)
library(readxl)
library(rlang)
library(cli)
library(writexl)
library(dplyr)


#1. Import
dataset1 <- read_excel("Dataset v. 1.xlsx")

#2. Change the str of relevant columns. Add new columns ln_ND and IUCN_grouped.
# Remove entries with <5 sequences.
dataset1$Nucleotide_Diversity <- as.numeric(dataset1$Nucleotide_Diversity)
dataset1$Haplotype_Diversity  <- as.numeric(dataset1$Haplotype_Diversity)
dataset1$IUCN_Status <- factor(dataset1$IUCN_Status,
                          levels = c("NE", "DD", "LC", "NT", "VU", "EN", "CR"))
dataset1$ln_ND <- log((dataset1$Nucleotide_Diversity) + 0.000001)
dataset1$IUCN_grouped <- case_when(
                          dataset1$IUCN_Status %in% c("EN", "CR") ~ "More Threatened",
                          dataset1$IUCN_Status %in% c("NT", "VU") ~ "Less Threatened",
                          dataset1$IUCN_Status == "LC" ~ "Least Concern")
dataset1 <- dataset1%>% filter(Number >= 5)

#3. Subsetting - logical vectors to a) retain HD, ND only, b) retain location data only, 
# c) remove DD, NE, and d) retain COI, cyt b, CR sequences only
has_div       <- !is.na(dataset1$Nucleotide_Diversity) | !is.na(dataset1$Haplotype_Diversity)
has_location  <- !is.na(dataset1$Population_Description) & dataset1$Population_Description == "Location"
no_dd_ne     <- !(dataset1$IUCN_Status %in% c("DD", "NE"))
target_region <- dataset1$Mitochondrial_Region %in% c("COI", "cyt b", "CR")

#4. Filtering hen combine all logical vectors to subset the imported excel file
dataset1 <- dataset1[has_div & has_location & no_dd_ne & target_region, ]
dataset1$Mitochondrial_Region <- factor(dataset1$Mitochondrial_Region,
              levels = c("COI", "cyt b", "CR"))
View(dataset1)

#segue: just to test count the number of the other mtDNA regions that are
#not COI, cyt b, CR
dataset2 <- dataset1[has_div & has_location & no_dd_ne, ]
table(dataset2$Mitochondrial_Region)
dataset2 %>% count(Mitochondrial_Region)
View(dataset2)

#5. Split filtered dataset into ND only and HD only
Nucleotide_Diversity <- dataset1[!is.na(dataset1$Nucleotide_Diversity), ]
Haplotype_Diversity  <- dataset1[!is.na(dataset1$Haplotype_Diversity),  ]

View(Nucleotide_Diversity)
View(Haplotype_Diversity)

table(Nucleotide_Diversity$IUCN_Status)
table(Haplotype_Diversity$IUCN_Status)

#6. Save split ND, HD to different .xlsx files. These files contain the same information
# as justND and justHD. The only difference is that justND and justHD both have the
# Basin column which is used for downstream analysis via GLMs.
write_xlsx(Nucleotide_Diversity, "justND1.xlsx")
write_xlsx(Haplotype_Diversity,  "justHD1.xlsx")
write_xlsx(dataset1, "NDHD.xlsx")

#7. Filter once more to include entries with BOTH ND & HD. Important for the ND
# vs. HD plot later.
both_div <- dataset1[!is.na(dataset1$Nucleotide_Diversity) & !is.na(dataset1$Haplotype_Diversity), ]
View(both_div)
write_xlsx(both_div, "NDHD_BOTH.xlsx")

#8. dS Processing: factorize IUCN, Basin, IUCN_grouped. Add ln_dS column.
dNdS_cytB2 <- read_excel("dNdS2.xlsx")
View(dNdS_cytB2)
str(dNdS_cytB2)

dNdS_cytB2$IUCN <- factor(dNdS_cytB2$IUCN,
                          levels = c("LC", "NT", "VU","EN","CR"))
dNdS_cytB2$Basin <- factor(dNdS_cytB2$Basin)
dNdS_cytB2$IUCN_grouped <- factor(dNdS_cytB2$IUCN_grouped,
                          levels = c("Least Concern", "Less Threatened", "More Threatened"))
dNdS_cytB2$ln_dS <- log(dNdS_cytB2$dS)
