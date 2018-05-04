################################################################################
#
# plot_properties.R
#


library(tidyverse)
library(ggplot2)

###############
## Import data

aod8 <- read.csv("properties/AOD8.csv") %>%
  select(SMILES, FractionCSP3, MolWt, RingCount, TotalRing, NumBridges, ringFusionDensity, ringComplexityIndex, NumRingSystems, NumChiralCenters) %>%
  mutate(library = 'Oncology Drugs')
cl <- read.csv("properties/chembridge_dvs_cl.csv") %>%
  select(SMILES, FractionCSP3, MolWt, RingCount, TotalRing, NumBridges, ringFusionDensity, ringComplexityIndex, NumRingSystems, NumChiralCenters) %>%
  mutate(library = 'Chembridge - CL')
exp <- read.csv("properties/chembridge_dvs_exp.csv") %>%
  select(SMILES, FractionCSP3, MolWt, RingCount, TotalRing, NumBridges, ringFusionDensity, ringComplexityIndex, NumRingSystems, NumChiralCenters) %>%
  mutate(library = 'Chembridge - EXP')
plm <- read.csv("properties/pleuromutilin.csv") %>%
  select(SMILES, FractionCSP3, MolWt, RingCount, TotalRing, NumBridges, ringFusionDensity, ringComplexityIndex, NumRingSystems, NumChiralCenters) %>%
  mutate(library = 'Pleuromutilin')

raw <- rbind(aod8, cl, exp, plm)

##############
## Histograms

# Fraction SP3
png(filename = "figures/hist_fraction_sp3.png",
    width=6, height=6,
    units = "in",
    res = 600)
ggplot(raw, aes(x = FractionCSP3, fill = library)) +
  geom_histogram(aes(y=..density..), binwidth = 0.1) +
  geom_vline(data = plyr::ddply(raw, "library", summarize, avg = mean(FractionCSP3)),
             aes(xintercept=avg),
             color="blue", linetype="dashed", size=1) +
  facet_wrap(~library)
dev.off()

# Num Chiral Centers
png(filename = "figures/hist_chiral_centers.png",
    width=6, height=6,
    units = "in",
    res = 600)
ggplot(raw, aes(x = NumChiralCenters, fill = library)) +
  geom_histogram(aes(y=..density..), binwidth = 1) +
  geom_vline(data = plyr::ddply(raw, "library", summarize, avg = mean(NumChiralCenters)),
             aes(xintercept=avg),
             color="blue", linetype="dashed", size=1) +
  xlim(0,12) +
  facet_wrap(~library)
dev.off()

# Ring Complexity Index
png(filename = "figures/hist_ring_complexity_index.png",
    width=6, height=6,
    units = "in",
    res = 600)
ggplot(filter(raw, !is.na(ringComplexityIndex)), aes(x = ringComplexityIndex, fill = library)) +
  geom_histogram(aes(y=..density..), binwidth = 0.1) +
  geom_vline(data = plyr::ddply(filter(raw, !is.na(ringComplexityIndex)), "library", summarize, avg = mean(ringComplexityIndex)),
             aes(xintercept=avg),
             color="blue", linetype="dashed", size=1) +
  xlim(1,2) +
  facet_wrap(~library)
dev.off()

