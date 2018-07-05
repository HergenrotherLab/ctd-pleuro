################################################################################
#
# plot_properties.R
#


library(tidyverse)
library(ggplot2)

###############
## Import data

properties <- c('SMILES', 'FractionCSP3', 'MolWt', 'RingCount', 'TotalRing', 'NumBridges', 'ringFusionDensity',
                'ringComplexityIndex', 'NumRingSystems', 'NumChiralCenters')

aod8 <- read.csv("properties/AOD8.csv") %>% select(properties) %>% mutate(library = 'Oncology Drugs')
oshea <- read.csv("properties/OSheaMoser_128.csv") %>% select(properties) %>% mutate(library = 'Antibiotics')
drugbank <- read.csv("properties/drugbank.csv") %>% select(properties) %>% mutate(library = 'Drugbank')
mlsmr <- read.csv("properties/MLSMR-NP.csv") %>% select(properties) %>% mutate(library = 'MLSMR - NP')
cl <- read.csv("properties/chembridge_dvs_cl.csv") %>% select(properties) %>% mutate(library = 'Chembridge - CL')
exp <- read.csv("properties/chembridge_dvs_exp.csv") %>% select(properties) %>% mutate(library = 'Chembridge - EXP')
microformat <- read.csv("properties/chembridge_microFormat.csv") %>% select(properties) %>% mutate(library = 'Microformat')
pnas_cc <- read.csv("properties/PNAS_CC.csv") %>% select(properties) %>% mutate(library = 'PNAS CC')
pnas_np <- read.csv("properties/PNAS_NP.csv") %>% select(properties) %>% mutate(library = 'PNAS NP')
pnas_dc <- read.csv("properties/PNAS_DC.csv") %>% select(properties) %>% mutate(library = 'PNAS DC')
plm <- read.csv("properties/pleuromutilin.csv") %>% select(properties) %>% mutate(library = 'Pleuromutilin')
lic <- read.csv("properties/lycorine.csv") %>% select(properties) %>% mutate(library = 'Lycorine')

df <- rbind(aod8, oshea, drugbank, mlsmr, cl, exp, microformat, pnas_cc, pnas_np, pnas_dc, plm, lic)

##############
## Histograms

# Fraction SP3
png(filename = "figures/hist_fraction_sp3.png",
    width=7, height=5,
    units = "in",
    res = 600)
ggplot(df, aes(x = FractionCSP3, fill = library)) +
  geom_histogram(color = "black",
                 binwidth = 0.1) +
  geom_vline(data = plyr::ddply(df, "library", summarize, avg = mean(FractionCSP3)),
             aes(xintercept=avg),
             color="blue", linetype="dashed", size=1) +
  facet_wrap(~library, scales = "free_y") +
  labs(title = "Distribution of Fraction sp3",
       x = "Fraction sp3",
       y = "Count") +
  theme(legend.position = "none")
dev.off()

# Num Chiral Centers
png(filename = "figures/hist_chiral_centers.png",
    width=7, height=5,
    units = "in",
    res = 600)
ggplot(df, aes(x = NumChiralCenters, fill = library)) +
  geom_histogram(binwidth = 1,
                 color = "black") +
  geom_vline(data = plyr::ddply(df, "library", summarize, avg = mean(NumChiralCenters)),
             aes(xintercept=avg),
             color="blue", linetype="dashed", size=1) +
  xlim(-1,12) +
  facet_wrap(~library, scales = "free_y") +
  labs(title = "Distribution of Number of Chiral Centers",
       x = "# Chiral Centers",
       y = "Count") +
  theme(legend.position = "none")
dev.off()

# Ring Complexity Index
png(filename = "figures/hist_ring_complexity_index.png",
    width=7, height=5,
    units = "in",
    res = 600)
ggplot(filter(df, !is.na(ringComplexityIndex)), aes(x = ringComplexityIndex, fill = library)) +
  geom_histogram(binwidth = 0.1, color = "black") +
  geom_vline(data = plyr::ddply(filter(df, !is.na(ringComplexityIndex)), "library", summarize, avg = mean(ringComplexityIndex)),
             aes(xintercept=avg),
             color="blue", linetype="dashed", size=1) +
  xlim(0.9,2.1) +
  facet_wrap(~library, scales = "free_y") +
  labs(title = "Distribution of Ring Complexity",
       x = "Ring Complexity Index",
       y = "Count") +
  theme(legend.position = "none")
dev.off()

####################
## Box-wisker plots

# Fraction SP3
png(filename = "figures/bw_fraction_sp3.png",
    width=7, height=5,
    units = "in",
    res = 600)
ggplot(df, aes(x = library, y = FractionCSP3)) +
  geom_boxplot(outlier.color = "NA") +
  labs(title = "Distribution of Fraction sp3",
       x = "Compound Library",
       y = "Fraction sp3") +
  coord_flip()
dev.off()

# Num Chiral Centers
png(filename = "figures/bw_chiral_centers.png",
    width=7, height=5,
    units = "in",
    res = 600)
ggplot(df, aes(x = library, y = NumChiralCenters)) +
  geom_boxplot(outlier.color = "NA") +
  ylim(0,12) +
  labs(title = "Distribution of Number of Chiral Centers",
       x = "Compound Library",
       y = "# Chiral Centers") +
  coord_flip()
dev.off()

# Ring Complexity Index
png(filename = "figures/bw_ring_complexity_index.png",
    width=7, height=5,
    units = "in",
    res = 600)
ggplot(df, aes(x = library, y = ringComplexityIndex)) +
  geom_boxplot(outlier.color = "NA") +
  ylim(1,2) +
  labs(title = "Distribution of Ring Complexity",
       x = "Compound Library",
       y = "Ring Complexity Index") +
  coord_flip()
dev.off()

################
## Violin plots
selection <- c('Oncology Drugs', 'Antibiotics', 'Drugbank', 'Chembridge - CL', 'Chembridge - EXP', 'Microformat', 'MLSMR - NP', 'PNAS CC', 'Pleuromutilin')
to_plot <- filter(df, library %in% selection)
to_plot$library <- factor(to_plot$library, levels = selection)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# Fraction SP3
png(filename = "figures/violin_fraction_sp3.png",
    width=7, height=5,
    units = "in",
    res = 600)
ggplot(to_plot, aes(x = library, y = FractionCSP3)) +
  geom_violin(scale = "width", adjust = 1.5) +
  stat_summary(fun.data=data_summary, color = "blue") +
  labs(title = "Distribution of Fraction sp3",
       x = "Compound Library",
       y = "Fraction sp3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Num Chiral Centers
png(filename = "figures/violin_chiral_centers.png",
    width=7, height=5,
    units = "in",
    res = 600)
ggplot(to_plot, aes(x = library, y = NumChiralCenters)) +
  geom_violin(scale = "width", adjust = 1.5) +
  stat_summary(fun.data=data_summary, color = "blue") +
  ylim(-2,20) +
  labs(title = "Distribution of Number of Chiral Centers",
       x = "Compound Library",
       y = "# Chiral Centers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Ring Complexity Index
png(filename = "figures/violin_ring_complexity_index.png",
    width=7, height=5,
    units = "in",
    res = 600)
ggplot(to_plot, aes(x = library, y = ringComplexityIndex)) +
  geom_violin(scale = "width", adjust = 1.5) +
  stat_summary(fun.data=data_summary, color = "blue") +
  ylim(0.9,3) +
  labs(title = "Distribution of Ring Complexity",
       x = "Compound Library",
       y = "Ring Complexity Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

###############
## 2D Hexbins

my_breaks = c(2,10,50,250,1250,6000,30000)

# Fraction sp3 vs Chiral Centers
png(filename = "figures/bin2d_chiral_fsp3.png",
    width=7, height=5,
    units = "in",
    res = 600)
ggplot(df, aes(x = NumChiralCenters, y = FractionCSP3)) +
  geom_bin2d(binwidth = c(1,0.1)) +
  scale_fill_distiller(palette= "Spectral", direction=-1, 
                       trans = "log",
                       breaks = my_breaks, labels = my_breaks) +
  xlim(NA,20) +
  facet_wrap(~library) +
  labs(title = "Fsp3 vs. # Chiral Centers",
       x = "# Chiral Centers",
       y = "Fsp3") +
  theme_bw()
dev.off()

# Fraction sp3 vs Ring Complexity
png(filename = "figures/bin2d_ring_fsp3.png",
    width=7, height=5,
    units = "in",
    res = 600)
ggplot(df, aes(x = ringComplexityIndex, y = FractionCSP3)) +
  geom_bin2d(binwidth = c(0.2,0.1)) +
  scale_fill_distiller(palette= "Spectral", direction=-1, 
                       trans = "log",
                       breaks = my_breaks, labels = my_breaks) +
  xlim(NA, 3) +
  facet_wrap(~library) +
  labs(title = "Fsp3 vs. Ring Complexity",
       x = "Ring Complexity Index",
       y = "Fsp3") +
  theme_bw()
dev.off()

# Num chiral centers vs Ring Complexity
png(filename = "figures/bin2d_chiral_ring.png",
    width=7, height=5,
    units = "in",
    res = 600)
ggplot(df, aes(x = ringComplexityIndex, y = NumChiralCenters)) +
  geom_bin2d(binwidth = c(0.2,1)) +
  scale_fill_distiller(palette= "Spectral", direction=-1, 
                       trans = "log",
                       breaks = my_breaks, labels = my_breaks) +
  lims(x = c(NA, 3), y = c(NA, 20)) +
  facet_wrap(~library) +
  labs(title = "Fsp3 vs. # Chiral Centers",
       x = "Ring Complexity Index",
       y = "# Chiral Centers") +
  theme_bw()
dev.off()
