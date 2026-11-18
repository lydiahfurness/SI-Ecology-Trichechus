###########################################################
################# Stable Isotope Analysis #################
############ Frontiers special issue - zooarch ############
###########################################################

library(readxl)
manatee_data <- read_excel("Dropbox (UiO)/4-OCEANS Sirenia DNA/African_M_Paper/African_manatee_metadata.xlsx", 
                                       sheet = "Table_3_Isotopes_metadata", 
                                       col_types = c("text", "numeric", "text", 
                                                     "numeric", "numeric", "numeric", 
                                                     "numeric", "numeric", "numeric", 
                                                     "text", "text", "text", "text", "text", 
                                                     "text", "numeric", "text", "text", 
                                                     "numeric", "numeric", "text", "text", 
                                                     "text", "text", "text"))

library(plotly)
library(dplyr)
library(readxl)
library(SIBER)
library(ggplot2)
library(tidyverse)
library(SuessR)
library(gridExtra)
library(patchwork)
library(dunn.test)
library(ggcorrplot)
library(car)
library(MVN)

###########################################################
################## DATA CALIBRATION #######################
###########################################################

# Make individual data frames per species
# African data:
df1 <- manatee_data
df1 <- df1 %>%
  filter(grepl("^SAF", .[[1]]))
# Format the data for SuessR:
# region is left empty for now but needs to be added
df1$region <- ""
df1_formatted <- df1[, c("NLD ID", "Year", "dC13", "region")]
colnames(df1_formatted) <- c("id", "year", "d13c", "region")
df1_formatted$year <- as.numeric(as.character(df1_formatted$year))

# Amazonian data:
df2 <- manatee_data
df2 <- df2 %>%
  filter(grepl("^SAM", .[[1]]))
# Format the data for SuessR:
# region is left empty for now but needs to be added
df2$region <- ""
df2_formatted <- df2[, c("NLD ID", "Year", "dC13", "region")]
colnames(df2_formatted) <- c("id", "year", "d13c", "region")
df2_formatted$year <- as.numeric(as.character(df2_formatted$year))

# American data:
df3 <- manatee_data
df3 <- df3 %>%
  filter(grepl("^STM", .[[1]]))
# Format the data for SuessR:
# region is left empty for now but needs to be added
df3$region <- ""
df3_formatted <- df3[, c("NLD ID", "Year", "dC13", "region")]
colnames(df3_formatted) <- c("id", "year", "d13c", "region")
df3_formatted$year <- as.numeric(as.character(df3_formatted$year))

# Note on r for the custom corrections #
# i have taken an estimate of 8 for freshwater and 5 for marine
# based on bigger phytoplankton radius in fresh systems
# cite Litchman et al 2008, Hellebrand et al 2022, Acevedo-Trejos et al 2015

### AFRICAN DATA CORRECTION
# Create custom data:
custom_region_data_africa <- df1_formatted
custom_region_data_africa$region <- "Africa"

# Custom correction: Freshwater
years <- 1850:2000
year1 <- 1850
year2 <- 2000
d13c.change <- 0.4
reg.up <- reg.uptake(year1 = year1,year2 = year2, d13c.change = d13c.change)
## taken from Eire et al 2017 - estimation for surface water (upper 100m) in Atlantic tropics/subtropics
## i do not have sufficient data to calculate using known years so have to use this instead
custom_africa_fresh_data <- data.frame(
  year = years,
  region = rep("Africa", length(years)),
  r = rep(8, length(years)),  # Phytoplankton cell radius (larger in rivers) - see references in exploratory analyses
  sst = rep(28.5, length(years)),  # Approximate temperature - compared across locations and looks similar across the range (Senegal slightly cooler)
  sss = rep(0, length(years)),  # Freshwater has salinity = 0 :)
  CO2atm = SuessR.reference.data$CO2atm[1:151], # atmospheric CO₂ values 
  up.con = rep(reg.up, length(years)),  # Uptake constant
  Cp = rep(0.41, length(years))  # Phytoplankton carbon content - just left as is in example
)

# correct all to 1950
corrected_africa_fresh <- SuessR.custom(
  data = custom_region_data_africa, 
  custom.region.data = custom_africa_fresh_data, 
  correct.to = 1950  # Adjust δ¹³C to a suitable baseline for the angle of my study
)

# Custom correction: Marine
# uptake will be the same so no need to change
custom_africa_marine_data <- data.frame(
  year = years,
  region = rep("Africa", length(years)),
  r = rep(5, length(years)),  # Phytoplankton cell radius
  sst = rep(28.5, length(years)),  # Approximate temperature - compared across locations and looks similar across the range (Senegal slightly cooler)
  sss = rep(35, length(years)),  # SSS estimate from https://www.mdpi.com/2072-4292/16/9/1578
  CO2atm = SuessR.reference.data$CO2atm[1:151], # atmospheric CO₂ values 
  up.con = rep(reg.up, length(years)),  # Uptake constant
  Cp = rep(0.41, length(years))  # Phytoplankton carbon content - just left as is in example
)

# correct to 1950
corrected_africa_marine <- SuessR.custom(
  data = custom_region_data_africa, 
  custom.region.data = custom_africa_marine_data, 
  correct.to = 1950)

# Freshwater corrected dataset:
fresh_df1 <- df1 %>%
  left_join(corrected_africa_fresh, by = c("NLD ID" = "id"))
# Marine corrected dataset:
marine_df1 <- df1 %>%
  left_join(corrected_africa_marine, by = c("NLD ID" = "id"))

# Now just plot C/N for each corrected dataset and the uncorrected dataset:
# Plot 1: Uncorrected
df1p1 <- ggplot(df1, aes(x = dC13, y = dN15)) +
  geom_point() +
  stat_ellipse(level = 0.95, color = "lightblue", linetype = "dashed") +
  labs(
    x = expression(delta^13 * "C"),
    y = expression(delta^15 * "N"),
    title = expression("Uncorrected " * delta^13 * "C data")
  ) +
  theme_minimal()

# Plot 2: Corrected for freshwater
df1p2 <- ggplot(fresh_df1, aes(x = d13c.cor, y = dN15)) +
  geom_point() +
  stat_ellipse(level = 0.95, color = "lightblue", linetype = "dashed") +
  labs(
    x = expression(delta^13 * "C"),
    y = expression(delta^15 * "N"),
    title = expression("Corrected for freshwater environment " * delta^13 * "C data")
  ) +
  theme_minimal()

# Plot 3: Corrected for marine
df1p3 <- ggplot(marine_df1, aes(x = d13c.cor, y = dN15)) +
  geom_point() +
  stat_ellipse(level = 0.95, color = "lightblue", linetype = "dashed") +
  labs(
    x = expression(delta^13 * "C"),
    y = expression(delta^15 * "N"),
    title = expression("Corrected for marine environment " * delta^13 * "C data")
  ) +
  theme_minimal()

# Compare plots:
df1p1 + df1p2 + df1p3

### AMAZONIAN DATA CORRECTION
# I will not correct for marine here as they are exclusively freshwater

# Create custom data:
custom_region_data_amazon <- df2_formatted
custom_region_data_amazon$region <- "Amazon"

# first calculate the uptake constant using known years
year1 <- 1917
year2 <- 1943
d13c.change <- 0.4725 # this is mean of 4 x 1886 d13C compared to 1910 sample (002 could be from another region/is captive, 008 not in MT tree)
reg.up <- reg.uptake(year1 = year1,year2 = year2, d13c.change = d13c.change)

# here i only need to correct for freshwater
years <- 1850:1970
custom_amazon_data <- data.frame(
  year = years,
  region = rep("Amazon", length(years)),
  r = rep(8, length(years)),  # Phytoplankton cell radius (larger in rivers)
  sst = rep(27, length(years)),  # Approximate river water temperature from dos Santos et al 2014 
  sss = rep(0, length(years)),  # Freshwater has salinity = 0 :)
  CO2atm = SuessR.reference.data$CO2atm[1:121],  # atmospheric CO₂ values
  up.con = rep(reg.up, length(years)),  # Uptake constant
  Cp = rep(0.41, length(years))  # Phytoplankton carbon content - just left as is in example
)

# correct to 1950
corrected_amazon <- SuessR.custom(
  data = custom_region_data_amazon, 
  custom.region.data = custom_amazon_data, 
  correct.to = 1950  # Adjust δ¹³C to a preindustrial baseline 
)

# Freshwater corrected dataset:
fresh_df2 <- df2 %>%
  left_join(corrected_amazon, by = c("NLD ID" = "id"))

# Now just plot C/N for the corrected dataset and the uncorrected dataset:
# Plot 1: Uncorrected
df2p1 <- ggplot(df2, aes(x = dC13, y = dN15)) +
  geom_point() +
  stat_ellipse(level = 0.95, color = "lightblue", linetype = "dashed") +
  labs(
    x = expression(delta^13 * "C"),
    y = expression(delta^15 * "N"),
    title = expression("Uncorrected " * delta^13 * "C data")
  ) +
  theme_minimal()

# Plot 2: Corrected for freshwater
df2p2 <- ggplot(fresh_df2, aes(x = d13c.cor, y = dN15)) +
  geom_point() +
  stat_ellipse(level = 0.95, color = "lightblue", linetype = "dashed") +
  labs(
    x = expression(delta^13 * "C"),
    y = expression(delta^15 * "N"),
    title = expression("Corrected for freshwater environment " * delta^13 * "C data")
  ) +
  theme_minimal()

# Compare plots:
df2p1 + df2p2

### AMERICAN DATA CORRECTION
# Create custom data:
custom_region_data_caribbean <- df3_formatted
custom_region_data_caribbean$region <- "Caribbean"
year1 <- 1876
year2 <- 1971
d13c.change <- 0.225
reg.up <- reg.uptake(year1 = year1,year2 = year2, d13c.change = d13c.change)

years <- 1850:2000
# here i correct for freshwater
custom_caribbean_fresh_data <- data.frame(
  year = years,
  region = rep("Caribbean", length(years)),
  r = rep(8, length(years)),  # Phytoplankton cell radius (larger in rivers)
  sst = rep(27, length(years)),  # Approximate temperature - seems like 24 lowest and 30 highest so split the difference
  sss = rep(0, length(years)),  # Freshwater has salinity = 0 :)
  CO2atm = SuessR.reference.data$CO2atm[1:151], # co2atm_values,  # atmospheric CO₂ values # fix this tho position correctly and use earliest value as constant before 1850
  up.con = rep(reg.up, length(years)),  # Uptake constant
  Cp = rep(0.41, length(years))  # Phytoplankton carbon content - just left as is in example
)

# correct to 1950
corrected_caribbean_fresh <- SuessR.custom(
  data = custom_region_data_caribbean, 
  custom.region.data = custom_caribbean_fresh_data, 
  correct.to = 1950  # Adjust δ¹³C to a preindustrial baseline
)

year1 <- 1850
year2 <- 2000
d13c.change <- 0.5
reg.up <- reg.uptake(year1 = year1,year2 = year2, d13c.change = d13c.change)
## taken from Eire et al 2017 - estimation for surface water (upper 100m) in Atlantic tropics/subtropics
# I increased it to 0.5 (difference between 2 zones N and S Atlantic)
# Caribbean is absent from exact estimations in this paper

years <- 1850:2000
# here i correct for freshwater
custom_caribbean_marine_data <- data.frame(
  year = years,
  region = rep("Caribbean", length(years)),
  r = rep(5, length(years)),  # Phytoplankton cell radius 
  sst = rep(27, length(years)),  # Approximate temperature - seems like 24 lowest and 30 highest so split the difference
  sss = rep(36, length(years)), 
  CO2atm = SuessR.reference.data$CO2atm[1:151], # co2atm_values,  # atmospheric CO₂ values # fix this tho position correctly and use earliest value as constant before 1850
  up.con = rep(reg.up, length(years)),  # Uptake constant
  Cp = rep(0.41, length(years))  # Phytoplankton carbon content - just left as is in example
)

# correct to 1950
corrected_caribbean_marine <- SuessR.custom(
  data = custom_region_data_caribbean, 
  custom.region.data = custom_caribbean_marine_data, 
  correct.to = 1950  # Adjust δ¹³C to a preindustrial baseline
)

# Freshwater corrected dataset:
fresh_df3 <- df3 %>%
  left_join(corrected_caribbean_fresh, by = c("NLD ID" = "id"))
# Marine corrected dataset:
marine_df3 <- df3 %>%
  left_join(corrected_caribbean_marine, by = c("NLD ID" = "id"))

# Now just plot C/N for each corrected dataset and the uncorrected dataset:
# Plot 1: Uncorrected
df3p1 <- ggplot(df3, aes(x = dC13, y = dN15)) +
  geom_point() +
  stat_ellipse(level = 0.95, color = "lightblue", linetype = "dashed") +
  labs(
    x = expression(delta^13 * "C"),
    y = expression(delta^15 * "N"),
    title = expression("Uncorrected " * delta^13 * "C data")
  ) +
  theme_minimal()

# Plot 2: Corrected for freshwater
df3p2 <- ggplot(fresh_df3, aes(x = d13c.cor, y = dN15)) +
  geom_point() +
  stat_ellipse(level = 0.95, color = "lightblue", linetype = "dashed") +
  labs(
    x = expression(delta^13 * "C"),
    y = expression(delta^15 * "N"),
    title = expression("Corrected for freshwater environment " * delta^13 * "C data")
  ) +
  theme_minimal()

# Plot 3: Corrected for marine
df3p3 <- ggplot(marine_df3, aes(x = d13c.cor, y = dN15)) +
  geom_point() +
  stat_ellipse(level = 0.95, color = "lightblue", linetype = "dashed") +
  labs(
    x = expression(delta^13 * "C"),
    y = expression(delta^15 * "N"),
    title = expression("Corrected for marine environment " * delta^13 * "C data")
  ) +
  theme_minimal()

# Compare plots:
df3p1 + df3p2 + df3p3

### DECISION:
# i will use fresh for amazon and marine corrections for the others - overall
# very little difference between the corrections...

###########################################################
################### DATA FOR SIBER ########################
###########################################################

updated_manatee_data <- manatee_data %>%
  left_join(select(corrected_amazon, id, temp1 = colnames(corrected_amazon)[7]), by = c("NLD ID" = "id")) %>%
  left_join(select(corrected_africa_marine, id, temp2 = colnames(corrected_africa_marine)[7]), by = c("NLD ID" = "id")) %>%
  left_join(select(corrected_caribbean_marine, id, temp3 = colnames(corrected_caribbean_marine)[7]), by = c("NLD ID" = "id")) %>%
  mutate(d13c.cor = coalesce(temp1, temp2, temp3)) %>%  # Combine into one column
  select(-temp1, -temp2, -temp3)

## plot with SIBER:
siber_data <- updated_manatee_data %>%
  select(d13c.cor, dN15, Species) %>%
  mutate(
    group = as.numeric(factor(Species)),
    community = 1 
  ) %>%
  rename(
    iso1 = d13c.cor, 
    iso2 = dN15
  ) %>%
  select(iso1, iso2, group, community) 

siber_object <- createSiberObject(siber_data)

plotSiberObject(
  siber_object,
  ax.pad = 2,
  hulls = FALSE,
  ellipses = TRUE,
  group.hulls = FALSE,
  bty = "n",
  xlab = expression(paste(delta^{13}, "C")),  
  ylab = expression(paste(delta^{15}, "N")),
  x.limits = c(-30, 0),
  y.limits = c(0, 14)
)

# Add extra tick marks every 1 unit (keep default labels)
axis(1, at = seq(-30, 0, by = 1), labels = FALSE, tck = -0.02)
axis(2, at = seq(0, 14, by = 1), labels = FALSE, tck = -0.02)

# Legend
legend(
  "topright",
  legend = c("T. senegalensis", "T. manatus", "T. inunguis"),
  col = c("cyan", "red", "black"),
  pch = 1,
  cex = 1.0,
  bty = "n",
  text.font = 3
)

# Title
title(main = expression(paste("", delta^{13}, "C / ", delta^{15}, "N by species")),
      font.main = 2,
      cex.main = 2.0)

###
my_cols <- c("black", "red", "cyan")

palette(my_cols)

plotSiberObject(
  siber_object,
  ax.pad = 2,
  hulls = FALSE,
  ellipses = TRUE,
  group.hulls = FALSE,
  bty = "n",
  xlab = expression(paste(delta^{13}, "C")),  
  ylab = expression(paste(delta^{15}, "N")),
  x.limits = c(-30, 0),
  y.limits = c(0, 14)
)

# Add extra tick marks every 1 unit (keep default labels)
axis(1, at = seq(-30, 0, by = 1), labels = FALSE, tck = -0.02)
axis(2, at = seq(0, 14, by = 1), labels = FALSE, tck = -0.02)

# Legend
legend(
  "topright",
  legend = c("T. senegalensis", "T. manatus", "T. inunguis"),
  col = c("cyan", "red", "black"),
  pch = 1,
  cex = 1.0,
  bty = "n",
  text.font = 3
)

# Title
title(main = expression(paste("", delta^{13}, "C / ", delta^{15}, "N by species")),
      font.main = 2,
      cex.main = 1.5)

###########################################################
############# REMAKE BOX AND SCATTERPLOTS #################
###########################################################

# Create plots on species level
updated_manatee_data_grouped <- updated_manatee_data %>%
  group_by(Species) %>%
  mutate(n = n())

my_cols <- c(
  adjustcolor("black", alpha.f = 0.75),
  "red",
  "cyan"
)

p_main_C <- ggplot(updated_manatee_data_grouped, aes(x = Species, y = d13c.cor, fill = Species)) +
  geom_boxplot(data = subset(updated_manatee_data_grouped, n > 3), outlier.shape = NA) +
  geom_jitter(width = 0.2, shape = 21, color = "black", size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = expression(paste(delta^{13}, "C per Species")),
    y = expression(delta^{13}~"C"),
    x = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_manual(values = my_cols)

p_main_N <- ggplot(updated_manatee_data_grouped, aes(x = Species, y = dN15, fill = Species)) +
  geom_boxplot(data = subset(updated_manatee_data_grouped, n > 3), outlier.shape = NA) +
  geom_jitter(width = 0.2, shape = 21, color = "black", size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = expression(paste(delta^{15}, "N per Species")),
    y = expression(delta^{15}~"N"),
    x = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_manual(values = my_cols)

p_main_S <- ggplot(updated_manatee_data_grouped, aes(x = Species, y = dS34, fill = Species)) +
  geom_boxplot(data = subset(updated_manatee_data_grouped, n > 3), outlier.shape = NA) +
  geom_jitter(width = 0.2, shape = 21, color = "black", size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = expression(paste(delta^{34}, "S per Species")),
    y = expression(delta^{34}~"S"),
    x = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_manual(values = my_cols) +
  guides(fill = guide_legend(
    label.theme = element_text(face = "italic")
  ))

p_main_C <- p_main_C + theme(legend.position = "none")  
p_main_N <- p_main_N + theme(legend.position = "none")   
p_main_S <- p_main_S + theme(legend.position = "right")
grid.arrange(p_main_C, p_main_N, p_main_S, ncol = 3, widths = c(1, 1, 1.5))

p_main_CN <- ggplot(updated_manatee_data_grouped, aes(x = d13c.cor, y = dN15, color = Species)) +
  geom_jitter(width = 0.2, height = 0.2, size = 3) +
  stat_ellipse(level = 0.95, size = 0.75, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = expression(delta^{13}~"C"),
    y = expression(delta^{15}~"N"),
    title = expression(paste(delta^{15}, "N/", delta^{13}, "C by Species"))
  ) +
  scale_color_manual(values = my_cols) +
  guides(color = guide_legend(label.theme = element_text(face = "italic")))

p_main_NS <- ggplot(updated_manatee_data_grouped, aes(x = dS34, y = dN15, color = Species)) +
  geom_jitter(width = 0.2, height = 0.2, size = 3) +
  stat_ellipse(level = 0.95, size = 0.75, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = expression(delta^{34}~"S"),
    y = expression(delta^{15}~"N"),
    title = expression(paste(delta^{15}, "N/", delta^{34}, "S by Species"))
  ) +
  scale_color_manual(values = my_cols) +
  guides(color = guide_legend(label.theme = element_text(face = "italic")))

p_main_CS <- ggplot(updated_manatee_data_grouped, aes(x = d13c.cor, y = dS34, color = Species)) +
  geom_jitter(width = 0.2, height = 0.2, size = 3) +
  stat_ellipse(level = 0.95, size = 0.75, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = expression(delta^{13}~"C"),
    y = expression(delta^{34}~"S"),
    title = expression(paste(delta^{34}, "S/", delta^{13}, "C by Species"))
  ) +
  scale_color_manual(values = my_cols) +
  guides(color = guide_legend(label.theme = element_text(face = "italic")))

p_main_CN <- p_main_CN + theme(legend.position = "none")
p_main_NS <- p_main_NS + theme(legend.position = "none")
p_main_CS <- p_main_CS + theme(legend.position = "right")

grid.arrange(p_main_CN, p_main_NS, p_main_CS, ncol = 3, widths = c(1, 1, 1.5))

###########################################################
############ LOCALITY BOX AND SCATTERPLOTS ################
###########################################################

amazon <- updated_manatee_data %>%
  filter(grepl("^SAM", .[[1]]))
caribbean <- updated_manatee_data %>%
  filter(grepl("^STM", .[[1]]))
africa <- updated_manatee_data %>%
  filter(grepl("^SAF", .[[1]]))

### PLOT AMAZON

# Count points per locality
amazon_counts <- amazon %>%
  group_by(Locality) %>%
  mutate(n_points = n())

AM_p_main_C <- ggplot(amazon_counts, aes(x = Locality, y = d13c.cor)) +
  geom_boxplot(
    data = subset(amazon_counts, n_points >= 3),
    fill = adjustcolor("black", alpha.f = 0.75),
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.2,
    shape = 21,
    color = "black",
    size = 2,
    alpha = 0.8
  ) +
  theme_minimal() +
  labs(
    title = expression(paste(delta^{13}, "C per Locality")),
    y = expression(delta^{13}~"C"),
    x = "Locality"
  )

AM_p_main_N <- ggplot(amazon_counts, aes(x = Locality, y = dN15)) +
  geom_boxplot(
    data = subset(amazon_counts, n_points >= 3),
    fill = adjustcolor("black", alpha.f = 0.75),
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.2,
    shape = 21,
    color = "black",
    size = 2,
    alpha = 0.8
  ) +
  theme_minimal() +
  labs(
    title = expression(paste(delta^{15}, "N per Locality")),
    y = expression(delta^{15}~"N"),
    x = "Locality"
  )

AM_p_main_S <- ggplot(amazon_counts, aes(x = Locality, y = dS34)) +
  geom_boxplot(
    data = subset(amazon_counts, n_points >= 3),
    fill = adjustcolor("black", alpha.f = 0.75),
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.2,
    shape = 21,
    color = "black",
    size = 2,
    alpha = 0.8
  ) +
  theme_minimal() +
  labs(
    title = expression(paste(delta^{34}, "S per Locality")),
    y = expression(delta^{34}~"S"),
    x = "Locality"
  )

AM_p_main_C <- AM_p_main_C + theme(legend.position = "none")  
AM_p_main_N <- AM_p_main_N + theme(legend.position = "none")   
AM_p_main_S <- AM_p_main_S + theme(legend.position = "right")
grid.arrange(AM_p_main_C, AM_p_main_N, AM_p_main_S, ncol = 3, widths = c(1, 1, 1))

AM_p_main_CN <- ggplot(amazon, aes(x = d13c.cor, y = dN15, shape = Locality)) +
  geom_jitter(width = 0.2, height = 0.2, size = 4) +
  geom_text(
    aes(label = paste(`NLD ID`, amazon$`Specific locality`, sep = "\n")),      
    size = 1.5,                                         
    vjust = -0.1,                                       
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(
    x = expression(delta^{13}~"C"),
    y = expression(delta^{15}~"N"),
    title = expression(paste(delta^{15}, "N/", delta^{13}, "C by Locality"))
  ) +
  scale_color_manual(adjustcolor("black", alpha.f = 0.75))

AM_p_main_NS <- ggplot(amazon, aes(x = dS34, y = dN15, shape = Locality)) +
  geom_jitter(width = 0.2, height = 0.2, size = 4) +
  geom_text(
    aes(label = paste(`NLD ID`, amazon$`Specific locality`, sep = "\n")),      
    size = 1.5,                                         
    vjust = -0.1,                                       
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(
    x = expression(delta^{34}~"S"),
    y = expression(delta^{15}~"N"),
    title = expression(paste(delta^{15}, "N/", delta^{34}, "S by Locality"))
  ) +
  scale_color_manual(adjustcolor("black", alpha.f = 0.75))

AM_p_main_CS <- ggplot(amazon, aes(x = d13c.cor, y = dS34, shape = Locality)) +
  geom_jitter(width = 0.2, height = 0.2, size = 4) +
  geom_text(
    aes(label = paste(`NLD ID`, amazon$`Specific locality`, sep = "\n")),      
    size = 1.5,                                         
    vjust = -0.1,                                       
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(
    x = expression(delta^{13}~"C"),
    y = expression(delta^{34}~"S"),
    title = expression(paste(delta^{34}, "S/", delta^{13}, "C by Locality"))
  ) +
  scale_color_manual(adjustcolor("black", alpha.f = 0.75))

AM_p_main_CN <- AM_p_main_CN + theme(legend.position = "none")
AM_p_main_NS <- AM_p_main_NS + theme(legend.position = "none")
AM_p_main_CS <- AM_p_main_CS + theme(legend.position = "right")

grid.arrange(AM_p_main_CN, AM_p_main_NS, AM_p_main_CS, ncol = 3, widths = c(1, 1, 1.2))

## PLOT CARIBBEAN

# Count points per locality
caribbean_counts <- caribbean %>%
  group_by(Locality) %>%
  mutate(n_points = n())

CA_p_main_C <- ggplot(caribbean_counts, aes(x = Locality, y = d13c.cor)) +
  geom_boxplot(
    data = subset(caribbean_counts, n_points >= 3),
    fill = "red",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.2,
    shape = 21,
    color = "black",
    fill = "red",
    size = 2,
    alpha = 0.8
  ) +
  theme_minimal() +
  labs(
    title = expression(paste(delta^{13}, "C per Locality")),
    y = expression(delta^{13}~"C"),
    x = "Locality"
  )

CA_p_main_N <- ggplot(caribbean_counts, aes(x = Locality, y = dN15)) +
  geom_boxplot(
    data = subset(caribbean_counts, n_points >= 3),
    fill = "red",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.2,
    shape = 21,
    color = "black",
    fill = "red",
    size = 2,
    alpha = 0.8
  ) +
  theme_minimal() +
  labs(
    title = expression(paste(delta^{15}, "N per Locality")),
    y = expression(delta^{15}~"N"),
    x = "Locality"
  )

CA_p_main_S <- ggplot(caribbean_counts, aes(x = Locality, y = dS34)) +
  geom_boxplot(
    data = subset(caribbean_counts, n_points >= 3),
    fill = "red",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.2,
    shape = 21,
    color = "black",
    fill = "red",
    size = 2,
    alpha = 0.8
  ) +
  theme_minimal() +
  labs(
    title = expression(paste(delta^{34}, "S per Locality")),
    y = expression(delta^{34}~"S"),
    x = "Locality"
  )

CA_p_main_C <- CA_p_main_C + theme(legend.position = "none")  
CA_p_main_N <- CA_p_main_N + theme(legend.position = "none")   
CA_p_main_S <- CA_p_main_S + theme(legend.position = "none")
grid.arrange(CA_p_main_C, CA_p_main_N, CA_p_main_S, ncol = 3, widths = c(1, 1, 1))

CA_p_main_CN <- ggplot(caribbean, aes(x = d13c.cor, y = dN15, shape = Locality)) +
  geom_jitter(width = 0.2, height = 0.2, size = 4, color = adjustcolor("red")) +
  geom_text(
    aes(label = paste(`NLD ID`, caribbean$`Specific locality`, sep = "\n")),      
    size = 1.5,                                         
    vjust = -0.1,                                       
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(
    x = expression(delta^{13}~"C"),
    y = expression(delta^{15}~"N"),
    title = expression(paste(delta^{15}, "N/", delta^{13}, "C by Locality"))
  ) +
  scale_color_manual(adjustcolor("red", alpha.f = 0.75))

CA_p_main_NS <- ggplot(caribbean, aes(x = dS34, y = dN15, shape = Locality)) +
  geom_jitter(width = 0.2, height = 0.2, size = 4, color = adjustcolor("red")) +
  geom_text(
    aes(label = paste(`NLD ID`, caribbean$`Specific locality`, sep = "\n")),      
    size = 1.5,                                         
    vjust = -0.1,                                       
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(
    x = expression(delta^{34}~"S"),
    y = expression(delta^{15}~"N"),
    title = expression(paste(delta^{15}, "N/", delta^{34}, "S by Locality"))
  ) +
  scale_color_manual(adjustcolor("red", alpha.f = 0.75))

CA_p_main_CS <- ggplot(caribbean, aes(x = d13c.cor, y = dS34, shape = Locality)) +
  geom_jitter(width = 0.2, height = 0.2, size = 4, color = adjustcolor("red")) +
  geom_text(
    aes(label = paste(`NLD ID`, caribbean$`Specific locality`, sep = "\n")),      
    size = 1.5,                                         
    vjust = -0.1,                                       
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(
    x = expression(delta^{13}~"C"),
    y = expression(delta^{34}~"S"),
    title = expression(paste(delta^{34}, "S/", delta^{13}, "C by Locality"))
  ) +
  scale_color_manual(adjustcolor("red", alpha.f = 0.75))

CA_p_main_CN <- CA_p_main_CN + theme(legend.position = "none")
CA_p_main_NS <- CA_p_main_NS + theme(legend.position = "none")
CA_p_main_CS <- CA_p_main_CS + theme(legend.position = "right")

grid.arrange(CA_p_main_CN, CA_p_main_NS, CA_p_main_CS, ncol = 3, widths = c(1, 1, 1.1))


### PLOT AFRICA

# Count points per locality
africa_counts <- africa %>%
  group_by(Locality) %>%
  mutate(n_points = n())

AF_p_main_C <- ggplot(africa_counts, aes(x = Locality, y = d13c.cor)) +
  geom_boxplot(
    data = subset(africa_counts, n_points >= 3),
    fill = "cyan",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.2,
    shape = 21,
    color = "black",
    fill = "cyan",
    size = 2,
    alpha = 0.8
  ) +
  theme_minimal() +
  labs(
    title = expression(paste(delta^{13}, "C per Locality")),
    y = expression(delta^{13}~"C"),
    x = "Locality"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

AF_p_main_N <- ggplot(africa_counts, aes(x = Locality, y = dN15)) +
  geom_boxplot(
    data = subset(africa_counts, n_points >= 3),
    fill = "cyan",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.2,
    shape = 21,
    color = "black",
    fill = "cyan",
    size = 2,
    alpha = 0.8
  ) +
  theme_minimal() +
  labs(
    title = expression(paste(delta^{15}, "N per Locality")),
    y = expression(delta^{15}~"N"),
    x = "Locality"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

AF_p_main_S <- ggplot(africa_counts, aes(x = Locality, y = dS34)) +
  geom_boxplot(
    data = subset(africa_counts, n_points >= 3),
    fill = "cyan",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.2,
    shape = 21,
    color = "black",
    fill = "cyan",
    size = 2,
    alpha = 0.8
  ) +
  theme_minimal() +
  labs(
    title = expression(paste(delta^{34}, "S per Locality")),
    y = expression(delta^{34}~"S"),
    x = "Locality"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

grid.arrange(AF_p_main_C, AF_p_main_N, AF_p_main_S,
             ncol = 3, widths = c(1, 1, 1))

AF_p_main_CN <- ggplot(africa, aes(x = d13c.cor, y = dN15, shape = Locality)) +
  geom_jitter(width = 0.2, height = 0.2, size = 3) +
  theme_minimal() +
  labs(
    x = expression(delta^{13}~"C"),
    y = expression(delta^{15}~"N"),
    title = expression(paste(delta^{15}, "N/", delta^{13}, "C by Locality"))
  ) +
  scale_color_manual(adjustcolor("cyan", alpha.f = 0.75))

AF_p_main_NS <- ggplot(africa, aes(x = dS34, y = dN15, shape = Locality)) +
  geom_jitter(width = 0.2, height = 0.2, size = 3) +
  theme_minimal() +
  labs(
    x = expression(delta^{34}~"S"),
    y = expression(delta^{15}~"N"),
    title = expression(paste(delta^{15}, "N/", delta^{34}, "S by Locality"))
  ) +
  scale_color_manual(adjustcolor("cyan", alpha.f = 0.75))

AF_p_main_CS <- ggplot(africa, aes(x = d13c.cor, y = dS34, shape = Locality)) +
  geom_jitter(width = 0.2, height = 0.2, size = 3) +
  theme_minimal() +
  labs(
    x = expression(delta^{13}~"C"),
    y = expression(delta^{34}~"S"),
    title = expression(paste(delta^{34}, "S/", delta^{13}, "C by Locality"))
  ) +
  scale_color_manual(adjustcolor("cyan", alpha.f = 0.75))

AF_p_main_CN <- AF_p_main_CN + theme(legend.position = "none")
AF_p_main_NS <- AF_p_main_NS + theme(legend.position = "none")
AF_p_main_CS <- AF_p_main_CS + theme(legend.position = "right")

grid.arrange(AF_p_main_CN, AF_p_main_NS, AF_p_main_CS, ncol = 3, widths = c(1, 1, 1.5))

###########################################################
################ COUNTS AND STATISTICS ####################
###########################################################

### My corrected carbon values ###
mean(updated_manatee_data$d13c.cor)
# -18.16128
sd(updated_manatee_data$d13c.cor)
# 4.614494
mean(corrected_amazon$d13c.cor)
# -21.177
sd(corrected_amazon$d13c.cor)
# 1.889216
mean(corrected_africa_marine$d13c.cor)
# -15.96811
sd(corrected_africa_marine$d13c.cor)
# 4.196999
mean(corrected_caribbean_marine$d13c.cor)
# -20.2173
sd(corrected_caribbean_marine$d13c.cor)
# 4.824344

### My nitrogen values ###
mean(updated_manatee_data$dN15)
# 8.258472
sd(updated_manatee_data$dN15)
# 1.828962
mean(amazon$dN15)
# 8.233571
sd(amazon$dN15)
# 0.6040607
mean(caribbean$dN15)
# 7.3625
sd(caribbean$dN15)
# 2.111377
mean(africa$dN15)
# 8.739211
sd(africa$dN15)
# 1.861077

### My sulphur values ###
mean(updated_manatee_data$dS34)
# 6.207959
sd(updated_manatee_data$dS34)
# 5.036675
mean(amazon$dS34)
# 3.202178
sd(amazon$dS34)
# 3.417543
mean(caribbean$dS34)
# 5.086274
sd(caribbean$dS34)
# 7.009274
mean(africa$dS34)
# 7.905713
sd(africa$dS34)
# 3.689676

###################
### DATA NORMALITY

# Shapiro-Wilk normality tests - do the data follow a normal distribution?
shap_test_N <- shapiro.test(updated_manatee_data$dN15)
# W = 0.95757, p-value = 0.1808
  ### YES
shap_test_C <- shapiro.test(updated_manatee_data$d13c.cor)
# W = 0.9103, p-value = 0.006607
  ### NO
shap_test_S <- shapiro.test(updated_manatee_data$dS34)
# W = 0.91067, p-value = 0.006769
  ### NO

# Plot...
par(mfrow = c(2, 3),          
    mar   = c(5, 4, 5, 2) + 0.1)  

## Q-Q plots with titles
qqnorm(updated_manatee_data$d13c.cor,
       main = expression(delta^{13} * C),
       cex.main = 1.5)
qqline(updated_manatee_data$d13c.cor, col = "red")

qqnorm(updated_manatee_data$dN15,
       main = expression(delta^{15} * N),
       cex.main = 1.5)
qqline(updated_manatee_data$dN15, col = "red")

qqnorm(updated_manatee_data$dS34,
       main = expression(delta^{34} * S),
       cex.main = 1.5)
qqline(updated_manatee_data$dS34, col = "red")

## Histograms
## Carbon: x-axis from -30 to 0, with label
hist(updated_manatee_data$d13c.cor,
     main = "",
     xlab = expression(delta^{13} * C),
     xlim = c(-25, -5))

## Nitrogen: own x-axis label
hist(updated_manatee_data$dN15,
     main = "",
     xlab = expression(delta^{15} * N))

## Sulfur: own x-axis label
hist(updated_manatee_data$dS34,
     main = "",
     xlab = expression(delta^{34} * S))

# MULTI-VARIATE NORMALITY TEST
mvnData <- subset(updated_manatee_data, select = c(dN15,dS34,d13c.cor))
mvn(mvnData)

#            Test Statistic p.value     Method          MVN
# Henze-Zirkler     1.208   0.001 asymptotic ✗ Not normal

# univariate_normality
# Test Variable Statistic p.value    Normality
# 1 Anderson-Darling     dN15     0.655    0.08     ✓ Normal
# 2 Anderson-Darling     dS34     0.893    0.02 ✗ Not normal
# 3 Anderson-Darling d13c.cor     1.016    0.01 ✗ Not normal

royston(mvnData)
# Test Statistic      p.value     Method
# 1 Royston  18.11464 0.0004208975 asymptotic

multivariate_diagnostic_plot(
  mvnData,
  type = c("qq", "persp", "contour"),
  tol = 1e-25,
  use_population = TRUE
)

### the data appear normal overall especially for nitrogen
### but with carbon and sulphur there are a few influential observations
### that are driving the non-normality...

# do Levene's test
leveneTest(d13c.cor~Species, data = updated_manatee_data)
leveneTest(dN15~Species, data = updated_manatee_data)
leveneTest(dS34~Species, data = updated_manatee_data)

# do royston within species
africa <- subset(africa, select = c(dN15,dS34,d13c.cor))
royston(africa)
caribbean <- subset(caribbean, select = c(dN15,dS34,d13c.cor))
royston(caribbean)
amazon <- subset(amazon, select = c(dN15,dS34,d13c.cor))
royston(amazon)

####################
### DATA CORRELATION

plot(updated_manatee_data$dN15,updated_manatee_data$d13c.cor)
abline(lm(d13c.cor ~ dN15, data = updated_manatee_data), col = "red")
cor.test(updated_manatee_data$dN15, updated_manatee_data$d13c.cor, method = "spearman")
# no correlation for N/C (ρ = 0.017, p = 0.921)

plot(updated_manatee_data$d13c.cor,updated_manatee_data$dS34)
abline(lm(dS34 ~ d13c.cor, data = updated_manatee_data), col = "red")
cor.test(updated_manatee_data$dS34, updated_manatee_data$d13c.cor, method = "spearman")
# no corelation S/C (ρ = 0.162, p = 0.343)

plot(updated_manatee_data$dN15,updated_manatee_data$dS34)
abline(lm(dS34 ~ dN15, data = updated_manatee_data), col = "red")
cor.test(updated_manatee_data$dN15, updated_manatee_data$dS34, method = "spearman")
# weak but non significant positive correlation for S/N (ρ = 0.220, p = 0.196)

#####################
### STATS PER SPECIES

###### Do nonparametric test for the 3 species per isotope:
# Nitrogen
kruskal.test(updated_manatee_data$dN15 ~ updated_manatee_data$Species)
# Kruskal-Wallis chi-squared = 2.261, df = 2, p-value = 0.3229

# Carbon
kruskal.test(updated_manatee_data$d13c.cor ~ updated_manatee_data$Species)
# Kruskal-Wallis chi-squared = 10.361, df = 2, p-value = 0.005627

# Sulphur
kruskal.test(updated_manatee_data$dS34 ~ updated_manatee_data$Species)
# Kruskal-Wallis chi-squared = 5.6573, df = 2, p-value = 0.05909

# Do post-hoc test for Carbon:
library(FSA)
dunn.test(updated_manatee_data$d13c.cor ~ updated_manatee_data$Species, method = "bonferroni")

##############################################
########### Niche size and overlap ###########
##############################################
##### Using: nicheROVER #####
library(nicheROVER)
NR_data <- updated_manatee_data
NR_data <- NR_data %>%
  select(Species, d13c.cor, dN15, dS34) %>%
  rename(
    D15N = dN15,
    D13C = d13c.cor,
    D34S = dS34)

aggregate(NR_data[2:4], NR_data[1], mean)
nsamples <- 1e3
system.time({
  fish.par <- tapply(1:nrow(NR_data), NR_data$Species,
                     function(ii) niw.post(nsamples = nsamples, X = NR_data[ii,2:4]))
})
clrs <- c("black", "red", "cyan") 

# mu1 (del15N), mu2 (del13C), and Sigma12
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
niche.par.plot(fish.par, col = clrs, plot.index = 1)
niche.par.plot(fish.par, col = clrs, plot.index = 2)
niche.par.plot(fish.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = names(fish.par), fill = clrs, cex = 0.55)

# all mu (del15N, del13C, del34S)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = FALSE)
legend("topleft", legend = names(fish.par), fill = clrs, cex = 0.55)

# all mu and Sigma
par(mar = c(4.2, 4.2, 2, 1)+.1)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = TRUE)
legend("topright", legend = names(fish.par), fill = clrs, cex = 0.6)

# 2-d projections of 10 niche regions
nsamples <- 10

# Generate niche parameters (just like in the example)
fish.par <- tapply(1:nrow(NR_data), NR_data$Species, 
                   function(ii) {niw.post(nsamples = nsamples, 
                                          X = NR_data[ii, 2:4])
                   })

# Format data for the plotting function
myNR.data <- tapply(1:nrow(NR_data), NR_data$Species, function(ii) {
  NR_data[ii, 2:4]
})
str(myNR.data) # check okay and then convert tibble to dataframe:
myNR.data <- lapply(myNR.data, as.data.frame)

niche.plot(niche.par = fish.par, niche.data = myNR.data, pfrac = .05,
           iso.names = expression(delta^{15}*N, delta^{13}*C, delta^{34}*S),
           col = clrs, xlab = expression("Isotope Ratio (per mille)"),
           legend = FALSE)

### Calculate Hypervolume niche size
# for each species, size of 95% niche region using sample variance
tapply(1:nrow(NR_data), NR_data$Species, function(ind) {
  X <- NR_data[ind,2:4] 
  Sighat <- var(X) # sample variance
  niche.size(Sigma = Sighat, alpha = .95)
})

# for each species, size of 50% niche region using sample variance
tapply(1:nrow(NR_data), NR_data$Species, function(ind) {
  X <- NR_data[ind,2:4] 
  Sighat <- var(X) # sample variance
  niche.size(Sigma = Sighat, alpha = .5)
})

### Overlap calculation
nsamples <- 500
over.stat <- overlap(fish.par, nreps = nsamples, nprob = nsamples, alpha = .95)
overlap.plot(over.stat, col = clrs, mean.cred.col = "darkblue",
             equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

# plotting the posterior distribution of niche size by locality
MY.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})
rbind(est = colMeans(MY.size),
      se = apply(MY.size, 2, sd))
transparent_clrs <- adjustcolor(clrs, alpha.f = 0.5)
boxplot(MY.size, col=transparent_clrs, pch=16, cex =0.35,
        ylab ="Niche size", xlab = "Species")

### modify the scale ###
boxplot(MY.size, col=transparent_clrs, pch=16, cex=0.35,
        ylab="Niche size", xlab="Species", yaxt="n")

# Adjust range to nearest 1000 for major ticks
y_range <- range(MY.size, na.rm=TRUE)
y_min <- floor(y_range[1] / 1000) * 1000
y_max <- ceiling(y_range[2] / 1000) * 1000

# Major ticks every 1000 with labels
major_ticks <- seq(y_min, y_max, by=1000)
axis(2, at=major_ticks, labels=major_ticks, las=1, cex.axis=0.8, tck=-0.02)

# Minor ticks every 100 without labels
minor_ticks <- seq(y_min, y_max, by=250)
minor_ticks <- setdiff(minor_ticks, major_ticks)
axis(2, at=minor_ticks, labels=FALSE, tck=-0.01)
