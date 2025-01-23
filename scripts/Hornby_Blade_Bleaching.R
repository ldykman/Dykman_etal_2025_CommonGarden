# BLADE BLEACHING AND AREA ANALYSIS

# Lauren Dykman and Alex Wiebe
# Created: February, 2024

rm(list=ls())

# INSTALLING PACKAGES

#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("lme4")
#install.packages("multcomp")

# ACCESSING LIBRARIES

library(ggplot2)
library(dplyr)
library(lme4)
library(multcomp)

# SETTING WORKING DIRECTORY

path <- "/Users/laurendykman/Desktop/github/Dykman_etal_2025_EcologicalApplications"

getwd()
setwd(path)
getwd()

# IMPORTING DATA

input_file = "hornby_blade_bleaching.csv"
data <- read.csv(paste(path, "data", input_file, sep = "/"), header=TRUE)

data$kelp_individual <- paste(data$population, data$line, data$individual, sep = ".")

# Create subset removing bad variables (remove NA)

data.simple <- data %>%
  filter(blade_location == "hole punch", !is.na(weight_of_hole_punch_mg),!is.na(mean_grey_scale)) %>%
  group_by(population)

# Rescale by the maximum and minimum darkness value that are possible for the kelp

# blank background mean = 166.9, bryozoan mean = 27.8
# 139.1 scale mean

# blank background max = 173.8, bryozoan min = 12.0
# 161.8 scale mean

# blank background max = 173.8, kelp min = 28.1
# 145.7 scale mean

# blank background mean = 166.9, kelp mean = 39.2
# 127.7 scale mean
Lmin = 39.2
Lmax = 166.9

data.simple$mean_darkness_normalized <- (1-(data.simple$mean_grey_scale - Lmin)/(Lmax - Lmin))/data.simple$weight_of_hole_punch_mg

# Count number of data points for each location (to use in standard error)
counts <- data.simple %>%
  group_by(population) %>%
  tally()

# Calculate the mean and standard deviation and standard error of the normalized mean greyscale for each population
data_msd <- data.simple %>%                          
  group_by(population) %>%
  summarise_at(vars(mean_darkness_normalized),
               list(mean = ~mean(as.numeric(.), na.rm = TRUE),
                    sd = ~sd(as.numeric(.), na.rm = TRUE))) 

# Merge the counts of population data with standard deviation and standard error for each population
data.error <- merge(counts, data_msd, by = "population")
data.error$se <-data.error$sd/sqrt(data.error$n)

# Order the population so plot shows locations closest to Hornby on left and furthest on right
data.error$population <- factor(data.error$population, levels = c("Oyster River", "Gabriola", "Browns Bay", "Vancouver", "Victoria"))

ggsave(paste(path, "Figures", paste0("Figure_Kelp_BladeBleaching_Hornby_Darkness_", format(Sys.Date(), "%Y-%m-%d"), ".pdf"), sep = "/"),
       ggplot(data.error, aes(fill = population, x=population, y=mean)) +
         geom_bar(stat = "identity", position = position_dodge(), color = "black", show.legend = FALSE, width = 0.6) +
         geom_errorbar(aes(x = population,
                           ymin = mean - se,
                           ymax = mean + se), linewidth = 0.5, width = 0.4, stat = "identity", position = position_dodge(0.9))+
         xlab("population") +
         ylab("darkness normalized by weight") +
         theme_bw() +
         theme(text=element_text(size = 10),
               aspect.ratio=0.6,
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"),
               panel.background = element_rect(fill = "transparent", colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA)), height=3, width=4) # Setting the theme

# Generalized linear mixed effect model for significance in darkness, normalized by hole punch

glm.results <- glmer(mean_darkness_normalized ~ population + (1|line), data=data.simple, family=Gamma)
summary(glm.results) # Nothing significant except almost Victoria

# PLOTTING AVERAGE BLADE WEIGHT

# Calculate the mean and standard deviation and standard error of hole punch weight for each population
data_msd <- data.simple %>%                          
  group_by(population) %>%
  summarise_at(vars(weight_of_hole_punch_mg),
               list(mean = ~mean(as.numeric(.), na.rm = TRUE),
                    sd = ~sd(as.numeric(.), na.rm = TRUE))) 

# Merge the counts of population data with standard deviation and standard error for each population
data.error <- merge(counts, data_msd, by = "population")
data.error$se <-data.error$sd/sqrt(data.error$n)

# Order the population so plot shows locations closest to Hornby on left and furthest on right
data.error$population <- factor(data.error$population, levels = c("Oyster River", "Gabriola", "Browns Bay", "Vancouver", "Victoria"))

ggsave(paste(path, "Figures", paste0("Figure_Kelp_BladeBleaching_Hornby_HolePunchWeight_", format(Sys.Date(), "%Y-%m-%d"), ".pdf"), sep = "/"),
       ggplot(data.error, aes(fill = population, x=population, y=mean)) +
         geom_bar(stat = "identity", position = position_dodge(), color = "black", show.legend = FALSE, width = 0.6) +
         geom_errorbar(aes(x = population,
                           ymin = ifelse(mean - se <= 0, 0, mean - se),
                           ymax = mean + se), linewidth = 0.5, width = 0.4, stat = "identity", position = position_dodge(0.9))+
         xlab("population") +
         ylab("hole punch weight (mg)") +
         theme_bw() +
         theme(text=element_text(size = 10),
               aspect.ratio=0.6,
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"),
               panel.background = element_rect(fill = "transparent", colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA)), height=3, width=4) # Setting the theme

# Generalized linear mixed effect model for significance in hole punch weight
glm.results <- glmer(weight_of_hole_punch_mg ~ population + (1|line), data=data.simple, family=Gamma)
summary(glm.results)
# Nothing significant

# ANALYZING TOTAL BLADE AREA

input_file = "hornby_blade_area.csv"
data.area <- read.csv(paste(path, "data", input_file, sep = "/"), header=TRUE)

# Count number of data points for each location (to use in standard error)
counts <- data.area %>%
  group_by(population) %>%
  tally()

# Calculate the mean and standard deviation and standard error of kelp blade area for each population
data_msd <- data.area %>%                          
  group_by(population) %>%
  summarise_at(vars(area_cm2),
               list(mean = ~mean(as.numeric(.), na.rm = TRUE),
                    sd = ~sd(as.numeric(.), na.rm = TRUE))) 

# Merge the counts of population data with standard deviation and standard error for each population
data.error <- merge(counts, data_msd, by = "population")
data.error$se <-data.error$sd/sqrt(data.error$n)

# Order the population so plot shows locations closest to Hornby on left and furthest on right
data.error$population <- factor(data.error$population, levels = c("Oyster River", "Gabriola", "Browns Bay", "Vancouver", "Victoria"))

ggsave(paste(path, "Figures", paste0("Figure_Kelp_BladeBleaching_Hornby_BladeArea_", format(Sys.Date(), "%Y-%m-%d"), ".pdf"), sep = "/"),
       ggplot(data.error, aes(fill = population, x=population, y=mean)) +
         geom_bar(stat = "identity", position = position_dodge(), color = "black", show.legend = FALSE, width = 0.6) +
         geom_errorbar(aes(x = population,
                           ymin = ifelse(mean - se <= 0, 0, mean - se),
                           ymax = mean + se), linewidth = 0.5, width = 0.4, stat = "identity", position = position_dodge(0.9))+
         labs(x = "population", y = bquote("blade area "(cm^2))) +
         theme_bw() +
         ylim(0,300) +
         theme(text=element_text(size = 10),
               aspect.ratio=0.6,
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"),
               panel.background = element_rect(fill = "transparent", colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA)), height=3, width=4) # Setting the theme

# Generalized linear mixed effect model for significance in hole punch weight
glm.results <- glmer(area_cm2/10 ~ population + (1|line), data=data.area, family=Gamma)
summary(glm.results)

# Oyster river significant
comp.test <- glht(glm.results, linfct = mcp(population="Tukey"))
summary(comp.test)
