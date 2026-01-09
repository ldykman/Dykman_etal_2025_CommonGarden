# BLADE AND STIPE GROWTH ON SEEDED LINES

# This script runs statistical analyses for Hornby Island Seeded line bull kelp blade growth.
# Lauren Dykman
# April 15, 2024

rm(list=ls())

# INSTALLING PACKAGES

#install.packages("tidyverse")
#install.packages("lme4")
#install.packages("ggfortify")
#install.packages("dplyr")
#install.packages("multcomp")
#install.packages("DHARMa")
#install.packages("pacman")
#pacman::p_load(Matrix)
#install.packages("TMB")
#install.packages("glmmTMB")
#install.packages("lmerTest")
#install.packages("rstatix")

library(tidyverse)
library(lme4)
library(ggfortify)
library(dplyr)
library(multcomp)
library(DHARMa)
library(glmmTMB)
library(TMB)
library(lmerTest)
library(rstatix)

# SETTING WORKING DIRECTORY

path <- "/Users/laurendykman/Desktop/github/Dykman_etal_2025_CommonGarden"

getwd()
setwd(path)
getwd()

# IMPORTING DATA

input_file = "hornby_blade_growth.csv"
data <- read.csv(paste(path, "data", input_file, sep = "/"), header=TRUE)

# CHANGING DATA TYPES

data$hole_punch_date <- as.Date(data$hole_punch_date)
data$hole_punch_measurement_date <- as.Date(data$hole_punch_measurement_date)
data$days <- as.numeric(data$hole_punch_measurement_date - data$hole_punch_date)

data <- data[!(data$line_number == 11 & data$kelp_individual == 2),] # Remove kelp 2 in line 11 from analyses because stipe was ripped in half

data.hole.punch <- gather(data, blade_position, hole_punch_distance_final_cm, hole_punch_distance_final_outer1_cm:hole_punch_distance_final_inner2_cm)

# RENAMING SOME VARIABLES

data.hole.punch$growth_total_cm <- data.hole.punch$hole_punch_distance_final_cm - data.hole.punch$hole_punch_distance_initial_cm

data.hole.punch$growth_per_day_cm <- data.hole.punch$growth_total_cm/data.hole.punch$days
data.hole.punch$growth_per_day_by_stipe_length_cm <- data.hole.punch$growth_per_day_cm/data.hole.punch$stipe_length_initial_cm
data.hole.punch$growth_stipe_cm <- data.hole.punch$stipe_length_final_cm - data.hole.punch$stipe_length_initial_cm
data.hole.punch$growth_stipe_per_day_cm <- data.hole.punch$growth_stipe_cm/data.hole.punch$days
data.hole.punch$growth_per_day_cm[data.hole.punch$growth_per_day_cm < 0] <- 0 # Setting some slight negative numbers to zero, since negatives were likely due to shrinkage or measurement error, indicates no growth

data.hole.punch$growth_per_day_cm_adjusted <- log(data.hole.punch$growth_per_day_cm + 1)

data.hole.punch <- data.hole.punch[!(is.na(data.hole.punch$hole_punch_date) | is.na(data.hole.punch$growth_total_cm)),]

ggplot(data.hole.punch, aes(x=stipe_length_initial_cm, y=growth_per_day_cm)) + geom_point()

# CREATING BARPLOT OF DAILY BLADE GROWTH

growth.mean <- aggregate(data.hole.punch$growth_per_day_cm, by = list(data.hole.punch$population), FUN=mean, na.rm = TRUE)
colnames(growth.mean) <- c("population", "growth_mean")

growth.sd <- aggregate(data.hole.punch$growth_per_day_cm, by = list(data.hole.punch$population), FUN=sd, na.rm = TRUE)
colnames(growth.sd) <- c("population", "growth_sd")

data.hole.punch$combine <- paste(data.hole.punch$line_number, data.hole.punch$kelp_individual, sep=".")

counts <- data.hole.punch %>%
  group_by(population, combine) %>%
  tally()

counts <- counts %>%
  group_by(population) %>%
  tally()

length2 <- merge(growth.mean, growth.sd, by = "population")

length3 <- merge(length2, counts, by = "population")
colnames(length3) <- c("population", "mean", "SD", "n")

length3$sqrtN <- sqrt(length3$n)
length3$SDerr <- length3$SD/length3$sqrtN

length3$population <- factor(length3$population, levels = c("Oyster River", "Gabriola", "Brown's Bay", "Vancouver", "Victoria"))

ggsave(paste(path, "Figures", paste0("Figure_Kelp_BladeGrowth_Hornby_", format(Sys.Date(), "%Y-%m-%d"), ".pdf"), sep = "/"),
  ggplot(length3, aes(fill = population, x=population, y=mean)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black", show.legend = FALSE, width = 0.6) +
    geom_errorbar(aes(x = population,
                  ymin = ifelse(mean - SDerr <= 0, 0, mean - SDerr),
                  ymax = mean + SDerr), linewidth = 0.5, width = 0.4, stat = "identity", position = position_dodge(0.9))+
    xlab("population") +
    ylab("blade elongation (cm/day)") +
    theme_bw() +
    ylim(0,0.8) +
    theme(text=element_text(size = 10),
          aspect.ratio=0.6,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA)), height=3, width=4) # Setting the theme

ggsave(paste(path, "Figures", paste0("Figure_Kelp_BladeGrowth_Hornby_Histogram_", format(Sys.Date(), "%Y-%m-%d"), ".pdf"), sep = "/"),
ggplot(data.hole.punch, aes(x = growth_per_day_cm_adjusted)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  facet_wrap(~population, ncol = 1) + # Creates a separate plot for each factor level, arranged in a single column
  labs(title = "Histogram of Blade Elongation",
       x = "Value",
       y = "Frequency") +
  theme_minimal(), height=5, width=3) # Setting the theme

# After normalizing data, exploring a linear mixed effects model with normal distribution
model.fit.1 <- lmer(growth_per_day_cm_adjusted ~ population + stipe_length_initial_cm + (1|line_number/kelp_individual), data = data.hole.punch)
summary(model.fit.1) # Interaction removed cause not significant

# Testing normality of groups with Shapiro Wilk test
results <- data.hole.punch %>%
  group_by(population) %>%
  shapiro_test(growth_per_day_cm_adjusted) # Shapiro wilk test shows that all are normally distributed after transformation except Oyster River

par(mar = c(1,1,1,1))

# Testing for homogeneity of variance and within-group uniformity
simulationOutput <- simulateResiduals(fittedModel = model.fit.1, plot = F, re.form = NULL)
plot(simulationOutput) # No significant values found

testUniformity(fittedModel = model.fit.1, plot = F, re.form = NULL)

comp.test <- glht(model.fit.1, linfct = mcp(population="Tukey"))
summary(comp.test)

# Testing generalized linear mixed effects model with Gamma distribution
glm.results <- glmer(growth_per_day_cm + 0.4 ~ population + stipe_length_initial_cm + (1|line_number/kelp_individual), family = Gamma, data = data.hole.punch) 
summary(glm.results)

# Testing for homogeneity of variance and within-group uniformity
simulationOutput <- simulateResiduals(fittedModel = glm.results, plot = F, re.form = NULL)
plot(simulationOutput) # No significant values found

# Calculating percent difference

bb <- length3[length3$population == "Brown's Bay",]$mean
or <- length3[length3$population == "Oyster River",]$mean
vic <- length3[length3$population == "Victoria",]$mean

percent.bb.or <- (or - bb)/bb
percent.bb.vic <- (or - vic)/vic

# Testing with a Tweedie distribution

glm.results = glmmTMB(growth_per_day_cm_adjusted ~ population + (1|line_number/kelp_individual), family=tweedie, data=data.hole.punch)
summary(glm.results)

# Testing for homogeneity of variance and within-group uniformity
simulationOutput <- simulateResiduals(fittedModel = glm.results, plot = F, re.form = NULL)
plot(simulationOutput) # No significant values found
