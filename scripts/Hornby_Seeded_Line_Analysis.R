# ANALYZING NEREOCYSTIS ON SEEDED LINES

# This script runs statistical analyses for Hornby Island Seeded line kelp lengths.
# Lauren Dykman
# Created: May 8, 2023
# Modified: Jan 9, 2026

rm(list=ls())

# INSTALLING PACKAGES

#install.packages("tidyverse")
#install.packages("lme4")
#install.packages("multcomp")
#install.packages("DHARMa")
#install.packages("glmmTMB")
#install.packages("TMB")
#install.packages("lmerTest")
#install.packages("rstatix")

library(tidyverse)
library(lme4)
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

input_file = "hornby_lines_measurements.csv"
data <- read.csv(paste(path, "data", input_file, sep = "/"), header=TRUE)

# Merging consecutive dates for analysis

data$date[data$date == "2023-05-05"] <- "2023-05-04"
data$date[data$date == "2023-05-24"] <- "2023-05-23"

# Formatting date data

data$date <- as.Date(data$date)
data$days_old <- as.character(data$date - as.Date("2023-01-03")) # Calculating dates since outplanting
data$days_old[data$days_old == "73"] <- "70" # Making some dates nicer for plot
data$days_old[data$days_old == "121"] <- "120"
data$days_old <- factor(data$days_old, levels = c("70", "120", "140")) # Check unique(data$days_old) and write the days as factors here
data$stipe_with_float_cm <- data$stipe_with_float/10 # Making length into cm rather than mm
data$line_ind <- paste(data$plot_no, data$plant)

data$population_days <- paste(data$population, data$days_old, sep = "_")

length.mean <- aggregate(data$stipe_with_float_cm, by = list(data$population_days, data$population, data$days_old), FUN=mean, na.rm = TRUE)
colnames(length.mean) <- c("population_days", "population", "days_old", "stipe_length_mean")

length.sd <- aggregate(data$stipe_with_float_cm, by = list(data$population_days, data$population, data$days_old), FUN=sd, na.rm = TRUE)
colnames(length.sd) <- c("population_days", "population", "days_old", "stipe_length_sd")

data.counts <- data[!is.na(data$stipe_with_float_cm),] # Getting counts of measured kelp by removing individuals that don't have stipe measurements

counts <- data.counts %>%
  group_by(population_days, population, days_old, line_ind) %>%
  tally()

counts.table <- aggregate(counts$n, by = list(counts$population_days), FUN=sum)
colnames(counts.table) <- c("population_days", "n")

length2 <- merge(length.mean, length.sd[c("population_days", "stipe_length_sd")], by = "population_days")

length3 <- merge(length2, counts.table, by = "population_days")

length3$STD_ERR <- length3$stipe_length_sd/sqrt(length3$n)

length3$population <- factor(length3$population, levels = c("Oyster River", "Gabriola", "Brown's Bay", "Vancouver", "Victoria"))

ggsave(paste(path, "Figures", paste0("Figure_Kelp_Length_SeededLines_Barplot_", format(Sys.Date(), "%Y-%m-%d"), ".pdf"), sep = "/"),
       ggplot(length3, aes(fill = population, x=days_old, y=stipe_length_mean)) +
         geom_bar(stat = "identity", position = position_dodge(), color = "black") +
         geom_errorbar(aes(x = days_old,
                           ymin = ifelse(stipe_length_mean - STD_ERR <= 0, 0, stipe_length_mean - STD_ERR),
                           ymax = stipe_length_mean + STD_ERR), linewidth = 0.5, width = 0.4, stat = "identity", position = position_dodge(0.9))+ #, size = 0.2, width=0.2, alpha=0.6) +# Setting the theme
         xlab("days since outplanting") +
         ylab("stipe length (cm)") +
         theme_bw() +
         labs(fill = "population") +
         labs(color = "population") +
         ylim(0,170) +
         theme(text=element_text(size = 20),
               aspect.ratio=0.6,
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"),
               legend.title = element_text(size=16),
               legend.text = element_text(size=14),
               panel.background = element_rect(fill = "transparent", colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA)), height=6, width=8) # Setting the theme

# Shifting data slightly to avoid zeros

data.counts$stipe_with_float_cm_adjusted <- log(data.counts$stipe_with_float_cm + 1)

# The first timepoint

data.length.first <- data.counts[data.counts$days_old == "70",]
glm.results.first <- lmer(stipe_with_float_cm_adjusted ~ population + (1|plot_no), data.length.first)
summary(glm.results.first)

results <- data.length.first %>%
  group_by(population) %>%
  shapiro_test(stipe_with_float_cm_adjusted) # 

ggsave(paste(path, "Figures", paste0("Figure_Kelp_Length_Hornby_Histogram_Time1_", format(Sys.Date(), "%Y-%m-%d"), ".pdf"), sep = "/"),
       ggplot(data.length.first, aes(x = stipe_with_float_cm_adjusted)) +
         geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
         facet_wrap(~population, ncol = 1) + # Creates a separate plot for each factor level, arranged in a single column
         labs(title = "Histogram of Stipe Length",
              x = "Value",
              y = "Frequency") +
         theme_minimal(), height=5, width=3) # Setting the theme

simulationOutput <- simulateResiduals(fittedModel = glm.results.first, plot = F, re.form = NULL)
plot(simulationOutput) # No significant values found

comp.test <- glht(glm.results.first, mcp(population="Tukey"))
summary(comp.test)

# The second timepoint

data.length.second <- data.counts[data.counts$days_old == "120",]
glm.results.second <- lmer(stipe_with_float_cm_adjusted ~ population + (1|plot_no), data.length.second)
summary(glm.results.second)

results <- data.length.second %>%
  group_by(population) %>%
  shapiro_test(stipe_with_float_cm_adjusted) # 

ggsave(paste(path, "Figures", paste0("Figure_Kelp_Length_Hornby_Histogram_Time2_", format(Sys.Date(), "%Y-%m-%d"), ".pdf"), sep = "/"),
       ggplot(data.length.second, aes(x = stipe_with_float_cm)) +
         geom_histogram(binwidth = 20, fill = "skyblue", color = "black") +
         facet_wrap(~population, ncol = 1) + # Creates a separate plot for each factor level, arranged in a single column
         labs(title = "Histogram of Stipe Length",
              x = "Value",
              y = "Frequency") +
         theme_minimal(), height=5, width=3) # Setting the theme

simulationOutput <- simulateResiduals(fittedModel = glm.results.second, plot = F, re.form = NULL)
plot(simulationOutput) # No significant values found

comp.test <- glht(glm.results.second, mcp(population="Tukey"))
summary(comp.test)

# The third timepoint

data.length.third <- data.counts[data.counts$days_old == "140",]
glm.results.third <- lmer(stipe_with_float_cm_adjusted ~ population + (1|plot_no), data.length.third) # Had to divide by 10 because dispersion was too high to run.
summary(glm.results.third)

results <- data.length.third %>%
  group_by(population) %>%
  shapiro_test(stipe_with_float_cm_adjusted) # 

ggsave(paste(path, "Figures", paste0("Figure_Kelp_Length_Hornby_Histogram_Time3_", format(Sys.Date(), "%Y-%m-%d"), ".pdf"), sep = "/"),
       ggplot(data.length.third, aes(x = stipe_with_float_cm_adjusted)) +
         geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
         facet_wrap(~population, ncol = 1) + # Creates a separate plot for each factor level, arranged in a single column
         labs(title = "Histogram of Stipe Length",
              x = "Value",
              y = "Frequency") +
         theme_minimal(), height=5, width=3) # Setting the theme

simulationOutput <- simulateResiduals(fittedModel = glm.results.third, plot = F, re.form = NULL)
plot(simulationOutput) # No significant values found

comp.test <- glht(glm.results.third, mcp(population="Tukey"))
summary(comp.test)

# CALCULATING PERCENTAGE DIFFERENCE

# 120 days Brown's Bay vs Oyster River and Victoria

bb.120 <- length3[length3$population_days == "Brown's Bay_120",]$stipe_length_mean
or.120 <- length3[length3$population_days == "Oyster River_120",]$stipe_length_mean
vic.120 <- length3[length3$population_days == "Victoria_120",]$stipe_length_mean

percent.bb.or.120 <- (bb.120 - or.120)/or.120
percent.bb.vic.120 <- (bb.120 - vic.120)/vic.120

# 140 days Brown's Bay vs Oyster River and Victoria

bb.140 <- length3[length3$population_days == "Brown's Bay_140",]$stipe_length_mean
or.140 <- length3[length3$population_days == "Oyster River_140",]$stipe_length_mean
vic.140 <- length3[length3$population_days == "Victoria_140",]$stipe_length_mean

percent.bb.or.140 <- (bb.140 - or.140)/or.140
percent.bb.vic.140 <- (bb.140 - vic.140)/vic.140
