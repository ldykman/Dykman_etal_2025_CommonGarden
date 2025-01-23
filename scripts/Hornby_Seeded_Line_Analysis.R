# ANALYZING NEREOCYSTIS ON SEEDED LINES

# This script runs statistical analyses for Hornby Island Seeded line kelp lengths.
# Lauren Dykman
# Created: May 8, 2023
# Modified: April 23, 2024

rm(list=ls())

# INSTALLING PACKAGES

#install.packages("tidyverse")
#install.packages("lme4")
#install.packages("multcomp")

library(tidyverse)
library(lme4)
library(multcomp)

# SETTING WORKING DIRECTORY

path <- "/Users/laurendykman/Desktop/github/Dykman_etal_2025_EcologicalApplications"

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

data.counts$stipe_with_float_cm <- data.counts$stipe_with_float_cm + 0.001

data.length.first <- data.counts[data.counts$days_old == "70",]
glm.results.first <- glmer(stipe_with_float_cm ~ population + (1|plot_no), data.length.first, family = Gamma)
summary(glm.results.first)

comp.test <- glht(glm.results.first, mcp(population="Tukey"))
summary(comp.test)

data.length.second <- data.counts[data.counts$days_old == "120",]
glm.results.second <- glmer(stipe_with_float_cm ~ population + (1|plot_no), data.length.second, family = Gamma)
summary(glm.results.second)

comp.test <- glht(glm.results.second, mcp(population="Tukey"))
summary(comp.test)

data.length.third <- data.counts[data.counts$days_old == "140",]
glm.results.third <- glmer(stipe_with_float_cm/10 ~ population + (1|plot_no), data.length.third, family = Gamma) # Had to divide by 10 because dispersion was too high to run.
summary(glm.results.third)

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
