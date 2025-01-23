# PLOTTING TEMPERATURE
# Lauren Dykman
# Oct 30, 2023

rm(list=ls())

# INSTALLING PACKAGES

install.packages("ggplot2")

# ACCESSING LIBRARIES

library(ggplot2)

# SETTING WORKING DIRECTORY

path <- "/Users/laurendykman/Desktop/github/Dykman_etal_2025_CommonGarden"

getwd()
setwd(path)
getwd()

# IMPORTING DATA

input_file = "hornby_bottom_temp_2023-05-23.csv"
data.1 <- read.csv(paste(path, "data", input_file, sep = "/"), header=TRUE)
data.1 <- data.1[50:dim(data.1)[1],]

input_file = "hornby_bottom_temp_2023-06-09.csv"
data.2 <- read.csv(paste(path, "data", input_file, sep = "/"), header=TRUE)

input_file = "hornby_bottom_temp_2023-07-17.csv"
data.3 <- read.csv(paste(path, "data", input_file, sep = "/"), header=TRUE)

data.2$number <- data.2$number + data.1[dim(data.1)[1], 1]
data.3$number <- data.3$number + data.2[dim(data.2)[1], 1]

data <- rbind(data.1, data.2, data.3)

ggsave(paste(path, "Figures", paste0("Figure_Temperature_Bottom_Hornby_", format(Sys.Date(), "%Y-%m-%d"), ".pdf"), sep = "/"),
    ggplot(data=data, aes(x=number, y=temp)) +
    theme_bw()+
    ylab("temperature (C)") +
    ylim(5,25) +
    geom_line(linewidth=1, color ="blue")+
      geom_hline(yintercept=11.9, color="black", linewidth = 0.8) +
      geom_hline(yintercept=18.0, color="red", linewidth = 0.8) +
      theme(text=element_text(size = 22),
        aspect.ratio=0.4,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)), height=4, width=10) # Setting the theme

data.first <- data[data$number <= 5260,]
data.third <- data[data$number >= 6317,]

data <- separate(data = data, col = date, into = c("date", "time"), sep=" ")
data.over.18 <- data[data$temp >= 18 & !is.na(data$temp),]
unique(data.over.18$date)
