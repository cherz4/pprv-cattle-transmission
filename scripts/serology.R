#############################################################################################################
# Empirical and model-based evidence for a negligible role of cattle in 
# peste des petits ruminants transmission and eradication

# Catherine M. Herzog#*, Fasil Aklilu#, Demeke Sibhatu, Dereje Shegu, 
# Redeat Belaineh, Abde Aliy Mohammed, Menbere Kidane,  Claudia Schulz,
# Brian J. Willett, Sarah Cleaveland, Dalan Bailey, Andrew R. Peters, 
# Isabella M. Cattadori, Peter J. Hudson, Hagos Asgedom, Joram Buza, 
# Mesfin Sahle Forza, Tesfaye Rufael Chibssa, Solomon Gebredufe, Nick Juleff, 
# Ottar N. Bjørnstad, Michael D. Baron, Vivek Kapur*

# Code for: analyze and visualize serology data for PPRV transmission trials
# Code written by: Catherine M. Herzog, PhD MPH

#############################################################################################################


######################
#  Load libraries ----
######################

library(ggplot2)   # for plotting figures
library(tidyverse) # for piping %>%

###############################################
#  Serology Data Import, Data Management  ----
###############################################

# Load the data by pointing to the location within the repo/project that this data resides (data folder)
trial1_sero <- read.csv("data/sero/trial1_sero.csv", header = TRUE)
trial2_sero <- read.csv("data/sero/trial2_sero.csv", header = TRUE)
trial3_sero <- read.csv("data/sero/trial3_sero.csv", header = TRUE)
trial3_chall_sero <- read.csv("data/sero/trial3_challenge_sero.csv", header = TRUE)
trial4_sero <- read.csv("data/sero/trial4_sero.csv", header = TRUE)
trial5_sero <- read.csv("data/sero/trial5_sero.csv", header = TRUE)

# Making sure each variable (column in data frame) is of the correct data type (factor, numeric, character, etc)
trial1_sero$eartag <- as.factor(trial1_sero$eartag)
trial1_sero$barn <- as.factor(trial1_sero$barn)
trial2_sero$eartag <- as.factor(trial2_sero$eartag)
trial2_sero$barn <- as.factor(trial2_sero$barn)
trial3_sero$eartag <- as.factor(trial3_sero$eartag)
trial3_sero$barn <- as.factor(trial3_sero$barn)
trial3_chall_sero$eartag <- as.factor(trial3_chall_sero$eartag)
trial3_chall_sero$barn <- as.factor(trial3_chall_sero$barn)
trial4_sero$eartag <- as.factor(trial4_sero$eartag)
trial4_sero$barn <- as.factor(trial4_sero$barn)
trial5_sero$eartag <- as.factor(trial5_sero$eartag)
trial5_sero$barn <- as.factor(trial5_sero$barn)


#####################
#  Trial 1 ----
#####################

trial1_cols <- c("1" = "black", "3" = "deeppink2", "5" = "cyan4") 

jpeg("output/sero/trial1_serology_FigS4C.jpeg", width = 6, height = 5, units = "in", quality = 100, res = 600)
ggplot(trial1_sero, aes(dpi, sn, color=barn)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 28)) +
  scale_y_reverse(limits = c(120, 0)) +
  scale_color_manual(name = "", 
                     labels = c("Control", "Sheep", "Goat"),
                     values = trial1_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1.6, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1.6, y = 53, label = "Negative", colour = "red") +
  theme_minimal() +
  theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 15))
invisible(dev.off())
# 5 missing values for control 808 that died
# 3 missing values for euthanized goat 25 after dpi 14


#####################
#  Trial 2 ----
#####################

trial2_g <- subset(trial2_sero, barn != 6) # goats 
trial2_s <- subset(trial2_sero, barn != 4) # sheep 
#controls are barn 5

trial2_cols <- c("4a" = "cyan4", "4b" = "cyan3", "5b" = "black", "6a" = "deeppink2", "6b" = "rosybrown2" )


jpeg("output/sero/trial2_serology_sheep_FigS5Cs.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ss <- ggplot(trial2_s, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn + Status", 
                     labels = c("All Controls", "Inoculated Sheep", "Sentinel Sheep"),
                     values = trial2_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1.6, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1.6, y = 53, label = "Negative", colour = "red") +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
ss
invisible(dev.off())


jpeg("output/sero/trial2_serology_goats_FigS5Cg.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
sg <- ggplot(trial2_g, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn + Status", 
                     labels = c("Inoculated Goats", "Sentinel Goats", "All Controls"),
                     values = trial2_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1.6, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1.6, y = 53, label = "Negative", colour = "red") +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
sg
invisible(dev.off())
# 4 missing rows are the 4 dpis that animal 21 was dead for 17, 21, 24, 28


#####################
#  Trial 3 ----
#####################

# split by barn
trial3_b3_sero <- subset(trial3_sero, barn == 3) # Box A Inoculated cattle to sentinel goats
trial3_b4_sero <- subset(trial3_sero, barn == 4) # Box C Inoculated goats to sentinel cattle and goats (confirmatory)
trial3_b5_sero <- subset(trial3_sero, barn == 5) # Box B positive control goats
trial3_b6_sero <- subset(trial3_sero, barn == 6) # negative controls - not sampled

trial3_cols <- c("5a" = "cyan4", "5b" = "cyan3", "6b" = "black", "6c" = "black", "3a" = "purple4", "3b" = "cyan3", 
                 "4a" = "cyan4", "4b" = "cyan3", "4c" = "mediumpurple1")

jpeg("output/sero/trial3_b3_serology.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial3_b3_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn A", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial3_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())


jpeg("output/sero/trial3_b5_poscontrols_serology_FigS10A.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial3_b5_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn B", 
                     labels = c("Goats (Inoculated)", "Goats (In Contact)"),
                     values = trial3_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())
# 19 rows missing, all dead

jpeg("output/sero/trial3_b4_serology_Fig1C.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial3_b4_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accomodate 1 point
  scale_color_manual(name = "Barn C", 
                     labels = c("Goats (Inoculated)", "Goats (In Contact)", "Cattle (In Contact)"),
                     values = trial3_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())
# 5 rows, 5 dead  

jpeg("output/sero/trial3_b6_negcontrols_serology_FigS10B.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial3_b6_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accomodate 1 point
  scale_color_manual(name = "Controls", 
                     labels = c("Goats", "Cattle"),
                     values = trial3_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())


#####################
#  Trial 3 Challenge ----
#####################

trial3_chall_cols <- c("3" = "cyan4", "4" = "black", "6" = "cyan2")

jpeg("output/sero/trial3_challenge_serology_FigS9C.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial3_chall_sero, aes(dpi, sn, color=barn)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,14), breaks = c(0, 4, 7, 10, 14)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "", 
                     labels = c("Cattle sentinel goats", "Seropositive goats", "Negative Controls"),
                     values = trial3_chall_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1.7, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 2.2, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
invisible(dev.off())


#####################
#  Trial 4 ----
#####################

# split by barn
trial4_b6_sero <- subset(trial4_sero, barn == 6) # box A Inoculated cattle to sentinel goats
trial4_b1_sero <- subset(trial4_sero, barn == 1) # box C Inoculated cattle to sentinel goats
trial4_b4_sero <- subset(trial4_sero, barn == 4) # Inoculated cattle to sentinel goats
trial4_b5_sero <- subset(trial4_sero, barn == 5) # positive controls 
trial4_b3_sero <- subset(trial4_sero, barn == 3) # negative controls 
trial4_b614_sero <- subset(trial4_sero, barn == 6 | barn == 1 | barn ==4) #all experimental together

trial4_cols <- c("6a" = "purple4", "6b" = "cyan3", "1a" = "purple4", "1b" = "cyan3", "4a" = "purple4", "4b" = "cyan3",
                 "5a" = "cyan4", "5b" = "cyan3", "3a" = "cyan3")
trial4_cols_sp <- c("Bovine" = "purple4", "Caprine" = "cyan3")


jpeg("output/sero/trial4_b6_serology.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial4_b6_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn A", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())


jpeg("output/sero/trial4_b1_serology.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial4_b1_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn B", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())

jpeg("output/sero/trial4_b4_serology.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial4_b4_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accommodate 1 point
  scale_color_manual(name = "Barn C", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())
# Eartag 181 died mechanical injury by dpi 21 so missing 3 rows

jpeg("output/sero/trial4_b5_poscontrols_serology_FigS11A.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial4_b5_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accommodate 1 point
  scale_color_manual(name = "Controls", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())
# Eartag 172 Inoculated died 12/6/2020 by dpi 10, 175 sentinal and 179 Inoculated died by dpi 21  

jpeg("output/sero/trial4_b3_negcontrols_serology_FigS11B.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial4_b3_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accommodate 1 point
  scale_color_manual(name = "Controls", 
                     labels = c("Goats (Sentinel)"),
                     values = trial4_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())

# All experimental together
jpeg("output/sero/trial4_b614_serology.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial4_b614_sero, aes(dpi, sn, color=species)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accommodate 1 point
  scale_color_manual(name = "Controls", 
                     #labels = c("Goats", "Cattle"),
                     values = trial4_cols_sp) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())
# Missing 3 rows due to goat 181



#####################
#  Trial 5 ----
#####################

trial5_cols <- c("1a" = "purple4", "2a" = "purple4", "5a" = "purple4", "6a" = "purple4",
                 "1b" = "cyan3", "2b" = "cyan3", "5b" = "cyan3", "6b" = "cyan3",
                 "3a" = "cyan4", "4b" = "cyan3") 

trial5_cols_sp <- c("Bovine" = "purple4", "Caprine" = "cyan3")

# split by barn
trial5_b1_sero <- subset(trial5_sero, barn == 1) # box D Inoculated cattle to sentinel goats
trial5_b2_sero <- subset(trial5_sero, barn == 2) # box C Inoculated cattle to sentinel goats
trial5_b5_sero <- subset(trial5_sero, barn == 5) # box A Inoculated cattle to sentinel goats
trial5_b6_sero <- subset(trial5_sero, barn == 6) # box B Inoculated cattle to sentinel goats
trial5_b3_sero <- subset(trial5_sero, barn == 3) # positive controls Box E
trial5_b4_sero <- subset(trial5_sero, barn == 4) # negative controls Box F
trial5_b1256_sero <- subset(trial5_sero, barn == 1 | barn == 2 | barn == 5 | barn == 6) # all experimental together

# Note although no sampling conducted on 24 dpi, this was in earlier trials to kept this on x axis for merging later.

jpeg("output/sero/trial5_b5_serology.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial5_b5_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn A", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())
#Warning removed 1 row - the insane SN value for goat 232 on dpi 4 due to od2 being very high. Demeke is retesting as of 4/7

jpeg("output/sero/trial5_b6_serology.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial5_b6_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn B", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())

jpeg("output/sero/trial5_b2_serology.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial5_b2_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accomodate 1 point
  scale_color_manual(name = "Barn C", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())

jpeg("output/sero/trial5_b1_serology.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial5_b1_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accomodate 1 point
  scale_color_manual(name = "Barn C", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())

jpeg("output/sero/trial5_b3_poscontrols_serology_FigS12A.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial5_b3_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accomodate 1 point
  scale_color_manual(name = "Controls", 
                     labels = c("Goats (Inoculated)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())
# Eartag died by dpi 8 so set value to blank in csv manually (so 7 sampling dates removed)

jpeg("output/sero/trial5_b4_negcontrols_serology_FigS12B.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial5_b4_sero, aes(dpi, sn, color=barni)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accomodate 1 point
  scale_color_manual(name = "Controls", 
                     labels = c("Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())
# One goat positive on dpi 10 and will have Demeke retest this as of 4/7/21

# All experimental together
jpeg("output/sero/trial5_b1256_serology.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(trial5_b1256_sero, aes(dpi, sn, color=species)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accomodate 1 point
  scale_color_manual(name = "Controls", 
                     #labels = c("Goats", "Cattle"),
                     values = trial5_cols_sp) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none")
invisible(dev.off())
#Warning removed 1 row - the insane SN value for goat 232 on dpi 4 due to od2 being very high. Demeke is retesting as of 4/7


#################################
#  Combined Trials 3-5 ----
#################################

#remove controls + any other barns put datasets together
trial3_sero_nc <- subset(trial3_sero, barn == 3) 
trial4_sero_nc <- subset(trial4_sero, barn != 5 & barn != 3) 
trial5_sero_nc <- subset(trial5_sero, barn != 3 & barn != 4) 

jointdat_sero <- rbind(trial3_sero_nc, trial4_sero_nc, trial5_sero_nc)
jointdat_sero <- jointdat_sero %>% mutate(barniv2 = case_when(barni == "3a" ~ "Ci1",
                                                              barni == "1a" ~ "Ci1",
                                                              barni == "4a" ~ "Ci1",
                                                              barni == "6a" ~ "Ci1",
                                                              barni == "2a" ~ "Ci1",
                                                              barni == "5a" ~ "Ci1",
                                                              barni == "3b" ~ "Gi0",
                                                              barni == "1b" ~ "Gi0",
                                                              barni == "4b" ~ "Gi0",
                                                              barni == "6b" ~ "Gi0",
                                                              barni == "2b" ~ "Gi0",
                                                              barni == "5b" ~ "Gi0"))
# Barns 3 (A) C to G experimental
# Barns 1, 4, 6 experimental
# Barns 1 (D), 2 (C), 5 (A), 6 (B) experimental
comb345_cols <- c("Ci1" = "purple4", "Gi0" = "cyan3")

# Manuscript Figure 2c ----
jpeg("output/sero/trials345_serology_Fig2C.jpeg", width = 8, height = 7, units = "in", quality = 100, res = 600)
ggplot(jointdat_sero, aes(dpi, sn, color=barniv2)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = comb345_cols) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
invisible(dev.off())



