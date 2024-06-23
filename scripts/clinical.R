#############################################################################################################
# Empirical and model-based evidence for a negligible role of cattle in 
# peste des petits ruminants transmission and eradication

# Catherine M. Herzog#*, Fasil Aklilu#, Demeke Sibhatu, Dereje Shegu, 
# Redeat Belaineh, Abde Aliy Mohammed, Menbere Kidane,  Claudia Schulz,
# Brian J. Willett, Sarah Cleaveland, Dalan Bailey, Andrew R. Peters, 
# Isabella M. Cattadori, Peter J. Hudson, Hagos Asgedom, Joram Buza, 
# Mesfin Sahle Forza, Tesfaye Rufael Chibssa, Solomon Gebredufe, Nick Juleff, 
# Ottar N. Bjørnstad, Michael D. Baron, Vivek Kapur*

# Code for: analyze and visualize clinical data for PPRV transmission trials
# Code written by: Catherine M. Herzog, PhD MPH

#############################################################################################################

######################
#  Load libraries ----
######################

library(ggplot2)   # for plotting figures
library(tidyverse) # for piping %>%
library(gridExtra) # for arranging several plots
library(grid) # for textGrob function in plots

###############################################
#  Clinical Data Import, Data Management  ----
###############################################

# Load the data by pointing to the location within the repo/project that this data resides (data folder)
trial1_clin <- read.csv("data/clin/trial1_clin.csv", header = TRUE)
trial2_clin <- read.csv("data/clin/trial2_clin.csv", header = TRUE)
trial3_clin <- read.csv("data/clin/trial3_clin.csv", header = TRUE)
trial3_chall_clin <- read.csv("data/clin/trial3_challenge_clin.csv", header = TRUE)
trial4_clin <- read.csv("data/clin/trial4_clin.csv", header = TRUE)
trial5_clin <- read.csv("data/clin/trial5_clin.csv", header = TRUE)

# Making sure each variable (column in data frame) is of the correct data type (factor, numeric, character, etc)
trial1_clin$eartag <- as.factor(trial1_clin$eartag)
trial1_clin$barn <- as.factor(trial1_clin$barn)
trial2_clin$eartag <- as.factor(trial2_clin$eartag)
trial2_clin$barn <- as.factor(trial2_clin$barn)
trial3_clin$eartag <- as.factor(trial3_clin$eartag)
trial3_clin$barn <- as.factor(trial3_clin$barn)
trial3_chall_clin$eartag <- as.factor(trial3_chall_clin$eartag)
trial3_chall_clin$barn <- as.factor(trial3_chall_clin$barn)
trial4_clin$eartag <- as.factor(trial4_clin$eartag)
trial4_clin$barn <- as.factor(trial4_clin$barn)
trial5_clin$eartag <- as.factor(trial5_clin$eartag)
trial5_clin$barn <- as.factor(trial5_clin$barn)

# Some trials have data before dpi 0, remove these days
trial1_clin <- subset(trial1_clin, !dpi <0)
trial1_clin <- subset(trial1_clin, eartag != 808) # Took out control animal that had severe anemia and died overnight into dpi 10
trial2_clin <- subset(trial2_clin, !dpi <0)
trial3_clin <- subset(trial3_clin, !dpi <0)
trial3_chall_clin <- subset(trial3_chall_clin, !dpi <0)
trial4_clin <- subset(trial4_clin, !dpi <0)
trial5_clin <- subset(trial5_clin, !dpi <0)

colstoconvert <- c("general", "feces", "odnd", "headmucos", "resp")
trial1_clin[colstoconvert] <- sapply(trial1_clin[colstoconvert], as.integer) # check warning (goat 808 and 28 die but there are a few singular missing scores)
trial2_clin[colstoconvert] <- sapply(trial2_clin[colstoconvert], as.integer)
trial3_clin[colstoconvert] <- sapply(trial3_clin[colstoconvert], as.integer)
trial3_chall_clin[colstoconvert] <- sapply(trial3_chall_clin[colstoconvert], as.integer)
trial4_clin[colstoconvert] <- sapply(trial4_clin[colstoconvert], as.integer)
trial5_clin[colstoconvert] <- sapply(trial5_clin[colstoconvert], as.integer)

# Create clinical score without temperature (Clinical Score - Modified)
trial1_clin <- trial1_clin %>% mutate(clinscore_mod = general + feces + odnd + headmucos + resp)
trial2_clin <- trial2_clin %>% mutate(clinscore_mod = general + feces + odnd + headmucos + resp)
trial3_clin <- trial3_clin %>% mutate(clinscore_mod = general + feces + odnd + headmucos + resp)
trial3_chall_clin <- trial3_chall_clin %>% mutate(clinscore_mod = general + feces + odnd + headmucos + resp)
trial4_clin <- trial4_clin %>% mutate(clinscore_mod = general + feces + odnd + headmucos + resp)
trial5_clin <- trial5_clin %>% mutate(clinscore_mod = general + feces + odnd + headmucos + resp)


#####################
#  Trial 1 ----
#####################

trial1_cols <- c("1" = "black", "3" = "deeppink2", "5" = "cyan4") 

# Calculating Peaks ----
# Rectal Temperature 
gp <- trial1_clin[trial1_clin$barn == 5 & trial1_clin$tempc>=36 & trial1_clin$tempc<=42,]
gp <- gp[!is.na(gp$eartag),]
gp_fit <- loess(tempc ~ dpi, gp)
gp_nd <- data.frame(dpi=seq(min(gp$dpi), max(gp$dpi), length=100))
gp_nd$fit <- predict(gp_fit, newdata=gp_nd) 
gp_tmax <- gp_nd$dpi[which.max(gp_nd$fit)]
gp[gp$barn == 5,]$color <- "deeppink2"
# plot(gp$tempc~gp$dpi)
# lines(gp_nd$fit~gp_nd$dpi, col="grey")

sp <- trial1_clin[trial1_clin$barn == 3 & trial1_clin$tempc>=36 & trial1_clin$tempc<=42,]
sp_fit <- loess(tempc ~ dpi, sp)
sp_nd <- data.frame(dpi=seq(min(sp$dpi), max(sp$dpi), length=100))
sp_nd$fit <- predict(sp_fit, newdata=sp_nd) 
sp_tmax <- sp_nd$dpi[which.max(sp_nd$fit)]
# plot(sp$tempc~sp$dpi)
# lines(sp_nd$fit~sp_nd$dpi, col="grey")

# Clinical Score - without temp (to do with temp in it use autosum)
gpc <- trial1_clin[trial1_clin$barn == 5,]
gpc_fit <- loess(clinscore_mod ~ dpi, gpc)
gpc_nd <- data.frame(dpi=seq(min(gpc$dpi), max(gpc$dpi), length=100))
gpc_nd$fit <- predict(gpc_fit, newdata=gpc_nd) 
gpc_cmax <- gpc_nd$dpi[which.max(gpc_nd$fit)]
# plot(gp$clinscore_mod~gpc$dpi)
# lines(gpc_nd$fit~gpc_nd$dpi, col="grey")

spc <- trial1_clin[trial1_clin$barn == 3,]
spc_fit <- loess(clinscore_mod ~ dpi, spc)
spc_nd <- data.frame(dpi=seq(min(spc$dpi), max(spc$dpi), length=100))
spc_nd$fit <- predict(spc_fit, newdata=spc_nd) 
spc_cmax <- spc_nd$dpi[which.max(spc_nd$fit)]
# plot(spc$clinscore_mod~spc$dpi)
# lines(spc_nd$fit~spc_nd$dpi, col="grey")

# Plots ---- 
# Rectal Temperature
rt <- ggplot(trial1_clin, aes(dpi, tempc, color=barn)) +
  geom_line(aes(group = eartag), alpha=0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Control", " Inoculated Sheep", " Inoculated Goat"),
                     values = trial1_cols) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0, 28, 2)) +
  scale_y_continuous(limits = c(37,41)) +
  geom_vline(aes(xintercept = gp_tmax), color = "deeppink2") +
  geom_vline(aes(xintercept = sp_tmax), color = "cyan4") +
  annotate(geom = "text", x = (gp_tmax - 1.75), y = 41, label = paste(round(gp_tmax,1)), color = "deeppink2") +
  annotate(geom = "text", x = (sp_tmax + 1.75), y = 41, label = paste(round(sp_tmax,1)), color = "cyan4") +
  theme_minimal() +
  theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 15))
# Warning about 16 rows is dpi 13-28 for the PM animal goat 28 barn 5
# 18 rows with non-finite values is:

# Clinical Score
cs <- ggplot(trial1_clin, aes(dpi, clinscore_mod, color=barn)) +
  geom_line(aes(group = eartag), alpha=0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Control", "Inoculated Sheep", " Inoculated Goat"),
                     values = trial1_cols) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0, 28, 2)) +
  geom_vline(aes(xintercept = gpc_cmax), color = "deeppink2") +
  geom_vline(aes(xintercept = spc_cmax), color = "cyan4") +
  annotate(geom = "text", x = (gpc_cmax - 1.75), y = 10, label = paste(round(gpc_cmax,1)), color = "deeppink2") +
  annotate(geom = "text", x = (spc_cmax + 1.75), y = 10, label = paste(round(spc_cmax,1)), color = "cyan4") +
  theme_minimal() +
  theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 15))

jpeg("output/clin/trial1_clin_FigS4AB.jpeg", width = 6.5, height = 4, units = "in", quality = 100, res = 600)
grid.arrange(rt, cs, ncol=2)
invisible(dev.off())


#####################
#  Trial 2 ----
#####################

trial2_cols <- c("4a" = "cyan4", "4b" = "cyan3", "5b" = "black", "6a" = "deeppink2", "6b" = "rosybrown3" )

trial2_clin_g <- subset(trial2_clin, barn != 6) #goats are 4 
trial2_clin_s <- subset(trial2_clin, barn != 4) #sheep are 6
#controls are barn 5

# Calculating Peaks ----
# Rectal Temperature 
gpi1 <- trial2_clin_g[trial2_clin_g$barn == 4 & trial2_clin_g$innoc == 1 & trial2_clin_g$tempc>=36 & trial2_clin_g$tempc<=42,]
gpi1 <- gpi1[!is.na(gpi1$eartag),]
gpi1_fit <- loess(tempc ~ dpi, gpi1)
gpi1_nd <- data.frame(dpi=seq(min(gpi1$dpi), max(gpi1$dpi), length=100))
gpi1_nd$fit <- predict(gpi1_fit, newdata=gpi1_nd) 
gpi1_tmax <- gpi1_nd$dpi[which.max(gpi1_nd$fit)]
gpi1[gpi1$barn == 4,]$color <- "cyan4"
# plot(gpi1$tempc~gpi1$dpi)
# lines(gpi1_nd$fit~gpi1_nd$dpi, col="grey")

gpi0 <- trial2_clin_g[trial2_clin_g$barn == 4 & trial2_clin_g$innoc == 0 & trial2_clin_g$tempc>=36 & trial2_clin_g$tempc<=42,]
gpi0 <- gpi0[!is.na(gpi0$eartag),]
gpi0_fit <- loess(tempc ~ dpi, gpi0)
gpi0_nd <- data.frame(dpi=seq(min(gpi0$dpi), max(gpi0$dpi), length=100))
gpi0_nd$fit <- predict(gpi0_fit, newdata=gpi0_nd) 
gpi0_tmax <- gpi0_nd$dpi[which.max(gpi0_nd$fit)]
gpi0[gpi0$barn == 4,]$color <- "cyan4"
# plot(gpi0$tempc~gpi0$dpi)
# lines(gpi0_nd$fit~gpi0_nd$dpi, col="grey")

spi1 <- trial2_clin_s[trial2_clin_s$barn == 6 & trial2_clin_s$innoc == 1 & trial2_clin_s$tempc>=36 & trial2_clin_s$tempc<=42,]
spi1_fit <- loess(tempc ~ dpi, spi1)
spi1_nd <- data.frame(dpi=seq(min(spi1$dpi), max(spi1$dpi), length=100))
spi1_nd$fit <- predict(spi1_fit, newdata=spi1_nd) 
spi1_tmax <- spi1_nd$dpi[which.max(spi1_nd$fit)]
spi1[spi1$barn == 6,]$color <- "deeppink2"
# plot(spi1$tempc~spi1$dpi)
# lines(spi1_nd$fit~spi1_nd$dpi, col="grey")

spi0 <- trial2_clin_s[trial2_clin_s$barn == 6 & trial2_clin_s$innoc == 0 & trial2_clin_s$tempc>=36 & trial2_clin_s$tempc<=42,]
spi0_fit <- loess(tempc ~ dpi, spi0)
spi0_nd <- data.frame(dpi=seq(min(spi0$dpi), max(spi0$dpi), length=100))
spi0_nd$fit <- predict(spi0_fit, newdata=spi0_nd) 
spi0_tmax <- spi0_nd$dpi[which.max(spi0_nd$fit)]
spi0[spi0$barn == 6,]$color <- "deeppink2"
# plot(spi0$tempc~spi0$dpi)
# lines(spi0_nd$fit~spi0_nd$dpi, col="grey")

# Clinical Score
gpi1c <- trial2_clin_g[trial2_clin_g$barn == 4 & trial2_clin_g$innoc == 1,]
gpi1c_fit <- loess(clinscore_mod ~ dpi, gpi1c)
gpi1c_nd <- data.frame(dpi=seq(min(gpi1c$dpi), max(gpi1c$dpi), length=100))
gpi1c_nd$fit <- predict(gpi1c_fit, newdata=gpi1c_nd) 
gpi1_cmax <- gpi1c_nd$dpi[which.max(gpi1c_nd$fit)]
gpi1c[gpi1c$barn == 4,]$color <- "cyan4"
# plot(gpi1$clinscore_mod~gpi1c$dpi)
# lines(gpi1c_nd$fit~gpi1c_nd$dpi, col="grey")

gpi0c <- trial2_clin_g[trial2_clin_g$barn == 4 & trial2_clin_g$innoc == 0,]
gpi0c_fit <- loess(clinscore_mod ~ dpi, gpi0c)
gpi0c_nd <- data.frame(dpi=seq(min(gpi0c$dpi), max(gpi0c$dpi), length=100))
gpi0c_nd$fit <- predict(gpi0c_fit, newdata=gpi0c_nd) 
gpi0_cmax <- gpi0c_nd$dpi[which.max(gpi0c_nd$fit)]
gpi0c[gpi0c$barn == 4,]$color <- "cyan4"
# plot(gpi1$clinscore_mod~gpi0c$dpi)
# lines(gpi0c_nd$fit~gpi0c_nd$dpi, col="grey")

spi1c <- trial2_clin_s[trial2_clin_s$barn == 6 & trial2_clin_s$innoc == 1,]
spi1c_fit <- loess(clinscore_mod ~ dpi, spi1c)
spi1c_nd <- data.frame(dpi=seq(min(spi1c$dpi), max(spi1c$dpi), length=100))
spi1c_nd$fit <- predict(spi1c_fit, newdata=spi1c_nd) 
spi1_cmax <- spi1c_nd$dpi[which.max(spi1c_nd$fit)]
spi1c[spi1c$barn == 6,]$color <- "deeppink2"
# plot(spi1c$clinscore_mod~spi1c$dpi)
# lines(spi1c_nd$fit~spic_nd$dpi, col="grey")

spi0c <- trial2_clin_s[trial2_clin_s$barn == 6 & trial2_clin_s$innoc == 0,]
spi0c_fit <- loess(clinscore_mod ~ dpi, spi0c)
spi0c_nd <- data.frame(dpi=seq(min(spi0c$dpi), max(spi0c$dpi), length=100))
spi0c_nd$fit <- predict(spi0c_fit, newdata=spi0c_nd) 
spi0_cmax <- spi0c_nd$dpi[which.max(spi0c_nd$fit)]
spi0c[spi0c$barn == 6,]$color <- "deeppink2"
# plot(spi0c$clinscore_mod~spi0c$dpi)
# lines(spi0c_nd$fit~spi0c_nd$dpi, col="grey")

# Plots ----
# SHEEP
# Rectal Temperature
rts <- ggplot(trial2_clin_s, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("All Controls", "Inoculated Sheep", "Sentinel Sheep"),
                     values = trial2_cols) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0, 28, 2)) +
  scale_y_continuous(limits = c(37,41)) +
  geom_vline(aes(xintercept = spi1_tmax), color = "deeppink2") +
  geom_vline(aes(xintercept = spi0_tmax), color = "rosybrown2") +
  annotate(geom = "text", x = (spi1_tmax - 1.75), y = 41, label = paste(round(spi1_tmax,1)), color = "deeppink2") +
  annotate(geom = "text", x = (spi0_tmax + 1.75), y = 41, label = paste(round(spi0_tmax,1)), color = "rosybrown2") +
  annotate(geom = "text", x = mean(c(spi1_tmax, spi0_tmax)), y = 41, label = paste(round(spi0_tmax - spi1_tmax,1)), color = "black") +
  geom_segment(x = spi1_tmax, y = 40.75, xend = spi0_tmax, yend = 40.75, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
# Warning about 16 rows is dpi13-28 for the PM animal goat 28 barn 5

# Clinical Score - Modified
css <- ggplot(trial2_clin_s, aes(dpi, clinscore_mod, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("All Controls", "Inoculated Sheep", "Sentinel Sheep"),
                     values = trial2_cols) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0, 28, 2)) +
  geom_vline(aes(xintercept = spi1_cmax), color = "deeppink2") +
  geom_vline(aes(xintercept = spi0_cmax), color = "rosybrown2") +
  annotate(geom = "text", x = (spi1_cmax - 1.75), y = 10, label = paste(round(spi1_cmax,1)), color = "deeppink2") +
  annotate(geom = "text", x = (spi0_cmax + 1.75), y = 10, label = paste(round(spi0_cmax,1)), color = "rosybrown2") +
  annotate(geom = "text", x = mean(c(spi1_cmax, spi0_cmax)), y = 10, label = paste(round(spi0_cmax - spi1_cmax,1)), color = "black") +
  geom_segment(x = spi1_cmax, y = 9.25, xend = spi0_cmax, yend = 9.25, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))


# GOAT
# Rectal Temperature
rtg <- ggplot(trial2_clin_g, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Inoculated Goats", "Sentinel Goats", "All Controls"),
                     values = trial2_cols) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0, 28, 2)) +
  scale_y_continuous(limits = c(37,41)) +
  geom_vline(aes(xintercept = gpi1_tmax), color = "cyan4") +
  geom_vline(aes(xintercept = gpi0_tmax), color = "cyan3") +
  annotate(geom = "text", x = (gpi1_tmax - 1.75), y = 41, label = paste(round(gpi1_tmax,1)), color = "cyan4") +
  annotate(geom = "text", x = (gpi0_tmax + 1.75), y = 41, label = paste(round(gpi0_tmax,1)), color = "cyan3") +
  annotate(geom = "text", x = mean(c(gpi1_tmax, gpi0_tmax)), y = 41, label = paste(round(gpi0_tmax - gpi1_tmax,1)), color = "black") +
  geom_segment(x = gpi1_tmax, y = 40.75, xend = gpi0_tmax, yend = 40.75, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
# Warning about 16 rows is dpi13-28 for the PM animal goat 28 barn 5

# Clinical Score - Modified
csg <- ggplot(trial2_clin_g, aes(dpi, clinscore_mod, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Inoculated Goats", "Sentinel Goats", "All Controls"),
                     values = trial2_cols) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0, 28, 2)) +
  geom_vline(aes(xintercept = gpi1_cmax), color = "cyan4") +
  geom_vline(aes(xintercept = gpi0_cmax), color = "cyan3") +
  annotate(geom = "text", x = (gpi1_cmax - 1.75), y = 10, label = paste(round(gpi1_cmax,1)), color = "cyan4") +
  annotate(geom = "text", x = (gpi0_cmax + 1.75), y = 10, label = paste(round(gpi0_cmax,1)), color = "cyan3") +
  annotate(geom = "text", x = mean(c(gpi1_cmax, gpi0_cmax)), y = 10, label = paste(round(gpi0_cmax - gpi1_cmax,1)), color = "black") +
  geom_segment(x = gpi1_cmax, y = 9.25, xend = gpi0_cmax, yend = 9.25, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))


jpeg("output/clin/trial2_clin_sheep_FigS5AB.jpeg", width = 8, height = 4, units = "in", quality = 100, res = 600)
grid.arrange(rts, css, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial2_clin_goats_FigS5AB.jpeg", width = 8, height = 4, units = "in", quality = 100, res = 600)
grid.arrange(rtg, csg, ncol=2)
invisible(dev.off())



#####################
#  Trial 3 ----
#####################

trial3_cols <- c("5a" = "cyan4", "5b" = "cyan3", "6b" = "black", "6c" = "black", "3a" = "purple4", "3b" = "cyan3", 
                 "4a" = "cyan4", "4b" = "cyan3", "4c" = "mediumpurple1")
trial3_cols_sp <- c("Bovine" = "mediumpurple1", "Caprine" = "cyan3")

trial3_b3 <- subset(trial3_clin, barn == 3) # box A Inoculated cattle to sentinel goats
trial3_b5 <- subset(trial3_clin, barn == 5) # box B positive control goats
trial3_b4 <- subset(trial3_clin, barn == 4) # box C Inoculated goats to sentinel cattle and goats (confirmatory)
trial3_b6 <- subset(trial3_clin, barn == 6) # negative controls - not sampled

# Barn 3 / Box A ----
# # Calculating Peaks ----
# # Rectal Temperature
# trial3_b3i1 <- trial3_b3[trial3_b3$innoc == 1 & trial3_b3$tempc>=36 & trial3_b3$tempc<=42,]
# trial3_b3i1_fit <- loess(tempc ~ dpi, trial3_b3i1)
# trial3_b3i1_nd <- data.frame(dpi=seq(min(trial3_b3i1$dpi), max(trial3_b3i1$dpi), length=100))
# trial3_b3i1_nd$fit <- predict(trial3_b3i1_fit, newdata=trial3_b3i1_nd) 
# trial3_b3i1tmax <- trial3_b3i1_nd$dpi[which.max(trial3_b3i1_nd$fit)]
# trial3_b3[trial3_b3$innoc == 1,]$color <- "cyan4"
# # plot(trial3_b3i1$tempc~trial3_b3i1$dpi)
# # lines(trial3_b3i1_nd$fit~trial3_b3i1_nd$dpi, col="grey")
# 
# trial3_b3i0 <- trial3_b3[trial3_b3$innoc == 0 & trial3_b3$tempc>=36 & trial3_b3$tempc<=42,]
# trial3_b3i0_fit <- loess(tempc ~ dpi, trial3_b3i0)
# trial3_b3i0_nd <- data.frame(dpi=seq(min(trial3_b3i0$dpi), max(trial3_b3i0$dpi), length=100))
# trial3_b3i0_nd$fit <- predict(trial3_b3i0_fit, newdata=trial3_b3i0_nd) 
# trial3_b3i0tmax <- trial3_b3i0_nd$dpi[which.max(trial3_b3i0_nd$fit)]
# # plot(trial3_b3i0$tempc~trial3_b3i0$dpi)
# # lines(trial3_b3i0_nd$fit~trial3_b3i0_nd$dpi, col="grey")
# 
# # Clinical Score
# trial3_b3i1c <- trial3_b3[trial3_b3$innoc == 1,]
# trial3_b3i1c_fit <- loess(autosum ~ dpi, trial3_b3i1c)
# trial3_b3i1c_nd <- data.frame(dpi=seq(min(trial3_b3i1c$dpi), max(trial3_b3i1c$dpi), length=100))
# trial3_b3i1c_nd$fit <- predict(trial3_b3i1c_fit, newdata=trial3_b3i1c_nd) 
# trial3_b3i1cmax <- trial3_b3i1c_nd$dpi[which.max(trial3_b3i1c_nd$fit)]
# # plot(trial3_b3i1c$autosum~trial3_b3i1c$dpi)
# # lines(trial3_b3i1c_nd$fit~trial3_b3i1c_nd$dpi, col="grey")
# 
# trial3_b3i0c <- trial3_b3[trial3_b3$innoc == 0,]
# trial3_b3i0c_fit <- loess(autosum ~ dpi, trial3_b3i0c)
# trial3_b3i0c_nd <- data.frame(dpi=seq(min(trial3_b3i0c$dpi), max(trial3_b3i0c$dpi), length=100))
# trial3_b3i0c_nd$fit <- predict(trial3_b3i0c_fit, newdata=trial3_b3i0c_nd) 
# trial3_b3i0cmax <- trial3_b3i0c_nd$dpi[which.max(trial3_b3i0c_nd$fit)]
# # plot(trial3_b3i0c$autosum~trial3_b3i0c$dpi)
# # lines(trial3_b3i0c_nd$fit~trial3_b3i0c_nd$dpi, col="grey")

# Plots ----
# Rectal Temperature
rtA <- ggplot(trial3_b3, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x, , k=7), size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial3_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(36,42)) +
  # geom_vline(aes(xintercept = trial3_b3i1tmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial3_b3i0tmax), color = "cyan3") +
  # annotate(geom = "text", x = (trial3_b3i1tmax + 1), y = 42, label = paste(round(trial3_b3i1tmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial3_b3i0tmax + 1), y = 42, label = paste(round(trial3_b3i0tmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = mean(c(trial3_b3i1tmax, trial3_b3i0tmax)), y = 41.5, label = paste(round(trial3_b3i0tmax - trial3_b3i1tmax,1)), color = "black") +
  # geom_segment(x = trial3_b3i1tmax, y = 41.25, xend = trial3_b3i0tmax, yend = 41.25, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))

# Clinical Score - Modified
csAm <- ggplot(trial3_b3, aes(dpi, clinscore_mod, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), size =1) +
  labs(x = "Days Post Infection", y = "Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial3_cols) +
  scale_y_continuous(limits = c(-1,10), breaks = seq(-1, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  # geom_vline(aes(xintercept = trial3_b3i1cmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial3_b3i0cmax), color = "cyan3") +
  # annotate(geom = "text", x = (trial3_b3i1cmax + 1), y = 10, label = paste(round(trial3_b3i1cmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial3_b3i0cmax + 1), y = 10, label = paste(round(trial3_b3i0cmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = mean(c(trial3_b3i1cmax, trial3_b3i0cmax)), y = 9.5, label = paste(round(trial3_b3i0cmax - trial3_b3i1cmax,1)), color = "black") +
  # geom_segment(x = trial3_b3i1cmax, y = 9.25, xend = trial3_b3i0cmax, yend = 9.25, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))


# Barn 5 / Box B / Positive Controls  ----
# Calculating Peaks ----
# Rectal Temperature
trial3_b5i1 <- trial3_b5[trial3_b5$innoc == 1 & trial3_b5$tempc>=36 & trial3_b5$tempc<=42,]
trial3_b5i1_fit <- loess(tempc ~ dpi, trial3_b5i1)
trial3_b5i1_nd <- data.frame(dpi=seq(min(trial3_b5i1$dpi), max(trial3_b5i1$dpi), length=100))
trial3_b5i1_nd$fit <- predict(trial3_b5i1_fit, newdata=trial3_b5i1_nd) 
trial3_b5i1tmax <- trial3_b5i1_nd$dpi[which.max(trial3_b5i1_nd$fit)]
trial3_b5[trial3_b5$innoc == 1,]$color <- "cyan4"
# plot(trial3_b5i1$tempc~trial3_b5i1$dpi)
# lines(trial3_b5i1_nd$fit~trial3_b5i1_nd$dpi, col="grey")

trial3_b5i0 <- trial3_b5[trial3_b5$innoc == 0 & trial3_b5$tempc>=36 & trial3_b5$tempc<=42,]
trial3_b5i0_fit <- loess(tempc ~ dpi, trial3_b5i0)
trial3_b5i0_nd <- data.frame(dpi=seq(min(trial3_b5i0$dpi), max(trial3_b5i0$dpi), length=100))
trial3_b5i0_nd$fit <- predict(trial3_b5i0_fit, newdata=trial3_b5i0_nd) 
trial3_b5i0tmax <- trial3_b5i0_nd$dpi[which.max(trial3_b5i0_nd$fit)]
# plot(trial3_b5i0$tempc~trial3_b5i0$dpi)
# lines(trial3_b5i0_nd$fit~trial3_b5i0_nd$dpi, col="grey")

# Clinical Score- Modified
trial3_b5i1cm <- trial3_b5[trial3_b5$innoc == 1,]
trial3_b5i1cm_fit <- loess(clinscore_mod ~ dpi, trial3_b5i1cm)
trial3_b5i1cm_nd <- data.frame(dpi=seq(min(trial3_b5i1cm$dpi), max(trial3_b5i1cm$dpi), length=100))
trial3_b5i1cm_nd$fit <- predict(trial3_b5i1cm_fit, newdata=trial3_b5i1cm_nd) 
trial3_b5i1cmmax <- trial3_b5i1cm_nd$dpi[which.max(trial3_b5i1cm_nd$fit)]
# plot(trial3_b5i1cm$clinscore_mod~trial3_b5i1cm$dpi)
# lines(trial3_b5i1cm_nd$fit~trial3_b5i1cm_nd$dpi, col="grey")

trial3_b5i0cm <- trial3_b5[trial3_b5$innoc == 0,]
trial3_b5i0cm_fit <- loess(clinscore_mod ~ dpi, trial3_b5i0cm)
trial3_b5i0cm_nd <- data.frame(dpi=seq(min(trial3_b5i0cm$dpi), max(trial3_b5i0cm$dpi), length=100))
trial3_b5i0cm_nd$fit <- predict(trial3_b5i0cm_fit, newdata=trial3_b5i0cm_nd) 
trial3_b5i0cmmax <- trial3_b5i0cm_nd$dpi[which.max(trial3_b5i0cm_nd$fit)]
# plot(trial3_b5i0cm$clinscore_mod~trial3_b5i0cm$dpi)
# lines(trial3_b5i0cm_nd$fit~trial3_b5i0cm_nd$dpi, col="grey")

# Plots ----
# Rectal Temperature
rtB <- ggplot(trial3_b5, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x, , k=7), size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Goats (Inoculated)", "Goats (In Contact)"),
                     values = trial3_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(37,41)) +
  geom_vline(aes(xintercept = trial3_b5i1tmax), color = "cyan4") +
  geom_vline(aes(xintercept = trial3_b5i0tmax), color = "cyan3") +
  annotate(geom = "text", x = (trial3_b5i1tmax + 1.75), y = 41, label = paste(round(trial3_b5i1tmax,1)), color = "cyan4") +
  annotate(geom = "text", x = (trial3_b5i0tmax + 1.75), y = 41, label = paste(round(trial3_b5i0tmax,1)), color = "cyan3") +
  annotate(geom = "text", x = mean(c(trial3_b5i1tmax, trial3_b5i0tmax)), y = 40.5, label = paste(round(trial3_b5i0tmax - trial3_b5i1tmax,1)), color = "black") +
  geom_segment(x = trial3_b5i1tmax, y = 40.35, xend = trial3_b5i0tmax, yend = 40.35, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
# 80 rows removed - 77 from dead animals, 3 from the animals going hypothermic before death (less than 36 y axis)

# Clinical Score - Modified
csBm <- ggplot(trial3_b5, aes(dpi, clinscore_mod, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), size =1) +
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Goats (Inoculated)", "Goats (In Contact)"),
                     values = trial3_cols) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  geom_vline(aes(xintercept = trial3_b5i1cmmax), color = "cyan4") +
  geom_vline(aes(xintercept = trial3_b5i0cmmax), color = "cyan3") +
  annotate(geom = "text", x = (trial3_b5i1cmmax + 1.75), y = 10, label = paste(round(trial3_b5i1cmmax,1)), color = "cyan4") +
  annotate(geom = "text", x = (trial3_b5i0cmmax + 1.75), y = 10, label = paste(round(trial3_b5i0cmmax,1)), color = "cyan3") +
  annotate(geom = "text", x = mean(c(trial3_b5i1cmmax, trial3_b5i0cmmax)), y = 9, label = paste(round(trial3_b5i0cmmax - trial3_b5i1cmmax,1)), color = "black") +
  geom_segment(x = trial3_b5i1cmmax, y = 8.25, xend = trial3_b5i0cmmax, yend = 8.25, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))


# Barn 4 / Box C ----
# Calculating Peaks ----
# Rectal Temperature
trial3_b4i1 <- trial3_b4[trial3_b4$innoc == 1 & trial3_b4$tempc>=36 & trial3_b4$tempc<=42,]
trial3_b4i1_fit <- loess(tempc ~ dpi, trial3_b4i1)
trial3_b4i1_nd <- data.frame(dpi=seq(min(trial3_b4i1$dpi), max(trial3_b4i1$dpi), length=100))
trial3_b4i1_nd$fit <- predict(trial3_b4i1_fit, newdata=trial3_b4i1_nd) 
trial3_b4i1tmax <- trial3_b4i1_nd$dpi[which.max(trial3_b4i1_nd$fit)]
trial3_b4[trial3_b4$innoc == 1,]$color <- "cyan4"
# plot(trial3_b4i1$tempc~trial3_b4i1$dpi)
# lines(trial3_b4i1_nd$fit~trial3_b4i1_nd$dpi, col="grey")

# Sentinels are two species
trial3_b4i0c <- trial3_b4[trial3_b4$innoc == 0 & trial3_b4$tempc>=36 & trial3_b4$tempc<=42 & trial3_b4$species == "Caprine",]
trial3_b4i0c_fit <- loess(tempc ~ dpi, trial3_b4i0c)
trial3_b4i0c_nd <- data.frame(dpi=seq(min(trial3_b4i0c$dpi), max(trial3_b4i0c$dpi), length=100))
trial3_b4i0c_nd$fit <- predict(trial3_b4i0c_fit, newdata=trial3_b4i0c_nd) 
trial3_b4i0ctmax <- trial3_b4i0c_nd$dpi[which.max(trial3_b4i0c_nd$fit)]
# plot(trial3_b4i0c$tempc~trial3_b4i0c$dpi)
# lines(trial3_b4i0c_nd$fit~trial3_b4i0c_nd$dpi, col="grey")

trial3_b4i0b <- trial3_b4[trial3_b4$innoc == 0 & trial3_b4$tempc>=36 & trial3_b4$tempc<=42 & trial3_b4$species == "Bovine",]
trial3_b4i0b_fit <- loess(tempc ~ dpi, trial3_b4i0b)
trial3_b4i0b_nd <- data.frame(dpi=seq(min(trial3_b4i0b$dpi), max(trial3_b4i0b$dpi), length=100))
trial3_b4i0b_nd$fit <- predict(trial3_b4i0b_fit, newdata=trial3_b4i0b_nd) 
trial3_b4i0btmax <- trial3_b4i0b_nd$dpi[which.max(trial3_b4i0b_nd$fit)]
# plot(trial3_b4i0b$tempc~trial3_b4i0b$dpi)
# lines(trial3_b4i0b_nd$fit~trial3_b4i0b_nd$dpi, col="grey")

# Clinical Score - Modified
trial3_b4i1cm <- trial3_b4[trial3_b4$innoc == 1,]
trial3_b4i1cm_fit <- loess(clinscore_mod ~ dpi, trial3_b4i1cm)
trial3_b4i1cm_nd <- data.frame(dpi=seq(min(trial3_b4i1cm$dpi), max(trial3_b4i1cm$dpi), length=100))
trial3_b4i1cm_nd$fit <- predict(trial3_b4i1cm_fit, newdata=trial3_b4i1cm_nd) 
trial3_b4i1cmmax <- trial3_b4i1cm_nd$dpi[which.max(trial3_b4i1cm_nd$fit)]
# plot(trial3_b4i1cm$clinscore_mod~trial3_b4i1cm$dpi)
# lines(trial3_b4i1cm_nd$fit~trial3_b4i1cm_nd$dpi, col="grey")

# Sentinels are two species
trial3_b4i0ccm <- trial3_b4[trial3_b4$innoc == 0 & trial3_b4$species == "Caprine",]
trial3_b4i0ccm_fit <- loess(clinscore_mod ~ dpi, trial3_b4i0ccm)
trial3_b4i0ccm_nd <- data.frame(dpi=seq(min(trial3_b4i0ccm$dpi), max(trial3_b4i0ccm$dpi), length=100))
trial3_b4i0ccm_nd$fit <- predict(trial3_b4i0ccm_fit, newdata=trial3_b4i0ccm_nd) 
trial3_b4i0ccmmax <- trial3_b4i0ccm_nd$dpi[which.max(trial3_b4i0ccm_nd$fit)]
# plot(trial3_b4i0ccm$clinscore_mod~trial3_b4i0ccm$dpi)
# lines(trial3_b4i0ccm_nd$fit~trial3_b4i0ccm_nd$dpi, col="grey")

trial3_b4i0cbm <- trial3_b4[trial3_b4$innoc == 0 & trial3_b4$species == "Bovine",]
trial3_b4i0cbm_fit <- loess(clinscore_mod ~ dpi, trial3_b4i0cbm)
trial3_b4i0cbm_nd <- data.frame(dpi=seq(min(trial3_b4i0cbm$dpi), max(trial3_b4i0cbm$dpi), length=100))
trial3_b4i0cbm_nd$fit <- predict(trial3_b4i0cbm_fit, newdata=trial3_b4i0cbm_nd) 
trial3_b4i0cbmmax <- trial3_b4i0cbm_nd$dpi[which.max(trial3_b4i0cbm_nd$fit)]
# plot(trial3_b4i0cbm$clinscore_mod~trial3_b4i0cbm$dpi)
# lines(trial3_b4i0cbm_nd$fit~trial3_b4i0cbm_nd$dpi, col="grey")

# Plots ----
# Rectal Temperature
rtC <- ggplot(trial3_b4, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x, , k=7), size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Goats (Inoculated)", "Goats (In Contact)", "Cattle (In Contact)"),
                     values = trial3_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(37,41)) +
  geom_vline(aes(xintercept = trial3_b4i1tmax), color = "cyan4") +
  geom_vline(aes(xintercept = trial3_b4i0ctmax), color = "cyan3") +
  geom_vline(aes(xintercept = trial3_b4i0btmax), color = "mediumpurple1") +
  annotate(geom = "text", x = (trial3_b4i1tmax + 1.75), y = 41, label = paste(round(trial3_b4i1tmax,1)), color = "cyan4") +
  annotate(geom = "text", x = (trial3_b4i0ctmax + 1.75), y = 41, label = paste(round(trial3_b4i0ctmax,1)), color = "cyan3") +
  annotate(geom = "text", x = (trial3_b4i0btmax + 1.75), y = 41, label = paste(round(trial3_b4i0btmax,1)), color = "mediumpurple1") +
  annotate(geom = "text", x = mean(c(trial3_b4i1tmax, trial3_b4i0ctmax)), y = 40.7, label = paste(round(trial3_b4i0ctmax - trial3_b4i1tmax,1)), color = "black") +
  annotate(geom = "text", x = mean(c(trial3_b4i1tmax, trial3_b4i0btmax)), y = 40.3, label = paste(round(trial3_b4i0btmax - trial3_b4i1tmax,1)), color = "black") +
  geom_segment(x = trial3_b4i1tmax, y = 40.5, xend = trial3_b4i0ctmax, yend = 40.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  geom_segment(x = trial3_b4i1tmax, y = 40.1, xend = trial3_b4i0btmax, yend = 40.1, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
# Removed 19 rows - the 1 dead Inoculated animal

# Clinical Score - Modified
csCm <- ggplot(trial3_b4, aes(dpi, clinscore_mod, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), size =1) +
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Goats (Inoculated)", "Goats (In Contact)", "Cattle (In Contact)"),
                     values = trial3_cols) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  geom_vline(aes(xintercept = trial3_b4i1cmmax), color = "cyan4") +
  geom_vline(aes(xintercept = trial3_b4i0ccmmax), color = "cyan3") +
  geom_vline(aes(xintercept = trial3_b4i0cbmmax), color = "mediumpurple1") +
  annotate(geom = "text", x = (trial3_b4i1cmmax + 1.75), y = 10, label = paste(round(trial3_b4i1cmmax,1)), color = "cyan4") +
  annotate(geom = "text", x = (trial3_b4i0ccmmax + 3.75), y = 10, label = paste(round(trial3_b4i0ccmmax,1)), color = "cyan3") +
  annotate(geom = "text", x = (trial3_b4i0cbmmax + 1.75), y = 10, label = paste(round(trial3_b4i0cbmmax,1)), color = "mediumpurple1") +
  annotate(geom = "text", x = mean(c(trial3_b4i1cmmax, trial3_b4i0ccmmax)), y = 9, label = paste(round(trial3_b4i0ccmmax - trial3_b4i1cmmax,1)), color = "black") +
  annotate(geom = "text", x = mean(c(trial3_b4i1cmmax, trial3_b4i0cbmmax)), y = 7.5, label = paste(round(trial3_b4i0cbmmax - trial3_b4i1cmmax,1)), color = "black") +
  geom_segment(x = trial3_b4i1cmmax, y = 8.5, xend = trial3_b4i0ccmmax, yend = 8.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  geom_segment(x = trial3_b4i1cmmax, y = 7, xend = trial3_b4i0cbmmax, yend = 7, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))


# Barn 6 / Negative Control ----
# Rectal Temperature
rtCont <- ggplot(trial3_b6, aes(dpi, tempc, color = species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x, , k=7), size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = " ", 
                     labels = c("Cattle", "Goats"),
                     values = trial3_cols_sp) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(37,41)) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))


# Clinical Score - Modified
csmCont <- ggplot(trial3_b6, aes(dpi, clinscore_mod, color=species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), size =1) +
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Goats", "Cattle"),
                     values = trial3_cols_sp) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

# Arranging Plots ----
jpeg("output/clin/trial3_clin_b3_modcs.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtA, csAm, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial3_clin_b5_poscontrols_modcs_FigS10A.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtB, csBm, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial3_clin_b4_modcs_Fig1AB.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtC, csCm, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial3_clin_b6_negcontrols_modcs_FigS10B.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtCont, csmCont, ncol=2)
invisible(dev.off())



#########################
#  Trial 3 Challenge ----
#########################

trial3_chall_cols <- c("3" = "cyan4", "4" = "black", "6" = "cyan2")

# Rectal Temperature
rt_chall <- ggplot(trial3_chall_clin, aes(dpi, tempc, color=barn)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Cattle sentinel goats", "Seropositive goats", "Negative Controls"),
                     values = trial3_chall_cols) +
  scale_x_continuous(limits = c(0,14), breaks = c(0, 4, 7, 10, 14)) +
  scale_y_continuous(limits = c(37,41)) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))


# Clinical Score - Modified
cs_chall <- ggplot(trial3_chall_clin, aes(dpi, clinscore_mod, color=barn)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Cattle sentinel goats", "Seropositive goats", "Negative Controls"),
                     values = trial3_chall_cols) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,14), breaks = c(0, 4, 7, 10, 14)) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))


jpeg("output/clin/trial3_challenge_clin_FigS9A.jpeg", width = 5, height = 5, units = "in", quality = 100, res = 600)
rt_chall
invisible(dev.off())
jpeg("output/clin/trial3_challenge_clin_FigS9B.jpeg", width = 5, height = 5, units = "in", quality = 100, res = 600)
cs_chall
invisible(dev.off())



#####################
#  Trial 4 ----
#####################

trial4_cols <- c("6a" = "purple4", "6b" = "cyan3", "1a" = "purple4", "1b" = "cyan3", "4a" = "purple4", "4b" = "cyan3",
          "5a" = "cyan4", "5b" = "cyan3", "3a" = "cyan3")

trial4_cols_sp <- c("Bovine" = "purple4", "Caprine" = "cyan3")

trial4_b6 <- subset(trial4_clin, barn == 6) # Inoculated cattle to sentinel goats
trial4_b1 <- subset(trial4_clin, barn == 1) # Inoculated cattle to sentinel goats
trial4_b4 <- subset(trial4_clin, barn == 4) # Inoculated cattle to sentinel goats
trial4_b5 <- subset(trial4_clin, barn == 5) # Positive Controls 
trial4_b3 <- subset(trial4_clin, barn == 3) # Negative Controls 
trial4_b614 <- subset(trial4_clin, barn == 6 | barn == 1 | barn ==4) #all experimental together


# Barn 6 ----
# # Calculating peaks ----
# # Rectal Temperature
# trial4_b6i1 <- trial4_b6[trial4_b6$innoc == 1 & trial4_b6$tempc>=36 & trial4_b6$tempc<=42,]
# trial4_b6i1_fit <- loess(tempc ~ dpi, trial4_b6i1)
# trial4_b6i1_nd <- data.frame(dpi=seq(min(trial4_b6i1$dpi), max(trial4_b6i1$dpi), length=100))
# trial4_b6i1_nd$fit <- predict(trial4_b6i1_fit, newdata=trial4_b6i1_nd) 
# trial4_b6i1tmax <- trial4_b6i1_nd$dpi[which.max(trial4_b6i1_nd$fit)]
# trial4_b6[trial4_b6$innoc == 1,]$color <- "cyan4"
# # plot(trial4_b6i1$tempc~trial4_b6i1$dpi)
# # lines(trial4_b6i1_nd$fit~trial4_b6i1_nd$dpi, col="grey")
# 
# trial4_b6i0 <- trial4_b6[trial4_b6$innoc == 0 & trial4_b6$tempc>=36 & trial4_b6$tempc<=42,]
# trial4_b6i0_fit <- loess(tempc ~ dpi, trial4_b6i0)
# trial4_b6i0_nd <- data.frame(dpi=seq(min(trial4_b6i0$dpi), max(trial4_b6i0$dpi), length=100))
# trial4_b6i0_nd$fit <- predict(trial4_b6i0_fit, newdata=trial4_b6i0_nd) 
# trial4_b6i0tmax <- trial4_b6i0_nd$dpi[which.max(trial4_b6i0_nd$fit)]
# # plot(trial4_b6i0$tempc~trial4_b6i0$dpi)
# # lines(trial4_b6i0_nd$fit~trial4_b6i0_nd$dpi, col="grey")
# 
# 
# # Clinical Score
# trial4_b6i1c <- trial4_b6[trial4_b6$innoc == 1,]
# trial4_b6i1c_fit <- loess(autosum ~ dpi, trial4_b6i1c)
# trial4_b6i1c_nd <- data.frame(dpi=seq(min(trial4_b6i1c$dpi), max(trial4_b6i1c$dpi), length=100))
# trial4_b6i1c_nd$fit <- predict(trial4_b6i1c_fit, newdata=trial4_b6i1c_nd) 
# trial4_b6i1cmax <- trial4_b6i1c_nd$dpi[which.max(trial4_b6i1c_nd$fit)]
# # plot(trial4_b6i1c$autosum~trial4_b6i1c$dpi)
# # lines(trial4_b6i1c_nd$fit~trial4_b6i1c_nd$dpi, col="grey")
# 
# trial4_b6i0c <- trial4_b6[trial4_b6$innoc == 0,]
# trial4_b6i0c_fit <- loess(autosum ~ dpi, trial4_b6i0c)
# trial4_b6i0c_nd <- data.frame(dpi=seq(min(trial4_b6i0c$dpi), max(trial4_b6i0c$dpi), length=100))
# trial4_b6i0c_nd$fit <- predict(trial4_b6i0c_fit, newdata=trial4_b6i0c_nd) 
# trial4_b6i0cmax <- trial4_b6i0c_nd$dpi[which.max(trial4_b6i0c_nd$fit)]
# # plot(trial4_b6i0c$autosum~trial4_b6i0c$dpi)
# # lines(trial4_b6i0c_nd$fit~trial4_b6i0c_nd$dpi, col="grey")

# Plots ---
# Rectal Temperature
rtA <- ggplot(trial4_b6, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial4_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(36,42)) +
  # geom_vline(aes(xintercept = trial4_b6i1tmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial4_b6i0tmax), color = "cyan3") +
  # annotate(geom = "text", x = (trial4_b6i1tmax + 1), y = 42, label = paste(round(trial4_b6i1tmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial4_b6i0tmax + 1), y = 42, label = paste(round(trial4_b6i0tmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = mean(c(trial4_b6i1tmax, trial4_b6i0tmax)), y = 41.5, label = paste(round(trial4_b6i0tmax - trial4_b6i1tmax,1)), color = "black") +
  # geom_segment(x = trial4_b6i1tmax, y = 41.25, xend = trial4_b6i0tmax, yend = 41.25, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))

# Clinical Score
csA <- ggplot(trial4_b6, aes(dpi, autosum, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial4_cols) +
  scale_y_continuous(limits = c(-1,10), breaks = seq(-1, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  # geom_vline(aes(xintercept = trial4_b6i1cmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial4_b6i0cmax), color = "cyan3") +
  # annotate(geom = "text", x = (trial4_b6i1cmax + 1), y = 10, label = paste(round(trial4_b6i1cmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial4_b6i0cmax + 1), y = 10, label = paste(round(trial4_b6i0cmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = mean(c(trial4_b6i1cmax, trial4_b6i0cmax)), y = 9.5, label = paste(round(trial4_b6i0cmax - trial4_b6i1cmax,1)), color = "black") +
  # geom_segment(x = trial4_b6i1cmax, y = 9.25, xend = trial4_b6i0cmax, yend = 9.25, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))

# Barn 1 ----
# Calculating peaks ----
# Rectal Temperature
# trial4_b1i1 <- trial4_b1[trial4_b1$innoc == 1 & trial4_b1$tempc>=36 & trial4_b1$tempc<=42,]
# trial4_b1i1_fit <- loess(tempc ~ dpi, trial4_b1i1)
# trial4_b1i1_nd <- data.frame(dpi=seq(min(trial4_b1i1$dpi), max(trial4_b1i1$dpi), length=100))
# trial4_b1i1_nd$fit <- predict(trial4_b1i1_fit, newdata=trial4_b1i1_nd) 
# trial4_b1i1tmax <- trial4_b1i1_nd$dpi[which.max(trial4_b1i1_nd$fit)]
# trial4_b1[trial4_b1$innoc == 1,]$color <- "cyan4"
# # plot(trial4_b1i1$tempc~trial4_b1i1$dpi)
# # lines(trial4_b1i1_nd$fit~trial4_b1i1_nd$dpi, col="grey")
# 
# trial4_b1i0 <- trial4_b1[trial4_b1$innoc == 0 & trial4_b1$tempc>=36 & trial4_b1$tempc<=42,]
# trial4_b1i0_fit <- loess(tempc ~ dpi, trial4_b1i0)
# trial4_b1i0_nd <- data.frame(dpi=seq(min(trial4_b1i0$dpi), max(trial4_b1i0$dpi), length=100))
# trial4_b1i0_nd$fit <- predict(trial4_b1i0_fit, newdata=trial4_b1i0_nd) 
# trial4_b1i0tmax <- trial4_b1i0_nd$dpi[which.max(trial4_b1i0_nd$fit)]
# # plot(trial4_b1i0$tempc~trial4_b1i0$dpi)
# # lines(trial4_b1i0_nd$fit~trial4_b1i0_nd$dpi, col="grey")

# Clinical Score
# trial4_b1i1c <- trial4_b1[trial4_b1$innoc == 1,]
# trial4_b1i1c_fit <- loess(autosum ~ dpi, trial4_b1i1c)
# trial4_b1i1c_nd <- data.frame(dpi=seq(min(trial4_b1i1c$dpi), max(trial4_b1i1c$dpi), length=100))
# trial4_b1i1c_nd$fit <- predict(trial4_b1i1c_fit, newdata=trial4_b1i1c_nd) 
# trial4_b1i1cmax <- trial4_b1i1c_nd$dpi[which.max(trial4_b1i1c_nd$fit)]
# # plot(trial4_b1i1c$autosum~trial4_b1i1c$dpi)
# # lines(trial4_b1i1c_nd$fit~trial4_b1i1c_nd$dpi, col="grey")
# 
# trial4_b1i0c <- trial4_b1[trial4_b1$innoc == 0,]
# trial4_b1i0c_fit <- loess(autosum ~ dpi, trial4_b1i0c)
# trial4_b1i0c_nd <- data.frame(dpi=seq(min(trial4_b1i0c$dpi), max(trial4_b1i0c$dpi), length=100))
# trial4_b1i0c_nd$fit <- predict(trial4_b1i0c_fit, newdata=trial4_b1i0c_nd) 
# trial4_b1i0cmax <- trial4_b1i0c_nd$dpi[which.max(trial4_b1i0c_nd$fit)]
# # plot(trial4_b1i0c$autosum~trial4_b1i0c$dpi)
# # lines(trial4_b1i0c_nd$fit~trial4_b1i0c_nd$dpi, col="grey")

# Plots ----
# Rectal Temperature
rtB <- ggplot(trial4_b1, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial4_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(36,42)) +
  # geom_vline(aes(xintercept = trial4_b1i1tmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial4_b1i0tmax), color = "cyan3") +
  # annotate(geom = "text", x = (trial4_b1i1tmax + 1.75), y = 42, label = paste(round(trial4_b1i1tmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial4_b1i0tmax + 1.75), y = 42, label = paste(round(trial4_b1i0tmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = mean(c(trial4_b1i1tmax, trial4_b1i0tmax)), y = 41.5, label = paste(round(trial4_b1i0tmax - trial4_b1i1tmax,1)), color = "black") +
  # geom_segment(x = trial4_b1i1tmax, y = 41.25, xend = trial4_b1i0tmax, yend = 41.25, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))

# Clinical Score
csB <- ggplot(trial4_b1, aes(dpi, autosum, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial4_cols) +
  scale_y_continuous(limits = c(-1,10), breaks = seq(-1, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  # geom_vline(aes(xintercept = trial4_b1i1cmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial4_b1i0cmax), color = "cyan3") +
  # annotate(geom = "text", x = (trial4_b1i1cmax + 1.75), y = 10, label = paste(round(trial4_b1i1cmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial4_b1i0cmax + 1.75), y = 10, label = paste(round(trial4_b1i0cmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = mean(c(trial4_b1i1cmax, trial4_b1i0cmax)), y = 9, label = paste(round(trial4_b1i0cmax - trial4_b1i1cmax,1)), color = "black") +
  # geom_segment(x = trial4_b1i1cmax, y = 8.25, xend = trial4_b1i0cmax, yend = 8.25, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))


# Barn 4 ----
# # Calculating peaks
# # Rectal Temperature
# trial4_b4i1 <- trial4_b4[trial4_b4$innoc == 1 & trial4_b4$tempc>=36 & trial4_b4$tempc<=42,]
# trial4_b4i1_fit <- loess(tempc ~ dpi, trial4_b4i1)
# trial4_b4i1_nd <- data.frame(dpi=seq(min(trial4_b4i1$dpi), max(trial4_b4i1$dpi), length=100))
# trial4_b4i1_nd$fit <- predict(trial4_b4i1_fit, newdata=trial4_b4i1_nd) 
# trial4_b4i1tmax <- trial4_b4i1_nd$dpi[which.max(trial4_b4i1_nd$fit)]
# trial4_b4[trial4_b4$innoc == 1,]$color <- "cyan4"
# # plot(trial4_b4i1$tempc~trial4_b4i1$dpi)
# # lines(trial4_b4i1_nd$fit~trial4_b4i1_nd$dpi, col="grey")
# 
# # Sentinels are two species
# trial4_b4i0c <- trial4_b4[trial4_b4$innoc == 0 & trial4_b4$tempc>=36 & trial4_b4$tempc<=42 & trial4_b4$species == "Caprine",]
# trial4_b4i0c_fit <- loess(tempc ~ dpi, trial4_b4i0c)
# trial4_b4i0c_nd <- data.frame(dpi=seq(min(trial4_b4i0c$dpi), max(trial4_b4i0c$dpi), length=100))
# trial4_b4i0c_nd$fit <- predict(trial4_b4i0c_fit, newdata=trial4_b4i0c_nd) 
# trial4_b4i0ctmax <- trial4_b4i0c_nd$dpi[which.max(trial4_b4i0c_nd$fit)]
# # plot(trial4_b4i0c$tempc~trial4_b4i0c$dpi)
# # lines(trial4_b4i0c_nd$fit~trial4_b4i0c_nd$dpi, col="grey")
# 
# trial4_b4i0b <- trial4_b4[trial4_b4$innoc == 0 & trial4_b4$tempc>=36 & trial4_b4$tempc<=42 & trial4_b4$species == "Bovine",]
# trial4_b4i0b_fit <- loess(tempc ~ dpi, trial4_b4i0b)
# trial4_b4i0b_nd <- data.frame(dpi=seq(min(trial4_b4i0b$dpi), max(trial4_b4i0b$dpi), length=100))
# trial4_b4i0b_nd$fit <- predict(trial4_b4i0b_fit, newdata=trial4_b4i0b_nd) 
# trial4_b4i0btmax <- trial4_b4i0b_nd$dpi[which.max(trial4_b4i0b_nd$fit)]
# # plot(trial4_b4i0b$tempc~trial4_b4i0b$dpi)
# # lines(trial4_b4i0b_nd$fit~trial4_b4i0b_nd$dpi, col="grey")
# 
# # Clinical Score
# trial4_b4i1c <- trial4_b4[trial4_b4$innoc == 1,]
# trial4_b4i1c_fit <- loess(autosum ~ dpi, trial4_b4i1c)
# trial4_b4i1c_nd <- data.frame(dpi=seq(min(trial4_b4i1c$dpi), max(trial4_b4i1c$dpi), length=100))
# trial4_b4i1c_nd$fit <- predict(trial4_b4i1c_fit, newdata=trial4_b4i1c_nd) 
# trial4_b4i1cmax <- trial4_b4i1c_nd$dpi[which.max(trial4_b4i1c_nd$fit)]
# # plot(trial4_b4i1c$autosum~trial4_b4i1c$dpi)
# # lines(trial4_b4i1c_nd$fit~trial4_b4i1c_nd$dpi, col="grey")
# 
# # Sentinels are two species
# trial4_b4i0cc <- trial4_b4[trial4_b4$innoc == 0 & trial4_b4$species == "Caprine",]
# trial4_b4i0cc_fit <- loess(autosum ~ dpi, trial4_b4i0cc)
# trial4_b4i0cc_nd <- data.frame(dpi=seq(min(trial4_b4i0cc$dpi), max(trial4_b4i0cc$dpi), length=100))
# trial4_b4i0cc_nd$fit <- predict(trial4_b4i0cc_fit, newdata=trial4_b4i0cc_nd) 
# trial4_b4i0ccmax <- trial4_b4i0cc_nd$dpi[which.max(trial4_b4i0cc_nd$fit)]
# # plot(trial4_b4i0cc$autosum~trial4_b4i0cc$dpi)
# # lines(trial4_b4i0cc_nd$fit~trial4_b4i0cc_nd$dpi, col="grey")
# 
# trial4_b4i0cb <- trial4_b4[trial4_b4$innoc == 0 & trial4_b4$species == "Bovine",]
# trial4_b4i0cb_fit <- loess(autosum ~ dpi, trial4_b4i0cb)
# trial4_b4i0cb_nd <- data.frame(dpi=seq(min(trial4_b4i0cb$dpi), max(trial4_b4i0cb$dpi), length=100))
# trial4_b4i0cb_nd$fit <- predict(trial4_b4i0cb_fit, newdata=trial4_b4i0cb_nd) 
# trial4_b4i0cbmax <- trial4_b4i0cb_nd$dpi[which.max(trial4_b4i0cb_nd$fit)]
# # plot(trial4_b4i0cb$autosum~trial4_b4i0cb$dpi)
# # lines(trial4_b4i0cb_nd$fit~trial4_b4i0cb_nd$dpi, col="grey")

# Plots ----
# Rectal Temperature
rtC <- ggplot(trial4_b4, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial4_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(36,42)) +
  # geom_vline(aes(xintercept = trial4_b4i1tmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial4_b4i0ctmax), color = "cyan3") +
  # geom_vline(aes(xintercept = trial4_b4i0btmax), color = "mediumpurple1") +
  # annotate(geom = "text", x = (trial4_b4i1tmax + 1.75), y = 42, label = paste(round(trial4_b4i1tmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial4_b4i0ctmax + 1.75), y = 42, label = paste(round(trial4_b4i0ctmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = (trial4_b4i0btmax + 1.75), y = 42, label = paste(round(trial4_b4i0btmax,1)), color = "mediumpurple1") +
  # annotate(geom = "text", x = mean(c(trial4_b4i1tmax, trial4_b4i0ctmax)), y = 41.35, label = paste(round(trial4_b4i0ctmax - trial4_b4i1tmax,1)), color = "black") +
  # annotate(geom = "text", x = mean(c(trial4_b4i1tmax, trial4_b4i0btmax)), y = 40.55, label = paste(round(trial4_b4i0btmax - trial4_b4i1tmax,1)), color = "black") +
  # geom_segment(x = trial4_b4i1tmax, y = 41, xend = trial4_b4i0ctmax, yend = 41, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  # geom_segment(x = trial4_b4i1tmax, y = 40.2, xend = trial4_b4i0btmax, yend = 40.2, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))
# Removed 16 rows - the 1 dead sentinel animal, eartag 181 - died by dpi 20 so has 16 missing temp measurements for row.

# Clinical Score
csC <- ggplot(trial4_b4, aes(dpi, autosum, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)", "Cattle (In Contact)"),
                     values = trial4_cols) +
  scale_y_continuous(limits = c(-1,10), breaks = seq(-1, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  # geom_vline(aes(xintercept = trial4_b4i1cmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial4_b4i0ccmax), color = "cyan3") +
  # geom_vline(aes(xintercept = trial4_b4i0cbmax), color = "mediumpurple1") +
  # annotate(geom = "text", x = (trial4_b4i1cmax + 1.75), y = 10, label = paste(round(trial4_b4i1cmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial4_b4i0ccmax + 3.75), y = 10, label = paste(round(trial4_b4i0ccmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = (trial4_b4i0cbmax + 1.75), y = 10, label = paste(round(trial4_b4i0cbmax,1)), color = "mediumpurple1") +
  # annotate(geom = "text", x = mean(c(trial4_b4i1cmax, trial4_b4i0ccmax)), y = 9, label = paste(round(trial4_b4i0ccmax - trial4_b4i1cmax,1)), color = "black") +
  # annotate(geom = "text", x = mean(c(trial4_b4i1cmax, trial4_b4i0cbmax)), y = 7.5, label = paste(round(trial4_b4i0cbmax - trial4_b4i1cmax,1)), color = "black") +
  # geom_segment(x = trial4_b4i1cmax, y = 8.5, xend = trial4_b4i0ccmax, yend = 8.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  # geom_segment(x = trial4_b4i1cmax, y = 7, xend = trial4_b4i0cbmax, yend = 7, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))
# removed 2 rows but all in range so unclear why


# Barn 5 / Positive Controls ----
# Plots ---
# Rectal Temperature
rtpcont <- ggplot(trial4_b5, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(37,41)) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
# Inoculated smooth goes down so far because 3 die
# Eartag 172 Inoculated died 12/6/2020 by dpi 10, 175 sentinel and 179 Inoculated died by dpi 21

# Clinical Score - Modified
cspcont <- ggplot(trial4_b5, aes(dpi, clinscore_mod, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
# Eartag 172 Inoculated died 12/6/2020 by dpi 10, 175 sentinal and 179 Inoculated died? by dpi 21


# Barn 3 / Negative Controls  ----
# Plots ----
# Rectal Temperature
rtncont <- ggplot(trial4_b3, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Goats (Negative Control)"),
                     values = trial4_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(37,41)) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
#Inoculated smooth goes down so far because 3 die

# Clinical Score
csncont <- ggplot(trial4_b3, aes(dpi, clinscore_mod, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Goats (Negative Control)"),
                     values = trial4_cols) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))


# All three experimental barns together ----
# Plots ----
# Rectal Temperature
rtABC <- ggplot(trial4_b614, aes(dpi, tempc, color=species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols_sp) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(36,42)) +
  theme_minimal() +
  theme(legend.position = "top")

# Clinical Score
csABC <- ggplot(trial4_b614, aes(dpi, autosum, color=species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols_sp) +
  scale_y_continuous(limits = c(-1,10), breaks = seq(-1, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  theme_minimal() +
  theme(legend.position = "top")


# Arranging Plots ----
jpeg("output/clin/trial4_clin_b6.jpeg", width = 7, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtA, csA, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial4_clin_b1.jpeg", width = 7, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtB, csB, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial4_clin_b4.jpeg", width = 7, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtC, csC, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial4_clin_b5_FigS11A.jpeg", width = 7, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtpcont, cspcont, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial4_clin_b3_FigS11B.jpeg", width = 7, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtncont, csncont, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial4_clin_b614.jpeg", width = 7, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtABC, csABC, ncol=2)
invisible(dev.off())


#####################
#  Trial 5 ----
#####################

trial5_cols <- c("1a" = "purple4", "2a" = "purple4", "5a" = "purple4", "6a" = "purple4",
          "1b" = "cyan3", "2b" = "cyan3", "5b" = "cyan3", "6b" = "cyan3",
          "3a" = "cyan4", "4b" = "cyan3") 

trial5_cols_sp <- c("Bovine" = "purple4", "Caprine" = "cyan3")


trial5_b5 <- subset(trial5_clin, barn == 5) # box A Inoculated cattle to sentinel goats
trial5_b6 <- subset(trial5_clin, barn == 6) # box B Inoculated cattle to sentinel goats
trial5_b2 <- subset(trial5_clin, barn == 2) # box C Inoculated cattle to sentinel goats
trial5_b1 <- subset(trial5_clin, barn == 1) # box D Inoculated cattle to sentinel goats
trial5_b3 <- subset(trial5_clin, barn == 3) # positive controls Box E
trial5_b4 <- subset(trial5_clin, barn == 4) # negative controls Box F
trial5_b1256 <- subset(trial5_clin, barn == 1 | barn == 2 | barn == 5 | barn == 6) #all experimental together


# Barn 5 ----
# # Calculating peaks ----
# # Rectal Temperature
# trial5_b5i1 <- trial5_b5[trial5_b5$innoc == 1 & trial5_b5$tempc>=36 & trial5_b5$tempc<=42,]
# trial5_b5i1_fit <- loess(tempc ~ dpi, trial5_b5i1)
# trial5_b5i1_nd <- data.frame(dpi=seq(min(trial5_b5i1$dpi), max(trial5_b5i1$dpi), length=100))
# trial5_b5i1_nd$fit <- predict(trial5_b5i1_fit, newdata=trial5_b5i1_nd) 
# trial5_b5i1tmax <- trial5_b5i1_nd$dpi[which.max(trial5_b5i1_nd$fit)]
# trial5_b5[trial5_b5$innoc == 1,]$color <- "cyan4"
# # plot(trial5_b5i1$tempc~trial5_b5i1$dpi)
# # lines(trial5_b5i1_nd$fit~trial5_b5i1_nd$dpi, col="grey")
# 
# trial5_b5i0 <- trial5_b5[trial5_b5$innoc == 0 & trial5_b5$tempc>=36 & trial5_b5$tempc<=42,]
# trial5_b5i0_fit <- loess(tempc ~ dpi, trial5_b5i0)
# trial5_b5i0_nd <- data.frame(dpi=seq(min(trial5_b5i0$dpi), max(trial5_b5i0$dpi), length=100))
# trial5_b5i0_nd$fit <- predict(trial5_b5i0_fit, newdata=trial5_b5i0_nd) 
# trial5_b5i0tmax <- trial5_b5i0_nd$dpi[which.max(trial5_b5i0_nd$fit)]
# # plot(trial5_b5i0$tempc~trial5_b5i0$dpi)
# # lines(trial5_b5i0_nd$fit~trial5_b5i0_nd$dpi, col="grey")
# 
# # Clinical Score
# trial5_b5i1c <- trial5_b5[trial5_b5$innoc == 1,]
# trial5_b5i1c_fit <- loess(autosum ~ dpi, trial5_b5i1c)
# trial5_b5i1c_nd <- data.frame(dpi=seq(min(trial5_b5i1c$dpi), max(trial5_b5i1c$dpi), length=100))
# trial5_b5i1c_nd$fit <- predict(trial5_b5i1c_fit, newdata=trial5_b5i1c_nd) 
# trial5_b5i1cmax <- trial5_b5i1c_nd$dpi[which.max(trial5_b5i1c_nd$fit)]
# # plot(trial5_b5i1c$autosum~trial5_b5i1c$dpi)
# # lines(trial5_b5i1c_nd$fit~trial5_b5i1c_nd$dpi, col="grey")
# 
# trial5_b5i0c <- trial5_b5[trial5_b5$innoc == 0,]
# trial5_b5i0c_fit <- loess(autosum ~ dpi, trial5_b5i0c)
# trial5_b5i0c_nd <- data.frame(dpi=seq(min(trial5_b5i0c$dpi), max(trial5_b5i0c$dpi), length=100))
# trial5_b5i0c_nd$fit <- predict(trial5_b5i0c_fit, newdata=trial5_b5i0c_nd) 
# trial5_b5i0cmax <- trial5_b5i0c_nd$dpi[which.max(trial5_b5i0c_nd$fit)]
# # plot(trial5_b5i0c$autosum~trial5_b5i0c$dpi)
# # lines(trial5_b5i0c_nd$fit~trial5_b5i0c_nd$dpi, col="grey")

# Plots ----
# Rectal Temperature
rtA <- ggplot(trial5_b5, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Rectal Temperature (Celsius)") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial5_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(36,42)) +
  # geom_vline(aes(xintercept = trial5_b5i1tmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial5_b5i0tmax), color = "cyan3") +
  # annotate(geom = "text", x = (trial5_b5i1tmax + 1), y = 42, label = paste(round(trial5_b5i1tmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial5_b5i0tmax + 1), y = 42, label = paste(round(trial5_b5i0tmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = mean(c(trial5_b5i1tmax, trial5_b5i0tmax)), y = 41.5, label = paste(round(trial5_b5i0tmax - trial5_b5i1tmax,1)), color = "black") +
  # geom_segment(x = trial5_b5i1tmax, y = 41.25, xend = trial5_b5i0tmax, yend = 41.25, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))


# Clinical Score
csA <- ggplot(trial5_b5, aes(dpi, autosum, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial5_cols) +
  scale_y_continuous(limits = c(-1,10), breaks = seq(-1, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  # geom_vline(aes(xintercept = trial5_b5i1cmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial5_b5i0cmax), color = "cyan3") +
  # annotate(geom = "text", x = (trial5_b5i1cmax + 1), y = 10, label = paste(round(trial5_b5i1cmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial5_b5i0cmax + 1), y = 10, label = paste(round(trial5_b5i0cmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = mean(c(trial5_b5i1cmax, trial5_b5i0cmax)), y = 9.5, label = paste(round(trial5_b5i0cmax - trial5_b5i1cmax,1)), color = "black") +
  # geom_segment(x = trial5_b5i1cmax, y = 9.25, xend = trial5_b5i0cmax, yend = 9.25, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))


# Barn 6 ----
# Calculating peaks
# Rectal Temperature
# trial5_b6i1 <- trial5_b6[trial5_b6$innoc == 1 & trial5_b6$tempc>=36 & trial5_b6$tempc<=42,]
# trial5_b6i1_fit <- loess(tempc ~ dpi, trial5_b6i1)
# trial5_b6i1_nd <- data.frame(dpi=seq(min(trial5_b6i1$dpi), max(trial5_b6i1$dpi), length=100))
# trial5_b6i1_nd$fit <- predict(trial5_b6i1_fit, newdata=trial5_b6i1_nd) 
# trial5_b6i1tmax <- trial5_b6i1_nd$dpi[which.max(trial5_b6i1_nd$fit)]
# trial5_b6[trial5_b6$innoc == 1,]$color <- "cyan4"
# # plot(trial5_b6i1$tempc~trial5_b6i1$dpi)
# # lines(trial5_b6i1_nd$fit~trial5_b6i1_nd$dpi, col="grey")
# 
# trial5_b6i0 <- trial5_b6[trial5_b6$innoc == 0 & trial5_b6$tempc>=36 & trial5_b6$tempc<=42,]
# trial5_b6i0_fit <- loess(tempc ~ dpi, trial5_b6i0)
# trial5_b6i0_nd <- data.frame(dpi=seq(min(trial5_b6i0$dpi), max(trial5_b6i0$dpi), length=100))
# trial5_b6i0_nd$fit <- predict(trial5_b6i0_fit, newdata=trial5_b6i0_nd) 
# trial5_b6i0tmax <- trial5_b6i0_nd$dpi[which.max(trial5_b6i0_nd$fit)]
# # plot(trial5_b6i0$tempc~trial5_b6i0$dpi)
# # lines(trial5_b6i0_nd$fit~trial5_b6i0_nd$dpi, col="grey")

# Clinical Score
# trial5_b6i1c <- trial5_b6[trial5_b6$innoc == 1,]
# trial5_b6i1c_fit <- loess(autosum ~ dpi, trial5_b6i1c)
# trial5_b6i1c_nd <- data.frame(dpi=seq(min(trial5_b6i1c$dpi), max(trial5_b6i1c$dpi), length=100))
# trial5_b6i1c_nd$fit <- predict(trial5_b6i1c_fit, newdata=trial5_b6i1c_nd) 
# trial5_b6i1cmax <- trial5_b6i1c_nd$dpi[which.max(trial5_b6i1c_nd$fit)]
# # plot(trial5_b6i1c$autosum~trial5_b6i1c$dpi)
# # lines(trial5_b6i1c_nd$fit~trial5_b6i1c_nd$dpi, col="grey")
# 
# trial5_b6i0c <- trial5_b6[trial5_b6$innoc == 0,]
# trial5_b6i0c_fit <- loess(autosum ~ dpi, trial5_b6i0c)
# trial5_b6i0c_nd <- data.frame(dpi=seq(min(trial5_b6i0c$dpi), max(trial5_b6i0c$dpi), length=100))
# trial5_b6i0c_nd$fit <- predict(trial5_b6i0c_fit, newdata=trial5_b6i0c_nd) 
# trial5_b6i0cmax <- trial5_b6i0c_nd$dpi[which.max(trial5_b6i0c_nd$fit)]
# # plot(trial5_b6i0c$autosum~trial5_b6i0c$dpi)
# # lines(trial5_b6i0c_nd$fit~trial5_b6i0c_nd$dpi, col="grey")

# Plots ----
# Rectal Temperature
rtB <- ggplot(trial5_b6, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Rectal Temperature (Celsius)") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial5_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(36,42)) +
  # geom_vline(aes(xintercept = trial5_b6i1tmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial5_b6i0tmax), color = "cyan3") +
  # annotate(geom = "text", x = (trial5_b6i1tmax + 1.75), y = 42, label = paste(round(trial5_b6i1tmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial5_b6i0tmax + 1.75), y = 42, label = paste(round(trial5_b6i0tmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = mean(c(trial5_b6i1tmax, trial5_b6i0tmax)), y = 41.5, label = paste(round(trial5_b6i0tmax - trial5_b6i1tmax,1)), color = "black") +
  # geom_segment(x = trial5_b6i1tmax, y = 41.25, xend = trial5_b6i0tmax, yend = 41.25, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))


# Clinical Score
csB <- ggplot(trial5_b6, aes(dpi, autosum, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial5_cols) +
  scale_y_continuous(limits = c(-1,10), breaks = seq(-1, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  # geom_vline(aes(xintercept = trial5_b6i1cmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial5_b6i0cmax), color = "cyan3") +
  # annotate(geom = "text", x = (trial5_b6i1cmax + 1.75), y = 10, label = paste(round(trial5_b6i1cmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial5_b6i0cmax + 1.75), y = 10, label = paste(round(trial5_b6i0cmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = mean(c(trial5_b6i1cmax, trial5_b6i0cmax)), y = 9, label = paste(round(trial5_b6i0cmax - trial5_b6i1cmax,1)), color = "black") +
  # geom_segment(x = trial5_b6i1cmax, y = 8.25, xend = trial5_b6i0cmax, yend = 8.25, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))


# Barn 2 ----
# # Calculating peaks
# # Rectal Temperature
# trial5_b2i1 <- trial5_b2[trial5_b2$innoc == 1 & trial5_b2$tempc>=36 & trial5_b2$tempc<=42,]
# trial5_b2i1_fit <- loess(tempc ~ dpi, trial5_b2i1)
# trial5_b2i1_nd <- data.frame(dpi=seq(min(trial5_b2i1$dpi), max(trial5_b2i1$dpi), length=100))
# trial5_b2i1_nd$fit <- predict(trial5_b2i1_fit, newdata=trial5_b2i1_nd) 
# trial5_b2i1tmax <- trial5_b2i1_nd$dpi[which.max(trial5_b2i1_nd$fit)]
# trial5_b2[trial5_b2$innoc == 1,]$color <- "cyan4"
# # plot(trial5_b2i1$tempc~trial5_b2i1$dpi)
# # lines(trial5_b2i1_nd$fit~trial5_b2i1_nd$dpi, col="grey")
# 
# # Sentinels are two species
# trial5_b2i0c <- trial5_b2[trial5_b2$innoc == 0 & trial5_b2$tempc>=36 & trial5_b2$tempc<=42 & trial5_b2$species == "Caprine",]
# trial5_b2i0c_fit <- loess(tempc ~ dpi, trial5_b2i0c)
# trial5_b2i0c_nd <- data.frame(dpi=seq(min(trial5_b2i0c$dpi), max(trial5_b2i0c$dpi), length=100))
# trial5_b2i0c_nd$fit <- predict(trial5_b2i0c_fit, newdata=trial5_b2i0c_nd) 
# trial5_b2i0ctmax <- trial5_b2i0c_nd$dpi[which.max(trial5_b2i0c_nd$fit)]
# # plot(trial5_b2i0c$tempc~trial5_b2i0c$dpi)
# # lines(trial5_b2i0c_nd$fit~trial5_b2i0c_nd$dpi, col="grey")
# 
# trial5_b2i0b <- trial5_b2[trial5_b2$innoc == 0 & trial5_b2$tempc>=36 & trial5_b2$tempc<=42 & trial5_b2$species == "Bovine",]
# trial5_b2i0b_fit <- loess(tempc ~ dpi, trial5_b2i0b)
# trial5_b2i0b_nd <- data.frame(dpi=seq(min(trial5_b2i0b$dpi), max(trial5_b2i0b$dpi), length=100))
# trial5_b2i0b_nd$fit <- predict(trial5_b2i0b_fit, newdata=trial5_b2i0b_nd) 
# trial5_b2i0btmax <- trial5_b2i0b_nd$dpi[which.max(trial5_b2i0b_nd$fit)]
# # plot(trial5_b2i0b$tempc~trial5_b2i0b$dpi)
# # lines(trial5_b2i0b_nd$fit~trial5_b2i0b_nd$dpi, col="grey")
# 
# # # Clinical Score
# trial5_b2i1c <- trial5_b2[trial5_b2$innoc == 1,]
# trial5_b2i1c_fit <- loess(autosum ~ dpi, trial5_b2i1c)
# trial5_b2i1c_nd <- data.frame(dpi=seq(min(trial5_b2i1c$dpi), max(trial5_b2i1c$dpi), length=100))
# trial5_b2i1c_nd$fit <- predict(trial5_b2i1c_fit, newdata=trial5_b2i1c_nd) 
# trial5_b2i1cmax <- trial5_b2i1c_nd$dpi[which.max(trial5_b2i1c_nd$fit)]
# # plot(trial5_b2i1c$autosum~trial5_b2i1c$dpi)
# # lines(trial5_b2i1c_nd$fit~trial5_b2i1c_nd$dpi, col="grey")
# 
# # Sentinels are two species
# trial5_b2i0cc <- trial5_b2[trial5_b2$innoc == 0 & trial5_b2$species == "Caprine",]
# trial5_b2i0cc_fit <- loess(autosum ~ dpi, trial5_b2i0cc)
# trial5_b2i0cc_nd <- data.frame(dpi=seq(min(trial5_b2i0cc$dpi), max(trial5_b2i0cc$dpi), length=100))
# trial5_b2i0cc_nd$fit <- predict(trial5_b2i0cc_fit, newdata=trial5_b2i0cc_nd) 
# trial5_b2i0ccmax <- trial5_b2i0cc_nd$dpi[which.max(trial5_b2i0cc_nd$fit)]
# # plot(trial5_b2i0cc$autosum~trial5_b2i0cc$dpi)
# # lines(trial5_b2i0cc_nd$fit~trial5_b2i0cc_nd$dpi, col="grey")
# 
# trial5_b2i0cb <- trial5_b2[trial5_b2$innoc == 0 & trial5_b2$species == "Bovine",]
# trial5_b2i0cb_fit <- loess(autosum ~ dpi, trial5_b2i0cb)
# trial5_b2i0cb_nd <- data.frame(dpi=seq(min(trial5_b2i0cb$dpi), max(trial5_b2i0cb$dpi), length=100))
# trial5_b2i0cb_nd$fit <- predict(trial5_b2i0cb_fit, newdata=trial5_b2i0cb_nd) 
# trial5_b2i0cbmax <- trial5_b2i0cb_nd$dpi[which.max(trial5_b2i0cb_nd$fit)]
# # plot(trial5_b2i0cb$autosum~trial5_b2i0cb$dpi)
# # lines(trial5_b2i0cb_nd$fit~trial5_b2i0cb_nd$dpi, col="grey")

# Plots ----
# Rectal Temperature
rtC <- ggplot(trial5_b2, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Rectal Temperature (Celsius)") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial5_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(36,42)) +
  # geom_vline(aes(xintercept = trial5_b2i1tmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial5_b2i0ctmax), color = "cyan3") +
  # geom_vline(aes(xintercept = trial5_b2i0btmax), color = "mediumpurple1") +
  # annotate(geom = "text", x = (trial5_b2i1tmax + 1.75), y = 42, label = paste(round(trial5_b2i1tmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial5_b2i0ctmax + 1.75), y = 42, label = paste(round(trial5_b2i0ctmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = (trial5_b2i0btmax + 1.75), y = 42, label = paste(round(trial5_b2i0btmax,1)), color = "mediumpurple1") +
  # annotate(geom = "text", x = mean(c(trial5_b2i1tmax, trial5_b2i0ctmax)), y = 41.35, label = paste(round(trial5_b2i0ctmax - trial5_b2i1tmax,1)), color = "black") +
  # annotate(geom = "text", x = mean(c(trial5_b2i1tmax, trial5_b2i0btmax)), y = 40.55, label = paste(round(trial5_b2i0btmax - trial5_b2i1tmax,1)), color = "black") +
  # geom_segment(x = trial5_b2i1tmax, y = 41, xend = trial5_b2i0ctmax, yend = 41, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  # geom_segment(x = trial5_b2i1tmax, y = 40.2, xend = trial5_b2i0btmax, yend = 40.2, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))
# Removed 16 rows - the 1 dead sentinel animal, eartag 181 - died by dpi 20 so has 16 missing temp measurements for row.
#What did this animal die of?

# Clinical Score
csC <- ggplot(trial5_b2, aes(dpi, autosum, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)", "Cattle (In Contact)"),
                     values = trial5_cols) +
  scale_y_continuous(limits = c(-1,10), breaks = seq(-1, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  # geom_vline(aes(xintercept = trial5_b2i1cmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial5_b2i0ccmax), color = "cyan3") +
  # geom_vline(aes(xintercept = trial5_b2i0cbmax), color = "mediumpurple1") +
  # annotate(geom = "text", x = (trial5_b2i1cmax + 1.75), y = 10, label = paste(round(trial5_b2i1cmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial5_b2i0ccmax + 3.75), y = 10, label = paste(round(trial5_b2i0ccmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = (trial5_b2i0cbmax + 1.75), y = 10, label = paste(round(trial5_b2i0cbmax,1)), color = "mediumpurple1") +
  # annotate(geom = "text", x = mean(c(trial5_b2i1cmax, trial5_b2i0ccmax)), y = 9, label = paste(round(trial5_b2i0ccmax - trial5_b2i1cmax,1)), color = "black") +
  # annotate(geom = "text", x = mean(c(trial5_b2i1cmax, trial5_b2i0cbmax)), y = 7.5, label = paste(round(trial5_b2i0cbmax - trial5_b2i1cmax,1)), color = "black") +
  # geom_segment(x = trial5_b2i1cmax, y = 8.5, xend = trial5_b2i0ccmax, yend = 8.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  # geom_segment(x = trial5_b2i1cmax, y = 7, xend = trial5_b2i0cbmax, yend = 7, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))


# Barn 1----
# # Calculating peaks
# # Rectal Temperature
# trial5_b1i1 <- trial5_b1[trial5_b1$innoc == 1 & trial5_b1$tempc>=36 & trial5_b1$tempc<=42,]
# trial5_b1i1_fit <- loess(tempc ~ dpi, trial5_b1i1)
# trial5_b1i1_nd <- data.frame(dpi=seq(min(trial5_b1i1$dpi), max(trial5_b1i1$dpi), length=100))
# trial5_b1i1_nd$fit <- predict(trial5_b1i1_fit, newdata=trial5_b1i1_nd) 
# trial5_b1i1tmax <- trial5_b1i1_nd$dpi[which.max(trial5_b1i1_nd$fit)]
# trial5_b1[trial5_b1$innoc == 1,]$color <- "cyan4"
# # plot(trial5_b1i1$tempc~trial5_b1i1$dpi)
# # lines(trial5_b1i1_nd$fit~trial5_b1i1_nd$dpi, col="grey")
# 
# # Sentinels are two species
# trial5_b1i0c <- trial5_b1[trial5_b1$innoc == 0 & trial5_b1$tempc>=36 & trial5_b1$tempc<=42 & trial5_b1$species == "Caprine",]
# trial5_b1i0c_fit <- loess(tempc ~ dpi, trial5_b1i0c)
# trial5_b1i0c_nd <- data.frame(dpi=seq(min(trial5_b1i0c$dpi), max(trial5_b1i0c$dpi), length=100))
# trial5_b1i0c_nd$fit <- predict(trial5_b1i0c_fit, newdata=trial5_b1i0c_nd) 
# trial5_b1i0ctmax <- trial5_b1i0c_nd$dpi[which.max(trial5_b1i0c_nd$fit)]
# # plot(trial5_b1i0c$tempc~trial5_b1i0c$dpi)
# # lines(trial5_b1i0c_nd$fit~trial5_b1i0c_nd$dpi, col="grey")
# 
# trial5_b1i0b <- trial5_b1[trial5_b1$innoc == 0 & trial5_b1$tempc>=36 & trial5_b1$tempc<=42 & trial5_b1$species == "Bovine",]
# trial5_b1i0b_fit <- loess(tempc ~ dpi, trial5_b1i0b)
# trial5_b1i0b_nd <- data.frame(dpi=seq(min(trial5_b1i0b$dpi), max(trial5_b1i0b$dpi), length=100))
# trial5_b1i0b_nd$fit <- predict(trial5_b1i0b_fit, newdata=trial5_b1i0b_nd) 
# trial5_b1i0btmax <- trial5_b1i0b_nd$dpi[which.max(trial5_b1i0b_nd$fit)]
# # plot(trial5_b1i0b$tempc~trial5_b1i0b$dpi)
# # lines(trial5_b1i0b_nd$fit~trial5_b1i0b_nd$dpi, col="grey")
# 
# # # Clinical Score
# trial5_b1i1c <- trial5_b1[trial5_b1$innoc == 1,]
# trial5_b1i1c_fit <- loess(autosum ~ dpi, trial5_b1i1c)
# trial5_b1i1c_nd <- data.frame(dpi=seq(min(trial5_b1i1c$dpi), max(trial5_b1i1c$dpi), length=100))
# trial5_b1i1c_nd$fit <- predict(trial5_b1i1c_fit, newdata=trial5_b1i1c_nd) 
# trial5_b1i1cmax <- trial5_b1i1c_nd$dpi[which.max(trial5_b1i1c_nd$fit)]
# # plot(trial5_b1i1c$autosum~trial5_b1i1c$dpi)
# # lines(trial5_b1i1c_nd$fit~trial5_b1i1c_nd$dpi, col="grey")
# 
# # Sentinels are two species
# trial5_b1i0cc <- trial5_b1[trial5_b1$innoc == 0 & trial5_b1$species == "Caprine",]
# trial5_b1i0cc_fit <- loess(autosum ~ dpi, trial5_b1i0cc)
# trial5_b1i0cc_nd <- data.frame(dpi=seq(min(trial5_b1i0cc$dpi), max(trial5_b1i0cc$dpi), length=100))
# trial5_b1i0cc_nd$fit <- predict(trial5_b1i0cc_fit, newdata=trial5_b1i0cc_nd) 
# trial5_b1i0ccmax <- trial5_b1i0cc_nd$dpi[which.max(trial5_b1i0cc_nd$fit)]
# # plot(trial5_b1i0cc$autosum~trial5_b1i0cc$dpi)
# # lines(trial5_b1i0cc_nd$fit~trial5_b1i0cc_nd$dpi, col="grey")
# 
# trial5_b1i0cb <- trial5_b1[trial5_b1$innoc == 0 & trial5_b1$species == "Bovine",]
# trial5_b1i0cb_fit <- loess(autosum ~ dpi, trial5_b1i0cb)
# trial5_b1i0cb_nd <- data.frame(dpi=seq(min(trial5_b1i0cb$dpi), max(trial5_b1i0cb$dpi), length=100))
# trial5_b1i0cb_nd$fit <- predict(trial5_b1i0cb_fit, newdata=trial5_b1i0cb_nd) 
# trial5_b1i0cbmax <- trial5_b1i0cb_nd$dpi[which.max(trial5_b1i0cb_nd$fit)]
# # plot(trial5_b1i0cb$autosum~trial5_b1i0cb$dpi)
# # lines(trial5_b1i0cb_nd$fit~trial5_b1i0cb_nd$dpi, col="grey")

# Plots ----
# Rectal Temperature
rtD <- ggplot(trial5_b1, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Rectal Temperature (Celsius)") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial5_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(36,42)) +
  # geom_vline(aes(xintercept = trial5_b1i1tmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial5_b1i0ctmax), color = "cyan3") +
  # geom_vline(aes(xintercept = trial5_b1i0btmax), color = "mediumpurple1") +
  # annotate(geom = "text", x = (trial5_b1i1tmax + 1.75), y = 42, label = paste(round(trial5_b1i1tmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial5_b1i0ctmax + 1.75), y = 42, label = paste(round(trial5_b1i0ctmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = (trial5_b1i0btmax + 1.75), y = 42, label = paste(round(trial5_b1i0btmax,1)), color = "mediumpurple1") +
  # annotate(geom = "text", x = mean(c(trial5_b1i1tmax, trial5_b1i0ctmax)), y = 41.35, label = paste(round(trial5_b1i0ctmax - trial5_b1i1tmax,1)), color = "black") +
  # annotate(geom = "text", x = mean(c(trial5_b1i1tmax, trial5_b1i0btmax)), y = 40.55, label = paste(round(trial5_b1i0btmax - trial5_b1i1tmax,1)), color = "black") +
  # geom_segment(x = trial5_b1i1tmax, y = 41, xend = trial5_b1i0ctmax, yend = 41, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  # geom_segment(x = trial5_b1i1tmax, y = 40.2, xend = trial5_b1i0btmax, yend = 40.2, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))

# Clinical Score
csD <- ggplot(trial5_b1, aes(dpi, autosum, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)", "Cattle (In Contact)"),
                     values = trial5_cols) +
  scale_y_continuous(limits = c(-1,10), breaks = seq(-1, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  # geom_vline(aes(xintercept = trial5_b1i1cmax), color = "cyan4") +
  # geom_vline(aes(xintercept = trial5_b1i0ccmax), color = "cyan3") +
  # geom_vline(aes(xintercept = trial5_b1i0cbmax), color = "mediumpurple1") +
  # annotate(geom = "text", x = (trial5_b1i1cmax + 1.75), y = 10, label = paste(round(trial5_b1i1cmax,1)), color = "cyan4") +
  # annotate(geom = "text", x = (trial5_b1i0ccmax + 3.75), y = 10, label = paste(round(trial5_b1i0ccmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = (trial5_b1i0cbmax + 1.75), y = 10, label = paste(round(trial5_b1i0cbmax,1)), color = "mediumpurple1") +
  # annotate(geom = "text", x = mean(c(trial5_b1i1cmax, trial5_b1i0ccmax)), y = 9, label = paste(round(trial5_b1i0ccmax - trial5_b1i1cmax,1)), color = "black") +
  # annotate(geom = "text", x = mean(c(trial5_b1i1cmax, trial5_b1i0cbmax)), y = 7.5, label = paste(round(trial5_b1i0cbmax - trial5_b1i1cmax,1)), color = "black") +
  # geom_segment(x = trial5_b1i1cmax, y = 8.5, xend = trial5_b1i0ccmax, yend = 8.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  # geom_segment(x = trial5_b1i1cmax, y = 7, xend = trial5_b1i0cbmax, yend = 7, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))



# Barn 3 ----
# Rectal Temperature
rtpcont <- ggplot(trial5_b3, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(37,41)) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
# 1 of the 2 died after 7 dpi so 28 rows missing for rest of it

# Clinical Score - Modified
cspcont <- ggplot(trial5_b3, aes(dpi, clinscore_mod, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
# 1 of the 2 died after 7 dpi so 28 rows of 0 rest of it - so manually replaced with blank in csv (NA)


# Barn 4 ----
# Plots ----
# Rectal Temperature
rtncont <- ggplot(trial5_b4, aes(dpi, tempc, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Goats (Negative Control)"),
                     values = trial5_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(37,41)) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

# Clinical Score - Modified
csncont <- ggplot(trial5_b4, aes(dpi, clinscore_mod, color=barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Goats (Negative Control)"),
                     values = trial5_cols) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))


# All three experimental barns together ----
# Plots ----
# Rectal Temperature
rtABCD <- ggplot(trial5_b1256, aes(dpi, tempc, color=species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Rectal Temperature (Celsius)") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols_sp) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(36,42)) +
  theme_minimal() +
  theme(legend.position = "top")

# Clinical Score
csABCD <- ggplot(trial5_b1256, aes(dpi, autosum, color=species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols_sp) +
  scale_y_continuous(limits = c(-1,10), breaks = seq(-1, 10, 1)) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  theme_minimal() +
  theme(legend.position = "top")


# Arranging Plots ----
jpeg("output/clin/trial5_clin_b5.jpeg", width = 7, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtA, csA, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial5_clin_b6.jpeg", width = 7, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtB, csB, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial5_clin_b2.jpeg", width = 7, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtC, csC, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial5_clin_b1.jpeg", width = 7, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtC, csC, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial5_clin_b3_FigS12A.jpeg", width = 7, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtpcont, cspcont, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial5_clin_b4_FigS12B.jpeg", width = 7, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtncont, csncont, ncol=2)
invisible(dev.off())

jpeg("output/clin/trial5_clin_b1256.jpeg", width = 7, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rtABCD, csABCD, ncol=2)
invisible(dev.off())



###########################
# Combined Trials 1-2 ----
###########################

#remove controls + put datasets together
trial1_clin_nc <- subset(trial1_clin, barn != 1) 
trial2_clin_nc <- subset(trial2_clin, barn != 5) %>% select(-X, -X.1, -X.2, -X.3)

jointdat12 <- rbind(trial1_clin_nc, trial2_clin_nc)
jointdat12 <- jointdat12 %>% mutate(barniv2 = case_when(barni == "5a" ~ "Gi1",
                                                    barni == "4a" ~ "Gi1",
                                                    barni == "4b" ~ "Gi0",
                                                    barni == "3a" ~ "Si1",
                                                    barni == "6a" ~ "Si1",
                                                    barni == "6b" ~ "Si0"))
comb12_cols <- c("Gi1" = "cyan4", "Gi0" = "cyan3",
                 "Si1" = "deeppink2", "Si0" = "rosybrown3")


# PLOTS ----
#getting peaks
#temp
gpi1 <- jointdat12[jointdat12$barniv2 == "Gi1" & jointdat12$tempc>=36 & jointdat12$tempc<=42,]
gpi1 <- gpi1[!is.na(gpi1$eartag),]
gpi1_fit <- loess(tempc ~ dpi, gpi1)
gpi1_nd <- data.frame(dpi=seq(min(gpi1$dpi), max(gpi1$dpi), length=100))
gpi1_nd$fit <- predict(gpi1_fit, newdata=gpi1_nd) 
gpi1_tmax <- gpi1_nd$dpi[which.max(gpi1_nd$fit)]
gpi1[gpi1$barniv2 == "Gi1",]$color <- "cyan4"
# plot(gpi1$tempc~gpi1$dpi)
# lines(gpi1_nd$fit~gpi1_nd$dpi, col="grey")

gpi0 <- jointdat12[jointdat12$barniv2 == "Gi0" & jointdat12$tempc>=36 & jointdat12$tempc<=42,]
gpi0 <- gpi0[!is.na(gpi0$eartag),]
gpi0_fit <- loess(tempc ~ dpi, gpi0)
gpi0_nd <- data.frame(dpi=seq(min(gpi0$dpi), max(gpi0$dpi), length=100))
gpi0_nd$fit <- predict(gpi0_fit, newdata=gpi0_nd) 
gpi0_tmax <- gpi0_nd$dpi[which.max(gpi0_nd$fit)]
gpi0[gpi0$barniv2 == "Gi0",]$color <- "cyan4"
# plot(gpi0$tempc~gpi0$dpi)
# lines(gpi0_nd$fit~gpi0_nd$dpi, col="grey")

spi1 <- jointdat12[jointdat12$barniv2 == "Si1" & jointdat12$tempc>=36 & jointdat12$tempc<=42,]
spi1_fit <- loess(tempc ~ dpi, spi1)
spi1_nd <- data.frame(dpi=seq(min(spi1$dpi), max(spi1$dpi), length=100))
spi1_nd$fit <- predict(spi1_fit, newdata=spi1_nd) 
spi1_tmax <- spi1_nd$dpi[which.max(spi1_nd$fit)]
spi1[spi1$barniv2 == "Si1",]$color <- "deeppink2"
# plot(spi1$tempc~spi1$dpi)
# lines(spi1_nd$fit~spi1_nd$dpi, col="grey")

spi0 <- jointdat12[jointdat12$barniv2 == "Si0" & jointdat12$tempc>=36 & jointdat12$tempc<=42,]
spi0_fit <- loess(tempc ~ dpi, spi0)
spi0_nd <- data.frame(dpi=seq(min(spi0$dpi), max(spi0$dpi), length=100))
spi0_nd$fit <- predict(spi0_fit, newdata=spi0_nd) 
spi0_tmax <- spi0_nd$dpi[which.max(spi0_nd$fit)]
spi0[spi0$barniv2 == "Si0",]$color <- "deeppink2"
# plot(spi0$tempc~spi0$dpi)
# lines(spi0_nd$fit~spi0_nd$dpi, col="grey")

#clinical score
gpi1c <- jointdat12[jointdat12$barniv2 == "Gi1" & jointdat12$innoc == 1,]
gpi1c_fit <- loess(autosum ~ dpi, gpi1c)
gpi1c_nd <- data.frame(dpi=seq(min(gpi1c$dpi), max(gpi1c$dpi), length=100))
gpi1c_nd$fit <- predict(gpi1c_fit, newdata=gpi1c_nd) 
gpi1_cmax <- gpi1c_nd$dpi[which.max(gpi1c_nd$fit)]
gpi1c[gpi1c$barniv2 == "Gi1",]$color <- "cyan4"
# plot(gpi1$autosum~gpi1c$dpi)
# lines(gpi1c_nd$fit~gpi1c_nd$dpi, col="grey")

gpi0c <- jointdat12[jointdat12$barniv2 == "Gi0" & jointdat12$innoc == 0,]
gpi0c_fit <- loess(autosum ~ dpi, gpi0c)
gpi0c_nd <- data.frame(dpi=seq(min(gpi0c$dpi), max(gpi0c$dpi), length=100))
gpi0c_nd$fit <- predict(gpi0c_fit, newdata=gpi0c_nd) 
gpi0_cmax <- gpi0c_nd$dpi[which.max(gpi0c_nd$fit)]
gpi0c[gpi0c$barniv2 == "Gi0",]$color <- "cyan4"
# plot(gpi0$autosum~gpi0c$dpi)
# lines(gpi0c_nd$fit~gpi0c_nd$dpi, col="grey")

spi1c <- jointdat12[jointdat12$barniv2 == "Si1" & jointdat12$innoc == 1,]
spi1c_fit <- loess(autosum ~ dpi, spi1c)
spi1c_nd <- data.frame(dpi=seq(min(spi1c$dpi), max(spi1c$dpi), length=100))
spi1c_nd$fit <- predict(spi1c_fit, newdata=spi1c_nd) 
spi1_cmax <- spi1c_nd$dpi[which.max(spi1c_nd$fit)]
spi1c[spi1c$barniv2 == "Si1",]$color <- "deeppink2"
# plot(spi1c$autosum~spi1c$dpi)
# lines(spi1c_nd$fit~spic_nd$dpi, col="grey")

spi0c <- jointdat12[jointdat12$barniv2 == "Si0"& jointdat12$innoc == 0,]
spi0c_fit <- loess(autosum ~ dpi, spi0c)
spi0c_nd <- data.frame(dpi=seq(min(spi0c$dpi), max(spi0c$dpi), length=100))
spi0c_nd$fit <- predict(spi0c_fit, newdata=spi0c_nd) 
spi0_cmax <- spi0c_nd$dpi[which.max(spi0c_nd$fit)]
spi0c[spi0c$barniv2 == "Si0",]$color <- "deeppink2"
# plot(spi0c$autosum~spi0c$dpi)
# lines(spi0c_nd$fit~spi0c_nd$dpi, col="grey")


#clinical score modified 
gpi1cm <- jointdat12[jointdat12$barniv2 == "Gi1" & jointdat12$innoc == 1,]
gpi1cm_fit <- loess(clinscore_mod ~ dpi, gpi1cm)
gpi1cm_nd <- data.frame(dpi=seq(min(gpi1cm$dpi), max(gpi1cm$dpi), length=100))
gpi1cm_nd$fit <- predict(gpi1cm_fit, newdata=gpi1cm_nd) 
gpi1_cmmax <- gpi1cm_nd$dpi[which.max(gpi1cm_nd$fit)]
gpi1cm[gpi1cm$barniv2 == "Gi1",]$color <- "cyan4"
# plot(gpi1$clinscore_mod~gpi1cm$dpi)
# lines(gpi1cm_nd$fit~gpi1cm_nd$dpi, col="grey")

gpi0cm <- jointdat12[jointdat12$barniv2 == "Gi0" & jointdat12$innoc == 0,]
gpi0cm_fit <- loess(clinscore_mod ~ dpi, gpi0cm)
gpi0cm_nd <- data.frame(dpi=seq(min(gpi0cm$dpi), max(gpi0cm$dpi), length=100))
gpi0cm_nd$fit <- predict(gpi0cm_fit, newdata=gpi0cm_nd) 
gpi0_cmmax <- gpi0cm_nd$dpi[which.max(gpi0cm_nd$fit)]
gpi0cm[gpi0cm$barniv2 == "Gi0",]$color <- "cyan4"
# plot(gpi0$clinscore_mod~gpi0cm$dpi)
# lines(gpi0cm_nd$fit~gpi0cm_nd$dpi, col="grey")

spi1cm <- jointdat12[jointdat12$barniv2 == "Si1" & jointdat12$innoc == 1,]
spi1cm_fit <- loess(clinscore_mod ~ dpi, spi1cm)
spi1cm_nd <- data.frame(dpi=seq(min(spi1cm$dpi), max(spi1cm$dpi), length=100))
spi1cm_nd$fit <- predict(spi1cm_fit, newdata=spi1cm_nd) 
spi1_cmmax <- spi1cm_nd$dpi[which.max(spi1cm_nd$fit)]
spi1cm[spi1cm$barniv2 == "Si1",]$color <- "deeppink2"
# plot(spi1cm$clinscore_mod~spi1cm$dpi)
# lines(spi1cm_nd$fit~spicm_nd$dpi, col="grey")

spi0cm <- jointdat12[jointdat12$barniv2 == "Si0"& jointdat12$innoc == 0,]
spi0cm_fit <- loess(clinscore_mod ~ dpi, spi0cm)
spi0cm_nd <- data.frame(dpi=seq(min(spi0cm$dpi), max(spi0cm$dpi), length=100))
spi0cm_nd$fit <- predict(spi0cm_fit, newdata=spi0cm_nd) 
spi0_cmmax <- spi0cm_nd$dpi[which.max(spi0cm_nd$fit)]
spi0cm[spi0cm$barniv2 == "Si0",]$color <- "deeppink2"
# plot(spi0cm$clinscore_mod~spi0cm$dpi)
# lines(spi0cm_nd$fit~spi0cm_nd$dpi, col="grey")


#PLOTS ----
#rectal temp
rt <- ggplot(jointdat12, aes(dpi, tempc, color=barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Goat (Sentinel)", "Goat (Inoculated)", "Sheep (Sentinel)", "Sheep (Inoculated)"),
                     values = comb12_cols) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0, 28, 2)) +
  scale_y_continuous(limits = c(37,41)) +
  geom_vline(aes(xintercept = spi1_tmax), color = "deeppink2") +
  geom_vline(aes(xintercept = spi0_tmax), color = "rosybrown2") +
  geom_vline(aes(xintercept = gpi1_tmax), color = "cyan4") +
  geom_vline(aes(xintercept = gpi0_tmax), color = "cyan3") +
  annotate(geom = "text", x = (spi1_tmax + 1.75), y = 41, label = paste(round(spi1_tmax,1)), color = "deeppink2") +
  annotate(geom = "text", x = (spi0_tmax - 1.75), y = 41, label = paste(round(spi0_tmax,1)), color = "rosybrown2") +
  annotate(geom = "text", x = (gpi1_tmax - 1.75), y = 41, label = paste(round(gpi1_tmax,1)), color = "cyan4") +
  annotate(geom = "text", x = (gpi0_tmax + 1.75), y = 41, label = paste(round(gpi0_tmax,1)), color = "cyan3") +
  annotate(geom = "text", x = mean(c(spi1_tmax, spi0_tmax)), y = 40.9, label = paste(round(spi0_tmax - spi1_tmax,1)), color = "black") +
  geom_segment(x = spi1_tmax, y = 40.8, xend = spi0_tmax, yend = 40.8, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  annotate(geom = "text", x = mean(c(gpi1_tmax, gpi0_tmax)), y = 40.6, label = paste(round(gpi0_tmax - gpi1_tmax,1)), color = "black") +
  geom_segment(x = gpi1_tmax, y = 40.5, xend = gpi0_tmax, yend = 40.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

#clinical score
cs <- ggplot(jointdat12, aes(dpi, autosum, color=barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Goat (Sentinel)", "Goat (Inoculated)", "Sheep (Sentinel)", "Sheep (Inoculated)"),
                     values = comb12_cols) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0, 28, 2)) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  geom_vline(aes(xintercept = spi1_cmax), color = "deeppink2") +
  geom_vline(aes(xintercept = spi0_cmax), color = "rosybrown2") +
  geom_vline(aes(xintercept = gpi1_cmax), color = "cyan4") +
  geom_vline(aes(xintercept = gpi0_cmax), color = "cyan3") +
  annotate(geom = "text", x = (spi1_cmax + 1.75), y = 10, label = paste(round(spi1_cmax,1)), color = "deeppink2") +
  annotate(geom = "text", x = (spi0_cmax - 1.75), y = 10, label = paste(round(spi0_cmax,1)), color = "rosybrown2") +
  annotate(geom = "text", x = (gpi1_cmax - 1.75), y = 10, label = paste(round(gpi1_cmax,1)), color = "cyan4") +
  annotate(geom = "text", x = (gpi0_cmax + 1.75), y = 10, label = paste(round(gpi0_cmax,1)), color = "cyan3") +
  annotate(geom = "text", x = mean(c(spi1_cmax, spi0_cmax)), y = 8.9, label = paste(round(spi0_cmax - spi1_cmax,1)), color = "black") +
  geom_segment(x = spi1_cmax, y = 9.15, xend = spi0_cmax, yend = 9.15, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  annotate(geom = "text", x = mean(c(gpi1_cmax, gpi0_cmax)), y = 9.65, label = paste(round(gpi0_cmax - gpi1_cmax,1)), color = "black") +
  geom_segment(x = gpi1_cmax, y = 9.35, xend = gpi0_cmax, yend = 9.35, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))


# modified clinical score - no temp
csm <- ggplot(jointdat12, aes(dpi, clinscore_mod, color=barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Goat (Sentinel)", "Goat (Inoculated)", "Sheep (Sentinel)", "Sheep (Inoculated)"),
                     values = comb12_cols) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0, 28, 2)) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  geom_vline(aes(xintercept = spi1_cmmax), color = "deeppink2") +
  geom_vline(aes(xintercept = spi0_cmmax), color = "rosybrown2") +
  geom_vline(aes(xintercept = gpi1_cmmax), color = "cyan4") +
  geom_vline(aes(xintercept = gpi0_cmmax), color = "cyan3") +
  annotate(geom = "text", x = (spi1_cmmax + 1.75), y = 10, label = paste(round(spi1_cmmax,1)), color = "deeppink2") +
  annotate(geom = "text", x = (spi0_cmmax + 1.75), y = 10, label = paste(round(spi0_cmmax,1)), color = "rosybrown2") +
  annotate(geom = "text", x = (gpi1_cmmax - 1.75), y = 10, label = paste(round(gpi1_cmmax,1)), color = "cyan4") +
  annotate(geom = "text", x = (gpi0_cmmax - 1.75), y = 10, label = paste(round(gpi0_cmmax,1)), color = "cyan3") +
  annotate(geom = "text", x = mean(c(spi1_cmmax, spi0_cmmax)), y = 8.9, label = paste(round(spi0_cmmax - spi1_cmmax,1)), color = "black") +
  geom_segment(x = spi1_cmmax, y = 9.15, xend = spi0_cmmax, yend = 9.15, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  annotate(geom = "text", x = mean(c(gpi1_cmmax, gpi0_cmmax)), y = 9.65, label = paste(round(gpi0_cmmax - gpi1_cmmax,1)), color = "black") +
  geom_segment(x = gpi1_cmmax, y = 9.35, xend = gpi0_cmmax, yend = 9.35, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))


jpeg("output/clin/trials12_clinicalscorecomparison_FigS6.jpeg", width = 12, height = 5, units = "in", quality = 100, res = 600)
rt2 <- rt + labs(x = "")
cs2 <- cs + labs(x = "")
csm2 <- csm + labs(x = "")
grid.arrange(rt2, cs2, csm2, ncol=3, bottom = textGrob("Days Post Infection", gp=gpar(fontsize=15)))
invisible(dev.off())


###########################
#  Combined Trials 3-5 ----
###########################

#remove controls + any other barns put datasets together
trial3_clin_nc <- subset(trial3_clin, barn == 3) %>% select(-X, -X.1, -X.2)
trial4_clin_nc <- subset(trial4_clin, barn != 5 & barn != 3) 
trial5_clin_nc <- subset(trial5_clin, barn != 3 & barn != 4) 

jointdat <- rbind(trial3_clin_nc, trial4_clin_nc, trial5_clin_nc)
jointdat <- jointdat %>% mutate(barniv2 = case_when(barni == "3a" ~ "Ci1",
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

comb345_cols <- c("Ci1" = "purple4", "Gi0" = "cyan3")

# Calculating peaks
# Rectal Temperature
ci1 <- jointdat[jointdat$barniv2 == "Ci1" & jointdat$tempc>=36 & jointdat$tempc<=42,]
ci1 <- ci1[!is.na(ci1$eartag),]
ci1_fit <- loess(tempc ~ dpi, ci1)
ci1_nd <- data.frame(dpi=seq(min(ci1$dpi), max(ci1$dpi), length=100))
ci1_nd$fit <- predict(ci1_fit, newdata=ci1_nd) 
ci1_tmax <- ci1_nd$dpi[which.max(ci1_nd$fit)]
ci1[ci1$barniv2 == "Ci1",]$color <- "purple4"
# plot(ci1$tempc~ci1$dpi)
# lines(ci1_nd$fit~ci1_nd$dpi, col="grey")

gpi0 <- jointdat[jointdat$barniv2 == "Gi0" & jointdat$tempc>=36 & jointdat$tempc<=42,]
gpi0 <- gpi0[!is.na(gpi0$eartag),]
gpi0_fit <- loess(tempc ~ dpi, gpi0)
gpi0_nd <- data.frame(dpi=seq(min(gpi0$dpi), max(gpi0$dpi), length=100))
gpi0_nd$fit <- predict(gpi0_fit, newdata=gpi0_nd) 
gpi0_tmax <- gpi0_nd$dpi[which.max(gpi0_nd$fit)]
gpi0[gpi0$barniv2 == "Gi0",]$color <- "cyan3"
# plot(gpi0$tempc~gpi0$dpi)
# lines(gpi0_nd$fit~gpi0_nd$dpi, col="grey")

# Clinical Score - Modified 
ci1cm <- jointdat[jointdat$barniv2 == "Ci1" & jointdat$innoc == 1,]
ci1cm_fit <- loess(clinscore_mod ~ dpi, ci1cm)
ci1cm_nd <- data.frame(dpi=seq(min(ci1cm$dpi), max(ci1cm$dpi), length=100))
ci1cm_nd$fit <- predict(ci1cm_fit, newdata=ci1cm_nd) 
ci1_cmmax <- ci1cm_nd$dpi[which.max(ci1cm_nd$fit)]
ci1cm[ci1cm$barniv2 == "Ci1",]$color <- "purple4"
# plot(gpi1$clinscore_mod~ci1cm$dpi)
# lines(ci1cm_nd$fit~ci1cm_nd$dpi, col="grey")

gpi0cm <- jointdat[jointdat$barniv2 == "Gi0" & jointdat$innoc == 0,]
gpi0cm_fit <- loess(clinscore_mod ~ dpi, gpi0cm)
gpi0cm_nd <- data.frame(dpi=seq(min(gpi0cm$dpi), max(gpi0cm$dpi), length=100))
gpi0cm_nd$fit <- predict(gpi0cm_fit, newdata=gpi0cm_nd) 
gpi0_cmmax <- gpi0cm_nd$dpi[which.max(gpi0cm_nd$fit)]
gpi0cm[gpi0cm$barniv2 == "Gi0",]$color <- "cyan3"
# plot(gpi0$clinscore_mod~gpi0cm$dpi)
# lines(gpi0cm_nd$fit~gpi0cm_nd$dpi, col="grey")

# Plots ----
# Rectal Temperature
rt <- ggplot(jointdat, aes(dpi, tempc, color=barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.1) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), size =1) +
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goat (Sentinel)"),
                     values = comb345_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(37,41)) +
  # geom_vline(aes(xintercept = ci1_tmax), color = "purple4") +
  # geom_vline(aes(xintercept = gpi0_tmax), color = "cyan3") +
  # annotate(geom = "text", x = (ci1_tmax + 1.75), y = 42, label = paste(round(ci1_tmax,1)), color = "purple4") +
  # annotate(geom = "text", x = (gpi0_tmax + 1.75), y = 42, label = paste(round(gpi0_tmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = mean(c(ci1_tmax, gpi0_tmax)), y = 41.05, label = paste(round(gpi0_tmax - ci1_tmax,1)), color = "black") +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
 

# Clinical Score - Modified
csm <- ggplot(jointdat, aes(dpi, clinscore_mod, color=barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), size =1) +
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goat (Sentinel)"),
                     values = comb345_cols) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  # geom_vline(aes(xintercept = ci1_cmmax), color = "purple4") +
  # geom_vline(aes(xintercept = gpi0_cmmax), color = "cyan3") +
  # annotate(geom = "text", x = (ci1_cmmax + 1.75), y = 10, label = paste(round(ci1_cmmax,1)), color = "purple4") +
  # annotate(geom = "text", x = (gpi0_cmmax - 1.75), y = 10, label = paste(round(gpi0_cmmax,1)), color = "cyan3") +
  # annotate(geom = "text", x = mean(c(ci1_cmmax, gpi0_cmmax)), y = 8.9, label = paste(round(gpi0_cmmax - ci1_cmmax,1)), color = "black") +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

# Get Legend Function
# get_legend<-function(myggplot){
#   tmp <- ggplot_gtable(ggplot_build(myggplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)
# }
# legend <- get_legend(rt)
# rt <- rt + theme(legend.position = "none")

# Arranging Plots ----
jpeg("output/clin/trials345_clin_Fig2AB.jpeg", width = 12, height = 5, units = "in", quality = 100, res = 600)
grid.arrange(rt, csm, ncol=2)
invisible(dev.off())
# jpeg("output/clin/trials345_clin_Fig2AB_legend.jpeg", width = 12, height = 5, units = "in", quality = 100, res = 600)
# grid.arrange(legend, ncol=1)
# invisible(dev.off())



