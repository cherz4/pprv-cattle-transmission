#############################################################################################################
# Empirical and model-based evidence for a negligible role of cattle in 
# peste des petits ruminants transmission and eradication

# Catherine M. Herzog#*, Fasil Aklilu#, Demeke Sibhatu, Dereje Shegu, 
# Redeat Belaineh, Abde Aliy Mohammed, Menbere Kidane,  Claudia Schulz,
# Brian J. Willett, Sarah Cleaveland, Dalan Bailey, Andrew R. Peters, 
# Isabella M. Cattadori, Peter J. Hudson, Hagos Asgedom, Joram Buza, 
# Mesfin Sahle Forza, Tesfaye Rufael Chibssa, Solomon Gebredufe, Nick Juleff, 
# Ottar N. Bjørnstad, Michael D. Baron, Vivek Kapur*

# Code for: analyze and visualize molecular data for PPRV transmission trials
# Code written by: Catherine M. Herzog, PhD MPH

#############################################################################################################


######################
#  Load libraries ----
######################

library(ggplot2)   # for plotting figures
library(tidyverse) # for piping %>%
library(gridExtra) # for arranging several plots
library(grid) # for textGrob function in plots

############################################################
#  Molecular & Isolation  Data Import, Data Management  ----
############################################################

# Load the data by pointing to the location within the repo/project that this data resides (data folder)
trial1_molec <- read.csv("data/molec/trial1_molec.csv", header = TRUE) %>% 
                select(-type) %>%
                mutate_at(c("eartag", "sex", "species" , "barn"), as.factor)
trial1_molec$o1[trial1_molec$o1 == "U"] <- 35
trial1_molec$o2[trial1_molec$o2 == "U"] <- 35
trial1_molec$n1[trial1_molec$n1 == "U"] <- 35
trial1_molec$n2[trial1_molec$n2 == "U"] <- 35
trial1_molec$r1[trial1_molec$r1 == "U"] <- 35
trial1_molec$r2[trial1_molec$r2 == "U"] <- 35
trial1_molec$wb1[trial1_molec$wb1 == "U"] <- 35
trial1_molec$wb2[trial1_molec$wb2 == "U"] <- 35
trial1_molec$omean[trial1_molec$omean == "#DIV/0!"] <- "35"
trial1_molec$nmean[trial1_molec$nmean == "#DIV/0!"] <- "35"
trial1_molec$rmean[trial1_molec$rmean == "#DIV/0!"] <- "35"
trial1_molec$wbmean[trial1_molec$wbmean == "#DIV/0!"] <- "35"
trial1_molec <- trial1_molec %>% mutate_at(c("o1","o2","n1","n2","r1","r2","wb1","wb2","omean","nmean","rmean","wbmean"), as.numeric)

trial2_molec <- read.csv("data/molec/trial2_molec.csv", header = TRUE) %>% 
                select(-type) %>%
                mutate_at(c("eartag", "sex", "species" , "barn", "innoc", "barni", "status"), as.factor)
trial2_molec$o1[trial2_molec$o1 == "U"] <- 35
trial2_molec$o2[trial2_molec$o2 == "U"] <- 35
trial2_molec$n1[trial2_molec$n1 == "U"] <- 35
trial2_molec$n2[trial2_molec$n2 == "U"] <- 35
trial2_molec$r1[trial2_molec$r1 == "U"] <- 35
trial2_molec$r2[trial2_molec$r2 == "U"] <- 35
trial2_molec$wb1[trial2_molec$wb1 == "U"] <- 35
trial2_molec$wb2[trial2_molec$wb2 == "U"] <- 35
trial2_molec$omean[trial2_molec$omean == "#DIV/0!"] <- "35"
trial2_molec$nmean[trial2_molec$nmean == "#DIV/0!"] <- "35"
trial2_molec$rmean[trial2_molec$rmean == "#DIV/0!"] <- "35"
trial2_molec$wbmean[trial2_molec$wbmean == "#DIV/0!"] <- "35"
trial2_molec <- trial2_molec %>% mutate_at(c("o1","o2","n1","n2","r1","r2","wb1","wb2","omean","nmean","rmean","wbmean"), as.numeric) # WARNING


trial3_molec <- read.csv("data/molec/trial3_molec.csv", header = TRUE)%>% 
                mutate_at(c("eartag", "species" , "barn"), as.factor)
trial3_molec$o1[trial3_molec$o1 == "U"] <- 35
trial3_molec$o2[trial3_molec$o2 == "U"] <- 35
trial3_molec$n1[trial3_molec$n1 == "U"] <- 35
trial3_molec$n2[trial3_molec$n2 == "U"] <- 35
trial3_molec$r1[trial3_molec$r1 == "U"] <- 35
trial3_molec$r2[trial3_molec$r2 == "U"] <- 35
trial3_molec$wb1[trial3_molec$wb1 == "U"] <- 35
trial3_molec$wb2[trial3_molec$wb2 == "U"] <- 35
trial3_molec$omean[trial3_molec$omean == "#DIV/0!"] <- "35"
trial3_molec$nmean[trial3_molec$nmean == "#DIV/0!"] <- "35"
trial3_molec$rmean[trial3_molec$rmean == "#DIV/0!"] <- "35"
trial3_molec$wbmean[trial3_molec$wbmean == "#DIV/0!"] <- "35"
trial3_molec <- trial3_molec %>% mutate_at(c("o1","o2","n1","n2","r1","r2","wb1","wb2","omean","nmean","rmean","wbmean"), as.numeric)

trial3_chall_molec <- read.csv("data/molec/trial3_challenge_molec.csv", header = TRUE) %>% 
                mutate_at(c("eartag", "species" , "barn"), as.factor)
trial3_chall_molec$o1[trial3_chall_molec$o1 == "U"] <- 35
trial3_chall_molec$o2[trial3_chall_molec$o2 == "U"] <- 35
trial3_chall_molec$n1[trial3_chall_molec$n1 == "U"] <- 35
trial3_chall_molec$n2[trial3_chall_molec$n2 == "U"] <- 35
trial3_chall_molec$r1[trial3_chall_molec$r1 == "U"] <- 35
trial3_chall_molec$r2[trial3_chall_molec$r2 == "U"] <- 35
trial3_chall_molec$wb1[trial3_chall_molec$wb1 == "U"] <- 35
trial3_chall_molec$wb2[trial3_chall_molec$wb2 == "U"] <- 35
trial3_chall_molec$omean[trial3_chall_molec$omean == "#DIV/0!"] <- "35"
trial3_chall_molec$nmean[trial3_chall_molec$nmean == "#DIV/0!"] <- "35"
trial3_chall_molec$rmean[trial3_chall_molec$rmean == "#DIV/0!"] <- "35"
trial3_chall_molec$wbmean[trial3_chall_molec$wbmean == "#DIV/0!"] <- "35"
trial3_chall_molec <- trial3_chall_molec %>% mutate_at(c("o1","o2","n1","n2","r1","r2","wb1","wb2","omean","nmean","rmean","wbmean"), as.numeric)


trial4_molec <- read.csv("data/molec/trial4_molec.csv", header = TRUE)%>% 
                mutate_at(c("eartag", "species" , "barn"), as.factor)
trial4_molec$o1[trial4_molec$o1 == "U"] <- 35
trial4_molec$o2[trial4_molec$o2 == "U"] <- 35
trial4_molec$n1[trial4_molec$n1 == "U"] <- 35
trial4_molec$n2[trial4_molec$n2 == "U"] <- 35
trial4_molec$r1[trial4_molec$r1 == "U"] <- 35
trial4_molec$r2[trial4_molec$r2 == "U"] <- 35
trial4_molec$wb1[trial4_molec$wb1 == "U"] <- 35
trial4_molec$wb2[trial4_molec$wb2 == "U"] <- 35
trial4_molec$omean[trial4_molec$omean == "#DIV/0!"] <- "35"
trial4_molec$nmean[trial4_molec$nmean == "#DIV/0!"] <- "35"
trial4_molec$rmean[trial4_molec$rmean == "#DIV/0!"] <- "35"
trial4_molec$wbmean[trial4_molec$wbmean == "#DIV/0!"] <- "35"
trial4_molec <- trial4_molec %>% mutate_at(c("o1","o2","n1","n2","r1","r2","wb1","wb2","omean","nmean","rmean","wbmean"), as.numeric)


trial5_molec <- read.csv("data/molec/trial5_molec.csv", header = TRUE)%>% 
                mutate_at(c("eartag", "species" , "barn"), as.factor)
trial5_molec$o1[trial5_molec$o1 == "U"] <- 35
trial5_molec$o2[trial5_molec$o2 == "U"] <- 35
trial5_molec$n1[trial5_molec$n1 == "U"] <- 35
trial5_molec$n2[trial5_molec$n2 == "U"] <- 35
trial5_molec$r1[trial5_molec$r1 == "U"] <- 35
trial5_molec$r2[trial5_molec$r2 == "U"] <- 35
trial5_molec$wb1[trial5_molec$wb1 == "U"] <- 35
trial5_molec$wb2[trial5_molec$wb2 == "U"] <- 35
trial5_molec$omean[trial5_molec$omean == "#DIV/0!"] <- "35"
trial5_molec$nmean[trial5_molec$nmean == "#DIV/0!"] <- "35"
trial5_molec$rmean[trial5_molec$rmean == "#DIV/0!"] <- "35"
trial5_molec$wbmean[trial5_molec$wbmean == "#DIV/0!"] <- "35"
trial5_molec <- trial5_molec %>% mutate_at(c("o1","o2","n1","n2","r1","r2","wb1","wb2","omean","nmean","rmean","wbmean"), as.numeric)


#####################
#  Trial 1 ----
#####################

t1_b1 <- subset(trial1_molec, barn == 1)
t1_b3 <- subset(trial1_molec, barn == 3)
t1_b5 <- subset(trial1_molec, barn == 5)

trial1_cols <- c("1" = "black", "3" = "deeppink2", "5" = "cyan4") 

# SHEEP
# getting peaks
sop_fit <- loess(35-omean ~ dpi, t1_b3)
sop_nd <- data.frame(dpi=seq(min(t1_b3$dpi), max(t1_b3$dpi), length=100))
sop_nd$fit <- predict(sop_fit, newdata=sop_nd) 
sop_max <- sop_nd$dpi[which.max(sop_nd$fit)]

snp_fit <- loess(35-nmean ~ dpi, t1_b3)
snp_nd <- data.frame(dpi=seq(min(t1_b3$dpi), max(t1_b3$dpi), length=100))
snp_nd$fit <- predict(snp_fit, newdata=snp_nd) 
snp_max <- snp_nd$dpi[which.max(snp_nd$fit)]

srp_fit <- loess(35-rmean ~ dpi, t1_b3)
srp_nd <- data.frame(dpi=seq(min(t1_b3$dpi), max(t1_b3$dpi), length=100))
srp_nd$fit <- predict(srp_fit, newdata=srp_nd) 
srp_max <- srp_nd$dpi[which.max(srp_nd$fit)]

swbp_fit <- loess(35-wbmean ~ dpi, t1_b3)
swbp_nd <- data.frame(dpi=seq(min(t1_b3$dpi), max(t1_b3$dpi), length=100))
swbp_nd$fit <- predict(swbp_fit, newdata=swbp_nd) 
swbp_max <- swbp_nd$dpi[which.max(swbp_nd$fit)]

so <- ggplot(t1_b3, aes(dpi, 35-omean, color=barn)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "Ocular", x = "", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial1_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = sop_max), color = "deeppink2") +
  annotate(geom = "text", x = (sop_max - 4), y = 20, label = paste(round(sop_max,1)), color = "deeppink2") +
  theme_minimal() +
  theme(legend.position = "none", axis.text = element_text(size = 12)) 

sn <- ggplot(t1_b3, aes(dpi, 35-nmean, color = barn)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "Nasal", x = "", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial1_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = snp_max), color = "deeppink2") +
  annotate(geom = "text", x = (snp_max - 4), y = 20, label = paste(round(snp_max,1)), color = "deeppink2") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 12),axis.text.y = element_blank()) 

sr <- ggplot(t1_b3, aes(dpi, 35-rmean, color=barn)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "Rectal", x = "", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial1_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = srp_max), color = "deeppink2") +
  annotate(geom = "text", x = (srp_max - 4), y = 20, label = paste(round(srp_max,1)), color = "deeppink2") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.text.x = element_text(size = 12)) 

swb <- ggplot(t1_b3, aes(dpi, 35-wbmean, color=barn)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "Whole Blood", x = "", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(name = "Sheep", values = trial1_cols, labels = c("Innoculated", "Sentinel")) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = swbp_max), color = "deeppink2") +
  annotate(geom = "text", x = (swbp_max - 4), y = 20, label = paste(round(swbp_max,1)), color = "deeppink2") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.text.x = element_text(size = 12)) 


# GOAT
# getting peaks
gop_fit <- loess(35-omean ~ dpi, t1_b5)
gop_nd <- data.frame(dpi=seq(min(t1_b5$dpi), max(t1_b5$dpi), length=100))
gop_nd$fit <- predict(gop_fit, newdata=gop_nd) 
gop_max <- gop_nd$dpi[which.max(gop_nd$fit)]

gnp_fit <- loess(35-nmean ~ dpi, t1_b5)
gnp_nd <- data.frame(dpi=seq(min(t1_b5$dpi), max(t1_b5$dpi), length=100))
gnp_nd$fit <- predict(gnp_fit, newdata=gnp_nd) 
gnp_max <- gnp_nd$dpi[which.max(gnp_nd$fit)]

grp_fit <- loess(35-rmean ~ dpi, t1_b5)
grp_nd <- data.frame(dpi=seq(min(t1_b5$dpi), max(t1_b5$dpi), length=100))
grp_nd$fit <- predict(grp_fit, newdata=grp_nd) 
grp_max <- grp_nd$dpi[which.max(grp_nd$fit)]

gwbp_fit <- loess(35-wbmean ~ dpi, t1_b5)
gwbp_nd <- data.frame(dpi=seq(min(t1_b5$dpi), max(t1_b5$dpi), length=100))
gwbp_nd$fit <- predict(gwbp_fit, newdata=gwbp_nd) 
gwbp_max <- gwbp_nd$dpi[which.max(gwbp_nd$fit)]

gn <- ggplot(t1_b5, aes(dpi, 35-nmean, color=barn)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "", x = "", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial1_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = gnp_max), color = "cyan4") +
  annotate(geom = "text", x = (gnp_max - 4), y = 20, label = paste(round(gnp_max,1)), color = "cyan4") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 12), axis.text.y = element_blank())

go <- ggplot(t1_b5, aes(dpi, 35-omean, color=barn)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "", x = "", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial1_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = gop_max), color = "cyan4") +
  annotate(geom = "text", x = (gop_max - 4), y = 20, label = paste(round(gop_max,1)), color = "cyan4") +
  theme_minimal() +
  theme(legend.position = "none", axis.text = element_text(size = 12))

gr <- ggplot(t1_b5, aes(dpi, 35-rmean, color=barn)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "", x = "", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial1_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = grp_max), color = "cyan4") +
  annotate(geom = "text", x = (grp_max - 4), y = 20, label = paste(round(grp_max,1)), color = "cyan4") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.text.x = element_text(size = 12))

gwb <- ggplot(t1_b5, aes(dpi, 35-wbmean, color=barn)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "", x = "", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(name = "Goat", values = trial1_cols, labels = c("Innoculated", "Sentinel")) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = gwbp_max), color = "cyan4") +
  annotate(geom = "text", x = (gwbp_max - 4), y = 20, label = paste(round(gwbp_max,1)), color = "cyan4") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.text.x = element_text(size = 12))

# Both species
jpeg("output/molec/trial1_molecular_experimental_FigS4D.jpeg", width = 10, height = 6, units = "in", quality = 100, res = 600)
grid.arrange(so, sn, sr, swb, go, gn, gr, gwb, nrow = 2, left = "35 - mean qRTPCR Ct", bottom = "Days Post Infection")
invisible(dev.off())

# Controls Barn 1
jpeg("output/molec/trial1_molecular_negcontrols.jpeg", width = 8, height = 5, units = "in", quality = 100, res = 600)
cn <-ggplot(t1_b1, aes(dpi, 35-nmean, color=barn)) +
  geom_line() +
  geom_point() +
  labs(title = "Nasal", x = "", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(name = "Goat", values = trial1_cols, labels = c("Innoculated", "Sentinel")) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() +
  theme(legend.position = "none")

co <- ggplot(t1_b1, aes(dpi, 35-omean, color=barn)) +
  geom_line() +
  geom_point() +
  labs(title = "Ocular", x = "", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(name = "Goat", values = trial1_cols, labels = c("Innoculated", "Sentinel")) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank())

cr <- ggplot(t1_b1, aes(dpi, 35-rmean, color=barn)) +
  geom_line() +
  geom_point() +
  labs(title = "Rectal", x = "", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  scale_color_manual(name = "Goat", values = trial1_cols, labels = c("Innoculated", "Sentinel")) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank())

cwb <- ggplot(t1_b1, aes(dpi, 35-wbmean, color=barn)) +
  geom_line() +
  geom_point() +
  labs(title = "Whole Blood", x = "", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  scale_color_manual(name = "Goat", values = trial1_cols, labels = c("Innoculated", "Sentinel")) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank())

grid.arrange(cn, co, cr, cwb, ncol = 4, left = "35 - mean qRTPCR Ct", bottom = "Days Post Infection")
invisible(dev.off())



#####################
#  Trial 2 ----
#####################

t2_b4 <- subset(trial2_molec, barn == 4) # goats
t2_b6 <- subset(trial2_molec, barn == 6) # sheep
# controls are barn 5

trial2_cols <- c("4a" = "cyan4", "4b" = "cyan3", "5b" = "black", "6a" = "deeppink2", "6b" = "rosybrown3" )

# getting peaks
# ocular
gpi1Qo_fit <- loess(35-omean ~ dpi, t2_b4[t2_b4$innoc == 1,])
gpi1Qo_nd <- data.frame(dpi=seq(min(t2_b4$dpi), max(t2_b4$dpi), length=100))
gpi1Qo_nd$fit <- predict(gpi1Qo_fit, newdata=gpi1Qo_nd) 
gpi1Qo_max <- gpi1Qo_nd$dpi[which.max(gpi1Qo_nd$fit)]

gpi0Qo_fit <- loess(35-omean ~ dpi, t2_b4[t2_b4$innoc == 0,])
gpi0Qo_nd <- data.frame(dpi=seq(min(t2_b4$dpi), max(t2_b4$dpi), length=100))
gpi0Qo_nd$fit <- predict(gpi0Qo_fit, newdata=gpi0Qo_nd) 
gpi0Qo_max <- gpi0Qo_nd$dpi[which.max(gpi0Qo_nd$fit)]

spi1Qo_fit <- loess(35-omean ~ dpi, t2_b6[t2_b6$innoc == 1,])
spi1Qo_nd <- data.frame(dpi=seq(min(t2_b6$dpi), max(t2_b6$dpi), length=100))
spi1Qo_nd$fit <- predict(spi1Qo_fit, newdata=spi1Qo_nd) 
spi1Qo_max <- spi1Qo_nd$dpi[which.max(spi1Qo_nd$fit)]

spi0Qo_fit <- loess(35-omean ~ dpi, t2_b6[t2_b6$innoc == 0,])
spi0Qo_nd <- data.frame(dpi=seq(min(t2_b6$dpi), max(t2_b6$dpi), length=100))
spi0Qo_nd$fit <- predict(spi0Qo_fit, newdata=spi0Qo_nd) 
spi0Qo_max <- spi0Qo_nd$dpi[which.max(spi0Qo_nd$fit)]

# nasal
gpi1Q_fit <- loess(35-nmean ~ dpi, t2_b4[t2_b4$innoc == 1,])
gpi1Q_nd <- data.frame(dpi=seq(min(t2_b4$dpi), max(t2_b4$dpi), length=100))
gpi1Q_nd$fit <- predict(gpi1Q_fit, newdata=gpi1Q_nd) 
gpi1Q_max <- gpi1Q_nd$dpi[which.max(gpi1Q_nd$fit)]

gpi0Q_fit <- loess(35-nmean ~ dpi, t2_b4[t2_b4$innoc == 0,])
gpi0Q_nd <- data.frame(dpi=seq(min(t2_b4$dpi), max(t2_b4$dpi), length=100))
gpi0Q_nd$fit <- predict(gpi0Q_fit, newdata=gpi0Q_nd) 
gpi0Q_max <- gpi0Q_nd$dpi[which.max(gpi0Q_nd$fit)]

spi1Q_fit <- loess(35-nmean ~ dpi, t2_b6[t2_b6$innoc == 1,])
spi1Q_nd <- data.frame(dpi=seq(min(t2_b6$dpi), max(t2_b6$dpi), length=100))
spi1Q_nd$fit <- predict(spi1Q_fit, newdata=spi1Q_nd) 
spi1Q_max <- spi1Q_nd$dpi[which.max(spi1Q_nd$fit)]

spi0Q_fit <- loess(35-nmean ~ dpi, t2_b6[t2_b6$innoc == 0,])
spi0Q_nd <- data.frame(dpi=seq(min(t2_b6$dpi), max(t2_b6$dpi), length=100))
spi0Q_nd$fit <- predict(spi0Q_fit, newdata=spi0Q_nd) 
spi0Q_max <- spi0Q_nd$dpi[which.max(spi0Q_nd$fit)]

# rectal
gpi1Qr_fit <- loess(35-rmean ~ dpi, t2_b4[t2_b4$innoc == 1,])
gpi1Qr_nd <- data.frame(dpi=seq(min(t2_b4$dpi), max(t2_b4$dpi), length=100))
gpi1Qr_nd$fit <- predict(gpi1Qr_fit, newdata=gpi1Qr_nd) 
gpi1Qr_max <- gpi1Qr_nd$dpi[which.max(gpi1Qr_nd$fit)]

gpi0Qr_fit <- loess(35-rmean ~ dpi, t2_b4[t2_b4$innoc == 0,])
gpi0Qr_nd <- data.frame(dpi=seq(min(t2_b4$dpi), max(t2_b4$dpi), length=100))
gpi0Qr_nd$fit <- predict(gpi0Qr_fit, newdata=gpi0Qr_nd) 
gpi0Qr_max <- gpi0Qr_nd$dpi[which.max(gpi0Qr_nd$fit)]

spi1Qr_fit <- loess(35-rmean ~ dpi, t2_b6[t2_b6$innoc == 1,])
spi1Qr_nd <- data.frame(dpi=seq(min(t2_b6$dpi), max(t2_b6$dpi), length=100))
spi1Qr_nd$fit <- predict(spi1Qr_fit, newdata=spi1Qr_nd) 
spi1Qr_max <- spi1Qr_nd$dpi[which.max(spi1Qr_nd$fit)]

spi0Qr_fit <- loess(35-rmean ~ dpi, t2_b6[t2_b6$innoc == 0,])
spi0Qr_nd <- data.frame(dpi=seq(min(t2_b6$dpi), max(t2_b6$dpi), length=100))
spi0Qr_nd$fit <- predict(spi0Qr_fit, newdata=spi0Qr_nd) 
spi0Qr_max <- spi0Qr_nd$dpi[which.max(spi0Qr_nd$fit)]

#whole blood
gpi1Qwb_fit <- loess(35-wbmean ~ dpi, t2_b4[t2_b4$innoc == 1,])
gpi1Qwb_nd <- data.frame(dpi=seq(min(t2_b4$dpi), max(t2_b4$dpi), length=100))
gpi1Qwb_nd$fit <- predict(gpi1Qwb_fit, newdata=gpi1Qwb_nd) 
gpi1Qwb_max <- gpi1Qwb_nd$dpi[which.max(gpi1Qwb_nd$fit)]

gpi0Qwb_fit <- loess(35-wbmean ~ dpi, t2_b4[t2_b4$innoc == 0,])
gpi0Qwb_nd <- data.frame(dpi=seq(min(t2_b4$dpi), max(t2_b4$dpi), length=100))
gpi0Qwb_nd$fit <- predict(gpi0Qwb_fit, newdata=gpi0Qwb_nd) 
gpi0Qwb_max <- gpi0Qwb_nd$dpi[which.max(gpi0Qwb_nd$fit)]

spi1Qwb_fit <- loess(35-wbmean ~ dpi, t2_b6[t2_b6$innoc == 1,])
spi1Qwb_nd <- data.frame(dpi=seq(min(t2_b6$dpi), max(t2_b6$dpi), length=100))
spi1Qwb_nd$fit <- predict(spi1Qwb_fit, newdata=spi1Qwb_nd) 
spi1Qwb_max <- spi1Qwb_nd$dpi[which.max(spi1Qwb_nd$fit)]

spi0Qwb_fit <- loess(35-wbmean ~ dpi, t2_b6[t2_b6$innoc == 0,])
spi0Qwb_nd <- data.frame(dpi=seq(min(t2_b6$dpi), max(t2_b6$dpi), length=100))
spi0Qwb_nd$fit <- predict(spi0Qwb_fit, newdata=spi0Qwb_nd) 
spi0Qwb_max <- spi0Qwb_nd$dpi[which.max(spi0Qwb_nd$fit)]


# plots
# ocular
so_comp <- ggplot(t2_b6, aes(dpi, 35-omean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "Ocular", x = "Days Post Infection", y = "RT-qPCR 35 - mean Ct") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial2_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = spi1Qo_max), color = "deeppink2") +
  geom_vline(aes(xintercept = spi0Qo_max), color = "rosybrown2") +
  annotate(geom = "text", x = (spi1Qo_max - 1.75), y = 20, label = paste(round(spi1Qo_max,1)), color = "deeppink2") +
  annotate(geom = "text", x = (spi0Qo_max + 1.75), y = 20, label = paste(round(spi0Qo_max,1)), color = "rosybrown2") +
  annotate(geom = "text", x = mean(c(spi1Qo_max, spi0Qo_max)), y = 19, label = paste(round(spi0Qo_max - spi1Qo_max,1)), color = "black") +
  geom_segment(x = spi1Qo_max, y = 17.5, xend = spi0Qo_max, yend = 17.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

go_comp <- ggplot(t2_b4, aes(dpi, 35-omean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "Ocular", x = "Days Post Infection", y = "RT-qPCR 35 - mean Ct") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial2_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = gpi1Qo_max), color = "cyan4") +
  geom_vline(aes(xintercept = gpi0Qo_max), color = "cyan3") +
  annotate(geom = "text", x = (gpi1Qo_max - 1.75), y = 20, label = paste(round(gpi1Qo_max,1)), color = "cyan4") +
  annotate(geom = "text", x = (gpi0Qo_max + 1.75), y = 20, label = paste(round(gpi0Qo_max,1)), color = "cyan3") +
  annotate(geom = "text", x = mean(c(gpi1Qo_max, gpi0Qo_max)), y = 19, label = paste(round(gpi0Qo_max - gpi1Qo_max,1)), color = "black") +
  geom_segment(x = gpi1Qo_max, y = 17.5, xend = gpi0Qo_max, yend = 17.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

# nasal
sn_comp <- ggplot(t2_b6, aes(dpi, 35-nmean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "Nasal", x = "Days Post Infection", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial2_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = spi1Q_max), color = "deeppink2") +
  geom_vline(aes(xintercept = spi0Q_max), color = "rosybrown2") +
  annotate(geom = "text", x = (spi1Q_max - 1.75), y = 20, label = paste(round(spi1Q_max,1)), color = "deeppink2") +
  annotate(geom = "text", x = (spi0Q_max + 1.75), y = 20, label = paste(round(spi0Q_max,1)), color = "rosybrown2") +
  annotate(geom = "text", x = mean(c(spi1Q_max, spi0Q_max)), y = 19, label = paste(round(spi0Q_max - spi1Q_max,1)), color = "black") +
  geom_segment(x = spi1Q_max, y = 17.5, xend = spi0Q_max, yend = 17.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12))

gn_comp <- ggplot(t2_b4, aes(dpi, 35-nmean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "Nasal", x = "Days Post Infection", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial2_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = gpi1Q_max), color = "cyan4") +
  geom_vline(aes(xintercept = gpi0Q_max), color = "cyan3") +
  annotate(geom = "text", x = (gpi1Q_max - 1.75), y = 20, label = paste(round(gpi1Q_max,1)), color = "cyan4") +
  annotate(geom = "text", x = (gpi0Q_max + 1.75), y = 20, label = paste(round(gpi0Q_max,1)), color = "cyan3") +
  annotate(geom = "text", x = mean(c(gpi1Q_max, gpi0Q_max)), y = 19, label = paste(round(gpi0Q_max - gpi1Q_max,1)), color = "black") +
  geom_segment(x = gpi1Q_max, y = 17.5, xend = gpi0Q_max, yend = 17.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12))

# rectal
sr_comp <- ggplot(t2_b6, aes(dpi, 35-rmean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "Rectal", x = "Days Post Infection", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial2_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = spi1Qr_max), color = "deeppink2") +
  geom_vline(aes(xintercept = spi0Qr_max), color = "rosybrown2") +
  annotate(geom = "text", x = (spi1Qr_max - 1.75), y = 20, label = paste(round(spi1Qr_max,1)), color = "deeppink2") +
  annotate(geom = "text", x = (spi0Qr_max + 1.75), y = 20, label = paste(round(spi0Qr_max,1)), color = "rosybrown2") +
  annotate(geom = "text", x = mean(c(spi1Qr_max, spi0Qr_max)), y = 19, label = paste(round(spi0Qr_max - spi1Qr_max,1)), color = "black") +
  geom_segment(x = spi1Qr_max, y = 17.5, xend = spi0Qr_max, yend = 17.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12))

gr_comp <- ggplot(t2_b4, aes(dpi, 35-rmean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "Rectal", x = "Days Post Infection", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial2_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = gpi1Qr_max), color = "cyan4") +
  geom_vline(aes(xintercept = gpi0Qr_max), color = "cyan3") +
  annotate(geom = "text", x = (gpi1Qr_max - 1.75), y = 20, label = paste(round(gpi1Qr_max,1)), color = "cyan4") +
  annotate(geom = "text", x = (gpi0Qr_max + 1.75), y = 20, label = paste(round(gpi0Qr_max,1)), color = "cyan3") +
  annotate(geom = "text", x = mean(c(gpi1Qr_max, gpi0Qr_max)), y = 19, label = paste(round(gpi0Qr_max - gpi1Qr_max,1)), color = "black") +
  geom_segment(x = gpi1Qr_max, y = 17.5, xend = gpi0Qr_max, yend = 17.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12))

# whole blood
swb_comp <- ggplot(t2_b6, aes(dpi, 35-wbmean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "Whole Blood", x = "Days Post Infection", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial2_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = spi1Qwb_max), color = "deeppink2") +
  geom_vline(aes(xintercept = spi0Qwb_max), color = "rosybrown2") +
  annotate(geom = "text", x = (spi1Qwb_max - 1.75), y = 20, label = paste(round(spi1Qwb_max,1)), color = "deeppink2") +
  annotate(geom = "text", x = (spi0Qwb_max + 1.75), y = 20, label = paste(round(spi0Qwb_max,1)), color = "rosybrown2") +
  annotate(geom = "text", x = mean(c(spi1Qwb_max, spi0Qwb_max)), y = 19, label = paste(round(spi0Qwb_max - spi1Qwb_max,1)), color = "black") +
  geom_segment(x = spi1Qwb_max, y = 17.5, xend = spi0Qwb_max, yend = 17.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12))

gwb_comp <- ggplot(t2_b4, aes(dpi, 35-wbmean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "Whole Blood", x = "Days Post Infection", y = "") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial2_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = gpi1Qwb_max), color = "cyan4") +
  geom_vline(aes(xintercept = gpi0Qwb_max), color = "cyan3") +
  annotate(geom = "text", x = (gpi1Qwb_max - 1.75), y = 20, label = paste(round(gpi1Qwb_max,1)), color = "cyan4") +
  annotate(geom = "text", x = (gpi0Qwb_max + 1.75), y = 20, label = paste(round(gpi0Qwb_max,1)), color = "cyan3") +
  annotate(geom = "text", x = mean(c(gpi1Qwb_max, gpi0Qwb_max)), y = 19, label = paste(round(gpi0Qwb_max - gpi1Qwb_max,1)), color = "black") +
  geom_segment(x = gpi1Qwb_max, y = 17.5, xend = gpi0Qwb_max, yend = 17.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12))

jpeg("output/molec/trial2_molecular_qRTPCR_sheep_FigS5D.jpeg", width = 12, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(so_comp,sn_comp,sr_comp, swb_comp, ncol=4)
invisible(dev.off())

jpeg("output/molec/trial2_molecular_qRTPCR_goats_FigS5D.jpeg", width = 12, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(go_comp,gn_comp,gr_comp, gwb_comp, ncol=4)
invisible(dev.off())

# IC-ELISA ----
# Comparing nasal icELISA to qRTPCR ----
# getting peaks
# icELISA
gpi1E_fit <- loess(sp ~ dpi, t2_b4[t2_b4$innoc == 1,])
gpi1E_nd <- data.frame(dpi=seq(min(t2_b4$dpi), max(t2_b4$dpi), length=100))
gpi1E_nd$fit <- predict(gpi1E_fit, newdata=gpi1E_nd) 
gpi1E_max <- gpi1E_nd$dpi[which.max(gpi1E_nd$fit)]

gpi0E_fit <- loess(sp ~ dpi, t2_b4[t2_b4$innoc == 0,])
gpi0E_nd <- data.frame(dpi=seq(min(t2_b4$dpi), max(t2_b4$dpi), length=100))
gpi0E_nd$fit <- predict(gpi0E_fit, newdata=gpi0E_nd) 
gpi0E_max <- gpi0E_nd$dpi[which.max(gpi0E_nd$fit)]

spi1E_fit <- loess(sp ~ dpi, t2_b6[t2_b6$innoc == 1,])
spi1E_nd <- data.frame(dpi=seq(min(t2_b6$dpi), max(t2_b6$dpi), length=100))
spi1E_nd$fit <- predict(spi1E_fit, newdata=spi1E_nd) 
spi1E_max <- spi1E_nd$dpi[which.max(spi1E_nd$fit)]

spi0E_fit <- loess(sp ~ dpi, t2_b6[t2_b6$innoc == 0,])
spi0E_nd <- data.frame(dpi=seq(min(t2_b6$dpi), max(t2_b6$dpi), length=100))
spi0E_nd$fit <- predict(spi0E_fit, newdata=spi0E_nd) 
spi0E_max <- spi0E_nd$dpi[which.max(spi0E_nd$fit)]

#qrtPCR
gpi1Q_fit <- loess(35-nmean ~ dpi, t2_b4[t2_b4$innoc == 1,])
gpi1Q_nd <- data.frame(dpi=seq(min(t2_b4$dpi), max(t2_b4$dpi), length=100))
gpi1Q_nd$fit <- predict(gpi1Q_fit, newdata=gpi1Q_nd) 
gpi1Q_max <- gpi1Q_nd$dpi[which.max(gpi1Q_nd$fit)]

gpi0Q_fit <- loess(35-nmean ~ dpi, t2_b4[t2_b4$innoc == 0,])
gpi0Q_nd <- data.frame(dpi=seq(min(t2_b4$dpi), max(t2_b4$dpi), length=100))
gpi0Q_nd$fit <- predict(gpi0Q_fit, newdata=gpi0Q_nd) 
gpi0Q_max <- gpi0Q_nd$dpi[which.max(gpi0Q_nd$fit)]

spi1Q_fit <- loess(35-nmean ~ dpi, t2_b6[t2_b6$innoc == 1,])
spi1Q_nd <- data.frame(dpi=seq(min(t2_b6$dpi), max(t2_b6$dpi), length=100))
spi1Q_nd$fit <- predict(spi1Q_fit, newdata=spi1Q_nd) 
spi1Q_max <- spi1Q_nd$dpi[which.max(spi1Q_nd$fit)]

spi0Q_fit <- loess(35-nmean ~ dpi, t2_b6[t2_b6$innoc == 0,])
spi0Q_nd <- data.frame(dpi=seq(min(t2_b6$dpi), max(t2_b6$dpi), length=100))
spi0Q_nd$fit <- predict(spi0Q_fit, newdata=spi0Q_nd) 
spi0Q_max <- spi0Q_nd$dpi[which.max(spi0Q_nd$fit)]

#Plots 
molecs_comp <- ggplot(t2_b6, aes(dpi, sp, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "AgELISA", x = "", y = "S/P % of Positive Control") +
  scale_color_manual(name = "Sheep", 
                     labels = c("Inoculated", "Sentinel"),
                     values = trial2_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = spi1E_max), color = "deeppink2") +
  geom_vline(aes(xintercept = spi0E_max), color = "rosybrown2") +
  annotate(geom = "text", x = (spi1E_max - 1.75), y = 500, label = paste(round(spi1E_max,1)), color = "deeppink2") +
  annotate(geom = "text", x = (spi0E_max + 1.75), y = 500, label = paste(round(spi0E_max,1)), color = "rosybrown2") +
  annotate(geom = "text", x = mean(c(spi1E_max, spi0E_max)), y = 485, label = paste(round(spi0E_max - spi1E_max,1)), color = "black") +
  geom_segment(x = spi1E_max, y = 440, xend = spi0E_max, yend = 440, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() + theme(axis.title = element_text(size = 15),
                          axis.text = element_text(size=12), legend.position = "none")

sn_comp <- ggplot(t2_b6, aes(dpi, 35-nmean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "qRTPCR", x = "Days Post Infection", y = "RT-qPCR 35 - mean Ct") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial2_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = spi1Q_max), color = "deeppink2") +
  geom_vline(aes(xintercept = spi0Q_max), color = "rosybrown2") +
  annotate(geom = "text", x = (spi1Q_max - 1.75), y = 20, label = paste(round(spi1Q_max,1)), color = "deeppink2") +
  annotate(geom = "text", x = (spi0Q_max + 1.75), y = 20, label = paste(round(spi0Q_max,1)), color = "rosybrown2") +
  annotate(geom = "text", x = mean(c(spi1Q_max, spi0Q_max)), y = 19, label = paste(round(spi0Q_max - spi1Q_max,1)), color = "black") +
  geom_segment(x = spi1Q_max, y = 17.5, xend = spi0Q_max, yend = 17.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme_minimal() + theme(axis.title = element_text(size = 15),
                          axis.text = element_text(size=12), legend.position = "none")

molecg_comp <- ggplot(t2_b4, aes(dpi, sp, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "AgELISA", x = "", y = "S/P % of Positive Control") +
  scale_color_manual(name = "Goat", 
                     labels = c("Inoculated", "Sentinel"),
                     values = trial2_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = gpi1E_max), color = "cyan4") +
  geom_vline(aes(xintercept = gpi0E_max), color = "cyan3") +
  annotate(geom = "text", x = (gpi1E_max - 1.75), y = 500, label = paste(round(gpi1E_max,1)), color = "cyan4") +
  annotate(geom = "text", x = (gpi0E_max + 1.75), y = 500, label = paste(round(gpi0E_max,1)), color = "cyan3") +
  annotate(geom = "text", x = mean(c(gpi1E_max, gpi0E_max)), y = 485, label = paste(round(gpi0E_max - gpi1E_max,1)), color = "black") +
  geom_segment(x = gpi1E_max, y = 440, xend = gpi0E_max, yend = 440, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() + theme(axis.title = element_text(size = 15),
                          axis.text = element_text(size=12), legend.position = "none")

gn_comp <- ggplot(t2_b4, aes(dpi, 35-nmean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "qRTPCR", x = "Days Post Infection", y = "RT-qPCR 35 - mean Ct") +
  scale_x_continuous(limits = c(0,28), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28)) +
  #scale_y_continuous(limits = c(0,35.1)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial2_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = gpi1Q_max), color = "cyan4") +
  geom_vline(aes(xintercept = gpi0Q_max), color = "cyan3") +
  annotate(geom = "text", x = (gpi1Q_max - 1.75), y = 20, label = paste(round(gpi1Q_max,1)), color = "cyan4") +
  annotate(geom = "text", x = (gpi0Q_max + 1.75), y = 20, label = paste(round(gpi0Q_max,1)), color = "cyan3") +
  annotate(geom = "text", x = mean(c(gpi1Q_max, gpi0Q_max)), y = 19, label = paste(round(gpi0Q_max - gpi1Q_max,1)), color = "black") +
  geom_segment(x = gpi1Q_max, y = 17.5, xend = gpi0Q_max, yend = 17.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme_minimal() + theme(axis.title = element_text(size = 15),
                          axis.text = element_text(size=12), legend.position = "none")

jpeg("output/molec/trial2_molecular_icELISA_FigS8.jpeg", width = 6.5, height = 6, units = "in", quality = 100, res = 600)
molecg_comp2 <- molecg_comp + labs(y ="", title = "") + theme(axis.text.y = element_blank())
gn_comp2 <- gn_comp + labs(x = "", y ="", title = "") + theme(axis.text.y = element_blank())
sn_comp2 <- sn_comp + labs(x = "")
grid.arrange(molecs_comp, molecg_comp2, sn_comp2, gn_comp2, ncol=2, nrow=2, 
             left = textGrob("Nasal Swab", rot = 90, gp=gpar(fontsize=15)),
             bottom = textGrob("Days Post Infection", gp=gpar(fontsize=15)))
invisible(dev.off())


#####################
#  Trial 3 ----
#####################

t3_b3 <- subset(trial3_molec, barn == 3) # Box A Inoculated cattle to sentinel goats
t3_b4 <- subset(trial3_molec, barn == 4) # Box C Inoculated goats to sentinel cattle and goats
t3_b5 <- subset(trial3_molec, barn == 5) # Box B Positive control goats
# barn 6 negative controls

trial3_cols <- c("5a" = "cyan4", "5b" = "cyan3", "6b" = "black", "6c" = "black", "3a" = "purple4", "3b" = "cyan3", 
                  "4a" = "cyan4", "4b" = "cyan3", "4c" = "mediumpurple1")


# Barn 3 / Box A ----
# Nasal
molecboxAn_ic <- ggplot(t3_b3, aes(dpi, sp, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "icELISA", x = "", y = "") +
  scale_color_manual(name = "Barn A: Nasal Status", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial3_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8),
                          axis.text=element_text(size=6))

molecboxAn_q <- ggplot(t3_b3, aes(dpi, 35-nmean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "qRTPCR", x = "Days Post Infection", y = "35 - mean Ct") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial3_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) 

grid.arrange(molecboxAn_ic,molecboxAn_q, nrow=2, left = "Nasal Swab")

# Fecal
molecboxAf_ic <- ggplot(t3_b3, aes(dpi, sp_f, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "icELISA", x = "", y = "") +
  scale_color_manual(name = "Barn A: Nasal Status", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial3_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8),
                          axis.text=element_text(size=6))

molecboxAf_q <- ggplot(t3_b3, aes(dpi, 35-rmean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "qRTPCR", x = "Days Post Infection", y = "35 - mean Ct") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial3_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) 

grid.arrange(molecboxAf_ic,molecboxAf_q, nrow=2, left = "Rectal Swab")


# Barn 4 / Box C ----
# getting peaks
# ocular
gpi1Qo_fit <- loess(35-omean ~ dpi, t3_b4[t3_b4$barni == "4a",])
gpi1Qo_nd <- data.frame(dpi=seq(min(t3_b4$dpi), max(t3_b4$dpi), length=100))
gpi1Qo_nd$fit <- predict(gpi1Qo_fit, newdata=gpi1Qo_nd) 
gpi1Qo_max <- gpi1Qo_nd$dpi[which.max(gpi1Qo_nd$fit)]

cpi0Qo_fit <- loess(35-omean ~ dpi, t3_b4[t3_b4$barni == "4c",])
cpi0Qo_nd <- data.frame(dpi=seq(min(t3_b4$dpi), max(t3_b4$dpi), length=100))
cpi0Qo_nd$fit <- predict(cpi0Qo_fit, newdata=cpi0Qo_nd) 
cpi0Qo_max <- cpi0Qo_nd$dpi[which.max(cpi0Qo_nd$fit)]

gpi0Qo_fit <- loess(35-omean ~ dpi, t3_b4[t3_b4$barni == "4b",])
gpi0Qo_nd <- data.frame(dpi=seq(min(t3_b4$dpi), max(t3_b4$dpi), length=100))
gpi0Qo_nd$fit <- predict(gpi0Qo_fit, newdata=gpi0Qo_nd) 
gpi0Qo_max <- gpi0Qo_nd$dpi[which.max(gpi0Qo_nd$fit)]

# nasal
gpi1Q_fit <- loess(35-nmean ~ dpi, t3_b4[t3_b4$barni == "4a",])
gpi1Q_nd <- data.frame(dpi=seq(min(t3_b4$dpi), max(t3_b4$dpi), length=100))
gpi1Q_nd$fit <- predict(gpi1Q_fit, newdata=gpi1Q_nd) 
gpi1Q_max <- gpi1Q_nd$dpi[which.max(gpi1Q_nd$fit)]

cpi0Q_fit <- loess(35-nmean ~ dpi, t3_b4[t3_b4$barni == "4c",])
cpi0Q_nd <- data.frame(dpi=seq(min(t3_b4$dpi), max(t3_b4$dpi), length=100))
cpi0Q_nd$fit <- predict(cpi0Q_fit, newdata=cpi0Q_nd) 
cpi0Q_max <- cpi0Q_nd$dpi[which.max(cpi0Q_nd$fit)]

gpi0Q_fit <- loess(35-nmean ~ dpi, t3_b4[t3_b4$barni == "4b",])
gpi0Q_nd <- data.frame(dpi=seq(min(t3_b4$dpi), max(t3_b4$dpi), length=100))
gpi0Q_nd$fit <- predict(gpi0Q_fit, newdata=gpi0Q_nd) 
gpi0Q_max <- gpi0Q_nd$dpi[which.max(gpi0Q_nd$fit)]

# rectal
gpi1Qr_fit <- loess(35-rmean ~ dpi, t3_b4[t3_b4$barni == "4a",])
gpi1Qr_nd <- data.frame(dpi=seq(min(t3_b4$dpi), max(t3_b4$dpi), length=100))
gpi1Qr_nd$fit <- predict(gpi1Qr_fit, newdata=gpi1Qr_nd) 
gpi1Qr_max <- gpi1Qr_nd$dpi[which.max(gpi1Qr_nd$fit)]

cpi0Qr_fit <- loess(35-rmean ~ dpi, t3_b4[t3_b4$barni == "4c",])
cpi0Qr_nd <- data.frame(dpi=seq(min(t3_b4$dpi), max(t3_b4$dpi), length=100))
cpi0Qr_nd$fit <- predict(cpi0Qr_fit, newdata=cpi0Qr_nd) 
cpi0Qr_max <- cpi0Qr_nd$dpi[which.max(cpi0Qr_nd$fit)]

gpi0Qr_fit <- loess(35-rmean ~ dpi, t3_b4[t3_b4$barni == "4b",])
gpi0Qr_nd <- data.frame(dpi=seq(min(t3_b4$dpi), max(t3_b4$dpi), length=100))
gpi0Qr_nd$fit <- predict(gpi0Qr_fit, newdata=gpi0Qr_nd) 
gpi0Qr_max <- gpi0Qr_nd$dpi[which.max(gpi0Qr_nd$fit)]

# Plots
# Ocular
molecboxCo_q <- ggplot(t3_b4, aes(dpi, 35-omean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(#title = "qRTPCR", x = "Days Post Infection", y = "35 - mean Ct"
    title = "Ocular", x = "", y = "RT-qPCR 35- mean Ct") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial3_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = gpi1Qo_max), color = "cyan4") +
  geom_vline(aes(xintercept = gpi0Qo_max), color = "cyan3") +
  geom_vline(aes(xintercept = cpi0Qo_max), color = "mediumpurple3") +
  annotate(geom = "text", x = (gpi1Qo_max - 1.75), y = 20, label = paste(round(gpi1Qo_max,1)), color = "cyan4") +
  annotate(geom = "text", x = (gpi0Qo_max + 1.75), y = 20, label = paste(round(gpi0Qo_max,1)), color = "cyan3") +
  annotate(geom = "text", x = (cpi0Qo_max - 1.75), y = 20, label = paste(round(cpi0Qo_max,1)), color = "mediumpurple3") +
  annotate(geom = "text", x = mean(c(gpi1Qo_max, gpi0Qo_max)), y = 19, label = paste(round(gpi0Qo_max - gpi1Qo_max,1)), color = "black") +
  geom_segment(x = gpi1Qo_max, y = 17.5, xend = gpi0Qo_max, yend = 17.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

# Nasal
molecboxCn_ic <- ggplot(t3_b4, aes(dpi, sp, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "icELISA", x = "", y = "") +
  scale_color_manual(name = "Barn C: Nasal Status", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial3_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12),
        axis.text.y = element_blank())

molecboxCn_q <- ggplot(t3_b4, aes(dpi, 35-nmean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(#title = "qRTPCR", x = "Days Post Infection", y = "35 - mean Ct"
    title = "Nasal", x = "", y = "") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial3_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = gpi1Q_max), color = "cyan4") +
  geom_vline(aes(xintercept = gpi0Q_max), color = "cyan3") +
  geom_vline(aes(xintercept = cpi0Q_max), color = "mediumpurple3") +
  annotate(geom = "text", x = (gpi1Q_max - 0.75), y = 20, label = paste(round(gpi1Q_max,1)), color = "cyan4") +
  annotate(geom = "text", x = (gpi0Q_max + 1.75), y = 20, label = paste(round(gpi0Q_max,1)), color = "cyan3") +
  annotate(geom = "text", x = (cpi0Q_max - 1.75), y = 20, label = paste(round(cpi0Q_max,1)), color = "mediumpurple3") +
  annotate(geom = "text", x = mean(c(gpi1Q_max, gpi0Q_max)), y = 19, label = paste(round(gpi0Q_max - gpi1Q_max,1)), color = "black") +
  geom_segment(x = gpi1Q_max, y = 17.5, xend = gpi0Q_max, yend = 17.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12),
        axis.text.y = element_blank())

grid.arrange(molecboxCn_ic,molecboxCn_q, nrow=2, left = "Nasal Swab")

# Fecal swab
molecboxCf_ic <- ggplot(t3_b4, aes(dpi, sp_f, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "icELISA", x = "", y = "") +
  scale_color_manual(name = "Barn C: Nasal Status", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial3_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12),
        axis.text.y = element_blank())

molecboxCf_q <- ggplot(t3_b4, aes(dpi, 35-rmean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(#title = "qRTPCR", x = "Days Post Infection", y = "35 - mean Ct"
    title = "Rectal", x = "", y = "") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial3_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  geom_vline(aes(xintercept = gpi1Qr_max), color = "cyan4") +
  geom_vline(aes(xintercept = gpi0Qr_max), color = "cyan3") +
  geom_vline(aes(xintercept = cpi0Qr_max), color = "mediumpurple3") +
  annotate(geom = "text", x = (gpi1Qr_max - 1.75), y = 20, label = paste(round(gpi1Qr_max,1)), color = "cyan4") +
  annotate(geom = "text", x = (gpi0Qr_max + 1.75), y = 20, label = paste(round(gpi0Qr_max,1)), color = "cyan3") +
  annotate(geom = "text", x = (cpi0Qr_max + 1.75), y = 20, label = paste(round(cpi0Qr_max,1)), color = "mediumpurple3") +
  annotate(geom = "text", x = mean(c(gpi1Qr_max, gpi0Qr_max)), y = 19, label = paste(round(gpi0Qr_max - gpi1Qr_max,1)), color = "black") +
  geom_segment(x = gpi1Qr_max, y = 17.5, xend = gpi0Qr_max, yend = 17.5, color = "black", arrow = arrow(length=unit(0.30,"cm"), ends = "both")) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12),
        axis.text.y = element_blank())
# Correct - there are only 2 points in the fecal in old one

grid.arrange(molecboxCf_ic,molecboxCf_q, nrow=2, left = "Rectal Swab")


jpeg("output/molec/trial3_molecular_qRTPCR_Fig1D.jpeg", width = 9, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(molecboxCo_q, molecboxCn_q, molecboxCf_q, nrow=1)
invisible(dev.off())



# Barn 5 / Box B Positive Controls  ----
# Nasal
molecboxBn_ic <- ggplot(t3_b5, aes(dpi, sp, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "Barn B: Nasal Status", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial3_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8),
                          axis.text=element_text(size=6))

molecboxBn_q <- ggplot(t3_b5, aes(dpi, 35-nmean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "qRTPCR", x = "Days Post Infection", y = "35 - mean Ct") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial3_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) 

grid.arrange(molecboxBn_ic,molecboxBn_q, nrow=2, left = "Nasal Swab")

# Fecal
molecboxBf_ic <- ggplot(t3_b5, aes(dpi, sp_f, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Rectal", x = "", y = "") +
  scale_color_manual(name = "Barn B: Nasal Status", 
                     labels = c("Cattle (Inoculated)", "Goats (In Contact)"),
                     values = trial3_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8),
                          axis.text=element_text(size=6))

molecboxBf_q <- ggplot(t3_b5, aes(dpi, 35-rmean, color=barni)) +
  geom_line(aes(group=eartag), alpha = 0.4) +
  geom_point(aes(group=eartag), alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  labs(title = "qRTPCR", x = "Days Post Infection", y = "35 - mean Ct") +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  scale_color_manual(values = trial3_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) 

grid.arrange(molecboxBf_ic,molecboxBf_q, nrow=2, left = "Rectal Swab")

jpeg("output/molec/trial3_molecular_icELISA_FigS10A.jpeg", width = 3.25, height = 3.25, units = "in", quality = 100, res = 600)
grid.arrange(molecboxBn_ic, molecboxBf_ic, nrow=2, left = "S/P Ratio")
invisible(dev.off())

#########################
#  Trial 3 Challenge ----
#########################

trial3_chall_cols <- c("3" = "cyan4", "4" = "black", "6" = "cyan2")

jpeg("output/molec/trial3_chall_molecular_FigS9D.jpeg", width = 3.5, height = 3.25, units = "in", quality = 100, res = 600)
molecN <- ggplot(trial3_chall_molec, aes(dpi, sp, color = barn)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,14), breaks = c(0, 4, 7, 10, 14)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "", 
                     labels = c("Cattle sentinel goats", "Seropositive goats", "Negative Controls"),
                     values = trial3_chall_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

molecF<- ggplot(trial3_chall_molec, aes(dpi, sp_f, color = barn)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,14), breaks = c(0, 4, 7, 10, 14)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Rectal", x = "Days Post Infection", y = "") +
  scale_color_manual(name = "", 
                     labels = c("Cattle sentinel goats", "Seropositive goats", "Negative Controls"),
                     values = trial3_chall_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))
grid.arrange(molecN, molecF, nrow=2, left = "S/P Ratio")
invisible(dev.off())


#####################
#  Trial 4 ----
#####################

t4_b6 <- subset(trial4_molec, barn == 6) # box A Inoculated cattle to sentinel goats
t4_b1 <- subset(trial4_molec, barn == 1) # box B Inoculated cattle to sentinel goats
t4_b4 <- subset(trial4_molec, barn == 4) # box C Inoculated cattle to sentinel goats
t4_b5 <- subset(trial4_molec, barn == 5) # positive controls
#t4_b3 <- subset(trial4_molec, barn == 3) # negative controls, not sampled
t4_b614 <- subset(trial4_molec, barn == 6 | barn == 1 | barn ==4) 


trial4_cols <- c("6a" = "purple4", "6b" = "cyan3", "1a" = "purple4", "1b" = "cyan3", "4a" = "purple4", "4b" = "cyan3",
                 "5a" = "cyan4", "5b" = "cyan3", "3a" = "cyan3")

trial4_cols_sp <- c("Bovine" = "purple4", "Caprine" = "cyan3")

# icELISA ----
# Barn 6 ----
molec_t4_b6n <- ggplot(t4_b6, aes(dpi, sp, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "Barn 6: Nasal Status", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8),
                          axis.text=element_text(size=6))

molec_t4_b6f <- ggplot(t4_b6, aes(dpi, sp_f, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Rectal", x = "Days Post Infection", y = "") +
  scale_color_manual(name = "Barn 6: Rectal Status", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8),
                          axis.text=element_text(size=6))

grid.arrange(molec_t4_b6n, molec_t4_b6f, nrow=2, left = "S/P Ratio")


# Barn 1 ----
molec_t4_b1n <-ggplot(t4_b1, aes(dpi, sp, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "Barn 1: Nasal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                          axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

molec_t4_b1f <-ggplot(t4_b1, aes(dpi, sp_f, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Rectal", x = "", y = "") +
  scale_color_manual(name = "Barn 1: Rectal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                          axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

grid.arrange(molec_t4_b1n, molec_t4_b1f, nrow=2, left = "S/P Ratio")

# Barn 4 ----
molec_t4_b4n <-ggplot(t4_b4, aes(dpi, sp, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "Barn 4: Nasal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                          axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

molec_t4_b4f <-ggplot(t4_b4, aes(dpi, sp_f, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Rectal", x = "Days Post Infection", y = "") +
  scale_color_manual(name = "Barn 4: Rectal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                          axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

grid.arrange(molec_t4_b4n, molec_t4_b4f, nrow=2, left = "S/P Ratio")

# Barn 5: Positive Controls ----
jpeg("output/molec/trial4_molecular_icELISA_FigS11A.jpeg", width = 3.5, height = 3.25, units = "in", quality = 100, res = 600)
molec_t4_b5n <-ggplot(t4_b5, aes(dpi, sp, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "Barn C: Nasal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

molec_t4_b5f <-ggplot(t4_b5, aes(dpi, sp_f, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  #stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Rectal", x = "Days Post Infection", y = "") +
  scale_color_manual(name = "Barn C: Rectal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

grid.arrange(molec_t4_b5n, molec_t4_b5f, nrow=2, left = "S/P Ratio")
invisible(dev.off())

#Barns 6, 1, 4: All Experimental ----
molec_t4_b614n <-ggplot(t4_b614, aes(dpi, sp, color = species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "Barn C: Nasal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols_sp) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                          axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

molec_t4_b614f <-ggplot(t4_b614, aes(dpi, sp_f, color = species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Rectal", x = "Days Post Infection", y = "") +
  scale_color_manual(name = "Barn C: Rectal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial4_cols_sp) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                          axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

grid.arrange(molec_t4_b614n, molec_t4_b614f, nrow=2, left = "S/P Ratio")


# qRTPCR ----
# Created in the combined section below

#####################
#  Trial 5 ----
#####################

t5_b5 <- subset(trial5_molec, barn == 5) # box A Inoculated cattle to sentinel goats
t5_b6 <- subset(trial5_molec, barn == 6) # box B Inoculated cattle to sentinel goats
t5_b2 <- subset(trial5_molec, barn == 2) # box C Inoculated cattle to sentinel goats
t5_b1 <- subset(trial5_molec, barn == 1) # box D Inoculated cattle to sentinel goats
t5_b3 <- subset(trial5_molec, barn == 3) # positive controls Box E
#t5_b4 <- subset(trial5_molec, barn == 4)  # negative controls Box F, not sampled
t5_b1256 <- subset(trial5_molec, barn == 5 | barn == 6 | barn == 2 | barn == 1) # all experimental together

trial5_cols <- c("1a" = "purple4", "2a" = "purple4", "5a" = "purple4", "6a" = "purple4",
          "1b" = "cyan3", "2b" = "cyan3", "5b" = "cyan3", "6b" = "cyan3",
          "3a" = "cyan4", "4b" = "cyan3") 

trial5_cols_sp <- c("Bovine" = "purple4", "Caprine" = "cyan3")

# icELISA ----
# Barn 5 ----
molec_t5_b5n <- ggplot(t5_b5, aes(dpi, sp, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "Barn 5: Nasal Status", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8),
                          axis.text=element_text(size=6))

molec_t5_b5f <- ggplot(t5_b5, aes(dpi, sp_f, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Rectal", x = "Days Post Infection", y = "") +
  scale_color_manual(name = "Barn 5: Rectal Status", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8),
                          axis.text=element_text(size=6))

grid.arrange(molec_t5_b5n, molec_t5_b5f, nrow=2, left = "S/P Ratio")


# Barn 6 ----
molec_t5_b6n <-ggplot(t5_b6, aes(dpi, sp, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "Barn 6: Nasal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                          axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

molec_t5_b6f <-ggplot(t5_b6, aes(dpi, sp_f, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Rectal", x = "", y = "") +
  scale_color_manual(name = "Barn 6: Rectal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                          axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

grid.arrange(molec_t5_b6n, molec_t5_b6f, nrow=2, left = "S/P Ratio")

#Barn 2 ----
molect5_b2n <-ggplot(t5_b2, aes(dpi, sp, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "Barn 2: Nasal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                          axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

molect5_b2f <-ggplot(t5_b2, aes(dpi, sp_f, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Rectal", x = "Days Post Infection", y = "") +
  scale_color_manual(name = "Barn 2: Rectal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                          axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

grid.arrange(molect5_b2n, molect5_b2f, nrow=2, left = "S/P Ratio")

#Barn 1 ----
molect5_b1n <-ggplot(t5_b1, aes(dpi, sp, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "Barn 1: Nasal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                          axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

molect5_b1f <-ggplot(t5_b1, aes(dpi, sp_f, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Rectal", x = "Days Post Infection", y = "") +
  scale_color_manual(name = "Barn 1: Rectal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                          axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

grid.arrange(molect5_b1n, molect5_b1f, nrow=2, left = "S/P Ratio")

#Barn 3: Positive Controls ----
jpeg("output/molec/trial5_molecular_icELISA_poscontrols_FigS12A.jpeg", width = 3.5, height = 3.25, units = "in", quality = 100, res = 600)
molec_t5_b3n <-ggplot(t5_b3, aes(dpi, sp, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  #stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "Barn 5: Nasal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

molecbox_t5_b3f <-ggplot(t5_b3, aes(dpi, sp_f, color = barni)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  #stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Rectal", x = "Days Post Infection", y = "") +
  scale_color_manual(name = "Barn 5: Rectal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

grid.arrange(molec_t5_b3n, molecbox_t5_b3f, nrow=2, left = "S/P Ratio")
invisible(dev.off())


#Barns 1, 2, 5, 6 : All Experimental ----
molec_t5_b1256n <-ggplot(t5_b1256, aes(dpi, sp, color = species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "Nasal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols_sp) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                          axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

molec_t5_b1256f <-ggplot(t5_b1256, aes(dpi, sp_f, color = species)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barni), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Rectal", x = "Days Post Infection", y = "") +
  scale_color_manual(name = "Rectal Status", 
                     labels = c("Goats (Inoculated)", "Goats (Sentinel)"),
                     values = trial5_cols_sp) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                          axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

grid.arrange(molec_t5_b1256n, molec_t5_b1256f, nrow=2, left = "S/P Ratio")
grid.arrange(molec_t5_b1256n, molec_t5_b1256f, ncol=2, left = "S/P % of Positive Control")

# qRTPCR ----
# Created in the combined section below

#################################
#  Combined Trials 1-2 ----
#################################

# Combined molecular results to look at biphasic whole blood response


trial1_molec_sub <- trial1_molec %>% select(-retest, -isolation1_dpbs, -i1o, -i1n, -i1r, -isolation2_dmem, -i2o, -i2n, -i2r)
trial1_molec_sub <- trial1_molec_sub %>% add_column(innoc =1, .before = "o1") 
trial1_molec_sub <- trial1_molec_sub %>% add_column(barni = NA, .before = "o1") 
trial1_molec_sub[trial1_molec_sub$barn == "1",]$barni <- "1a"
trial1_molec_sub[trial1_molec_sub$barn == "3",]$barni <- "3a"
trial1_molec_sub[trial1_molec_sub$barn == "5",]$barni <- "5a"

trial2_molec_sub <- trial2_molec %>% select(-od, -sp, -status, -note, -retest, -isolation1_dpbs, -i1o, -i1n, -i1r, -notes)


#remove controls + any other barns put datasets together
trial1_molec_subt_nc <- subset(trial1_molec_sub, barn != 1) 
trial2_molec_sub_nc <- subset(trial2_molec_sub, barn != 5) 

jointdat12_mol <- rbind(trial1_molec_subt_nc, trial2_molec_sub_nc)
jointdat12_mol <- jointdat12_mol %>% mutate(barniv2 = case_when(barni == "3a" ~ "Si1",
                                                            barni == "6a" ~ "Si1",
                                                            barni == "6b" ~ "Si0",
                                                            barni == "5a" ~ "Gi1",
                                                            barni == "4a" ~ "Gi1",
                                                            barni == "4b" ~ "Gi0"))
comb12_cols <- c("Gi1" = "cyan4", "Gi0" = "cyan3",
                 "Si1" = "deeppink2", "Si0" = "rosybrown3")

molec_o <- ggplot(jointdat12_mol, aes(dpi, 35-omean, color = barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barniv2), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Ocular", x = "", y = "qRT-PCR 35 - mean Ct") +
  scale_color_manual(name = "",
                     values = comb12_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

molec_n <- ggplot(jointdat12_mol, aes(dpi, 35-nmean, color = barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barniv2), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "",
                     values = comb12_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

molec_f <- ggplot(jointdat12_mol, aes(dpi, 35-rmean, color = barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barniv2), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Rectal", x = "", y = "") +
  scale_color_manual(name = "",
                     values = comb12_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

molec_wb <- ggplot(jointdat12_mol, aes(dpi, 35-wbmean, color = barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barniv2), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Whole Blood", x = "", y = "") +
  scale_color_manual(name = "",
                     values = comb12_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "top", axis.title = element_text(size = 15), axis.text = element_text(size = 12))


# by species
molec_wb_g <- ggplot(jointdat12_mol[jointdat12_mol$species == "G",], aes(dpi, 35-wbmean, color = barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barniv2), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Whole Blood", x = "", y = "") +
  scale_color_manual(name = "",
                     values = comb12_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

molec_wb_s <- ggplot(jointdat12_mol[jointdat12_mol$species == "S",], aes(dpi, 35-wbmean, color = barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barniv2), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Whole Blood", x = "", y = "") +
  scale_color_manual(name = "",
                     values = comb12_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12), 
        axis.text.y = element_blank())

jpeg("output/molec/trials12_molecular_qRTPCR_biphasic_FigS7.jpeg", width = 6, height = 4, units = "in", quality = 100, res = 600)
molec_wb_s2 <- molec_wb_s + labs(y = "", title = "")
grid.arrange(molec_wb_g, molec_wb_s2, ncol=2, 
             left = textGrob("35 - mean qRTPCR Ct",rot = 90, gp=gpar(fontsize=15)), 
             bottom = textGrob("Days Post Infection", gp=gpar(fontsize=15)))
invisible(dev.off())


#################################
#  Combined Trials 3-5 ----
#################################

#remove controls + any other barns put datasets together
trial3_molec_nc <- subset(trial3_molec, barn == 3) 
trial4_molec_nc <- subset(trial4_molec, barn != 5 & barn != 3) 
trial5_molec_nc <- subset(trial5_molec, barn != 3 & barn != 4) 

jointdat_mol <- rbind(trial3_molec_nc, trial4_molec_nc, trial5_molec_nc)
jointdat_mol <- jointdat_mol %>% mutate(barniv2 = case_when(barni == "3a" ~ "Ci1",
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

# icELISA ----
molec_n_ag <- ggplot(jointdat_mol, aes(dpi, sp, color = barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barniv2), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = comb345_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

molec_f_ag <- ggplot(jointdat_mol, aes(dpi, sp_f, color = barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barniv2), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(-5, 515)) +
  labs(title = "Rectal", x = "", y = "") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = comb345_cols) +
  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

jpeg("output/molec/trials345_molecular_icELISA_FigS13.jpeg", width = 8, height = 3.25, units = "in", quality = 100, res = 600)
grid.arrange(molec_n_ag, molec_f_ag, ncol=2, left = "S/P Ratio", bottom = "Days Post Infection")
invisible(dev.off())


# qRTPCR ----
molec_o <- ggplot(jointdat_mol, aes(dpi, 35-omean, color = barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barniv2), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Ocular", x = "", y = "qRT-PCR 35 - mean Ct") +
  scale_color_manual(values = comb345_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

molec_n <- ggplot(jointdat_mol, aes(dpi, 35-nmean, color = barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barniv2), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = comb345_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))

molec_f <- ggplot(jointdat_mol, aes(dpi, 35-rmean, color = barniv2)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(aes(group=barniv2), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,35), breaks = c(0, 4, 7, 10, 14, 17, 21, 24, 28, 32, 35)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Rectal", x = "", y = "") +
  scale_color_manual(name = "", 
                     labels = c("Cattle (Inoculated)", "Goats (Sentinel)"),
                     values = comb345_cols) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 12))


jpeg("output/molec/trials345_molecular_qRTPCR_Fig2D.jpeg", width = 9, height = 3.25, units = "in", quality = 100, res = 600)
grid.arrange(molec_o, molec_n, molec_f, ncol=3)
invisible(dev.off())

