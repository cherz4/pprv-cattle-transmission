# Empirical and model-based evidence for a negligible role of cattle in 
# peste des petits ruminants transmission and eradication

# Catherine M. Herzog#*, Fasil Aklilu#, Demeke Sibhatu, Dereje Shegu, 
# Redeat Belaineh, Abde Aliy Mohammed, Menbere Kidane,  Claudia Schulz,
# Brian J. Willett, Sarah Cleaveland, Dalan Bailey, Andrew R. Peters, 
# Isabella M. Cattadori, Peter J. Hudson, Hagos Asgedom, Joram Buza, 
# Mesfin Sahle Forza, Tesfaye Rufael Chibssa, Solomon Gebredufe, Nick Juleff, 
# Ottar N. Bjørnstad, Michael D. Baron, Vivek Kapur*

# Code for: 2 species SIR model with vaccination + transmission asymmetries
# Code authors: Ottar N. Bjørnstad + Catherine Herzog

#######################################################
# Next Generation Matrix function & example usage
######################################################
require(deSolve)

nextgenR0=function(Istates, Flist, Vlist, params, dfe){
  paras = as.list(c(dfe, params)) 
  
  k=0
  vl=fl=list(NULL)
  for(i in 1:length(Istates)){
    assign(paste("f", i, sep = "."), lapply(lapply(Flist,deriv, Istates[i]), eval, paras))
    assign(paste("v", i, sep = "."), lapply(lapply(Vlist,deriv, Istates[i]), eval, paras))
    for(j in 1:length(Istates)){
      k=k+1
      fl[[k]]=attr(eval(as.name(paste("f", i, sep=".")))[[j]], "gradient")[1,]
      vl[[k]]=attr(eval(as.name(paste("v", i, sep=".")))[[j]], "gradient")[1,]
    }
  }
  
  f=matrix(as.numeric(as.matrix(fl)[,1]), ncol=length(Istates))
  v=matrix(as.numeric(as.matrix(vl)[,1]), ncol=length(Istates))
  R0=max(eigen(f%*%solve(v))$values)
  return(R0)
}


# Recipe for use 
# Step 1: All states
istates=c("I1", "I2")

# Step 2: All new infections: 
flist=c(dI1dt=quote(beta11 * S1 *(1-p1) * I1 / N1 + beta21 * S1 *(1-p1)* I2 / N1), dI2dt=quote(beta12 * S2*(1-p2) * I1 / N2 + beta22 * S2*(1-p2)* I2 / N2))

#Step 3-5
#All losses 
Vm1=quote(mu1 * I1 + gamma1 * I1)
Vm2=quote(mu2 * I2 + gamma2 * I2)
#All gained transfers$
Vp1=0
Vp2=0
#Subtract Vp from Vm
V1=substitute(a-b, list(a=Vm1, b=Vp1))
V2=substitute(a-b, list(a=Vm2, b=Vp2))
#Make Vlist
vlist = c(V1,V2)

#Define list of parameter vectors:
# This is done for each figure below

#Specify disease-free equilibrium
df = list(S1 = 1, S2 = 1, I1 = 0, I2 = 0)

#Call nextgen function


########################################
# Figure 3 usage 
########################################

# Figure Margin Setup
# set working directory
setwd("C:/Users/Catherine Herzog/Box Sync/BMGF_PPRVTransmissionTrial/Manuscripts/")
jpeg("Fig4_ModelingSpillback_03142023_notitle_xaxislimited_10d.jpeg", width = 12, height = 10, units = "in", quality = 100, res = 600)
par(mfrow = c(2,2),
    oma = c(4,6.5,0,1) + 0.1,
    mar = c(1,2,1,1) + 0.1)

# # Transmission rates - include rates from only 1 infectious period per run of this script
# # 14 day infectious period
# beta_rate_ctg = 0.002132355  # (n=32 goats)
# beta_rate_ctc = 0.002132355  # no empirical data, have assumed same as C->G
# beta_rate_gtg = 0.1885041    # (n=12 goats)
# beta_rate_gtc = 0.09902103    # (n=2 cattle)

# 10 day infectious period
beta_rate_ctg = 0.002985296  # (n=32 goats)
beta_rate_ctc = 0.002985296  # no empirical data, have assumed same as C->G
beta_rate_gtg = 0.2639057    # (n=12 goats)
beta_rate_gtc = 0.1386294    # (n=2 cattle)
# 
# # 8 day infectious period
# beta_rate_ctg = 0.00373162  # (n=32 goats)
# beta_rate_ctc = 0.00373162  # no empirical data, have assumed same as C->G
# beta_rate_gtg = 0.3298822    # (n=12 goats)
# beta_rate_gtc = 0.1732868    # (n=2 cattle)


#Mortality rates (converted to daily rates) 
# Small ruminant rate from literature. Inverse of rate gives life span ~ 2.15 years.
# Yitagesu et al 2022 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9514490
# Cattle mortality rate order of magnitude smaller than small ruminants; inverse gives 10 year life span
mu_s = ((0.629+.302)/2)/365 # mortality rate (mean of kids and adults) per animal year 
mu_c = .10/365

# Recovery rate
# Recover usually occurs by ~ 14 days
gamma_s = 1/14
gamma_c = 1/14

########################################
# Symmetric spill forward and spillback 
########################################

parms  = list(mu = c(mu_s, mu_c), 
              N = c(1,1), 
              beta =  matrix(c(beta_rate_gtg, beta_rate_gtg, beta_rate_gtg, beta_rate_gtg), ncol=2, byrow=TRUE), 
              gamma = c(gamma_s, gamma_c),
              p=c(0, 0))
#Extract to work with nextgenR0:
para = list(mu1 = parms$mu[1],
            mu2 = parms$mu[2],
            beta11 = parms$beta[1,1],
            beta12 = parms$beta[1,2],
            beta21 = parms$beta[2,1],
            beta22 = parms$beta[2,2],
            gamma1 = parms$gamma[1],
            gamma2 = parms$gamma[2],
            p1 = parms$p[1],
            p2 = parms$p[2],
            N1 = parms$N[1],
            N2 = parms$N[2])

beta2=seq(0, round(beta_rate_gtg, 2), length.out = 101)
p=seq(0,1, by=0.01)  
RE=matrix(NA, ncol=length(beta2), nrow=length(p))
for(i in 1:length(beta2)){
  for(j in 1:length(p)){
    para$p1=p[j]
    para$beta21=beta2[i]
    RE[i,j]=nextgenR0(Istates=istates, Flist=flist, Vlist=vlist, params=para, dfe=df)
  }
}

# plot
contour(beta2, p, RE, levels = seq(0,9,1), lwd = 2, vfont = c("sans serif", "bold"), labcex = 1.25, 
        #main = "Symmetric", 
        xaxt = "n", yaxt = "n", cex.axis = 1.5)
axis(side = 1, at = seq(0,max(beta2), by =0.2), labels = FALSE)
axis(side = 2, at = seq(0,1,by =0.2), labels = c(0, 20, 40, 60, 80, 100), cex.axis = 1.5)
abline(v=beta_rate_ctg, lty = "dashed")
abline(v=beta_rate_gtg, lty = "dotdash")


########################################################
# Species Specific Symmetric spill forward and spillback 
########################################################

# with 3x more from shoats to cattle than back
parms  = list(mu = c(mu_s, mu_c), 
              N = c(1,1), 
              beta =  matrix(c(beta_rate_gtg, beta_rate_gtg, beta_rate_ctg, beta_rate_ctg), ncol=2, byrow=TRUE), 
              gamma = c(gamma_s, gamma_c),
              p=c(0, 0))
#Extract to work with nextgenR0:
para = list(mu1 = parms$mu[1],
            mu2 = parms$mu[2],
            beta11 = parms$beta[1,1],
            beta12 = parms$beta[1,2],
            beta21 = parms$beta[2,1],
            beta22 = parms$beta[2,2],
            gamma1 = parms$gamma[1],
            gamma2 = parms$gamma[2],
            p1 = parms$p[1],
            p2 = parms$p[2],
            N1 = parms$N[1],
            N2 = parms$N[2])

beta2=seq(0, round(beta_rate_gtg, 2), length.out = 101)
p=seq(0,1, by=0.01)  
RE=matrix(NA, ncol=length(beta2), nrow=length(p))
for(i in 1:length(beta2)){
  for(j in 1:length(p)){
    para$p1=p[j]
    para$beta21=beta2[i]
    RE[i,j]=nextgenR0(Istates=istates, Flist=flist, Vlist=vlist, params=para, dfe=df)
  }
}

# plot
contour(beta2, p, RE, levels = seq(0,9,1), lwd = 2, vfont = c("sans serif", "bold"), labcex = 1.25, 
        #main = "Species-Specific Symmetric",
        xaxt = "n", yaxt = "n", cex.axis = 1.5)
axis(side = 1, at = seq(0, max(beta2),by =0.2), labels = FALSE)
axis(side = 2, at = seq(0,1,by =0.2), labels = FALSE, cex.axis = 1.5)
abline(v=beta_rate_ctg, lty = "dashed")
abline(v=beta_rate_gtg, lty = "dotdash")


#######################################
#Asymmetric spill forward 
#######################################
parms  = list(mu = c(mu_s, mu_c),
              N = c(1,1),
              beta =  matrix(c(beta_rate_gtg, beta_rate_ctg, beta_rate_ctg, beta_rate_ctg), ncol=2, byrow=TRUE),
              gamma = c(gamma_s, gamma_c),
              p=c(0, 0))
#Extract to work with nextgenR0:
para = list(mu1 = parms$mu[1],
            mu2 = parms$mu[2],
            beta11 = parms$beta[1,1],
            beta12 = parms$beta[1,2],
            beta21 = parms$beta[2,1],
            beta22 = parms$beta[2,2],
            gamma1 = parms$gamma[1],
            gamma2 = parms$gamma[2],
            p1 = parms$p[1],
            p2 = parms$p[2],
            N1 = parms$N[1],
            N2 = parms$N[2])

beta2=seq(0, round(beta_rate_gtg, 2), length.out = 101)
p=seq(0,1, by=0.01)  
RE=matrix(NA, ncol=length(beta2), nrow=length(p))
for(i in 1:length(beta2)){
  for(j in 1:length(p)){
    para$p1=p[j]
    para$beta21=beta2[i]
    RE[i,j]=nextgenR0(Istates=istates, Flist=flist, Vlist=vlist, params=para, dfe=df)
  }
}

#plot
contour(beta2, p, RE, levels = seq(0,9,1), lwd = 2, vfont = c("sans serif", "bold"), labcex = 1.25, 
        #main = "Asymmetric Spill Forward", 
        yaxt = "n", cex.axis = 1.5)
axis(side = 2, at = seq(0,1,by =0.2), labels = c(0, 20, 40, 60, 80, 100), cex.axis = 1.5)
abline(v=beta_rate_ctg, lty = "dashed")
abline(v=beta_rate_gtg, lty = "dotdash")


#######################################
#Asymmetric spill forward and spillback 
#######################################
parms  = list(mu = c(mu_s, mu_c), 
              N = c(1,1), 
              beta =  matrix(c(beta_rate_gtg, beta_rate_gtc, beta_rate_ctg, 0), ncol=2, byrow=TRUE), 
              gamma = c(gamma_s, gamma_c),
              p=c(0, 0))
#Extract to work with nextgenR0:
para = list(mu1 = parms$mu[1],
            mu2 = parms$mu[2],
            beta11 = parms$beta[1,1],
            beta12 = parms$beta[1,2],
            beta21 = parms$beta[2,1],
            beta22 = parms$beta[2,2],
            gamma1 = parms$gamma[1],
            gamma2 = parms$gamma[2],
            p1 = parms$p[1],
            p2 = parms$p[2],
            N1 = parms$N[1],
            N2 = parms$N[2])

beta2=seq(0, round(beta_rate_gtg, 2), length.out = 101)
p=seq(0,1, by=0.01)  
RE=matrix(NA, ncol=length(beta2), nrow=length(p))
for(i in 1:length(beta2)){
  for(j in 1:length(p)){
    para$p1=p[j]
    para$beta21=beta2[i]
    RE[i,j]=nextgenR0(Istates=istates, Flist=flist, Vlist=vlist, params=para, dfe=df)
  }
}


#plot
contour(beta2, p, RE, levels = seq(0,9,1), lwd = 2, vfont = c("sans serif", "bold"), labcex = 1.25, 
        #main = "Asymmetric Spill Forward & Spillback", 
        yaxt = "n", cex.axis = 1.5)
axis(side = 2, at = seq(0,1,by =0.2), labels = FALSE, cex.axis = 1.5)
abline(v=beta_rate_ctg, lty = "dashed")
abline(v=beta_rate_gtg, lty = "dotdash")


# Overall plot title and labels
title(xlab = expression(paste("Cattle-to-Goat Transmission Rate (", beta[CS], ")")),
      ylab = "Small Ruminant Population Immunity (%)",
      outer = TRUE,
      line = 3,
      cex.lab = 2)

invisible(dev.off())
