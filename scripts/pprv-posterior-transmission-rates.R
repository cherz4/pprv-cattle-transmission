# Empirical and model-based evidence for a negligible role of cattle in 
# peste des petits ruminants transmission and eradication

# Catherine M. Herzog#*, Fasil Aklilu#, Demeke Sibhatu, Dereje Shegu, 
# Redeat Belaineh, Abde Aliy Mohammed, Menbere Kidane,  Claudia Schulz,
# Brian J. Willett, Sarah Cleaveland, Dalan Bailey, Andrew R. Peters, 
# Isabella M. Cattadori, Peter J. Hudson, Hagos Asgedom, Joram Buza, 
# Mesfin Sahle Forza, Tesfaye Rufael Chibssa, Solomon Gebredufe, Nick Juleff, 
# Ottar N. Bjørnstad, Michael D. Baron, Vivek Kapur*

# Code for: calculating transmission probabilities and rates and credible intervals
# Code authors: Catherine M. Herzog



library(bayesrules)
library(dplyr)
library(ggplot2)

set.seed(84735)

# transmission probability C->G 
# the prob that an exponentially distributed waiting time <  trial length (35 days)
# convert posterior mean from probability to rate using p(t) = 1-e^-rt or -ln(1-p(t))/t = r
plot_beta(1, 1) # visualize uninformative prior (uniform - special case of beta)
plot_beta_binomial(alpha = 1, beta = 1, y = 0, n = 32) # look at prior and likelihood based on experimental results
summarize_beta_binomial(alpha = 1, beta = 1, y = 0, n = 32)
post_mean_ctg <- summarize_beta_binomial(alpha = 1, beta = 1, y = 0, n = 32)[2,4]  # pull posterior mean from summaries of prior and posterior distribution
# rbeta(alpha + y, beta + n - y), here rbeta(1 + 0, 1 + 32 - 0) or rbeta(1,33) - the distribution the beta binomial is summarizing
qbeta(c(.025, 0.975), 1, 33) # 95% credible interval for C->G transmission probability
# 0.0007669121 0.1057628101

# Calculating rates
# 14d
# beta_rate_ctg = -log(1-post_mean_ctg)/14 # convert mean probability to a rate during the 35 day period of the C->G trials 
# beta_rate_ctg # 0.002132355
# sample_ctg <- rbeta(10000, 1, 33)
# rates_ctg <- -log(1-sample_ctg)/14
# quantile(rates_ctg, probs = c(.025, 0.975))
# # 2.5%         97.5% 
# # 0.0000491082 0.0081290883 

#10d
beta_rate_ctg = -log(1-post_mean_ctg)/10 # convert mean probability to a rate during the 35 day period of the C->G trials 
beta_rate_ctg # 0.002985296
sample_ctg <- rbeta(10000, 1, 33)
rates_ctg <- -log(1-sample_ctg)/10 
quantile(rates_ctg, probs = c(.025, 0.975))
# 2.5%       97.5% 
# 8.59845e-05 1.12721e-02 

#8d
# beta_rate_ctg = -log(1-post_mean_ctg)/8 # convert mean probability to a rate during the 35 day period of the C->G trials 
# beta_rate_ctg # 0.00373162
# sample_ctg <- rbeta(10000, 1, 33)
# rates_ctg <- -log(1-sample_ctg)/8 
# quantile(rates_ctg, probs = c(.025, 0.975))
# # 2.5%          97.5% 
# # 9.777366e-05  1.461837e-02 



# transmission probability G->G
plot_beta(1, 1) 
plot_beta_binomial(alpha = 1, beta = 1, y = 12, n = 12) 
summarize_beta_binomial(alpha = 1, beta = 1, y = 12, n = 12) 
post_mean_gtg <- summarize_beta_binomial(alpha = 1, beta = 1, y = 12, n = 12)[2,4]
# rbeta(alpha + y, beta + n - y), here rbeta(1 + 12, 1 + 12 - 12) or rbeta(13,1) - the distribution the beta binomial is summarizing
qbeta(c(.025, 0.975), 13, 1) # 95% credible interval for C->G transmission probability
# 0.7529474 0.9980544

# Calculating rates
# # 14d
# beta_rate_gtg = -log(1-post_mean_gtg)/14
# beta_rate_gtg #0.1885041
# sample_gtg <- rbeta(10000, 13, 1)
# rates_gtg <- -log(1-sample_gtg)/14
# quantile(rates_gtg, probs = c(.025, 0.975))
# # 2.5%      97.5% 
# # 0.1005720 0.4435633  

# 10d
beta_rate_gtg = -log(1-post_mean_gtg)/10
beta_rate_gtg #0.2639057
sample_gtg <- rbeta(10000, 13, 1)
rates_gtg <- -log(1-sample_gtg)/10
quantile(rates_gtg, probs = c(.025, 0.975))
# 2.5%      97.5% 
# 0.1437850 0.6335271 

# # 8d
# beta_rate_gtg = -log(1-post_mean_gtg)/8
# beta_rate_gtg #0.3298822
# sample_gtg <- rbeta(10000, 13, 1)
# rates_gtg <- -log(1-sample_gtg)/8
# quantile(rates_gtg, probs = c(.025, 0.975))
# # 2.5%      97.5% 
# # 0.1771663 0.7814859



# transmission probability G->C
plot_beta(1, 1)
plot_beta_binomial(alpha = 1, beta = 1, y = 2, n = 2)
summarize_beta_binomial(alpha = 1, beta = 1, y = 2, n = 2)
post_mean_gtc <- summarize_beta_binomial(alpha = 1, beta = 1, y = 2, n = 2)[2,4]
# rbeta(alpha + y, beta + n - y), here rbeta(1 + 2, 1 + 2 - 2) or rbeta(3,1) - the distribution the beta binomial is summarizing
qbeta(c(.025, 0.975), 3, 1) # 95% credible interval for C->G transmission probability
# 0.2924018 0.9915962

# Calculating rates
# # 14d
# beta_rate_gtc = -log(1-post_mean_gtc)/14
# beta_rate_gtc # 0.09902103
# sample_gtc <- rbeta(10000, 3, 1)
# rates_gtc <- -log(1-sample_gtc)/14
# quantile(rates_gtc, probs = c(.025, 0.975))
# # 2.5%       97.5% 
# # 0.02467719 0.33366332  

# 10d
beta_rate_gtc = -log(1-post_mean_gtc)/10
beta_rate_gtc # 0.1386294
sample_gtc <- rbeta(10000, 3, 1)
rates_gtc <- -log(1-sample_gtc)/10
quantile(rates_gtc, probs = c(.025, 0.975))
# 2.5%     97.5% 
# 0.034547 0.487640 

# # 8d
# beta_rate_gtc = -log(1-post_mean_gtc)/8
# beta_rate_gtc # 0.1732868
# sample_gtc <- rbeta(10000, 3, 1)
# rates_gtc <- -log(1-sample_gtc)/8
# quantile(rates_gtc, probs = c(.025, 0.975))
# # 2.5%       97.5% 
# # 0.04376008 0.58085127


# Comparison of G->G vs C->G transmission rate
beta_rate_gtg/beta_rate_ctg
# 88.40186  Goats have a ~88x greater transmission rate to other goats than cattle to goats

# Comparison of G->G vs G->C
beta_rate_gtg/beta_rate_gtc
# 1.903677