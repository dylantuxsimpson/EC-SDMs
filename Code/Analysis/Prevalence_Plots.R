library(ggplot2)
library(bayesplot)

dat.sum <- function(x){
  m <- mean(x)
  hdi <- HDIofMCMC(x)
  low <- hdi[1]
  high <- hdi[2]
  return(c(y = m, ymin = low, ymax = high))
}

theme_dts <- function (base_size = 14, base_family = "serif"){
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.text = element_text(size = 12), 
          axis.ticks = element_line(colour = "black"), 
          legend.key = element_rect(colour = "grey80"), 
          panel.background = element_blank(),
          panel.border = element_blank(), axis.line.x.bottom = element_line(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          strip.background = element_rect(fill = "grey80", colour = "grey50", 
                                          size = 0.2),
          plot.background = element_blank()
    )
}

pdf(file = "Prevalence_modSel_HabFrag.pdf", width = 7, height = 3)
barplot(table(prev.modSel.dat$I1)/320, ylim = c(0, 40), xlab = "Indicator value", ylab = "Percent of inclusions", col = "grey90")
dev.off()
pdf(file = "Prevalence_modSel_EG.pdf", width = 3, height = 3)
barplot(table(prev.modSel.dat$I2)/320, ylim = c(0,80), xlab = "Indicator value", ylab = "Percent of inclusions", col = "grey90")
dev.off()

# Global model:
pdf(file = "Prevalence_Betas.pdf", width= 7, height = 3.5)
ggplot(data = modDat.global) + 
  geom_violin(aes(x = as.factor(5), y = beta.1.), fill = 'grey90') +
  geom_violin(aes(x = as.factor(4), y = beta.2.), fill = 'grey90') +
  geom_violin(aes(x = as.factor(3), y = beta.3.), fill = 'grey90') +
  geom_violin(aes(x = as.factor(2), y = beta.4.), fill = 'grey90') + 
  geom_violin(aes(x = as.factor(1), y = beta.5.), fill = 'grey90') +
  stat_summary(aes(x = as.factor(5), y = beta.1.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor(4), y = beta.2.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor(3), y = beta.3.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor(2), y = beta.4.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor(1), y = beta.5.), fun.data = dat.sum) +
  labs(x = "Variable", y = "Coefficient") +
  scale_x_discrete(labels = c("5" = "Forest cover", 
                              "4" = "Forest edge", 
                              "3" = expression("Forest cover" %*% "Edge"),
                              "2" = "Evergreen cover",
                              "1" = "Distance to water")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  coord_flip() +
  theme_dts() +
  theme(panel.border = element_blank(), axis.line.x.bottom = element_line())
dev.off()


#Variance
pdf(file = "Prevalence_Variance.pdf", width = 7, height = 2)
ggplot(data = modDat.global) +
  geom_violin(aes(x = factor(1), y = sigma), fill = 'grey90') +
  stat_summary(aes(x = factor(1), y = sigma), fun.data = dat.sum) +
  geom_violin(aes(x = factor(2), y = sigma0), fill = 'grey90') +
  stat_summary(aes(x = factor(2), y = sigma0), fun.data = dat.sum) +
  labs(x = "Variance \n Parameter", y = "Standard deviation") +
  scale_x_discrete(labels = c("1" = expression(sigma["plot"]),
                              "2" = expression(sigma["site"]))) +
  coord_flip() +
  theme_dts()
dev.off()


# Annual prevalence:
pdf(file = "Prevalence_Intercepts.pdf", width = 7, height = 4)
ggplot() + 
  geom_violin(data = modDat.global, aes(x = as.factor(6), y = plogis(mu)*100), fill = 'grey90') +
  geom_violin(data = annual.prev, aes(x = as.factor(5), y = prev.2012), fill = 'grey90') +
  geom_violin(data = annual.prev, aes(x = as.factor(4), y = prev.2013), fill = 'grey90') +
  geom_violin(data = annual.prev, aes(x = as.factor(3), y = prev.2015), fill = 'grey90') +
  geom_violin(data = annual.prev, aes(x = as.factor(2), y = prev.2016), fill = 'grey90') + 
  geom_violin(data = annual.prev, aes(x = as.factor(1), y = prev.2017), fill = 'grey90') +
  stat_summary(data = modDat.global, aes(x = as.factor(6), y = plogis(mu)*100), fun.data = dat.sum) +
  stat_summary(data = annual.prev, aes(x = as.factor(5), y = prev.2012), fun.data = dat.sum) +
  stat_summary(data = annual.prev, aes(x = as.factor(4), y = prev.2013), fun.data = dat.sum) +
  stat_summary(data = annual.prev, aes(x = as.factor(3), y = prev.2015), fun.data = dat.sum) +
  stat_summary(data = annual.prev, aes(x = as.factor(2), y = prev.2016), fun.data = dat.sum) +
  stat_summary(data = annual.prev, aes(x = as.factor(1), y = prev.2017), fun.data = dat.sum) +
  labs(x = "Year", y = "Prevalence (%)") +
  scale_x_discrete(labels = c("6" = "Mean",
                              "5" = "2012", 
                              "4" = "2013", 
                              "3" = "2015",
                              "2" = "2016",
                              "1" = "2017")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  coord_flip() +
  theme_dts()
dev.off()

# Wireframe surface plot to display interaction effect:

library(lattice)
library(geoR)

prevPred <- function(x){
  -.84*x[1] - .81*x[2]
}

prevPred.dat <- expand.grid(list(x = seq(-1, 1, 0.05), y = seq(-1, 1, 0.05)))

prevPred.dat$z <- apply(prevPred.dat[,1:2], 1, prevPred)

pdf(file = "Prevalence_predSurface.pdf", width =  4, height = 4)
wireframe(z ~ x*y, data = prevPred.dat,
          screen = list(z = -115, x = -70, y = 0),
          xlab = "Forest cover",
          ylab = "Forest edge"
)
dev.off()