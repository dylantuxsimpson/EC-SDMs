## Set up:
library(ggplot2)
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
          panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          strip.background = element_rect(fill = "grey80", colour = "grey50", 
                                          size = 0.2),
          plot.background = element_blank(),
          axis.line.x.bottom = element_line()
    )
}

# Betas:
pdf(file = "Occupancy_Betas.pdf", width = 7, height = 5)
ggplot(data = OccDat.global) + 
  geom_violin(aes(x = as.factor("12"), y = beta.eps.1.), fill = 'grey90') +
  geom_violin(aes(x = as.factor("11"), y = beta.eps.2.), fill = 'grey90') +
  geom_violin(aes(x = as.factor("10"), y = beta.eps.3.), fill = 'grey90') +
  geom_violin(aes(x = as.factor("09"), y = beta.eps.4.), fill = 'grey90') + 
  geom_violin(aes(x = as.factor("08"), y = beta.eps.5.), fill = 'grey90') +
  geom_violin(aes(x = as.factor("07"), y = beta.gam.1.), fill = 'grey90') +
  geom_violin(aes(x = as.factor("06"), y = beta.gam.2.), fill = 'grey90') +
  geom_violin(aes(x = as.factor("05"), y = beta.gam.3.), fill = 'grey90') +
  geom_violin(aes(x = as.factor("04"), y = beta.gam.4.), fill = 'grey90') + 
  geom_violin(aes(x = as.factor("03"), y = beta.gam.5.), fill = 'grey90') +
  geom_violin(aes(x = as.factor("02"), y = beta.p.1.), fill = 'grey90') +
  geom_violin(aes(x = as.factor("01"), y = beta.p.2.), fill = 'grey90') +
  stat_summary(aes(x = as.factor("12"), y = beta.eps.1.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor('11'), y = beta.eps.2.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor('10'), y = beta.eps.3.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor('09'), y = beta.eps.4.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor('08'), y = beta.eps.5.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor('07'), y = beta.gam.1.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor('06'), y = beta.gam.2.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor('05'), y = beta.gam.3.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor('04'), y = beta.gam.4.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor('03'), y = beta.gam.5.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor('02'), y = beta.p.1.), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor('01'), y = beta.p.2.), fun.data = dat.sum) +
  labs(x = "Variable", y = "Coefficient") +
  scale_x_discrete(labels = c("12" = "Forest cover", 
                              "11" = "Forest edge", 
                              "10" = expression("Forest cover" %*% "Edge"),
                              "09" = "Evergreen cover",
                              "08" = "Distance to water",
                              "07" = "Forest cover", 
                              "06" = "Proximity", 
                              "05" = expression("Forest cover" %*% "Proximity"),
                              "04" = "Evergreen cover",
                              "03" = "Distance to water",
                              "02" = "Temperature",
                              "01" = "Slope"
  )) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  coord_flip() +
  theme_dts()
dev.off()

# Intercepts
## Prepare summary data:
mean.eps <- plogis(OccDat.global$alpha.e)
mean.gam <- plogis(OccDat.global$alpha.g)
mean.p <- plogis(OccDat.global$alpha.p)

means <- data.frame("eps" = mean.eps, "gam" <- mean.gam, "p" = mean.p)

pdf(file = "Occupancy_Intercepts.pdf", width = 7, height = 2.5)
ggplot(data = means) + 
  geom_violin(aes(x = as.factor(3), y = eps), fill = 'grey90') +
  geom_violin(aes(x = as.factor(2), y = gam), fill = 'grey90') + 
  geom_violin(aes(x = as.factor(1), y = p), fill = 'grey90') +
  stat_summary(aes(x = as.factor(3), y = eps), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor(2), y = gam), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor(1), y = p), fun.data = dat.sum) +
  labs(x = "Parameter", y = "Intercept") +
  scale_x_discrete(labels = c("3" = "Extinction",
                              "2" = "Colonization",
                              "1" = "Detection")) +
  geom_hline(yintercept = c(0,1), linetype = "dotted") +
  coord_flip() +
  theme_dts()
dev.off()

pdf(file = "Occupancy_Variance.pdf", width = 7, height = 3)
ggplot(data = OccDat.global) + 
  geom_violin(aes(x = as.factor(4), y = sigma.e), fill = 'grey90') +
  geom_violin(aes(x = as.factor(3), y = sigma.g), fill = 'grey90') + 
  geom_violin(aes(x = as.factor(2), y = sigma.p), fill = 'grey90') +
  geom_violin(aes(x = as.factor(1), y = det.error), fill = 'grey90') +
  stat_summary(aes(x = as.factor(4), y = sigma.e), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor(3), y = sigma.g), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor(2), y = sigma.p), fun.data = dat.sum) +
  stat_summary(aes(x = as.factor(1), y = det.error), fun.data = dat.sum) +
  labs(x = "Parameter", y = "Estimate") +
  scale_x_discrete(labels = c("4" = expression(sigma[epsilon]),
                              "3" = expression(sigma[gamma]),
                              "2" = expression(sigma[p]),
                              "1" = expression(varsigma)
                              )
                    ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  coord_flip() +
  theme_dts()
dev.off()

# Wireframe surface plot to display interaction effect:

library(lattice)
library(geoR)

eps.pred <- function(x){
  -.042*x[1] + 1.04*x[2] - 1.49*x[1]*x[2]
}

epsPred.dat <- expand.grid(list(x = seq(-1, 1, 0.05), y = seq(-1, 1, 0.05)))

epsPred.dat$z <- apply(epsPred.dat[,1:2], 1, eps.pred)

pdf(file = "Occupancy_PredictionSurface.pdf", width = 4, height = 4)
wireframe(z ~ x*y, data = epsPred.dat)
dev.off()

# Model selection:


pdf(file = "Occupancy_modSel_HabFrag_eps.pdf", width = 7, height = 4)
barplot(table(OccDat.modSel$ie1)/360, main = "Forest cover and fragmentation - epsilon", col = "grey90", ylim = c(0, 70))
dev.off()
pdf(file = "Occupancy_modSel_EG_eps.pdf", width = 4, height = 4)
barplot(table(OccDat.modSel$ie2)/360, main = "Proportion evergreen - epsilon", col = "grey90", ylim = c(0,70))
dev.off()

pdf(file = "Occupancy_modSel_HabFrag_gam.pdf", width = 7, height = 4)
barplot(table(OccDat.modSel$ig1)/360, main = "Forest cover and fragmentation - gamma", col = "grey90", ylim = c(0,70))
dev.off()
pdf(file = "Occupancy_modSel_EG_gam.pdf", width = 4, height = 4)
barplot(table(OccDat.modSel$ig2)/360, main = "Proportion evergreen - gamma", col = "grey90", ylim =c(0,70))
dev.off()
