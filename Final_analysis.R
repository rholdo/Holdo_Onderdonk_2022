rm(list = ls(all = TRUE))
library(raster)
library(rgdal)
library(dplyr)
library(ggplot2)
library(rgdal)
library(SDMTools)
library(nlme)
library(quantreg)
library(reshape2)
library(wvtool)
library(segmented)

df <- read.csv('Final_30m_data.csv')
# Cubic interpolation when converting to 30-m pixels resulted in some negative TC values
# Restrict minimum TC to 0
df$TC <- ifelse(df$TC < 0, 0, df$TC)
# 131 plots have trees
dmin <- aggregate(driver ~ ID, df, min)
IDsub <- subset(dmin, driver < 0.05)
dfriv <- df[df$ID %in% IDsub$ID, ]
# 113 plots are within 50 m (0.05 km) of a river
# Remove distance to river of 0 (i.e., pixels that may contain surface water)
dfriv <- subset(dfriv, driver > 0)

# Fig 2a
fig2a <- ggplot(data = dfriv) +
  geom_point(aes(x = driver, y = TC), size = 4, shape = 21) +
  labs(y = 'TC',
       x = expression(paste('d'[river],' (km)'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = rel(1)),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(vjust = 0, size = rel(2)),
        axis.text = element_text(size = rel(2)),
        axis.line = element_line(colour="black", size = 1.2),
        panel.border = element_rect(colour="black", size = 1.2),
        axis.ticks = element_line(colour="black", size = 1.2),
        axis.ticks.length = unit(.25, "cm"))

# Supplementary figure
figS1 <- ggplot(data = df) +
  geom_point(aes(x = driver, y = TC), size = 4, shape = 21) +
  labs(y = 'TC',
       x = expression(paste('d'[river],' (km)'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = rel(1)),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(vjust = 0, size = rel(2)),
        axis.text = element_text(size = rel(2)),
        axis.line = element_line(colour="black", size = 1.2),
        panel.border = element_rect(colour="black", size = 1.2),
        axis.ticks = element_line(colour="black", size = 1.2),
        axis.ticks.length = unit(.25, "cm"))

# Put pixels into distance to river bins of 30 m
dfriv$drbin <- cut(dfriv$driver, breaks = seq(from = 0, to = 1.2, by = 0.03),
                   labels = seq(from = 0.015, to = 1.185, by = 0.03))
dfriv.ag <- aggregate(log(TC + 0.01) ~ drbin + ID, dfriv, mean)
names(dfriv.ag)[3] <- 'logTC'
dfriv.mn <- aggregate(logTC ~ drbin, dfriv.ag, mean)
dfriv.sd <- aggregate(logTC ~ drbin, dfriv.ag, sd)
dfriv.ag2 <- merge(dfriv.mn, dfriv.sd, by = 'drbin')
names(dfriv.ag2) <- c('driver', 'Mean', 'SD')
dfriv.ag2$driver <- as.numeric(as.character(dfriv.ag2$driver))

# Fig 2b
fig2b <- ggplot(data = dfriv.ag2) +
  geom_errorbar(aes(x = driver, ymin = Mean - SD, ymax = Mean + SD), width = 0.02, size = 1) +
  geom_point(aes(x = driver, y = Mean), size = 4, shape = 21) +
  labs(y = 'logTC',
       x = expression(paste('d'[river],' (km)'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = rel(1)),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(vjust = 0, size = rel(2)),
        axis.text = element_text(size = rel(2)),
        axis.line = element_line(colour="black", size = 1.2),
        panel.border = element_rect(colour="black", size = 1.2),
        axis.ticks = element_line(colour="black", size = 1.2),
        axis.ticks.length = unit(.25, "cm"))

# Examine overall patterns of tree cover vs. distance to river
# Now test for overall effect of distance to river on TC using lme and gls
dfriv$Xcen <- (dfriv$X - mean(dfriv$X)) / sd(dfriv$X)
dfriv$Ycen <- (dfriv$Y - mean(dfriv$Y)) / sd(dfriv$Y)
dfriv$logTC <- log(dfriv$TC + 0.01)
# Note: this is very slow
riv.mod <- lme(logTC ~ dr, data = dfriv,
                       random = ~ 1 | ID,
                       corr = corSpatial(form = ~ Xcen + Ycen, type ="exponential", nugget = F), 
                       method = "ML", na.action = 'na.exclude')

# Test for more complex nonlinear patterns
# First, the cubic interpolation resulted in additional plots with TC values of 0
# Discard these to enable model convergence
df.ag <- aggregate(TC ~ ID, dfriv, sum)
df.ag <- subset(df.ag, TC > 0)
dfriv <- dfriv[dfriv$ID %in% df.ag$ID, ]

IDs <- unique(dfriv$ID)
for (i in 1:length(ID)){
  sub <- subset(dfriv, ID == IDs[i])
  # Center X and Y values
  sub$Xcen <- (sub$X - mean(sub$X)) / sd(sub$X)
  sub$Ycen <- (sub$Y - mean(sub$Y)) / sd(sub$Y)
  sub$dr2 <- sub$driver ^ 2
  sub$dr3 <- sub$driver ^ 3
  gls0 <- gls(logTC ~ 1, correlation=corExp(1, form=~Xcen+Ycen), data=sub, method="ML")
  gls1 <- gls(logTC ~ driver, correlation=corExp(1, form=~Xcen+Ycen), data=sub, method="ML")
  gls2 <- gls(logTC ~ driver + dr2, correlation=corExp(1, form=~Xcen+Ycen), data=sub, method="ML")
  gls3 <- gls(logTC ~ driver + dr2 + dr3, correlation=corExp(1, form=~Xcen+Ycen), data=sub, method="ML")
  dfr <- data.frame(IDs[i])
  dfr$AIC.null <- AIC(gls0)
  dfr$AIC.lin <- AIC(gls1)
  dfr$AIC.quad <- AIC(gls2)
  dfr$AIC.cub <- AIC(gls3)
  names(dfr)[1] <- 'ID'
  if (i == 1) dfspat <- dfr
  else dfspat <- rbind(dfspat, dfr)
  cat(i, '\n')
}
# The above loop is slow
# To skip, see results dataframe below
write.csv(dfspat, 'Dist_to_riv_spatial_regressions_poly.csv', 
          row.names = FALSE)
dfspat <- read.csv('Dist_to_riv_spatial_regressions_poly.csv')

# Now conduct segmented regressions
for (i in 1:length(IDs)){
  sub <- subset(dfriv, ID == IDs[i])
  if (sum(sub$TC) > 0){
    out.lm <- lm(TC ~ driver, data = sub)
    seg <- segmented.lm(out.lm, seg.Z = ~ driver, 
                        control=seg.control(stop.if.error=FALSE, n.boot=0, it.max=1000))
    sumseg <- summary(seg)
    psi <- sumseg$psi["psi1.driver", "Est."] # First break point
    dfr <- data.frame(IDs[i])
    if (!is.null(psi)){
      dfr$psi <- psi
      # Subset again to consider data before first breakpoint
      sub2 <- subset(sub, driver <= psi)
      # Center X and Y values
      sub2$Xcen <- (sub2$X - mean(sub2$X)) / sd(sub2$X)
      sub2$Ycen <- (sub2$Y - mean(sub2$Y)) / sd(sub2$Y)
      gls.exp <- gls(logTC ~ driver, correlation=corExp(1, form=~Xcen+Ycen), 
                     data=sub2, method="ML")
      sumexp <- summary(gls.exp)
      dfr$Slope <- coef(sumexp)[2]
      dfr$p.val <- coef(sumexp)["driver", "p-value"]
    }
    else {
      dfr$psi <- NA
      dfr$Slope <- NA
      dfr$p.val <- NA
    }
    names(dfr)[1] <- 'ID'
    if (i == 1) dfspatseg <- dfr
    else dfspatseg <- rbind(dfspatseg, dfr)
  }
  cat(i, '\n')
}
# Save results
write.csv(dfspatseg, 'Distance_to_river_cutoff.csv', row.names = FALSE)
dfspatseg <- read.csv('Distance_to_river_cutoff.csv')

# Plot histogram of psi values for cases with negative slopes
fig4b <- ggplot(subset(dfspatseg, Slope < 0), aes(psi)) + geom_histogram() +
  labs(y = 'Frequency',
       x = expression(paste(psi,' (km)'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = rel(1)),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(vjust = 0, size = rel(2)),
        axis.text = element_text(size = rel(2)),
        axis.line = element_line(colour="black", size = 1.2),
        panel.border = element_rect(colour="black", size = 1.2),
        axis.ticks = element_line(colour="black", size = 1.2),
        axis.ticks.length = unit(.25, "cm"))

# What proportion of plots have river effects?
dfreff1 <- subset(dfspatseg, p.val < 0.05 & Slope < 0)
dfreff2 <- subset(dfspatseg, Slope < 0)
# Psi estimate for all negative slopes
mean(dfreff2$psi, na.rm = TRUE) # Mean psi = 0.182 km
median(dfreff2$psi, na.rm = TRUE) # Median psi = 0.099 km
# Psi estimate for all significant negative slopes
mean(dfreff1$psi, na.rm = TRUE) # Mean psi = 140 km
median(dfreff1$psi, na.rm = TRUE) # Median psi = 0.077 km
# 46 of 104 plots have clear river effects with significant negative slopes
# 82 of 104 plots have estimated negative slopes in the first segment
# Only 2 plots of 104 have significant positive slopes
# Also mean tree cover near and far from river
# Measure tree cover influenced by river vs not
# psi for significant cases = 93.0 m; for all cases = 110.0 m
# Find the Psi cutoff for every plot at 30 m resolution
# Assign a type (river vs. non-river) to each 30 x 30 m pixel
dfreff2 <- dfreff2[, 1:2]
df <- merge(df, dfreff2, by = 'ID', all.x = TRUE)
df$Type <- ifelse(df$dist_to_river < df$psi, 'River', 'NonRiver')

# Case 1: plots with a psi estimate
dfsub <- subset(df, !is.na(df$psi))
ag.area1 <- aggregate(TC ~ ID + Type, dfsub, function(x){length(x)})
names(ag.area1)[3] <- 'Pixels'
ag.area1.w <- dcast(ag.area1, ID ~ Type, value.var = "Pixels")
ag.area1.w$Priver <- ag.area1.w$River / (ag.area1.w$River + ag.area1.w$NonRiver)
ag.area1.w$Priver <- ifelse(is.na(ag.area1.w$Priver), 0, ag.area1.w$Priver)
# Aggregate TC by distance to river
ag.tc1 <- aggregate(TC ~ ID + Type, dfsub, mean)
ag.tc1$TC <- ifelse(is.na(ag.tc1$TC), 0, ag.tc1$TC)
ag.tc1.w <- dcast(ag.tc1, ID ~ Type, value.var = "TC")
ag.tc1.w$NonRiver <- ifelse(is.na(ag.tc1.w$NonRiver), 0, ag.tc1.w$NonRiver)
ag.tc1.w$River <- ifelse(is.na(ag.tc1.w$River), 0, ag.tc1.w$River)
ag.comb1 <- merge(ag.area1.w, ag.tc1.w, by ='ID')
names(ag.comb1) <- c("ID", "NRpix", "Rpix", "Priver", "NRTC", "RTC")
ag.comb1$NRweightedTC <- (1 - ag.comb1$Priver) * ag.comb1$NRTC
ag.comb1$RweightedTC <- ag.comb1$Priver * ag.comb1$RTC
# Summary values
mean(ag.comb1$Priver) # 41.5 % of SNP is influenced by rivers
mean(ag.comb1$RweightedTC) / (mean(ag.comb1$RweightedTC) + mean(ag.comb1$NRweightedTC)) # 53.2 % of all trees are 'river trees'

df2 <- df
# Case 2: all plots (if no psi estimate, give median value of 111 m)
df$psi <- ifelse(is.na(df$psi), 110, df$psi)
df$Type <- ifelse(df$dist_to_river < df$psi, 'River', 'NonRiver')
ag.area2 <- aggregate(TC ~ ID + Type, df, function(x){length(x)})
names(ag.area2)[3] <- 'Pixels'
ag.area2.w <- dcast(ag.area2, ID ~ Type, value.var = "Pixels")
ag.area2.w$Priver <- ag.area2.w$River / (ag.area2.w$River + ag.area2.w$NonRiver)
ag.area2.w$Priver <- ifelse(is.na(ag.area2.w$Priver), 0, ag.area2.w$Priver)
# Aggregate TC by distance to river
ag.tc2 <- aggregate(TC ~ ID + Type, df, mean)
ag.tc2$TC <- ifelse(is.na(ag.tc2$TC), 0, ag.tc2$TC)
ag.tc2.w <- dcast(ag.tc2, ID ~ Type, value.var = "TC")
ag.tc2.w$NonRiver <- ifelse(is.na(ag.tc2.w$NonRiver), 0, ag.tc2.w$NonRiver)
ag.tc2.w$River <- ifelse(is.na(ag.tc2.w$River), 0, ag.tc2.w$River)
ag.comb2 <- merge(ag.area2.w, ag.tc2.w, by ='ID')
names(ag.comb2) <- c("ID", "NRpix", "Rpix", "Priver", "NRTC", "RTC")
ag.comb2$NRweightedTC <- (1 - ag.comb2$Priver) * ag.comb2$NRTC
ag.comb2$RweightedTC <- ag.comb2$Priver * ag.comb2$RTC
# Summary values
mean(ag.comb2$Priver) # 25.0 % of SNP is influenced by rivers
mean(ag.comb2$RweightedTC) / (mean(ag.comb2$RweightedTC) + mean(ag.comb2$NRweightedTC)) # 40.0 % of all trees are 'river trees'

# Case 3: all plots (if no psi estimate, give median value of 0 m, i.e., assume no river effect)
df2$psi <- ifelse(is.na(df2$psi), 0, df2$psi)
df2$Type <- ifelse(df2$dist_to_river < df2$psi, 'River', 'NonRiver')
ag.area3 <- aggregate(TC ~ ID + Type, df2, function(x){length(x)})
names(ag.area3)[3] <- 'Pixels'
ag.area3.w <- dcast(ag.area3, ID ~ Type, value.var = "Pixels")
ag.area3.w$Priver <- ag.area3.w$River / (ag.area3.w$River + ag.area3.w$NonRiver)
ag.area3.w$Priver <- ifelse(is.na(ag.area3.w$Priver), 0, ag.area3.w$Priver)
# Aggregate TC by distance to river
ag.tc3 <- aggregate(TC ~ ID + Type, df2, mean)
ag.tc3$TC <- ifelse(is.na(ag.tc3$TC), 0, ag.tc3$TC)
ag.tc3.w <- dcast(ag.tc3, ID ~ Type, value.var = "TC")
ag.tc3.w$NonRiver <- ifelse(is.na(ag.tc3.w$NonRiver), 0, ag.tc3.w$NonRiver)
ag.tc3.w$River <- ifelse(is.na(ag.tc3.w$River), 0, ag.tc3.w$River)
ag.comb3 <- merge(ag.area3.w, ag.tc3.w, by ='ID')
names(ag.comb3) <- c("ID", "NRpix", "Rpix", "Priver", "NRTC", "RTC")
ag.comb3$NRweightedTC <- (1 - ag.comb3$Priver) * ag.comb3$NRTC
ag.comb3$RweightedTC <- ag.comb3$Priver * ag.comb3$RTC
# Summary values
mean(ag.comb3$Priver) # 21.1 % of SNP is influenced by rivers
mean(ag.comb3$RweightedTC) / (mean(ag.comb3$RweightedTC) + mean(ag.comb3$NRweightedTC)) # 37.4 % of all trees are 'river trees'

# We do not include the raw high-resolution rasters for lacunarity calculations
# but provide all the lacunarity metrics at the three spatial scales used
ag.lac <- read.csv('Aggregated_lacunarity_data.csv')

# Find the functional form describing lacunearity with distance to river
# 10 m
gls0.10 <- gls(Lac10 ~ 1, correlation = corAR1(form = ~ 1 | ID), 
            data = ag.lac, method="ML", na.action = na.omit)
gls1.10 <- gls(Lac10 ~ dr, correlation = corAR1(form = ~ 1 | ID),  
            data = ag.lac, method="ML", na.action = na.omit)
gls2.10 <- gls(Lac10 ~ dr + dr2, correlation = corAR1(form = ~ 1 | ID),  
            data = ag.lac, method="ML", na.action = na.omit)
gls3.10 <- gls(Lac10 ~ dr + dr2 + dr3, correlation = corAR1(form = ~ 1 | ID),  
            data = ag.lac, method="ML", na.action = na.omit)
# 25 m
gls0.25 <- gls(Lac25 ~ 1, correlation = corAR1(form = ~ 1 | ID), 
               data = ag.lac, method="ML", na.action = na.omit)
gls1.25 <- gls(Lac25 ~ dr, correlation = corAR1(form = ~ 1 | ID),  
               data = ag.lac, method="ML", na.action = na.omit)
gls2.25 <- gls(Lac25 ~ dr + dr2, correlation = corAR1(form = ~ 1 | ID),  
               data = ag.lac, method="ML", na.action = na.omit)
gls3.25 <- gls(Lac25 ~ dr + dr2 + dr3, correlation = corAR1(form = ~ 1 | ID),  
               data = ag.lac, method="ML", na.action = na.omit)
# 50 m
gls0.50 <- gls(Lac50 ~ 1, correlation = corAR1(form = ~ 1 | ID), 
               data = ag.lac, method="ML", na.action = na.omit)
gls1.50 <- gls(Lac50 ~ dr, correlation = corAR1(form = ~ 1 | ID),  
               data = ag.lac, method="ML", na.action = na.omit)
gls2.50 <- gls(Lac50 ~ dr + dr2, correlation = corAR1(form = ~ 1 | ID),  
               data = ag.lac, method="ML", na.action = na.omit)
gls3.50 <- gls(Lac50 ~ dr + dr2 + dr3, correlation = corAR1(form = ~ 1 | ID),  
               data = ag.lac, method="ML", na.action = na.omit)

AIC(gls0.10, gls1.10, gls2.10, gls3.10)
AIC(gls0.25, gls1.25, gls2.25, gls3.25)
AIC(gls0.50, gls1.50, gls2.50, gls3.50)

# Aggregate for plotting
ag.mn <- aggregate(cbind(Lac10, Lac25, Lac50) ~ dr, ag.lac, mean)
ag.sd <- aggregate(cbind(Lac10, Lac25, Lac50) ~ dr, ag.lac, sd)
ag <- merge(ag.mn, ag.sd, by = 'dr')
names(ag) <- c('dr', 'Mean10', 'Mean25', 'Mean50', 'SD10', 'SD25', 'SD50')
ag$Pred10 <- coef(gls3.10)[1] + coef(gls3.10)[2] * ag$dr + coef(gls3.10)[3] * ag$dr ^ 2 + 
  coef(gls3.10)[4] * ag$dr ^ 3
ag$Pred25 <- coef(gls3.25)[1] + coef(gls3.25)[2] * ag$dr + coef(gls3.25)[3] * ag$dr ^ 2 + 
  coef(gls3.25)[4] * ag$dr ^ 3
ag$Pred50 <- coef(gls2.50)[1] + coef(gls2.50)[2] * ag$dr + coef(gls2.50)[3] * ag$dr ^ 2

# Plot lacunarity versus distance to river

fig3a <- ggplot(data = ag) +
  geom_errorbar(aes(x = dr, ymin = Mean10 - SD10, ymax = Mean10 + SD10), width = 0.05, size = 1) +
  geom_point(aes(x = dr, y = Mean10), size = 4, shape = 21) +
  geom_line(aes(x = dr, y = Pred10), size = 2, col = 'red') +
  labs(y = expression(paste('L'[norm],'(10)')),
       x = expression(paste('d'[river],' (km)'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = rel(1)),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(vjust = 0, size = rel(2)),
        axis.text = element_text(size = rel(2)),
        axis.line = element_line(colour="black", size = 1.2),
        panel.border = element_rect(colour="black", size = 1.2),
        axis.ticks = element_line(colour="black", size = 1.2),
        axis.ticks.length = unit(.25, "cm"))

fig3b <- ggplot(data = ag) +
  geom_errorbar(aes(x = dr, ymin = Mean25 - SD25, ymax = Mean25 + SD25), width = 0.05, size = 1) +
  geom_point(aes(x = dr, y = Mean25), size = 4, shape = 21) +
  geom_line(aes(x = dr, y = Pred25), size = 2, col = 'red') +
  labs(y = expression(paste('L'[norm],'(25)')),
       x = expression(paste('d'[river],' (km)'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = rel(1)),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(vjust = 0, size = rel(2)),
        axis.text = element_text(size = rel(2)),
        axis.line = element_line(colour="black", size = 1.2),
        panel.border = element_rect(colour="black", size = 1.2),
        axis.ticks = element_line(colour="black", size = 1.2),
        axis.ticks.length = unit(.25, "cm"))

fig3c <- ggplot(data = ag) +
  geom_errorbar(aes(x = dr, ymin = Mean50 - SD50, ymax = Mean50 + SD50), width = 0.05, size = 1) +
  geom_point(aes(x = dr, y = Mean50), size = 4, shape = 21) +
  geom_line(aes(x = dr, y = Pred50), size = 2, col = 'red') +
  labs(y = expression(paste('L'[norm],'(50)')),
       x = expression(paste('d'[river],' (km)'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = rel(1)),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(vjust = 0, size = rel(2)),
        axis.text = element_text(size = rel(2)),
        axis.line = element_line(colour="black", size = 1.2),
        panel.border = element_rect(colour="black", size = 1.2),
        axis.ticks = element_line(colour="black", size = 1.2),
        axis.ticks.length = unit(.25, "cm"))