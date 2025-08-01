# LOAD PACKAGES #
library(tidyverse)
library(lme4)
library(DHARMa)
library(olsrr)
library(car)
library(MASS)
library(lmerTest)
library(MuMIn)
library(emmeans)

# FUNCTION FOR BAR GRAPHS #
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean, na.rm=T)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd, na.rm=T)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
} 


BiomassNoduleAMF <- read.csv("Field_ANPP_NodNum.csv")

Nfix <- read.csv("Field_GClog.csv") %>%
  group_by(Site, Treatment, Plant_ID) %>%
  mutate(Ethylene_Area = ifelse(Ethylene_Area <1, 0, Ethylene_Area)) %>%
  summarise(Area = mean(Ethylene_Area, na.rm=T), AvgEthPPM = (100*Area)/204.6585)

LeafNuts <- read.csv("Field_LeafNutContent.csv")

data1 <- merge(Nfix, BiomassNoduleAMF, by=c("Site", "Treatment", "Plant_ID"), all=T)
data2 <- merge(data1, LeafNuts, by=c("Site", "Treatment", "Plant_ID", "Group"), all=T) %>%
  group_by(Site, Treatment, Plant_ID) %>%
  mutate(Rate = ((AvgEthPPM*0.01*10^6)/(24.45*0.75)))
  # rate is nmol C2H4 per hour

Summarized <- data2 %>%
  group_by(Treatment) %>%
  summarise(mANPP = mean(ANPP, na.rm=T),seANPP = sd(ANPP, na.rm = TRUE) / sqrt(sum(!is.na(ANPP))),
            mNod = mean(NodNum, na.rm=T),seNod = sd(NodNum, na.rm = TRUE) / sqrt(sum(!is.na(NodNum))),
            mAMF = mean(PercMC, na.rm=T),seAMF = sd(PercMC, na.rm = TRUE) / sqrt(sum(!is.na(PercMC))),
            mRate = mean(Rate, na.rm=T),seRate = sd(Rate, na.rm = TRUE) / sqrt(sum(!is.na(Rate))),
            mANPP = mean(ANPP, na.rm=T),seANPP = sd(ANPP, na.rm = TRUE) / sqrt(sum(!is.na(ANPP))),
            mN= mean(N, na.rm=T),seN = sd(N, na.rm = TRUE) / sqrt(sum(!is.na(N))),
            mP = mean(P, na.rm=T),seP = sd(P, na.rm = TRUE) / sqrt(sum(!is.na(P))))

# CHANGE SITE, REPLICATE, AND TREATMENT TO FACTOR #
data2$Site <- as.factor(data2$Site)
data2$Group <- as.factor(data2$Group)
data2$Treatment <- as.factor(data2$Treatment)



#################################
# Shoot Biomass #

## CHECK FOR NORMALITY, ETC ####
hist(data$ANPP)
res_ANPP <- lmer(ANPP ~ Treatment + (1|Site) + (1|Site:Group), data = data2)
resanpp <- residuals(res_ANPP, type="pearson")
plot(resanpp)
qqnorm(resanpp)
simulateResiduals(fittedModel = res_ANPP, plot = TRUE)
ols_test_normality(resanpp)
# Using leveneTest() to test for homoscedasticity
result = leveneTest(ANPP ~ interaction(Treatment), data = data2)
print(result)

# Shapiro-Wilk, Kolmogorov-Smirnov, and Anderson-Darling tests all improved after transformation #
res_ANPP <- lmer(sqrt(ANPP) ~ Treatment + (1|Site) + (1|Site:Group), data = data2)
resanpp <- residuals(res_ANPP, type="pearson")
plot(resanpp)
qqnorm(resanpp)
simulateResiduals(fittedModel = res_ANPP, plot = TRUE)
ols_test_normality(resanpp)
# Using leveneTest() to test for homoscedasticity
result = leveneTest(sqrt(ANPP)~ interaction(Treatment), data = data2)
print(result)


# RUN MODEL, GET R2 AND PAIRWISE SIGNIFICANT DIFFERENCES #
a <- lmer(sqrt(ANPP) ~ Treatment + (1|Site) + (1|Site:Group), data = data2)
anova(a)
r.squaredGLMM(a)
emmeans(a, pairwise ~ Treatment, adjust="BH") 


ggplot(data=barGraphStats(data=data2, variable="ANPP",byFactorNames=c("Treatment")),aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', position="dodge", width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Shoot Biomass (g)") +
  xlab("Treatment") +
  annotate("text", x= 1, y = 28, label= "a", size = 15)+ 
  annotate("text", x= 2, y = 26, label= "a", size = 15)+ 
  annotate("text", x= 3, y = 43, label= "b", size = 15)+ 
  scale_fill_manual(values = c("#C9B793", "#6E6C81", "#93AD90")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),axis.text.x=element_text(size = 30), 
        legend.position="none",axis.text.y=element_text(size = 40))
# EXPORT 1400 x 1400 #


##########################
# Nodule Number #

## CHECK FOR NORMALITY, ETC ####
hist(data$NodNum)
res_Nod <- lmer(NodNum ~ Treatment + (1|Site) + (1|Site:Group), data = data2)
resnod <- residuals(res_Nod, type="pearson")
plot(resnod)
qqnorm(resnod)
simulateResiduals(fittedModel = res_Nod, plot = TRUE)
ols_test_normality(resnod)
# Using leveneTest() to test for homoscedasticity
result = leveneTest(NodNum ~ interaction(Treatment), data = anpp2)
print(result)

# Shapiro-Wilk, Kolmogorov-Smirnov, and Anderson-Darling tests all improved after transformation #
res_Nod <- lmer(log1p(NodNum) ~ Treatment + (1|Site) + (1|Site:Group), data = data2)
resnod <- residuals(res_Nod, type="pearson")
plot(resnod)
qqnorm(resnod)
shapiro.test(resnod)
simulateResiduals(fittedModel = res_Nod, plot = TRUE)
ols_test_normality(resnod)
# Using leveneTest() to test for homoscedasticity
result = leveneTest(log1p(NodNum)~ interaction(Treatment), data = data2)
print(result)


# RUN MODEL, GET R2 AND PAIRWISE SIGNIFICANT DIFFERENCES #
b <- lmer(log1p(NodNum) ~ Treatment + (1|Site) + (1|Site:Group), data = data2)
anova(b)
r.squaredGLMM(b)
emmeans(b, pairwise ~ Treatment, adjust="BH") 


ggplot(data=barGraphStats(data=data2, variable="NodNum",byFactorNames=c("Treatment")),aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', position="dodge", width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Nodule Number") +
  xlab("Treatment") +
  annotate("text", x= 1, y = 9, label= "ab", size = 15)+ 
  annotate("text", x= 2, y = 6, label= "a", size = 15)+ 
  annotate("text", x= 3, y = 21, label= "b", size = 15)+ 
  scale_fill_manual(values = c("#C9B793", "#6E6C81", "#93AD90")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),axis.text.x=element_text(size = 30), 
        legend.position="none",axis.text.y=element_text(size = 40))
# 1400 x 1400 #



##########################
# AMF Root Colonization #

## CHECK FOR NORMALITY, ETC ####
hist(data$PercMC)
res_MC <- lmer(PercMC ~ Treatment + (1|Site) + (1|Site:Group), data = data2)
resmc <- residuals(res_MC, type="pearson")
plot(resmc)
qqnorm(resmc)
ols_test_normality(resmc)
simulateResiduals(fittedModel = res_MC, plot = TRUE)
# Using leveneTest() to test for homoscedasticity
result = leveneTest(PercMC ~ interaction(Treatment), data = data2)
print(result)


# RUN MODEL, GET R2 AND PAIRWISE SIGNIFICANT DIFFERENCES #
c <- lmer(PercMC ~ Treatment + (1|Site) + (1|Site:Group), data = data2)
anova(c)
r.squaredGLMM(c)
emmeans(c, pairwise ~ Treatment, adjust="BH")


ggplot(data=barGraphStats(data=data2, variable="PercMC",byFactorNames=c("Treatment")),aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', position="dodge", width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("% AMF Root Colonization") +
  xlab("Treatment") +
  annotate("text", x= 1, y = 56, label= "b", size = 15)+ 
  annotate("text", x= 2, y = 39, label= "a", size = 15)+ 
  annotate("text", x= 3, y = 49, label= "b", size = 15)+ 
  scale_fill_manual(values = c("#C9B793", "#6E6C81", "#93AD90")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),axis.text.x=element_text(size = 30), 
        legend.position="none",axis.text.y=element_text(size = 40))
# 1400 x 1400 #


#################################
# Leaf Nutrient Content #

# Nitrogen #
## CHECK FOR NORMALITY, ETC ####
hist(data2$N)
res_N <- lmer(N ~ Treatment + (1|Site) + (1|Site:Group), data = data2)
resN <- residuals(res_N, type="pearson")
plot(resN)
qqnorm(resN)
simulateResiduals(fittedModel = res_N, plot = TRUE)
ols_test_normality(resN)
# Using leveneTest() to test for homoscedasticity
result = leveneTest(N ~ interaction(Treatment), data = data2)
print(result)


# RUN MODEL, GET R2 AND PAIRWISE SIGNIFICANT DIFFERENCES #
d <- lmer(N ~ Treatment + (1|Site) + (1|Site:Group), data = data2)
anova(d)
r.squaredGLMM(d)
emmeans(d, pairwise ~ Treatment, adjust="BH") 


ggplot(data=barGraphStats(data=data2, variable="N",byFactorNames=c("Treatment")),aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', position="dodge", width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Leaf Nitrogen Content (%)") +
  xlab("Treatment") +
  annotate("text", x= 1, y = 1.9, label= "a", size = 15)+ 
  annotate("text", x= 2, y = 1.8, label= "a", size = 15)+ 
  annotate("text", x= 3, y = 2.2, label= "b", size = 15)+ 
  scale_fill_manual(values = c("#C9B793", "#6E6C81", "#93AD90")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),axis.text.x=element_text(size = 30), 
        legend.position="none",axis.text.y=element_text(size = 40))
# EXPORT 1400 x 1400 #


# Phosphorus #
## CHECK FOR NORMALITY, ETC ####
hist(data2$P)
res_P <- lmer(P ~ Treatment + (1|Site) + (1|Site:Group), data = data2)
resP <- residuals(res_N, type="pearson")
plot(resP)
qqnorm(resP)
simulateResiduals(fittedModel = res_P, plot = TRUE)
ols_test_normality(resP)
# Using leveneTest() to test for homoscedasticity
result = leveneTest(P ~ interaction(Treatment), data = data2)
print(result)

# RUN MODEL, GET R2 AND PAIRWISE SIGNIFICANT DIFFERENCES #
e <- lmer(P ~ Treatment + (1|Site) + (1|Site:Group), data = data2)
anova(e)
r.squaredGLMM(e)
emmeans(e, pairwise ~ Treatment, adjust="BH") 


ggplot(data=barGraphStats(data=data2, variable="P",byFactorNames=c("Treatment")),aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', position="dodge", width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Leaf Phosphorus Content (%)") +
  xlab("Treatment") +
  annotate("text", x= 1, y = 0.12, label= "a", size = 15)+ 
  annotate("text", x= 2, y = 0.115, label= "a", size = 15)+ 
  annotate("text", x= 3, y = 0.18, label= "b", size = 15)+ 
  scale_fill_manual(values = c("#C9B793", "#6E6C81", "#93AD90")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),axis.text.x=element_text(size = 30), 
        legend.position="none",axis.text.y=element_text(size = 40))
# EXPORT 1400 x 1400 #

##########################
# N-Fixation #

## CHECK FOR NORMALITY, ETC ####
hist(data2$Rate)
res_fix <- lmer(Rate ~ Treatment + (1|Site) + (1|Site:Group), data = data2)
resfix <- residuals(res_fix, type="pearson")
plot(resfix)
qqnorm(resfix)
simulateResiduals(fittedModel = res_fix, plot = TRUE)
ols_test_normality(resfix)
# Using leveneTest() to test for homoscedasticity
result = leveneTest(Rate ~ interaction(Treatment), data = data2)
print(result)


e <- lmer(Rate ~ Treatment + (1|Site) + (1|Site:Group), data = data2)
anova(e)
r.squaredGLMM(e)

ggplot(data=barGraphStats(data=data2, variable="Rate",byFactorNames=c("Treatment")),aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', position="dodge", width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Ethylene Produced (nmol per hour)") +
  xlab("Treatment") +
  scale_fill_manual(values = c("#C9B793", "#6E6C81", "#93AD90")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),axis.text.x=element_text(size = 30), 
        legend.position="none",axis.text.y=element_text(size = 40))
# 1400 x 1400 #

plot(e)
