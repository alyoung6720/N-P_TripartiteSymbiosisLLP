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


BiomassNoduleAMF <- read.csv("Field_ANPP_NodNum.csv") %>%
  # 2-N-15 was dead so exclude it #
  filter(Treatment!="N" | ID!=15)

# CHANGE SITE, REPLICATE, AND TREATMENT TO FACTOR #
BiomassNoduleAMF$Site <- as.factor(BiomassNoduleAMF$Site)
BiomassNoduleAMF$Replicate <- as.factor(BiomassNoduleAMF$Replicate)
BiomassNoduleAMF$Treatment <- as.factor(BiomassNoduleAMF$Treatment)


## CHECK FOR NORMALITY, ETC ####
hist(BiomassNoduleAMF$ANPP)
res_ANPP <- lmer(ANPP ~ Treatment + (1|Site/Replicate), data = BiomassNoduleAMF)
resanpp <- residuals(res_ANPP, type="pearson")
plot(resanpp)
qqnorm(resanpp)
simulateResiduals(fittedModel = res_ANPP, plot = TRUE)
ols_test_normality(resanpp)
# Using leveneTest() to test for homoscedasticity
result = leveneTest(ANPP ~ interaction(Treatment), data = BiomassNoduleAMF)
print(result)

# Shapiro-Wilk, Kolmogorov-Smirnov, and Anderson-Darling tests all improved after transformation #
res_ANPP <- lmer(sqrt(ANPP) ~ Treatment + (1|Site/Replicate), data = BiomassNoduleAMF)
resanpp <- residuals(res_ANPP, type="pearson")
plot(resanpp)
qqnorm(resanpp)
simulateResiduals(fittedModel = res_ANPP, plot = TRUE)
ols_test_normality(resanpp)
# Using leveneTest() to test for homoscedasticity
result = leveneTest(sqrt(ANPP)~ interaction(Treatment), data = BiomassNoduleAMF)
print(result)


# RUN MODEL, GET R2 AND PAIRWISE SIGNIFICANT DIFFERENCES #
a <- lmer(sqrt(ANPP) ~ Treatment + (1|Site/Replicate), data = anpp)
anova(a)
r.squaredGLMM(a)
emmeans(a, pairwise ~ Treatment, adjust="BH") 


ggplot(data=barGraphStats(data=anpp, variable="ANPP",byFactorNames=c("Treatment")),aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', position="dodge", width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Aboveground Biomass (g)") +
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
