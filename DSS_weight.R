setwd("~/PSU/DSSMice")

###Load Packages
library(readxl)
library(dplyr)
library(ggpubr)
library(agricolae)


#Load Data on weights for mice
dss <- read_excel("DSS_data.xlsx")
View(dss)
dss <- na.omit(dss)

#Set Factors
dss$Cage <- as.factor(dss$Cage)
dss$Mouse <- as.factor(dss$Mouse)
dss$Treatment <- as.factor(dss$Treatment)

str(dss)


###Compute ANOVA
#Colon length by Treatment
Colon.aov <- aov(Colon ~ Treatment, data = dss)
#Spleen weight by Treatment
Spleen.aov <- aov(Spleen ~ Treatment, data = dss)
#Final Weight by Treatment
Final.aov <- aov(Final ~ Treatment, data = dss)


# Summary of the analysis#
summary(Colon.aov)
summary(Spleen.aov)
summary(Final.aov)

#Tukey Post Hoc Test
TukeyHSD(Colon.aov)
TukeyHSD(Spleen.aov)
TukeyHSD(Final.aov)

#Group based on significance (assigns letters based on sig)#
HSD.test(Colon.aov, "Treatment", console=TRUE)
HSD.test(Spleen.aov, "Treatment", console=TRUE)
HSD.test(Final.aov, "Treatment", console=TRUE)


#Compute ANOVA 
#percent body weight loss by treament
perc.aov <- aov(perc ~ Treatment, data = dss)

# Summary of the analysis
summary(perc.aov)

#Tukey Post Hoc Test
TukeyHSD(perc.aov)

#Group based on significance (assigns letters based on sig)
HSD.test(perc.aov, "Treatment", console=TRUE)


###ANOVA on Gut Permeability
#Load data for gut permeability
gp <- read_excel("gp_ANOVA.xlsx")
View(gp)

#Set Factors
gp$Sample <- as.factor(gp$Sample)
gp$Treatment <- as.factor(gp$Treatment)

str(gp)

#ANOVA by FITC concentration by Treatment
conc.aov <- aov(conc ~ Treatment, data = gp)

#Summary of the analysis
summary(conc.aov)

#Tukey Post Hoc Test
TukeyHSD(conc.aov)

#Group based on significance (assigns letters based on sig)
HSD.test(conc.aov, "Treatment", console=TRUE)


