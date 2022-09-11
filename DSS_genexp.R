setwd("~/PSU/DSSMice")

###Load Packages
library(readxl)
library(dplyr)
library(ggpubr)
library(agricolae)
library(ggplot2)


#Load Data for gene expression (normalized 2^ddCT)
gene <- read_excel("ddCT_Stats.xlsx")
View(gene)


#1 is negative control (no DSS or cocoa)
#2 is positive control (only 2.5% DSS treatment)
#3 is DSS treatment plus defatted cocoa powder (CP)
#4 is DSS treatment plus extractable polyphenols (ExPP)
#5 is residual leftover after polyphenol extraction (NonExPP)

#Set Treatment as a factor
gene$Treatment <- as.factor(gene$Treatment)

####DOTPLOTS
#tnfa
ggplot(gene, aes(x=Treatment, y=tnfa)) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center')
#il1b
ggplot(gene, aes(x=Treatment, y=il1B)) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center')
#il6
ggplot(gene, aes(x=Treatment, y=il6)) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center')
#il10
ggplot(gene, aes(x=Treatment, y=il10)) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center')

####ANOVA
#Set Factors
gene$Sample <- as.factor(gene$Sample)
gene$Treatment <- as.factor(gene$Treatment)

str(gene)

#TNF-a ANOVA by Treatment
#Compute ANOVA
tnf.aov <- aov(tnfa ~ Treatment, data = gene)

# Summary of the analysis#
summary(tnf.aov)

#Tukey Post Hoc Test
TukeyHSD(tnf.aov)

#Group based on significance (assigns letters based on sig)#
HSD.test(tnf.aov, "Treatment", console=TRUE)


#il6 ANOVA by Treatment
il6.aov <- aov(il6 ~ Treatment, data = gene)
summary(il6.aov)
TukeyHSD(il6.aov)
HSD.test(il6.aov, "Treatment", console=TRUE)

#il1B-a ANOVA by Treatment
il1B.aov <- aov(il1B ~ Treatment, data = gene)
summary(il1B.aov)
TukeyHSD(il1B.aov)
HSD.test(il1B.aov, "Treatment", console=TRUE)

#il10 ANOVA by Treatment
il10.aov <- aov(il10 ~ Treatment, data = gene)
summary(il10.aov)
TukeyHSD(il10.aov)
HSD.test(il10.aov, "Treatment", console=TRUE)



