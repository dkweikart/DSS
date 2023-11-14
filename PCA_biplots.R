setwd("~/PSU/DSSMice")

#Load library
library(readxl)

#Load Data for gene expression (normalized 2^ddCT)
gene <- read_excel("gene_express.xlsx")
View(gene)


#1 is negative control (no DSS or cocoa)
#2 is positive control (only 2.5% DSS treatment)
#3 is DSS treatment plus defatted cocoa powder (CP)
#4 is DSS treatment plus extractable polyphenols (ExPP)
#5 is residual leftover after polyphenol extraction (NonExPP)

#convert into dataframe and set factors
gene$Subject <- as.factor(gene$Subject)
gene$Treatment <- as.factor(gene$Treatment)
gene <- as.data.frame(gene)
str(gene)
rownames(gene) <- gene$Subject
gene

# run PCA using FactoMineR and factoextra packages
# follow http://www.sthda.com/english/articles/31-principal-component-methods-in-rpractical-guide/112-pca-principal-component-analysis-essentials/
library(FactoMineR)
library(factoextra)
# run PCA on variables, using treatment as grouping variable
# run withouth ID and Treatment columns
# run as correlation matrix
data.pca <- PCA(gene[,-c(1:2)], quali.sup = 1:2,
                 scale.unit = TRUE,
                 graph = TRUE)

#Eigenvalues and scree plot
fviz_eig(data.pca, addlabels = TRUE) # keep 4 dimensions for 59% expl. variance

# which variables are sign. correlated to PC 1 through PC4
geneimdesc <- dimdesc(data.pca, axes = 1:4)
56
# PC1
geneimdesc$Dim.1
#Percent weight los, IL1B and spleen weight all sign. pos. corr
#Colon and final weight all sign. neg corr
# PC2
geneimdesc$Dim.2
#Percent weight loss and IL10 all sign. pos. corr
#IL1B and spleen weight all negative


#contribution of each var to PCs for first 5 PCs
library("corrplot")
library(DescTools)
corrplot(data.pca$var$contrib, is.corr=FALSE)
# color individual mice by treatment
# PC 1 and PC 2
PCA_biplot <- fviz_pca_ind(data.pca,
                           axes = c(1,2),
                           geom.ind = c('text', 'point'),
                           labelsize = 3,
                           col.ind = gene$Treatment,
                           palette = 'jco',
                           mean.point = FALSE,
                           addEllipses = TRUE,
                           ellipse.type = 'confidence',
                           legend.title = "Groups")


PCA_biplot


## HH .. played only with lines up to here !!! Aug-18-2021 HH
# and did the ANOVAs (lines 125 ff)
# PC 3 and PC 4
fviz_pca_ind(data.pca,
             axes = c(3,4),
             geom.ind = 'point',
             col.ind = gene$Treatment,
             palette = 'jco',
             addEllipses = TRUE,
             ellipse.type = 'confidence',
             legend.title = "Groups")
# biplot
library(ggpubr)
pc12 <- fviz_pca_biplot(data.pca,
                        axes = c(1,2),
                        geom.ind = 'point',
                        pointshape = 21,
                        fill.ind = gene$Treatment,
                        col.ind = 'black',
                        addEllipses = TRUE,
                        ellipse.type = 'confidence')


pc12



pc34 <- fviz_pca_biplot(data.pca, axes = c(3,4),geom.ind = 'point',
                        pointshape = 21, fill.ind = gene$Treatment, col.ind = 'black',
                        addEllipses = TRUE, ellipse.type = 'confidence', 
                        repel = TRUE, ggtheme = theme_classic2(), title = NULL) +
  ggpubr::fill_palette('jco') + ggpubr::color_palette('lancet') 
ggexport(plotlist = list(pc12, pc34), filename = "data_PCA_biplot.pdf")

##ANOVA for Kianas data
library(readxl)
library(dplyr)
library(ggpubr)
library(agricolae)
library(ggplot2)

ge <- read_excel("Kiana_PCR.xlsx")
View(ge)

#convert into dataframe and set factors
ge$Tretment <- as.factor(ge$Tretment)
ge$Mouse <- as.factor(ge$Mouse)

#Convert to data frame
ge <- as.data.frame(ge)
#Check that all values that will non-numerical values are a factor 
str(ge)
#Shows you the new file, but this is also updated in the tab 
ge
#Compute ANOVA for tnfa by Treatment
tnf.aov <- aov(tnf ~ Tretment, data = ge)
summary(tnf.aov)
TukeyHSD(tnf.aov)
HSD.test(tnf.aov, "Tretment", console=TRUE)
#ns

#Compute ANOVA il1b by Treatment
il1b.aov <- aov(il1b ~ Tretment, data = ge)
summary(il1b.aov)
TukeyHSD(il1b.aov)
HSD.test(il1b.aov, "Tretment", console=TRUE)
#ns

#ANOVA for my data on dCT
library(readxl)
library(dplyr)
library(ggpubr)
library(agricolae)
library(ggplot2)

genet10 <- read_excel("gene_express.xlsx", sheet="no71")
View(genet10)

genet10$Subject <- as.factor(genet10$Subject)
genet10$Treatment <- as.factor(genet10$Treatment)

genet10 <- as.data.frame(genet10)
str(genet10)
genet10

tnf.aov <- aov(tnf ~ Treatment, data = genet10)
summary(tnf.aov)
TukeyHSD(tnf.aov)
HSD.test(tnf.aov, "Treatment", console=TRUE)
#neg and NonExPP a; CP and ExPP ab; ctrl b

il10.aov <- aov(il10 ~ Treatment, data = genet10)
summary(il10.aov)
TukeyHSD(il10.aov)
HSD.test(il10.aov, "Treatment", console=TRUE)
#ExPP a; neg and ctrl ab; CP NonExPP b

il1b6 <- read_excel("gene_express.xlsx", sheet="no43")
View(il1b6)

il1b6$Subject <- as.factor(il1b6$Subject)
il1b6$Treatment <- as.factor(il1b6$Treatment)

il1b6 <- as.data.frame(il1b6)
str(il1b6)
il1b6

il1b.aov <- aov(il1b ~ Treatment, data = il1b6)
summary(il1b.aov)
TukeyHSD(il1b.aov)
HSD.test(il1b.aov, "Treatment", console=TRUE)
#neg a; NonExPP, ExPP, ctrl b; CP c

il6.aov <- aov(il6 ~ Treatment, data = il1b6)
summary(il6.aov)
TukeyHSD(il6.aov)
HSD.test(il6.aov, "Treatment", console=TRUE)
#neg a; ExPP b; ctrl, CP bc; NonExPP c

