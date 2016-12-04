#### T -Test - The difference between two groups

# Cohen's d is determined by calculating the mean difference between your two groups
# and then dividing the result by the pooled standard deviation. Cohen's d is the appropriate effect size 
# measure if two groups have similar standard deviations and are of similar size

## TRADITIONAL POWER
# assume effect size is .90, want power of .85
# use R to determine approproate sample size for the study

library(pwr)
pwr.t.test(d=.90,power=.85,sig.level=.05,type="two.sample",alternative="greater")

## SAFEGUARD POWER
# calculate d value if one isn't given
library(MBESS)
smd(Mean.1=3.4, Mean.2=2.499, s.1=.95,s.2=1.05, n.1=15,n.2=15)
#d value (or standardized mean difference, smd) = .90

# calculating d value from raw data using CIs
ci.smd(smd=0.90, n.1=15, n.2=15)
# CIs for d value [.1396, 1.659]
# we see d value of .90 could plausibky be caused by pop'n d value as low as .14 or as high as 1.65

# use lower bound to find CI
pwr.t.test(d=0.1396,power=0.85,sig.level=0.05,type="two.sample",alternative="greater")
#safeguard analysis shows we need N of 738 per cell 


# __________________________________________________________________

# CI interval width when you KNOW the population effect size
ss.aipe.smd.full(delta=.90, conf.level=.95,width=.90)
# need N of 42 PER GROUP to ensure CI width doesn't exceed effect size 

# CI interval width when you are basing your effect size on a STUDY RESULT (don't know pop'n effect size)
ci.smd(smd=.90, n.1=15, n.2=15)
# get CI of [0.1395, 1.6459]
# CI interval width when you KNOW the population effect size
ss.aipe.smd.full(delta=.1396, conf.level=.95,width=.1396)
# need N of 1581 PER GROUP   

# __________________________________________________________________

# if original study had better design (i.e., more people per cell, CI would be narrower)
ci.smd(smd=.90, n.1=100, n.2=100)
# CI for larger sample size [.60, 1.18]

#sageguard analysis therefore would be:
ss.aipe.smd.full(delta=.60, conf.level=.95,width=.60)
# N would be 90 per cell needed 

# _____________________________________________________________________

### CONDUCTING THE T-TEST
library(tidyverse)
my.data <- read_csv("drugData.csv")
head(my.data)


# ensure you convert categorical data to factors!!!!!!
my.data$Group <- as.factor(my.data$Group)

glimpse(my.data)

#assign factor levels 
levels(my.data$Group) <- list("Drug"=0, "Control"=1)

# get descriptive statistics
psych::describeBy(my.data$Arousal,group=my.data$Group)

#determine if variances are equal (will use different ttest if equal/not equal)
car::leveneTest(my.data$Arousal, group=my.data$Group,center="median")
#sig = 0.6146 so DO have homogeneity of variance - proceed accordingly 

library(tidyverse)

# select scores for experimental and control group
#Arousal scores for each group
exp.group.rows <- my.data %>% filter(Group == "Drug")
control.group.rows <- my.data %>% filter(Group == "Control")

## EQUAL VARIANCES ASSUMED
t.test(x=exp.group.rows$Arousal,
       y=control.group.rows$Arousal,var.equal = TRUE)

#then calculate effect size with CI
library(MBESS)

d.value <- smd(Group.1=exp.group.rows$Arousal,
               Group.2=control.group.rows$Arousal, Unbiased=TRUE)
print(d.value)
#d value = .62

ci.smd(smd=d.value,
       n.1=length(exp.group.rows$Arousal),
       n.2=length(control.group.rows$Arousal))
# CI [.11, 1.11]

## EQUAL VARIANCES NOT ASSUMED
t.test(x=exp.group.rows$Arousal,
       y=control.group.rows$Arousal,var.equal = FALSE)

d.value <- smd.c(Group.T=exp.group.rows$Arousal,
               Group.C=control.group.rows$Arousal)
print(d.value)
#  d value = 0.6354861

ci.smd.c(smd.c=d.value,n.E=length(exp.group.rows$Arousal),
       n.C=length(control.group.rows$Arousal))
# # CI [.11, 1.14]


# _______________________________________________________________________________
# _______________________________________________________________________________

#### ANOVA

# load data
#library(tidyverse)
oneway.data <- read.csv("Viagra.csv")

# set factors correctly
oneway.data$dose <- as.factor(oneway.data$dose)
levels(oneway.data$dose) <- list("Placebo"=1, "Low Dose"=2, "High Dose"=3)

# display descriptive statistics
psych::describeBy(oneway.data$libido, group=oneway.data$dose)
# will gets Means, SDs
# libido is DV

# test for homogeneity of variance
# levenes test
car::leveneTest(oneway.data$libido, group=oneway.data$dose, center="median")
# not significant so assumption not violated 

## MAIN ANALYSES, EQUAL VARIANCES ASSUMED (i.e., Levene's not significant)
# set contrasts
options(contrasts = c("contr.sum", "contr.poly"))

# run ANOVa
summary(oneway.results)

car::Anova(oneway.results,type=3)

# use APA Tables
library(apaTables)
apa.aov.table(oneway.results,table.number = 1, filename = "Table1.doc")

# APA Mean/SD table
apa.1way.table(iv=dose, dv=libido,
               data=oneway.data,
               table.number = 2,
               filename = "Table2.doc",
               show.conf.interval = TRUE)

# graphing ANOVA 
myBar<-ggplot(oneway.data,aes(dose,libido))
myBar<-myBar + stat_summary(fun.y=mean, geom="bar",fill="Grey",colour="Black")
myBar<-myBar + stat_summary(fun.data=mean_cl_normal, geom="errorbar",width=0.2)
myBar<-myBar + coord_cartesian(ylim=c(1, 7))
myBar<-myBar + scale_y_continuous(breaks=seq(1, 7))
myBar<-myBar + labs(x="Dose", y="Libido")
myBar<-myBar + theme_classic(12)
print(myBar)

## MAIN ANALYSES, EQUAL NOT VARIANCES ASSUMED (i.e., Levene's IS significant)
# see sheet


### POST HOC COMPARISONS 
## when you don't have a specific hypothesis

## Bonferroni correction
pairwise.t.test(oneway.data$libido,oneway.data$dose,p.adjust.method="bonferroni")

## Tukey
library(multcomp)
posthocs<-glht(oneway.results,linfct=mcp(dose="Tukey"))
summary(posthocs)

##use apa tables to to calculate all paired comparisons (w/o controlling for type I error)
library(apaTables)
apa.d.table(iv = dose, dv = libido,
            data =  oneway.data,
            filename="dValueTable.doc"",
            show.conf.interval = TRUE)
              