# Heavy Metals in American Robins: Impacts & Implications
# Makayla Ohrberg
# 3/14/24
#### Hypotheses ##############

#' Question 1
#' What is the most effective/accurate way of testing birds for lead 
#' when comparing RBCs and Whole Blood for ICP-MS?

#' Question 2
#' How have heavy metal concentrations in American Robins shifted 
#' over time, both short-term and long-term?

#' Question 3
#' How much lead is too much?

#' Question 4
#' Are there any correlations between physical attributes and heavy 
#' metal burden? Is heavy metal burden impacting parasite
#' susceptibility or other physical attributes? Are there any 
#' correlations between heavy metals?

#' Question 5
#' Are there any continental differences in heavy metal concentrations?

#' Question 6
#' Does estrogen correlate with heavy metal concentrations, 
#' particularly lead? Estrogen leads to an increase in medullary
#' bone deposition, so calcium content in brooding females should
#' be higher in bone and lower in blood. In theory lead deposition
#' would follow the same pattern since it often replaces calcium.
#' Thus, free lead in the blood may be lower.Hypotheses ####################
#### Libraries/Files ####################
## load in libraries
library(tidyverse)
library(here)
library(readxl)
library(ggpubr)
library(vtable)
library(chatgpt)
library(car)
library(qqplotr)
library(EnvStats)
library(bestNormalize)
library(PerformanceAnalytics)
library(DescTools)
library(factoextra)
library(ggbiplot)
library(vegan)
library(XICOR)
library(rstatix)


## read in all data
Whole_Blood <- read_excel(here("Research Data/Master Whole Blood Bird Data - Final.xlsx"),
                          col_types = c(rep("text", 6),"numeric", "numeric", rep("text", 5), rep("numeric", 12), "text", 
                                        "text", "text"), na = "<DL"
)
RBCs <- read_excel(here("Research Data/Master RBC Bird Data - Final.xlsx"),
                   col_types = c(rep("text", 5),rep("numeric", 3), "text", rep("numeric", 4), "text", rep("numeric", 12), "text", 
                                 "text", "text"), na = "<DL"
)
Epiphysis <- read_excel(here("Research Data/Master XRF Bone Data.xlsx"),
                        col_types = c(rep ("text", 12), rep ("numeric", 10), "text")
                        , na = "NA",
                        sheet = "Epiphysis"
)                  
Diaphysis <- read_excel(here("Research Data/Master XRF Bone Data.xlsx"),
                        col_types = c(rep ("text", 12), rep ("numeric", 10), "text")
                        , na = "NA",
                        sheet = "Diaphysis"
)  
Keel <- read_excel(here("Research Data/Master XRF Bone Data.xlsx"),
                   col_types = c(rep ("text", 12), rep ("numeric", 10), "text")
                   , na = "NA",
                   sheet = "Keel"
)
All_Bone <- read_excel(here("Research Data/Master XRF Bone Data.xlsx"),
                       col_types = c(rep ("text", 12), rep ("numeric", 10), "text")
                       , na = "NA",
                       sheet = "XRF Individual Specimens"
)
Correlated_Pb_Hg <- read_excel(here("Research Data/Master Combined Pb Blood Data.xlsx"),
                               col_types = c(rep ("text", 5),rep("numeric",3), "text", "numeric", "numeric", rep ("text", 6), "numeric",
                                             "numeric", "numeric", "numeric", "text", "text"))
#### Lead in Bone ##################
###### Diaphysis vs Epiphysis vs Keel ####
#' I'm comparing diaphysis to epiphysis to explore the changes in Pb with
#' trabecular bone content. It's established that there is more trabeculae
#' in the epiphysis or end of a long bone and less in the diaphysis,
#' or mid-shaft of the long bone. The keel has even more trabeculae than
#' than long bones.

All_Bone %>%
  ggplot(aes(y = Pb, x = `Bone Type`)) +
  geom_boxplot() +
  ggtitle("Bone Type vs Pb (Historic & Modern)") +
  xlab("Bone Type") +
  ylab("Pb Conc. (mcg/dL)") +
  stat_summary(fun="mean", color = "red", shape = 13)

All_Bone %>%
  filter(Set == "Modern") %>%
  ggplot(aes(y = Pb, x = `Bone Type`)) +
  geom_boxplot() +
  ggtitle("Bone Type vs Pb (Modern)") +
  xlab("Bone Type") +
  ylab("Pb Conc. (mcg/dL)") +
  stat_summary(fun="mean", color = "red", shape = 13)

All_Bone %>%
  filter(Set == "Historic") %>%
  ggplot(aes(y = Pb, x = `Bone Type`)) +
  geom_boxplot() +
  ggtitle("Bone Type vs Pb (Historic)") +
  xlab("Bone Type") +
  ylab("Pb Conc. (mcg/dL)") +
  stat_summary(fun="mean", color = "red", shape = 13)

All_Bone_Modern <- All_Bone %>%
  filter(Set == "Modern")
wilcox.test(All_Bone_Modern$Pb ~ All_Bone_Modern$`Bone Type`) # sig p-value = 0.03672

group_by(All_Bone_Modern, `Bone Type`) %>%
  summarise(
    count = n(),
    median = median(Pb, na.rm = TRUE),
    IQR = IQR(Pb, na.rm = TRUE)
  )
# DIA = 13.0 EPI = 17.3

All_Bone_Historic <- All_Bone %>%
  filter(Set == "Historic")
kruskal.test(All_Bone_Historic$Pb ~ All_Bone_Historic$`Bone Type`) # sig p-value = 0.01855

pairwise.wilcox.test(All_Bone_Historic$Pb, All_Bone_Historic$`Bone Type`,
                     p.adjust.method = "BH")
# sig difference between 
# EPI and DIA = 0.026
# KEEL and DIA = 0.026

group_by(All_Bone_Historic, `Bone Type`) %>%
  summarise(
    count = n(),
    median = median(Pb, na.rm = TRUE),
    IQR = IQR(Pb, na.rm = TRUE)
  )
# DIA = 10.5 EPI = 28.2 KEEL = 34.9

###### Long-Term Changes in Bone ######
######## Diaphysis ########

hist(Diaphysis$Pb)
shapiro.test(Diaphysis$Pb)
qPb <- quantile(Diaphysis$Pb, probs = c(.25, .75), na.rm = TRUE)
IQR_Pb <- IQR(Diaphysis$Pb)
Lower_Pb <- qPb[1] - 1.5*IQR_Pb
Upper_Pb <- qPb[2] + 1.5*IQR_Pb
Diaphysis_No_Outliers <- subset(Diaphysis, Diaphysis$Pb > Lower_Pb & Diaphysis$Pb < Upper_Pb)
hist(Diaphysis_No_Outliers$Pb)
shapiro.test(Diaphysis_No_Outliers$Pb) # non-normal p-value = 0.001684 needs transformed
hist(1/(Diaphysis_No_Outliers$Pb)^2)
shapiro.test(1/(Diaphysis_No_Outliers$Pb)^2) # this is the closest I can get to normal p-value = 0.0188

wilcox.test(data = Diaphysis_No_Outliers, Pb ~ Set)
# accept null hypothesis bc non-sig p-value = 0.6888
# Diaphysis has no sig difference in Pb between historic and modern.

######## Epiphysis ########
ggplot(Epiphysis, aes(x = Year, y = Pb)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, geom = "point", shape = 21, size = 2, fill = 'red') +
  ggtitle("Bone vs Pb") +
  xlab("Year") +
  ylab("Pb (ppm)")
mean(Diaphysis$Pb) # 16.60244
sd(Diaphysis$Pb) # 14.38707
geoMean(Epiphysis$Pb) # 22.238
geoSD(Epiphysis$Pb) # 2.429632
geoMean(Keel$Pb) # 38.57986
geoSD(Keel$Pb) # 2.526579
#' Epiphysis and diaphysis show a decrease in Pb over time based on
#' range. The mean of 2023 may be artificially high due to an 
#' outlier. This is likely due to the phasing out of leaded gasoline.
#' The keel has the highest geometric mean followed by the epiphysis.

# Let me remove the 2023 outlier and regraph it.
Epiphysis_No_Pb_Outlier <- Epiphysis %>%
  filter(Pb < 55)
ggplot(Epiphysis_No_Pb_Outlier, aes(x = Year, y = Pb)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, geom = "point", shape = 21, size = 2, fill = 'red') +
  ggtitle("Bone vs Pb") +
  xlab("Year") +
  ylab("Pb (ppm)")
#' So this shows more of a decrease in 2023 without the outlier but
#' 1997 throws it off a bit because I only have one bird for that year.

hist(Epiphysis_No_Pb_Outlier$Pb)
shapiro.test(Epiphysis_No_Pb_Outlier$Pb) 
# normal p-value = 0.2895 no need to transform

leveneTest(data = Epiphysis_No_Pb_Outlier, Pb ~ Set)
# no sig p-value so equal variances

t.test(data = Epiphysis_No_Pb_Outlier, Pb ~ Set)
# non-sig p-value = 0.485 so no sig difference in mean Pb between
# modern and historic.

p1 <- ggplot(Epiphysis_No_Pb_Outlier, aes(x = Set, y = Pb)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, geom = "point", shape = 21, size = 2, fill = 'red') +
  ggtitle("Pb in Modern vs Historic Birds") +
  xlab("Set") +
  ylab("Pb (ppm)")

###################### Lead in Blood
###### Diagnosing Lead Poisoning/Toxicity ####

# Measures of Central Tendency
Diaphysis_No_Pb_Outlier <- Diaphysis %>%
  filter(Pb < 50)

Diaphysis_Modern <- Diaphysis %>%
  filter(Set == 'Modern')

Diaphysis_Historic <- Diaphysis %>%
  filter(Set == 'Historic')

mean(Diaphysis$Pb) # 16.60244 ppm
median(Diaphysis$Pb) # 12.91964 ppm
var(Diaphysis$Pb) # 206.9877

mean(Diaphysis_No_Pb_Outlier$Pb) # 13.655 ppm
median(Diaphysis_No_Pb_Outlier$Pb) # 12.86044 ppm
var(Diaphysis_No_Pb_Outlier$Pb) # 16.65865

mean(Diaphysis_Modern$Pb) # 17.94629 ppm
median(Diaphysis_Modern$Pb) # 12.97884 ppm
var(Diaphysis_Modern$Pb) # 294.2536

Diaphysis_Modern_No_Outlier <- Diaphysis_Modern %>%
  filter(Pb < 50)

mean(Diaphysis_Modern_No_Outlier$Pb) # 13.62113 ppm
median(Diaphysis_Modern_No_Outlier$Pb) # 12.91964 ppm
var(Diaphysis_Modern_No_Outlier$Pb) # 14.69775

mean(Diaphysis_Historic$Pb) # 13.723 ppm
median(Diaphysis_Historic$Pb) # 10.538 ppm
var(Diaphysis_Historic$Pb) # 23.676

Epiphysis_Modern <- Epiphysis %>%
  dplyr::filter(Set == 'Modern')

Epiphysis_Historic <- Epiphysis %>%
  dplyr::filter(Set == 'Historic')

mean(Epiphysis_Modern$Pb) # 29.2132 ppm (n = 15)
median(Epiphysis_Modern$Pb) # 17.27194 ppm
var(Epiphysis_Modern$Pb) # 1470.282

mean(Epiphysis_Historic$Pb) # 38.930 ppm (n = 7)
median(Epiphysis_Historic$Pb) # 28.229 ppm
var(Epiphysis_Historic$Pb) # 919.46

mean(Keel$Pb) # 56.231 ppm
median(Keel$Pb) # 34.873 ppm
var(Keel$Pb) # 3366.636

### Now I want to find percentiles. 

## Diaphysis
quantile(Diaphysis$Pb, probs = 0.9, na.rm = TRUE) # 22.223 full data set

quantile(Diaphysis_Modern$Pb, probs = 0.9, na.rm = TRUE) # 20.956 modern data set with outlier

quantile(Diaphysis_Modern_No_Outlier$Pb, probs = .9, na.rm = TRUE) # 17.591 modern data set without outlier

quantile(Diaphysis_Historic$Pb, probs = 0.9, na.rm = TRUE) # 18.917 historic data set

## Epiphysis
quantile(Epiphysis$Pb, probs = 0.9, na.rm = TRUE) # 54.937 full data set

quantile(Epiphysis_Modern$Pb, probs = 0.9, na.rm = TRUE) # 37.753 modern data set

quantile(Epiphysis_Historic$Pb, probs = 0.9, na.rm = TRUE) # 73.73667 historic data set

## Keel
quantile(Keel$Pb, probs = 0.9, na.rm = TRUE) # 117.314 full data set (only historic)

### Now I need to find what percentage of birds meet the 20 ppm guideline

## Diaphysis
Diaphysis2 <- Diaphysis %>%
  filter(Pb > 20)
3/22 * 100 # 3 birds out of 22 = 13.636%

Diaphysis3<- Diaphysis_Modern %>%
  filter(Pb >20)
2/15 * 100 # 2 birds out of 15 = 13.33%

Diaphysis4<- Diaphysis_Historic %>%
  filter(Pb >20)
1/7 * 100 # 1 out of 7 birds = 14.29%

## Epiphysis
Epiphysis2 <- Epiphysis %>%
  dplyr::filter(Pb > 20)
11/22 * 100 # 11 birds out of 22 = 50%

Epiphysis3 <- Epiphysis_Modern %>%
  filter(Pb > 20)
6/15 * 100 # 6 birds out of 15 = 40%

Epiphysis4 <- Epiphysis_Historic %>%
  filter(Pb > 20)
5/7 * 100 # 5 birds out of 7 = 71.429%

## Keel
Keel2 <- Keel %>%
  filter(Pb > 20)
5/7 * 100 # 5 birds out of 7 = 71.429%

### Now for subclinical toxicity (10-19.9)

## Diaphysis
Diaphysis2 <- Diaphysis %>%
  filter(between(Pb, 10.0, 19.9)) 
18/22 * 100 # 18 birds out of 22 = 81.818%

Diaphysis3<- Diaphysis_Modern %>%
  filter(between(Pb, 10.0, 19.9))
13/15 * 100 # 13 birds out of 15 = 86.667%

Diaphysis4<- Diaphysis_Historic %>%
  filter(between(Pb, 10.0, 19.9))
5/7 * 100 # 5 birds out of 7 = 71.429

## Epiphysis
Epiphysis2 <- Epiphysis %>%
  filter(between(Pb, 10.0, 19.9))
10/22 * 100 # 10 birds out of 22 = 45.455%

Epiphysis3 <- Epiphysis_Modern %>%
  filter(between(Pb, 10.0, 19.9))
8/15 * 100 # 8 birds out of 15 = 53.333%

Epiphysis4 <- Epiphysis_Historic %>%
  filter(between(Pb, 10.0, 19.9))
2/7 * 100 # 2 birds out of 7 = 28.571%

## Keel

Keel2 <- Keel %>%
  filter(between(Pb, 10.0, 19.9))
1/7* 100 # 1 bird out of 7 = 14.286%

#### Lead in Blood #################
###### Correlating Pb between Whole Blood and RBCs ####

# First I'm filtering out non-correlated data points.
Corr_RBCs <- RBCs %>%
  filter(Correlated == "Y")
Corr_Whole <- Whole_Blood %>%
  filter(Correlated == "Y")
# Now I'm plotting correlated data points against each other.
p2 <- ggplot(mapping = aes(x = Corr_RBCs$Pb, y = Corr_Whole$Pb)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_regline_equation(label.x = 5.6, label.y = 7.5) +
  stat_cor(label.x = 5.2, label.y = 6.8) +
  ggtitle("Correlated RBCs vs Whole Blood for Pb") +
  xlab("RBCs (ppm)") +
  ylab("Whole Blood (ppm)")
#' Based on the graph, the regression explains 85% of the variation 
#' between our variables. Using this formula of 'y = 1.315 + 0.396x'
#' we can convert our RBCs to be equivalent to whole blood. To do this,
#' we can create a new dataframe by applying the formula to the RBCs.

#' First I created a function that will be applied to the RBCs$Pb
#' column. Then I created a clone of the RBCs dataset to apply the
#' function to. I used sapply to apply the function then viewed the
#' dataset to ensure it worked.
Pb_Correlation <- function(Pb){
  0.396*Pb + 1.315
}
Correlated_All <- RBCs
Correlated_All$Pb <- sapply(Correlated_All$Pb, Pb_Correlation)

#' Next I'll create a new table by binding portions of the Whole_Blood
#' table with the Correlated_Pb table. First, I'll have to remove
#' the extra elements from Whole_Blood and Correlated_Pb. Then I can
#' bind Whole_Pb and Correlated_Pb.
Whole_Pb <- Whole_Blood[c(1:13, 22, 24, 25)]
Correlated_Pb_Manual <- Correlated_All[c(1:13, 22, 24, 25)]
All_Blood <- rbind(Whole_Pb, Correlated_Pb)
#' Finally we have a dataset showing whole blood and RBCs that have
#' been converted to be comparable to whole blood.This gives us one
#' big dataset that is good for comparing things like age, sex, and
#' capture site to Pb concentration. 

###### Basic Statistics ######

mean(Correlated_Pb_Hg$Pb)
# 4.117
min(Correlated_Pb_Hg$Pb)
# 0.8497
max(Correlated_Pb_Hg$Pb)
# 30.873

Correlated_Pb_Hg %>% 
  group_by(State) %>%
  summarise(
    count = n(),
    mean = mean(Pb, na.rm = TRUE),
    median = median(Pb, na.rm = TRUE),
    IQR = IQR(Pb, na.rm = TRUE)
  )

sumtable(Correlated_Pb_Hg, 
         vars = c('T0', 'T2', 'E2', 'Wing', 'Mass', 'Hg', 'Pb'), 
         digits = 2, 
         summ = c("mean", "sd", "min", "max", "n"))

sumtable(Whole_Blood,
         vars = c('Cu', 'Zn', 'As'),
         digits = 2)
####
###### Seasonality in Blood ######
#' Now I want to think about seasonality in my blood data.To do this,
#' I'm going to work with my Correlated_Pb_Hg file and my Whole_Blood 
#' file for the other elements.

View(Correlated_Pb_Hg)

ggplot(Correlated_Pb_Hg, mapping = aes(x = Month, y = Pb)) +
  geom_boxplot()
# I want to reorder my months to be in year sequence.
Correlated_Pb_Hg$Month <- factor(Correlated_Pb_Hg$Month, levels = c("1", "2", "3", "4", "5", "6", "9", "10", "11", "12"))
ggplot(Correlated_Pb_Hg, mapping = aes(x = Month, y = Pb)) +
  geom_boxplot() +
  ggtitle("Month by Month Pb Concentrations") +
  xlab("Month") +
  ylab("Pb Conc. (mcg/dL)") +
  stat_summary(fun="mean", color = "red", shape = 13)

# I think I want to remove outliers just to clear up my month to month variation.

Q <- quantile(Correlated_Pb_Hg$Pb, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Correlated_Pb_Hg$Pb)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
No_Outliers_All_Pb <- subset(Correlated_Pb_Hg, Correlated_Pb_Hg$Pb > Lower & Correlated_Pb_Hg$Pb < Upper)

p3 <- ggplot(No_Outliers_All_Pb, mapping = aes(x = Month, y = Pb)) +
  geom_boxplot() +
  ggtitle("Month by Month Pb Concentrations") +
  xlab("Month") +
  ylab("Pb Conc. (mcg/dL)") +
  stat_summary(fun="mean", color = "red", shape = 13)

#' Averages seem to peak in the summer months, which would make sense
#' as the birds diet shifts to earthworms in the summer months.

# Now I need to run some kind of statistical test to verify the 
# difference I'm seeing between winter and summer.

hist(No_Outliers_All_Pb$Pb)
shapiro.test(No_Outliers_All_Pb$Pb) 
# non-normal p-value = 6.102e-13 so needs transformed

hist(log(No_Outliers_All_Pb$Pb))
shapiro.test(log(No_Outliers_All_Pb$Pb))

# ok none of the basic transformations are working so I'm going to use
# a wilcox test

wilcox.test(data = No_Outliers_All_Pb, Pb ~ Season)
# sig p-value = 0.006899

group_by(No_Outliers_All_Pb, Season) %>%
  summarise(
    count = n(),
    median = median(Pb, na.rm = TRUE),
    IQR = IQR(Pb, na.rm = TRUE)
  )
# median summer = 3.12 median winter = 2.59

###### Breeding Season vs Pb ######

Correlated_Pb_Hg_2023 <- Correlated_Pb_Hg %>%
  filter(Year == '2023')

Correlated_Pb_Hg_2023 %>%
  filter(Sex  == 'F') %>%
ggplot(mapping = aes(x = Breeding, y = Pb)) +
  geom_boxplot() +
  ggtitle("Breeding Season vs Pb (2023)") +
  xlab("Breeding Season") +
  ylab("Pb Conc. (mcg/dL)") +
  stat_summary(fun="mean", color = "red", shape = 13)
#' Female birds in their breeding season seem to have a lower blood lead level
#' but there are a few outliers. Let me remove outliers then retry this and 
#' look for statistical significance.

Q <- quantile(Correlated_Pb_Hg_2023$Pb, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Correlated_Pb_Hg_2023$Pb)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
No_Outliers <- subset(Correlated_Pb_Hg_2023, Correlated_Pb_Hg_2023$Pb > Lower & Correlated_Pb_Hg_2023$Pb < Upper)

No_Outliers %>%
  filter(Sex  == 'F') %>%
  ggplot(mapping = aes(x = Breeding, y = Pb)) +
  geom_boxplot() +
  ggtitle("Breeding Season vs Pb (2023) No Outliers") +
  xlab("Breeding Season") +
  ylab("Pb Conc. (mcg/dL)") +
  stat_summary(fun="mean", color = "red", shape = 13)
#' With outliers removed the difference between non-breeding and breeding is 
#' more obvious. With this, I can test for statistical significance.

No_Outliers_F <- No_Outliers %>%
  filter(Sex == 'F')

hist(No_Outliers_F$Pb)
shapiro.test(No_Outliers_F$Pb) # non-normal
hist(log(No_Outliers_F$Pb))
shapiro.test(log(No_Outliers_F$Pb)) # normal p-value = 0.0596

leveneTest(No_Outliers_F$Pb ~ No_Outliers_F$Breeding) # non-sig so equal variances

my.aov_Breeding_Pb <- aov(No_Outliers_F$Pb ~ No_Outliers_F$Breeding)
summary(my.aov_Breeding_Pb) # slightly significant difference p-value = 0.0352

group_by(No_Outliers_F, Breeding) %>%
  summarise(
    count = n(),
    mean = mean(Pb, na.rm = TRUE),
    IQR = IQR(Pb, na.rm = TRUE)
  )
# Non-breeding = 1.68, breeding = 1.63
# supports my theory that estradiol is resulting in translocation of Pb

#' Just occurred to me that population (AK vs IN) may be having an effect on these
#' findings so I'll repeat this on just Indiana. 

Correlated_Pb_Hg_2023 %>%
  filter(State = 'IN') %>%
  filter(Sex  == 'F') %>%
  ggplot(mapping = aes(x = Breeding, y = Pb)) +
  geom_boxplot() +
  ggtitle("Breeding Season vs Pb (2023)") +
  xlab("Breeding Season") +
  ylab("Pb Conc. (mcg/dL)") +
  stat_summary(fun="mean", color = "red", shape = 13)
# appears to be significant difference in mean Pb between non and breeding. Some
# outliers that need stripped out.

Correlated_Pb_Hg_2023_IN <- Correlated_Pb_Hg_2023 %>%
  filter(State == 'IN') %>%
  filter(sex == 'F')

Q <- quantile(Correlated_Pb_Hg_2023_IN$Pb, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Correlated_Pb_Hg_2023_IN$Pb)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
No_Outliers <- subset(Correlated_Pb_Hg_2023_IN, Correlated_Pb_Hg_2023_IN$Pb > Lower & Correlated_Pb_Hg_2023_IN$Pb < Upper)

hist(No_Outliers$Pb)
shapiro.test(No_Outliers$Pb) # non-normal
hist((No_Outliers$Pb)^(1/3))
shapiro.test((No_Outliers$Pb)^(1/3))

###### Location Impacts on Blood Conc. ######
#' First I'm going to look at my Pb (No Outliers) to check for normality,
#' run a wilcox.test because it is non-normal and hard to transform then,
#' I'll group by state (AK vs IN) to calculate the median for each one.

hist((No_Outliers_All_Pb$Pb))
shapiro.test((No_Outliers_All_Pb$Pb))

wilcox.test(data = No_Outliers_All_Pb, Pb ~ State)
# sig p-value = 0.0007707

group_by(No_Outliers_All_Pb, State) %>%
  summarise(
    count = n(),
    median = median(Pb, na.rm = TRUE),
    IQR = IQR(Pb, na.rm = TRUE)
  )
# AK median = 2.28 IN median = 2.99

#' Now I want to look at just Indiana Pb median because the lower Alaska
#'  values may be skewing them lower.

Indiana <- No_Outliers_All_Pb %>%
  filter(State == 'IN')

ggplot(Indiana, mapping = aes(x = Month, y = Pb)) +
  geom_boxplot() +
  ggtitle("Month by Month Pb Concentrations No Outliers - Indiana") +
  xlab("Month") +
  ylab("Pb Conc. (mcg/dL)") +
  stat_summary(fun="mean", color = "red", shape = 13)

hist(Indiana$Pb)
shapiro.test(Indiana$Pb)
# non-normal p-value = 1.956e-10
hist(1/sqrt(Indiana$Pb))
shapiro.test(1/sqrt(Indiana$Pb))
# can't be forced normal so using wilcox test
wilcox.test(data = Indiana, Pb ~ Season)
# sig diff p-value = 4.981e-06
group_by(Indiana, Season) %>%
  summarise(
    count = n(),
    median = median(Pb, na.rm = TRUE),
    IQR = IQR(Pb, na.rm = TRUE)
  )
# Summer = 3.59 Winter = 2.59

###### Sex vs Pb ######

ggboxplot(Correlated_Pb_Hg, x = "Sex", y = "Pb", 
          color = "Sex", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("M", "F", "U"),
          ylab = "Pb Conc. (mcg/dL)", xlab = "Sex")
# shows females have a slighly higher median and bigger range but males
# seem to have more outliers. Unknown range is very large but has less
# outliers.

hist(Correlated_Pb_Hg$Pb) # shows outliers
Q <- quantile(Correlated_Pb_Hg$Pb, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Correlated_Pb_Hg$Pb, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Correlated_Pb_No_Outliers <- subset(Correlated_Pb_Hg, Correlated_Pb_Hg$Pb > Lower & Correlated_Pb_Hg$Pb < Upper)
hist(Correlated_Pb_No_Outliers$Pb)
shapiro.test(Correlated_Pb_No_Outliers$Pb) # non-normal p-value must transform
hist(log(Correlated_Pb_No_Outliers$Pb))
shapiro.test(1/(Correlated_Pb_No_Outliers$Pb)^2)
# cant force normal so will use non-parametric tests

kruskal.test(Pb ~ Sex, data = Correlated_Pb_No_Outliers)
# non-sig p-value - 0.2377 meaning there is no significant variation
# in distribution of Pb concentrations according to sex.

###### Age vs Pb ######

ggboxplot(Correlated_Pb_Hg, x = "Age", y = "Pb", 
          color = "Age", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("HY", "AHY", "NA"),
          ylab = "Pb Conc. (mcg/dL)", xlab = "Age")
# Shows HY has a slightly higher median than AHY

hist(Correlated_Pb_Hg$Pb) # potential outliers
Q <- quantile(Correlated_Pb_Hg$Pb, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Correlated_Pb_Hg$Pb, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Correlated_Pb_No_Outliers <- subset(Correlated_Pb_Hg, Correlated_Pb_Hg$Pb > Lower & Correlated_Pb_Hg$Pb < Upper)
hist(Correlated_Pb_No_Outliers$Pb)
shapiro.test(Correlated_Pb_No_Outliers$Pb) 
# non-normal p-value = 6.102e-13 should try to transform
hist(log10(Correlated_Pb_No_Outliers$Pb))
shapiro.test(log10(Correlated_Pb_No_Outliers$Pb)) 
# non-normal p-value = 0.0008555 let me try a more powerful trans.

bn_Age_Pb <- bestNormalize(Correlated_Pb_No_Outliers$Pb)
orderNorm_Age_Pb <- orderNorm(Correlated_Pb_No_Outliers$Pb)
hist(orderNorm_Age_Pb$x.t)
shapiro.test(orderNorm_Age_Pb$x.t) # normal p-value = 1
leveneTest(orderNorm_Age_Pb$x.t ~ Correlated_Pb_No_Outliers$Age)
# non-sig p-value = 0.0814 means equal variance
my.aov_Age_Pb <- aov(orderNorm_Age_Pb$x.t ~ Correlated_Pb_No_Outliers$Age)
summary(my.aov_Age_Pb)
# non-sig p-value = 0.51 means no sig difference in Pb means between
# age groups.

# non-transformed data non-parametric
# Kruskal-Wallis
kruskal.test(Pb ~ Age, data = Correlated_Pb_No_Outliers)
# non-sig p-value = 0.4952 means there is no significant difference in
# distribution between age groups for Pb.

###### Fat vs Pb ######

hist(Correlated_Pb_Hg$Pb) # potential outliers
Q <- quantile(Correlated_Pb_Hg$Pb, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Correlated_Pb_Hg$Pb, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Correlated_Pb_No_Outliers <- subset(Correlated_Pb_Hg, Correlated_Pb_Hg$Pb > Lower & Correlated_Pb_Hg$Pb < Upper)
hist(Correlated_Pb_No_Outliers$Pb)
shapiro.test(Correlated_Pb_No_Outliers$Pb) # non-normal p-value
# cant force normal so must use non-parametric tests

kruskal.test(data = Correlated_Pb_No_Outliers, Pb ~ Fat)
# non-sig p-value = 0.37

kruskal.test(data = Correlated_Pb_Hg, Pb ~ condition_resid)
# non-sig p-value = 0.1817

###### Estrogen vs Pb ######
hist(RBCs$E2)
shapiro.test(RBCs$E2) # non-normal p-value = 0.007244

RBCs %>%
  dplyr::filter(Year == 2023) %>%
  ggplot(mapping = aes(x = E2, y = Pb)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("E2 vs Pb in RBCs") +
  xlab("E2") +
  ylab("Pb") +
  stat_regline_equation(label.x = 75, label.y = 40) +
  stat_cor(label.x = 70, label.y = 37)
# linear regression does not require normalization so this graph
# is valid.
# y = 11 - 0.032x R = -0.32 p = 0.045

######## 2023 ########

RBCs_2023 <- RBCs %>%
  filter(Year == 2023)

hist(RBCs_2023$E2)
shapiro.test(RBCs_2023$E2)

RBCs_2023 %>%
  ggplot(mapping = aes(x = E2, y = Pb)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("E2 vs Pb in RBCs") +
  xlab("E2") +
  ylab("Pb") +
  stat_regline_equation(label.x = 75, label.y = 40) +
  stat_cor(label.x = 70, label.y = 37)
# y = 11-0.032x R = -0.32 p = 0.045

hist((RBCs$E2)^(1/3))
shapiro.test((RBCs$E2)^(1/3)) # normal p-value = 0.05151

hist((RBCs_2023$E2)^(1/3))
shapiro.test((RBCs_2023$E2)^(1/3))
# can't force normal so need to try non-parametric stuff

cor.test(RBCs_2023$E2, RBCs_2023$Pb, method=("spearman"))
# rho = -0.32 p-value = 0.04037

RBCs_2023 %>%
  ggplot(mapping = aes(x = E2, y = Pb)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("E2 vs Pb in RBCs (2023)") +
  xlab("E2") +
  ylab("Pb Conc. (mcg/dL)") +
  stat_regline_equation(label.x = 75, label.y = 40) +
  stat_cor(method = "spearman", label.x = 65, label.y = 37, digits = 3)


RBCs_Numeric_E2 <- RBCs %>%
  dplyr::select(E2, Cu, Zn, As, Hg, Pb)

chart.Correlation(RBCs_Numeric_E2)
# showing a slight correlation between E2 and Pb.I should probably
# use my correlated_Pb_Hg worksheet for this.

######## Correlated to Whole Blood ########

Correlated_Pb_Hg_2023 <- Correlated_Pb_Hg %>%
  filter(Year == 2023)

Correlated_Pb_Hg_2023 %>%
  ggplot(mapping = aes(x = E2, y = Pb)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("E2 vs Pb in RBCs") +
  xlab("E2") +
  ylab("Pb") +
  stat_regline_equation(label.x = 75, label.y = 17) +
  stat_cor(label.x = 70, label.y = 15)
# y = 5.5 - 0.013x R = -0.32 p = 0.045

Correlated_Pb_Hg_2023 %>%
  dplyr::select(E2, Hg, Pb) %>%
  chart.Correlation()

p1 <- Correlated_Pb_Hg_2023 %>%
  ggplot(mapping = aes(x = E2, y = Pb)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("E2 vs Pb, Correlated Blood (2023)") +
  xlab("E2 (pg/mL)") +
  ylab("Pb Conc. (mcg/dL)") +
  ylim(0,15) +
  stat_regline_equation(label.x = 75, label.y = 14.5) +
  stat_cor(method = "spearman", label.x = 65, label.y = 13.5, cor.coef.name = "rho", digits = 3)
# y=5.5-0.013x rho = -0.322 p = 0.0404

# need to filter by state and check for difference in variance between states

Correlated_Pb_Hg_2023 %>%
  filter(State == "IN") %>%
  ggplot(mapping = aes(x = E2, y = Pb)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("E2 vs Pb, Correlated Whole Blood (2023, IN)") +
  xlab("E2 (pg/mL)") +
  ylab("Pb (mcg/dL)") +
  stat_regline_equation(label.x = 70, label.y = 15.2) +
  stat_cor(method = "spearman", label.x = 65, label.y = 14, cor.coef.name = "rho", digits = 3)
# y = 5.9 - 0.018x R = -0.184 p = 0.366
# non-sig correlation

Correlated_Pb_Hg_2023 %>%
  filter(State == "AK") %>%
  ggplot(mapping = aes(x = E2, y = Pb)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("E2 vs Pb, Correlated Whole Blood (2023, AK)") +
  xlab("E2 (pg/mL)") +
  ylab("Pb (mcg/dL)") +
  stat_regline_equation(label.x = 75, label.y = 9.5) +
  stat_cor(method = "spearman", label.x = 65, label.y = 9, cor.coef.name = "rho", digits = 3)
# y = 5.2 - 0.01x R = -0.404 p = 0.137 
# non-sig correlation

Correlated_Pb_Hg_2023 %>%
  dplyr::filter(State == "IN") %>%
  dplyr::select(E2, Hg, Pb) %>%
  chart.Correlation(method = "spearman")
# shows non-sig correlation between E2 and Pb for IN

Correlated_Pb_Hg_2023 %>%
  dplyr::filter(State == "AK") %>%
  dplyr::select(E2, Hg, Pb) %>%
  chart.Correlation(method = "spearman")
# no correlation between E2 and Pb for AK

hist(Correlated_Pb_Hg_2023$E2)
shapiro.test(Correlated_Pb_Hg_2023$E2) # non-normal p-value = 0.005457 need to normalize
hist(log(Correlated_Pb_Hg_2023$E2))
shapiro.test(log(Correlated_Pb_Hg_2023$E2)) # non-normal p-value = 0.03884
hist(1/sqrt(Correlated_Pb_Hg_2023$E2))
shapiro.test(1/sqrt(Correlated_Pb_Hg_2023$E2)) # non-normal p-value = 0.007854
hist((Correlated_Pb_Hg_2023$E2)^(1/3))
shapiro.test((Correlated_Pb_Hg_2023$E2)^(1/3)) # non-normal p-value = 0.04393
hist((Correlated_Pb_Hg_2023$E2)^(1/2))
shapiro.test((Correlated_Pb_Hg_2023$E2)^(1/2)) # non-normal p-value = 0.03424
hist((Correlated_Pb_Hg_2023$E2)^(1/4))
hist((Correlated_Pb_Hg_2023$E2) + 1)
shapiro.test((Correlated_Pb_Hg_2023$E2) + 1) # non-normal p-value = 0.005457
# can't force normal so have to use non-parametric tests

ks.test(Correlated_Pb_Hg_2023$E2, "pnorm")
# non-normal data

leveneTest(Correlated_Pb_Hg_2023$E2 ~ Correlated_Pb_Hg_2023$State)
# non-sig p-value = 0.1605 so equal variances

kruskal.test(data = Correlated_Pb_Hg_2023, E2 ~ State)
# sig p-value = 0.00079 so sig difference in E2 between states

group_by(Correlated_Pb_Hg_2023, State) %>%
  summarise(
    count = n(),
    median = median(E2, na.rm = TRUE),
    IQR = IQR(E2, na.rm = TRUE)
  )
# median AK = 158, IN = 66.4

######## Estrogen filtered by breeding season vs non-breeding season and state ########

Correlated_Pb_Hg_2023 <- Correlated_Pb_Hg %>%
  filter(Year == '2023')

boxplot(Correlated_Pb_Hg_2023$E2 ~ Correlated_Pb_Hg_2023$State)
# obvious difference in median E2 between states, with AK being much higher

# now to look for statistical significance
# E2 is non-normal and can't successfully be transformed so must use non-parametric tests

stat.test <- Correlated_Pb_Hg_2023 %>%
  wilcox_test(E2 ~ State) %>%
  add_significance()
stat.test
# highly significant difference in E2 between states p-value = 0.000504

stat.test <- Correlated_Pb_Hg_2023 %>%
  wilcox_test(E2 ~ Breeding) %>%
  add_significance()
stat.test

# only have one non-breeding bird in 2023 so this doesn't hold much weight.

Correlated_Pb_Hg_2023 %>%
  group_by(State) %>%
  get_summary_stats(E2, type = "median_iqr")
# medians AK = 158, IN = 66.4

Correlated_Pb_Hg_2023 %>%
  filter(State == "IN") %>%
  filter(Month %in% c('4', '5', '6', '7')) %>%
  ggplot(mapping = aes(x = E2, y = Pb)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("E2 vs Pb, Correlated Whole Blood (2023, IN)") +
  xlab("E2 (pg/mL)") +
  ylab("Pb Conc. (mcg/dL)") +
  stat_regline_equation(label.x = 75, label.y = 17) +
  stat_cor(method = "spearman", label.x = 65, label.y = 13.5, cor.coef.name = "rho", digits = 3)

Correlated_Pb_Hg_2023 %>%
  filter(State == "AK") %>%
  filter(Month %in% c('4', '5', '6', '7')) %>%
  ggplot(mapping = aes(x = E2, y = Pb)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("E2 vs Pb, Correlated Whole Blood (2023, AK)") +
  xlab("E2 (pg/mL)") +
  ylab("Pb Conc. (mcg/dL)") +
  stat_regline_equation(label.x = 75, label.y = 17) +
  stat_cor(method = "spearman", label.x = 65, label.y = 13.5, cor.coef.name = "rho", digits = 3)

calculateXI(Correlated_Pb_Hg_2023$E2, Correlated_Pb_Hg_2023$Pb)
# 0.1453298

# no correlation between E2 and Pb when filtered by state and year.

# now for IN only

Correlated_Pb_Hg_2023_IN <- Correlated_Pb_Hg %>%
  filter(Year == '2023') %>%
  filter(State == 'IN')
hist(Correlated_Pb_Hg_2023_IN$E2) # non-normal data

boxplot(Correlated_Pb_Hg_2023_IN$E2 ~ Correlated_Pb_Hg_2023_IN$Breeding)

stat.test <- Correlated_Pb_Hg_2023_IN %>%
  wilcox_test(E2 ~ Breeding) %>%
  add_significance()
stat.test
# no significant difference in E2 between breeding and non-breeding
# but only have one non-breeding from IN in 2023 so this doesn't hold
# much weight.
  
###### Testosterone vs Pb ######
hist(RBCs$T0)
shapiro.test(RBCs$T0) # non-normal 

leveneTest(Correlated_Pb_Hg_2023$T0 ~ Correlated_Pb_Hg_2023$State)
# non-sig p-value = 0.7498 so equal variances

kruskal.test(data = Correlated_Pb_Hg_2023, T0 ~ State)
# non-sig p-value = 0.1966 so non-sig difference in T0 between states

group_by(Correlated_Pb_Hg_2023, State) %>%
  summarise(
    count = n(),
    median = median(T0, na.rm = TRUE),
    IQR = IQR(T0, na.rm = TRUE)
  )

Correlated_Pb_Hg_2023 %>%
  dplyr::filter(State == "IN") %>%
  dplyr::select(T0, T2, Hg, Pb) %>%
  chart.Correlation(method = "spearman")

p2 <- Correlated_Pb_Hg_2023 %>%
  ggplot(mapping = aes(x = T0, y = Pb)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("T0 vs Pb, Correlated Blood (2023)") +
  xlab("T0 (ng/mL)") +
  ylab("Pb Conc. (mcg/dL)") +
  ylim(0,15) +
  stat_regline_equation(label.x = 2.3, label.y = 14.5) +
  stat_cor(method = "spearman", label.x = 2, label.y = 13.5, cor.coef.name = "rho", digits = 3)
# y=6.2-0.83x R = -0.323 p = 0.00876
# sig correlation between T0 and Pb

hist(Correlated_Pb_Hg_2023$T2)
shapiro.test(Correlated_Pb_Hg_2023$T2) # non-normal 

leveneTest(Correlated_Pb_Hg_2023$T2 ~ Correlated_Pb_Hg_2023$State)
# non-sig p-value = 0.08042 so equal variances

kruskal.test(data = Correlated_Pb_Hg_2023, T2 ~ State)
# non-sig p-value = 0.1001 so non-sig difference in T2 between states

group_by(Correlated_Pb_Hg_2023, State) %>%
  summarise(
    count = n(),
    median = median(T2, na.rm = TRUE),
    IQR = IQR(T2, na.rm = TRUE)
  )

p3 <- Correlated_Pb_Hg_2023 %>%
  ggplot(mapping = aes(x = T2, y = Pb)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("T2 vs Pb, Correlated Blood (2023)") +
  xlab("T2 (ng/mL)") +
  ylab("Pb Conc. (mcg/dL)") +
  ylim(0,15)+
  stat_regline_equation(label.x = 3, label.y = 14.5) +
  stat_cor(method = "spearman", label.x = 2.5, label.y = 13.5, cor.coef.name = "rho", digits = 3)
# y=3.4-0.064x R = -0.112 p = 0.413
# non-sig correlation

ggarrange(p1, p2, p3, ncol = 3, nrow = 1)

#### Mercury ####
###### Correlating Hg between Whole Blood and RBCs ############

ggplot(mapping = aes(x = Corr_RBCs$Hg, y = Corr_Whole$Hg)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_regline_equation(label.x = 150, label.y = 820) +
  stat_cor(label.x = 180, label.y = 750) +
  ggtitle("Correlated RBCs vs Whole Blood for Hg") +
  xlab("RBCs (ppm)") +
  ylab("Whole Blood (ppm)")
# Hg says it has an R2 of 1. Which is great. But doesn't seem real.
# I'm going to fact-check this.
rsq <- function(x, y) cor(x, y) ^ 2
rsq(Corr_RBCs$Hg, Corr_Whole$Hg)
#' Creating an rsquare function and running my data through returned 
#' an r2 of 0.9940806. I'm thinking it has something to do with the
#' outlier at the very top of the graph. So I think what I'll do is
#' cut out outliers and re-graph it. 
quartiles_Hg <- quantile(Corr_Whole$Hg, probs = c(.25, .75), na.rm = TRUE)
IQR_Hg <- IQR(Corr_Whole$Hg)
Lower_Hg <- quartiles_Hg[1] - 1.5*IQR_Hg
Upper_Hg <- quartiles_Hg[2] + 1.5*IQR_Hg
Corr_Whole_2 <- subset(Corr_Whole, Hg > Lower_Hg & Hg < Upper_Hg)

quartiles_Hg_2 <- quantile(Corr_RBCs$Hg, probs = c(.25, .75), na.rm = TRUE)
IQR_Hg_2 <- IQR(Corr_RBCs$Hg)
Lower_Hg_2 <- quartiles_Hg_2[1] - 1.5*IQR_Hg_2
Upper_Hg_2 <- quartiles_Hg_2[2] + 1.5*IQR_Hg_2
Corr_RBCs_2 <- subset(Corr_RBCs, Hg > Lower_Hg_2 & Hg < Upper_Hg_2)
# ok now that I've cut out outliers let me regraph it
ggplot(mapping = aes(x = Corr_RBCs_2$Hg, y = Corr_Whole_2$Hg)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_regline_equation(label.x = 72, label.y = 57) +
  stat_cor(label.x = 70, label.y = 54) +
  ggtitle("Correlated RBCs vs Whole Blood for Hg No Outliers") +
  xlab("RBCs (ppm)") +
  ylab("Whole Blood (ppm)")
#' Okay so with outliers removed, the r2 is 0.82 which is pretty good.
#' Now I'm thinking I should correlate for Hg. So I'll copy what I 
#' did for Pb.
Hg_Correlation <- function(Hg){
  0.456*Hg + 0.831
}
Correlated_All$Hg <- sapply(Correlated_All$Hg, Hg_Correlation)

Whole_Hg_Pb <- Whole_Blood[c(1:13, 20, 22, 24, 25)]
Correlated_Hg_Pb <- Correlated_All[c(1:13, 20, 22, 24, 25)]
All_Blood_Hg_Pb <- rbind(Whole_Hg_Pb, Correlated_Hg_Pb)
#' Now I have a file of all blood with Hg and Pb.This dataset,
#' "All_Blood_Hg_Pb" is comparable to the manually created excel
#' version titled,"Correlated_Pb_Hg".
# real quick for organization, I'll rename All_Blood to All_Blood_Pb
All_Blood_Pb <- All_Blood
#' Now that I have my correlations completed, I can start to think
#' about why you may want to use RBCs over whole blood. In my case,
#' it was a matter of sample availability. Whole blood hosts a huge
#' variety of uses and oftentimes, plasma and red blood cells are 
#' separated for analysis. Plasma on its own can be analyzed for
#' all sorts of things but hormones are a really common one. This 
#' leaves RBCs as a sort of waste material that aren't often utilized.
#' My correlations between RBCs and whole blood show that with some
#' forethought, RBCs can be used to measure for Hg and Pb when whole
#' blood is mostly unavailable. By setting aside some samples with
#' matching whole blood and RBCs and analyzing both, linear 
#' regressions can be generated. This allows for the conversion of 
#' RBC results to be equivalent to whole blood and render a larger,
#' more powerful dataset.

###### Seasonality in Blood  ######

ggplot(Correlated_Pb_Hg, mapping = aes(x = Month, y = Hg)) +
  geom_boxplot() +
  ggtitle("Month by Month Hg Concentrations") +
  xlab("Month") +
  ylab("Hg Conc. (mcg/dL)")
#' For Hg, the spring and summer months have higher averages,
#' meaning it could also be diet-related.

hist(Correlated_Pb_Hg$Hg) # need to remove some outliers.
Q <- quantile(Correlated_Pb_Hg$Hg, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Correlated_Pb_Hg$Hg)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
No_Outliers_Hg <- subset(Correlated_Pb_Hg, Correlated_Pb_Hg$Hg > Lower & Correlated_Pb_Hg$Hg < Upper)
hist(No_Outliers_Hg$Hg)
shapiro.test(No_Outliers_Hg$Hg) 
# non-normal p-value = 3.248e-15 need to transform
hist(log(No_Outliers_Hg$Hg))
shapiro.test(log(No_Outliers_Hg$Hg)) 

wilcox.test(data = No_Outliers_Hg, Hg ~ Season)
# sig p-value = 1.66e-08 means median Hg is significantly different
# between seasons.

group_by(No_Outliers_Hg, Season) %>%
  summarise(
    count = n(),
    median = median(Hg, na.rm = TRUE),
    IQR = IQR(Hg, na.rm = TRUE)
  )
# summer median = 233 winter median = 60.3

###### Locational Differences in Blood ######

boxplot(Correlated_Pb_Hg$Hg) # need to remove outliers
Q <- quantile(Correlated_Pb_Hg$Hg, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Correlated_Pb_Hg$Hg)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
No_Hg_Outliers <- subset(Correlated_Pb_Hg, Correlated_Pb_Hg$Hg > Lower & Correlated_Pb_Hg$Hg < Upper)
hist((No_Hg_Outliers$Hg))
shapiro.test((No_Hg_Outliers$Hg))
hist(1/(No_Hg_Outliers$Hg)^(1/3))
shapiro.test(1/(No_Hg_Outliers$Hg)^(1/3)) # can't force normal so must use non-parametric

wilcox.test(data = No_Hg_Outliers, Hg ~ State)
# sig p-value = 0.1551 so no sig difference in Hg between states

group_by(Correlated_Pb_Hg, State) %>%
  summarise(
    count = n(),
    median = median(Pb, na.rm = TRUE),
    IQR = IQR(Pb, na.rm = TRUE)
  )
# AK median = 2.28 IN median = 2.99

#' Now I want to look at just Indiana Pb median because the lower Alaska
#'  values may be skewing them lower.

Indiana <- No_Outliers_All_Pb %>%
  filter(State == 'IN')

ggplot(Indiana, mapping = aes(x = Month, y = Pb)) +
  geom_boxplot() +
  ggtitle("Month by Month Pb Concentrations No Outliers - Indiana") +
  xlab("Month") +
  ylab("Pb Conc. (mcg/dL)") +
  stat_summary(fun="mean", color = "red", shape = 13)

hist(Indiana$Pb)
shapiro.test(Indiana$Pb)
# non-normal p-value = 1.956e-10
hist(1/sqrt(Indiana$Pb))
shapiro.test(1/sqrt(Indiana$Pb))
# can't be forced normal so using wilcox test
wilcox.test(data = Indiana, Pb ~ Season)
# sig diff p-value = 4.981e-06
group_by(Indiana, Season) %>%
  summarise(
    count = n(),
    median = median(Pb, na.rm = TRUE),
    IQR = IQR(Pb, na.rm = TRUE)
  )
# Summer = 3.59 Winter = 2.59
###### Long-Term Changes in Bone ######

ggplot(Epiphysis, aes(x = Year, y = Hg)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, geom = "point", shape = 21, size = 2, fill = 'red') +
  ggtitle("Bone vs Hg") +
  xlab("Year") +
  ylab("Hg (ppm)")
mean(Diaphysis$Hg, na.rm = TRUE) # 2.518307
sd(Diaphysis$Hg, na.rm = TRUE) # 1.087374
mean(Epiphysis$Hg, na.rm = TRUE) # 4.234681
sd(Epiphysis$Hg, na.rm = TRUE) # 2.412718
mean(Keel$Hg) # 2.586711
sd(Keel$Hg) #0.5066992

# this plot shows a relatively steady Hg over time.

hist(Epiphysis$Hg)
shapiro.test(Epiphysis$Hg) # non-normal p-value = 2.567e-05 needs transformed
hist(1/(Epiphysis$Hg))
shapiro.test(1/(Epiphysis$Hg)) # normal p-value = 0.06891
Epiphysis_3 <- Epiphysis %>% 
  mutate(Hg = (1/(Epiphysis$Hg)))

leveneTest(data = Epiphysis_3, Hg ~ Set)
# non-sig p-value so equal variances

t.test(data = Epiphysis_3, Hg ~ Set)
# non-sig p-value = 0.9317 so no sig difference between Hg modern vs historic

###### Sex vs Hg #######

ggboxplot(Correlated_Pb_Hg, x = "Sex", y = "Hg", 
          color = "Sex", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("M", "F", "U"),
          ylab = "Hg Conc. (mcg/dL)", xlab = "Sex")
# Shows males have a higher median than females as well as more extreme
# outliers.

hist(Correlated_Pb_Hg$Hg)
# shows outliers so need to remove those

Q <- quantile(Correlated_Pb_Hg$Hg, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Correlated_Pb_Hg$Hg)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Correlated_Hg_No_Outliers <- subset(Correlated_Pb_Hg, Correlated_Pb_Hg$Hg > Lower & Correlated_Pb_Hg$Hg < Upper)
hist(Correlated_Hg_No_Outliers$Hg)
shapiro.test(Correlated_Hg_No_Outliers$Hg)
# non-normal p-value so must be transformed.

hist(log10(Correlated_Hg_No_Outliers$Hg))
shapiro.test(log(Correlated_Hg_No_Outliers$Hg))
# can't force to be normal so must use non-parametric tests

kruskal.test(data = Correlated_Hg_No_Outliers, Hg ~ Sex)
# sig p-value = 7.824e-07 means significant difference in medians

pairwise.wilcox.test(Correlated_Hg_No_Outliers$Hg, Correlated_Hg_No_Outliers$Sex,
                     p.adjust.method = "BH")
# sig difference between all groups 

group_by(Correlated_Hg_No_Outliers, Sex) %>%
  summarise(
    count = n(),
    median = median(Pb, na.rm = TRUE),
    IQR = IQR(Pb, na.rm = TRUE)
  )
# medians M = 2.99, F = 2.62, U = 3.74

ggboxplot(Correlated_Hg_No_Outliers, x = "Sex", y = "Hg", 
          color = "Sex", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("M", "F", "U"),
          ylab = "Hg Conc. (mcg/dL)", xlab = "Sex")

###### Age vs Hg ######

ggboxplot(Correlated_Pb_Hg, x = "Age", y = "Hg", 
          color = "Age", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("HY", "AHY", "NA"),
          ylab = "Hg Conc. (mcg/dL)", xlab = "Age")
# Shows AHY has a higher median than HY

hist(Correlated_Pb_Hg$Hg)
shapiro.test(Correlated_Pb_Hg$Hg) # non-normal p-value = <2.2e-16

# first I should get rid of some outliers
Q <- quantile(Correlated_Pb_Hg$Hg, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Correlated_Pb_Hg$Hg)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Correlated_Hg_No_Outliers <- subset(Correlated_Pb_Hg, Correlated_Pb_Hg$Hg > Lower & Correlated_Pb_Hg$Hg < Upper)

hist(Correlated_Hg_No_Outliers$Hg)
shapiro.test(Correlated_Hg_No_Outliers$Hg) # non-normal p-value = 3.248e-15
# need to transform
hist(sqrt(Correlated_Hg_No_Outliers$Hg))
shapiro.test(log(Correlated_Hg_No_Outliers$Hg))

# can't normalize any closer than p-value =  1.4e-09 so must use 
# non-parametric tests

# I'm using the Kruskal-Wallis test because my data is non-normal.

kruskal.test(Hg ~ Age, data = Correlated_Hg_No_Outliers)
# sig p-value = 0.002712 means there is a significant difference in
# distribution between mean Hg by age group.

pairwise.wilcox.test(Correlated_Hg_No_Outliers$Hg, Correlated_Hg_No_Outliers$Age,
                     p.adjust.method = "BH")
#    AHY   HY   
# HY 0.010 -    
# U  0.072 0.609
# sig difference between HY and AHY

group_by(Correlated_Hg_No_Outliers, Age) %>%
  summarise(
    count = n(),
    median = median(Hg, na.rm = TRUE),
    IQR = IQR(Hg, na.rm = TRUE)
  )
# medians AHY = 162, HY = 43.2 U = 30.6

###### Fat vs Hg ######

hist(No_Outliers_All_Hg$Hg)
shapiro.test(No_Outliers_All_Hg$Hg) # non-normal p-value
kruskal.test(data = No_Outliers_All_Hg, Hg ~ Fat)
# non-sig p-value = 0.2999

kruskal.test(data = Correlated_Pb_Hg, Hg ~ condition_resid)
# non-sig p-value = 0.09484

###### Estrogen vs Hg ######
# noticed a strong correlation in RBCs between Hg and E2 so wanted
# to dig into that a little more.

Hg_Filtered <- Correlated_Pb_Hg %>%
  filter(Hg < 1000)

boxplot(Hg_Filtered$Hg ~ Hg_Filtered$State)

hist(Hg_Filtered$E2)
shapiro.test(Hg_Filtered$E2) # non-normal p-value = 0.01241

leveneTest(Hg_Filtered$E2 ~ Hg_Filtered$State)
# non-sig p-value = 0.1851 so equal variances

kruskal.test(data = Hg_Filtered, E2 ~ State)
# sig p-value = 0.001302 so sig difference in E2 between states

group_by(Hg_Filtered, State) %>%
  summarise(
    count = n(),
    median = median(E2, na.rm = TRUE),
    IQR = IQR(E2, na.rm = TRUE)
  )
# median AK = 158, IN = 80.0

Hg_Filtered_2023 <- Hg_Filtered %>%
  filter(Year == '2023')

Hg_Filtered_2023 %>%
  dplyr::select(E2, Hg, Pb) %>%
  chart.Correlation(method = "spearman", exact = FALSE)
# shows significant correlation *** between Hg and E2 rho = -0.49

p1 <- Hg_Filtered_2023 %>%
  ggplot(mapping = aes(x = E2, y = Hg)) +
  geom_point(aes(color = State)) +
  geom_smooth(method = "lm") +
  ggtitle("E2 vs Hg, Correlated Blood (2023)") +
  ylim(0, 1000) +
  xlab("E2 (pg/mL)") +
  ylab("Hg (mcg/dL)") +
  stat_regline_equation(label.x = 75, label.y = 1000) +
  stat_cor(method = "spearman", label.x = 60, label.y = 900, cor.coef.name = "rho", digits = 3)
# y =550-1.8x R = -0.487 p = 0.00164

# is Hg significantly different by state?
leveneTest(Hg_Filtered_2023$Hg ~ Hg_Filtered_2023$State)
# sig p-value = 1.155e-08 so unequal variances must run non-parametric test

kruskal.test(data = Hg_Filtered_2023, Hg ~ State)
# sig p-value = 0.0005086 so sig difference in Hg between states

group_by(Hg_Filtered_2023, State) %>%
  summarise(
    count = n(),
    median = median(Hg, na.rm = TRUE),
    IQR = IQR(Hg, na.rm = TRUE)
  )
# median AK = 74.1, IN = 327

p2 <- Hg_Filtered_2023 %>%
  filter(State == "AK") %>%
  ggplot(mapping = aes(x = E2, y = Hg)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("E2 vs Hg, Correlated Blood (2023, AK)") +
  ylim(0, 1000) +
  xlab("E2 (pg/mL)") +
  ylab("Hg (mcg/dL)") +
  stat_regline_equation(label.x = 75, label.y = 1000) +
  stat_cor(method = "spearman", label.x = 60, label.y = 900, cor.coef.name = "rho", digits = 3)
# y =76+0.25x R = 0.339 p = 0.216

p3 <- Hg_Filtered_2023 %>%
  filter(State == "IN") %>%
  ggplot(mapping = aes(x = E2, y = Hg)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("E2 vs Hg, Correlated Blood (2023, IN)") +
  ylim(0, 1000) +
  xlab("E2 (pg/mL)") +
  ylab("Hg (mcg/dL)") +
  stat_regline_equation(label.x = 75, label.y = 1000) +
  stat_cor(method = "spearman", label.x = 60, label.y = 900, cor.coef.name = "rho", digits = 3)
# y =590-1.4x R = -0.308 p = 0.134

ggarrange(p1, p2, p3, ncol = 3, nrow = 1)

# interesting that each state has opposite Hg trends (non-significant) associated
# with E2.

###### Testosterone vs Hg ###### 

Correlated_Pb_Hg_2023 %>%
  dplyr::select(T0, T2, Hg, Pb) %>%
  chart.Correlation(method = "spearman")
# no significance between testosterone and Hg

#### Arsenic ####
###### Seasonality in Blood ######

ggplot(Whole_Blood, mapping = aes(x = Month, y = As)) +
  geom_boxplot() +
  ggtitle("Month by Month As Concentrations") +
  xlab("Month") +
  ylab("As Conc. (mcg/dL)")
#' As concentrations seem to peak in June and otherwise trend 
#' downwards through the year.

hist(Whole_Blood$As)
Whole_Blood_No_Outliers <- Whole_Blood_No_Outliers %>%
  filter(As < 1.8)
hist(Whole_Blood_No_Outliers$As)
shapiro.test(Whole_Blood_No_Outliers$As)
# normal p-value = 0.1319 no transformation needed

t.test(data = Whole_Blood_No_Outliers, As ~ Season)
# sig p-value = 9.659e-05 means sig difference is As by season.
# mean summer = 0.9942620 winter = 0.6216037

###### Long-Term Changes in Bone ######

ggplot(Epiphysis, aes(x = Year, y = As)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, geom = "point", shape = 21, size = 2, fill = 'red') +
  ggtitle("Bone vs As") +
  xlab("Year") +
  ylab("As (ppm)")
geoMean(Diaphysis$As) # 2.788084
geoSD(Diaphysis$As) # 3.14769
geoMean(Epiphysis$As, na.rm = TRUE) # 3.129889
geoSD(Epiphysis$As, na.rm = TRUE) # 3.397847
geoMean(Keel$As) # 7.269598
geoSD(Keel$As) # 2.39412
#' Epiphysis and Diaphysis don't show much change over time. Keel shows
#' an increase in mean between 1980 and 1985. The Keel has the 
#' highest concentration of As, then epiphysis then diaphysis. The
#' keel has the smallest SD.

hist(Epiphysis$As)
# need to remove outlier
Q <- quantile(Epiphysis$As, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Epiphysis$As, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Epiphysis_No_As_Outliers <- subset(Epiphysis, Epiphysis$As > Lower & Epiphysis$As < Upper)
hist(Epiphysis_No_As_Outliers$As)
shapiro.test(Epiphysis_No_As_Outliers$As) # non-normal p-value = 0.04925 need to transform
hist((Epiphysis_No_As_Outliers$As)^(1/3))
shapiro.test((Epiphysis_No_As_Outliers$As)^(1/3)) # normal p-value = 0.9724

Epiphysis_2 <- Epiphysis_No_As_Outliers %>%
  mutate(As = (Epiphysis_No_As_Outliers$As)^(1/3))

leveneTest(data = Epiphysis_2, Zn ~ Set)
# non-sig p-value so equal variances

t.test(data = Epiphysis_2, As ~ Set)
# sig p-value = 0.04346
# transformed means historic = 1.579533 modern = 1.202658
# non-transformed means historic = 3.940816 modern = 1.739508

# let me repeat this process removing negative values
Epiphysis_2$As <- replace(Epiphysis_2$As, which(Epiphysis_2$As < 0), NA)
hist(Epiphysis_2$As)
shapiro.test(Epiphysis_2$As) # non-normal p-value = 0.0002263 need to transform
hist(log(Epiphysis_2$As))
shapiro.test(log(Epiphysis_2$As)) # normal p-value = 0.9999

leveneTest(data = Epiphysis_2, As ~ Set) # non-sig p-value so equal variances
t.test(data = Epiphysis_2, As ~ Set) # non-sig p-value = 0.2217
# historic mean = 5.957 modern = 3.078

###### Sex vs As ######

hist(Whole_Blood$As) # may have outliers
Q <- quantile(Whole_Blood$As, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$As, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_No_As_Outliers <- subset(Whole_Blood, Whole_Blood$As > Lower & Whole_Blood$As < Upper)

hist(Whole_Blood_No_As_Outliers$As)
shapiro.test(Whole_Blood_No_As_Outliers$As) #
# normal p-value = 0.06791 so can run levenes test
leveneTest(data = Whole_Blood_No_As_Outliers, As ~ Sex)
# non-sig p-value = 0.8382 so no sig difference in mean As by sex.

my.aov_Sex_As <- aov(data = Whole_Blood_No_As_Outliers, As ~ Sex)
summary(my.aov_Sex_As)
# non-sig p-value = 0.817 so no sig difference in mean by sex


###### Age vs As ######

hist(Whole_Blood$As) # may have outliers
Q <- quantile(Whole_Blood$As, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$As, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_No_As_Outliers <- subset(Whole_Blood, Whole_Blood$As > Lower & Whole_Blood$As < Upper)

hist(Whole_Blood_No_As_Outliers$As)
shapiro.test(Whole_Blood_No_As_Outliers$As) #
# normal p-value = 0.06791 so can run levenes test
leveneTest(data = Whole_Blood_No_As_Outliers, As ~ Age)
# non-sig p-value = 0.1506 so equal variances and can run anova

my.aov_Age_As <- aov(data = Whole_Blood_No_As_Outliers, As ~ Age)
summary(my.aov_Age_As)
# non-sig p-value = 0.188 so no sig difference in means between ages.

###### Fat vs As ######
hist(Whole_Blood$As)
# need to remove outliers
Q <- quantile(Whole_Blood$As, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$As, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_No_As_Outliers <- subset(Whole_Blood, Whole_Blood$As > Lower & Whole_Blood$As < Upper)

hist(Whole_Blood_No_As_Outliers$As)
shapiro.test(Whole_Blood_No_As_Outliers$As) # normal p-value = 0.06791
leveneTest(data = Whole_Blood_No_As_Outliers, As ~ Fat)
#      Df F value Pr(>F)
#group  4  0.7594 0.5548
#      78               
# Non-sig p-value means equal variances and can run ANOVA
my.aov_Fat_As <- aov(data = Whole_Blood_No_As_Outliers, As ~ Fat)
summary(my.aov_Fat_As)
# non-sig p-value so no sig difference between fat classes for As

mod=lm(Mass ~ Wing,data=Whole_Blood_No_As_Outliers,na.action=na.exclude)
Whole_Blood_No_As_Outliers$condition_resid=resid(mod)
kruskal.test(data = Whole_Blood_No_As_Outliers, As ~ condition_resid)
# non-sig p-value = 0.4443

#### Zinc ####
###### Seasonality in Blood ######

ggplot(Whole_Blood, mapping = aes(x = Month, y = Zn)) +
  geom_boxplot() +
  ggtitle("Month by Month Zn Concentrations") +
  xlab("Month") +
  ylab("Zn Conc. (mcg/dL)")
# Zn concentrations are very steady through the year with the exception of outliers.

hist(Whole_Blood$Zn)
Whole_Blood_No_Outliers <- Whole_Blood_No_Outliers %>%
  filter(Zn < 550)
hist(Whole_Blood_Zn_No_Outliers$Zn)
shapiro.test(Whole_Blood_Zn_No_Outliers$Zn) 
# non-normal p-value = 0.0111 needs transformed

hist(sqrt(max(Whole_Blood_Zn_No_Outliers$Zn+1) - Whole_Blood_Zn_No_Outliers$Zn))
shapiro.test(sqrt(max(Whole_Blood_Zn_No_Outliers$Zn+1) - Whole_Blood_Zn_No_Outliers$Zn))
# pretty close to normal p-value = 0.04952

Whole_Blood_Zn_No_Outliers$Zn <- sqrt(max(Whole_Blood_Zn_No_Outliers$Zn+1) - Whole_Blood_Zn_No_Outliers$Zn)

leveneTest(data = Whole_Blood_Zn_No_Outliers, Zn ~ Season)
# no sig p-value so equal variances

t.test(data = Whole_Blood_Zn_No_Outliers, Zn ~ Season)
# non-sig p-value = 0.1496 no sig difference in Zn between
# seasons.

###### Long-Term Changes in Bone ######

ggplot(Epiphysis, aes(x = Year, y = Zn)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, geom = "point", shape = 21, size = 2, fill = 'red') +
  ggtitle("Bone vs Zn") +
  xlab("Year") +
  ylab("Zn (ppm)")
geoMean(Diaphysis$Zn, na.rm = TRUE) # 35.19094
geoSD(Diaphysis$Zn, na.rm = TRUE) # 5.329661
geoMean(Epiphysis$Zn, na.rm = TRUE) # 103.0343
geoSD(Epiphysis$Zn, na.rm = TRUE) # 4.60849
geoMean(Keel$Zn) # 247.4957
geoSD(Keel$Zn) # 3.406634
#' Zn epiphysis doesn't reveal anything interesting.Diaphysis somewhat
#' shows a downturn.Interestingly, the keel shows a dramatic decrease
#' from 1980 to 1985 in mean. The means of Zn go in order from highest
#' to lowest: Keel, Epiphysis, Diaphysis while SD is opposite.

hist(Epiphysis$Zn)
shapiro.test(Epiphysis$Zn)
# non-normal p-value 4.217e-07 needs outliers removed and transformed
Epiphysis_2 <- Epiphysis %>%
  filter(Zn < 1000)
hist(Epiphysis_2$Zn)
shapiro.test(Epiphysis_2$Zn) # non-normal p-value = 0.001035
hist(sqrt(Epiphysis_2$Zn))
shapiro.test(sqrt(Epiphysis_2$Zn)) # normal p-value = 0.1802
Epiphysis_2$Zn <- sqrt(Epiphysis_2$Zn)

# now I can do levenes test
leveneTest(data = Epiphysis_2, Zn ~ Set)
# non-significant therefore equal variances and now I can run t-test

t.test(data = Epiphysis_2, Zn ~ Set)
#' So this shows a non-sig p-value = 0.3336 so there is no sig 
#' difference between means of historic and modern datasets for Zn.

###### Sex vs Zn ######

hist(Whole_Blood_Zn_No_Outliers$Zn)
shapiro.test(Whole_Blood_Zn_No_Outliers$Zn) 
# non-sig p-value = 0.5798 so can run levenes test
leveneTest(data = Whole_Blood_Zn_No_Outliers, Zn ~ Sex)
# non-sig p-value = 0.4042 so can run ANOVA
my.aov_Sex_Zn <- aov(data = Whole_Blood_Zn_No_Outliers, Zn ~ Sex)
summary(my.aov_Sex_Zn)
# non-sig p-value = 0.519 so no sig difference in mean Zn by sex.

###### Age vs Zn ######

hist(Whole_Blood$Zn)
# need to remove outliers
Whole_Blood_Zn_No_Outliers <- Whole_Blood %>%
  filter(Zn < 550)
hist(Whole_Blood_Zn_No_Outliers$Zn)
shapiro.test(Whole_Blood_Zn_No_Outliers$Zn) 

# non-normal p-value = 0.0111 need to transform
hist(sqrt(max(Whole_Blood_Zn_No_Outliers$Zn+1) - Whole_Blood_Zn_No_Outliers$Zn))
shapiro.test(sqrt(max(Whole_Blood_Zn_No_Outliers$Zn+1) - Whole_Blood_Zn_No_Outliers$Zn))

hist((Whole_Blood_Zn_No_Outliers$Zn)^2)
shapiro.test((Whole_Blood_Zn_No_Outliers$Zn)^2)
# normal p-value = 0.5798

Whole_Blood_Zn_No_Outliers$Zn <- (Whole_Blood_Zn_No_Outliers$Zn)^2
leveneTest(data = Whole_Blood_Zn_No_Outliers, Zn ~ Age)
# slight significance p-value = 0.01186 so unequal variances between
# groups. Need to run non-parametric test, Kruskal-Wallis.
kruskal.test(data = Whole_Blood_Zn_No_Outliers, Zn ~ Age)
# non-sig p-value = .8212 so no sig difference in distributions

###### Fat vs Zn ######

hist(Whole_Blood$Zn) # outliers
Q <- quantile(Whole_Blood$Zn, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$Zn, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_No_Zn_Outliers <- subset(Whole_Blood, Whole_Blood$Zn > Lower & Whole_Blood$Zn < Upper)
hist(Whole_Blood_No_Zn_Outliers$Zn)
shapiro.test(Whole_Blood_No_Zn_Outliers$Zn) # normal p-value = 0.7053
leveneTest(data = Whole_Blood_No_Zn_Outliers, Zn ~ Fat) # non-sig p-value so equal variances
my.aov_Fat_Zn <- aov(data = Whole_Blood_No_Zn_Outliers, Zn ~ Fat)
summary(my.aov_Fat_Zn)
# non-sig p-value = 0.618 so no sig differences in mean Zn by fat class.

mod=lm(Mass ~ Wing,data=Whole_Blood_No_Zn_Outliers,na.action=na.exclude)
Whole_Blood_No_Zn_Outliers$condition_resid=resid(mod)
kruskal.test(data = Whole_Blood_No_Zn_Outliers, Zn ~ condition_resid)
# non-sig p-value = 0.4301

#### Copper ####
###### Seasonality in Blood ######

Whole_Blood$Month <- factor(Whole_Blood$Month, levels = c("1", "2", "3", "4", "5", "6", "9", "10", "11", "12"))
ggplot(Whole_Blood, mapping = aes(x = Month, y = Cu)) +
  geom_boxplot() +
  ggtitle("Month by Month Cu Concentrations") +
  xlab("Month") +
  ylab("Cu Conc. (mcg/dL)")
#' Cu concentrations stay relatively steady throughout the year.

hist(Whole_Blood$Cu)
Whole_Blood_No_Cu_Outliers <- Whole_Blood %>%
  filter(between(Cu, 1, 50))
hist(Whole_Blood_No_Cu_Outliers$Cu)
shapiro.test(Whole_Blood_No_Cu_Outliers$Cu)
# non-normal p-value = 0.0008619 needs transformed

hist((Whole_Blood_No_Cu_Outliers$Cu)^3)
shapiro.test(log(Whole_Blood_No_Cu_Outliers$Cu)^3)
# can't get it normal so I will try a wilcox test

wilcox.test(data = Whole_Blood_No_Cu_Outliers, Cu ~ Season)
# non-sig p-value so Cu mean doesn't vary by season.

###### Long-Term Changes in Bone ######

ggplot(Diaphysis, aes(x = Year, y = Cu)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, geom = "point", shape = 21, size = 2, fill = 'red') +
  ggtitle("Bone vs Cu") +
  xlab("Year") +
  ylab("Cu (ppm)")
geoMean(Diaphysis$Cu)# 16.6498
geoSD(Diaphysis$Cu) # 1.63777
geoMean(Epiphysis$Cu)# 24.33961
geoSD(Epiphysis$Cu) # 3.605915
geoMean(Keel$Cu) # 162.2879
geoSD(Keel$Cu) # 3.707069
#' Epiphysis doesn't show any clear trends but diaphysis shows a clear
#' downward trend of the mean disregarding 1997. Keel doesn't have 
#' enough data points to be very meaningful. On average the epiphysis
#' has more Cu than the diaphysis, and the keel has drastically more
#' than the entire femur.The SD of the diaphysis is smallest indicating
#' there is less spread there. The keel has the greatest Cu SD.

# I'm thinking I should check how normal it is and run some other kind
# of statistical test.
hist(Epiphysis$Cu)
shapiro.test(Epiphysis$Cu) 
# non-normal p-value = 0.006892 so needs transformed

hist(sqrt(Epiphysis$Cu))
shapiro.test(sqrt(Epiphysis$Cu)) # normal p-value = 0.5055

Epiphysis_2 <- Epiphysis %>%
  mutate(Cu = sqrt(Cu))

# I need to check for equal variances as well.
leveneTest(data = Epiphysis_2, Cu ~ Year)
#Levene's Test for Homogeneity of Variance (center = median)
#      Df F value Pr(>F)
#group  3  0.4548 0.7173
#      17
# likely equal variances due to no significance

t.test(Cu ~ Set, data = Epiphysis_2) 
# non-sig p-value = 0.1001 means no significant difference in mean Cu
# by group, Modern/Historic.

###### Sex vs Cu ######

hist(Whole_Blood$Cu)
No_Outliers_No_Outliers <- Whole_Blood %>%
  filter(between(Cu, 1, 50))
hist(Whole_Blood_Cu_No_Outliers$Cu)
shapiro.test(Whole_Blood_Cu_No_Outliers$Cu)
# normal p-value = 0.1083 can run levenes next
leveneTest(data = Whole_Blood_Cu_No_Outliers, Cu ~ Sex)
# non-sig p-value so equal variances
my.aov_Sex_Cu <- aov(data = Whole_Blood_Cu_No_Outliers, Cu ~ Sex)
summary(my.aov_Sex_Cu)
# non-sig p-value = 0.751 so no sig difference in mean Cu by sex.

###### Age vs Cu ######

hist(Whole_Blood$Cu)
shapiro.test(Whole_Blood$Cu) #non-normal p-value = 8.942e-06

# need to remove outliers and normalize
Q <- quantile(Whole_Blood$Cu, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$Cu, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_No_Outliers <- subset(Whole_Blood, Whole_Blood$Cu > Lower & Whole_Blood$Cu < Upper)
hist(Whole_Blood_No_Outliers$Cu)
shapiro.test(Whole_Blood_No_Outliers$Cu) # normal p-value = 0.1287
# no transformations needed and can use parametric tests
leveneTest(data = Whole_Blood_No_Outliers, Cu ~ Age)
# non-sig p-value = 0.4529 so can run ANOVA
my.aov_Age_Cu <- aov(data = Whole_Blood_No_Outliers, Cu ~ Age)
summary(my.aov_Age_Cu)
# non-sig p-value = 0.122 so no sig difference in mean Cu by age.

###### Fat vs Cu ######

hist(Whole_Blood_No_Outliers$Cu)
shapiro.test(Whole_Blood_No_Outliers$Cu) 
# normal p-value = 0.1083
leveneTest(data = Whole_Blood_No_Outliers, Cu ~ Fat)
# non-sig p-value equal variances
Whole_Blood_No_Outliers %>%
  dplyr::filter(Fat == 'N' | Fat == 'T' | Fat == 'L') %>%
  dplyr::group_by(Fat) %>%
  summarise(`W Statistic` = shapiro.test(Cu)$statistic,
            `p-value` = shapiro.test(Cu)$p.value)
# non-sig p-values so I can run an ANOVA

my.aov_Cu <- aov(data = Whole_Blood_No_Outliers, Cu ~ Fat)
summary(my.aov_Cu)
# non-sig p-value so no sig differences in mean Cu by fat class.

mod=lm(Mass ~ Wing,data=Whole_Blood_No_Outliers,na.action=na.exclude)
Whole_Blood_No_Outliers$condition_resid=resid(mod)
kruskal.test(data = Whole_Blood_No_Outliers, Cu ~ condition_resid)
# non-sig p-value = 0.4744

#### Inter-Elemental Relationships in Blood ####

# to start out I'm going to run histograms for each element in blood
# to see if things need transformed or outliers removed.

# Cu
hist(Whole_Blood$Cu)
Q <- quantile(Whole_Blood$Cu, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$Cu, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_2 <- subset(Whole_Blood, Whole_Blood$Cu > Lower & Whole_Blood$Cu < Upper)
hist(Whole_Blood_2$Cu)
shapiro.test(Whole_Blood_2$Cu) # normal p-value = 0.1083

# Zn
hist(Whole_Blood$Zn)
Q <- quantile(Whole_Blood$Zn, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$Zn, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_2 <- subset(Whole_Blood_2, Whole_Blood_2$Zn > Lower & Whole_Blood_2$Zn < Upper)
hist(Whole_Blood_2$Zn)
shapiro.test(Whole_Blood_2$Zn) # normal p-value = 0.9708

# As
hist(Whole_Blood$As)
Q <- quantile(Whole_Blood$As, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$As, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_2 <- subset(Whole_Blood_2, Whole_Blood_2$As > Lower & Whole_Blood_2$As < Upper)
hist(Whole_Blood_2$As)
shapiro.test(Whole_Blood_2$As) # normal p-value = 0.1319

# Hg
hist(Whole_Blood$Hg)
Q <- quantile(Whole_Blood$Hg, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$Hg, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_2 <- subset(Whole_Blood_2, Whole_Blood_2$Hg > Lower & Whole_Blood_2$Hg < Upper)
hist(Whole_Blood_2$Hg)
shapiro.test(Whole_Blood_2$Hg) # non-normal p-value = 8.199e-12 needs transformed
hist(-1/(Whole_Blood_2$Hg))
shapiro.test(-1/(Whole_Blood_2$Hg))

# trying boxcox transformation
y = Whole_Blood_2$Hg
hist(y,breaks = 12)
result = MASS::boxcox(y~1, lambda = seq(-5,5,0.5))
mylambda = result$x[which.max(result$y)]
mylambda
y2 = (y^mylambda-1)/mylambda
hist(y2)
shapiro.test(y2) # still non-normal p-value = 0.018

# trying bestNormalize
bestNormalize(y) # picked orderNorm
orderNorm_obj <- orderNorm(y)
hist(orderNorm_obj$x.t)
shapiro.test(orderNorm_obj$x.t) # normal p-value = 1

# now to replace my Hg in Whole_Blood_2 with the transformed version
Whole_Blood_3 <- Whole_Blood_2 %>%
  dplyr::mutate(Hg = orderNorm_obj$x.t)

# now Hg is normal but any results will need to be backconverted.
hist(Whole_Blood_2$Hg)
shapiro.test(Whole_Blood_2$Hg)

# Pb
hist(Whole_Blood$Pb)
Q <- quantile(Whole_Blood$Pb, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$Pb, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_2 <- subset(Whole_Blood_2, Whole_Blood_2$Pb > Lower & Whole_Blood_2$Pb < Upper)
hist(Whole_Blood_2$Pb)
shapiro.test(Whole_Blood_2$Pb) # non-normal p-value = 5.725e-06 needs transformed
hist(log(Whole_Blood_2$Pb))
shapiro.test(log(Whole_Blood_2$Pb)) # normal p-value = 0.06522
# need to replace Pb in Whole_Blood_2 with transformed values
Whole_Blood_3 <- Whole_Blood_3 %>%
  mutate(Pb = log(Pb))
# now Pb is normal but any results will need backconverted
hist(Whole_Blood_2$Pb)
shapiro.test(Whole_Blood_2$Pb)

Whole_Blood_4 <- Whole_Blood_3 %>%
  dplyr::select(Cu, Zn, As, Hg, Pb) # this is no outliers and transformed
p1 <- chart.Correlation(Whole_Blood_4)

Whole_Blood_5 <- Whole_Blood_2 %>%
  dplyr::select(Cu, Zn, As, Hg, Pb) # this is no outliers non-transformed
p2 <- chart.Correlation((Whole_Blood_5))

# now I'm going to run pairs to look for potential correlation.
pairs(Whole_Blood_4)
chart.Correlation(Whole_Blood_4)
colnames(Whole_Blood_4) = c('Cu', 'Zn', 'As', 'orderNorm(Hg)', 'log(Pb)')
chart.Correlation(Whole_Blood_4)
# seems to show some potential correlations between As and Hg, Hg and Pb,
# and maybe As and Pb. 

chart.Correlation(Whole_Blood_4)

# now to run cor.test to check for correlation.
cor.test(Whole_Blood_3$As, Whole_Blood_3$Hg)
# sig p-value = 1.943e-05 cor = 0.5213629
# strong correlation between As and orderNorm(Hg)
cor.test(Whole_Blood_3$Hg, Whole_Blood_3$Pb)
# sig p-value = 0.03019, cor = 0.2800972
# weak to moderate correlation between log(Pb) and orderNorm(Hg)
cor.test(Whole_Blood_3$As, Whole_Blood_3$Pb)
# sig p-value = 0.0002456, cor = 0.4565793
# moderate correlation between log(Pb) and As

# I'm going to repeat this on non-transformed but outlier stripped
# data using a more robust correlation test, spearman's rank correlation.

pairs(Whole_Blood_5)
p2 <- chart.Correlation(Whole_Blood_5, method = "spearman")
# shows potential correlation between As/Hg, Hg/Pb, As/Pb, Zn/As.
cor.test(Whole_Blood_3$As, Whole_Blood_3$Hg,
                  alternative = ("two.sided"),
                  method = ("spearman"),
                  exact = NULL, conf.level = 0.95, continuity = FALSE)
# sig p-value = 3.138e-05 rho = 0.5164768
cor.test(Whole_Blood_3$Hg, Whole_Blood_3$Pb,
         alternative = ("two.sided"),
         method = ("spearman"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)
# non-sig p-value = 0.06447 rho = 0.2404001
cor.test(Whole_Blood_3$As, Whole_Blood_3$Pb,
         alternative = ("two.sided"),
         method = ("spearman"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)
# sig p-value = 0.0003038 rho = 0.4547374
cor.test(Whole_Blood_3$Zn, Whole_Blood_3$As,
         alternative = ("two.sided"),
         method = ("spearman"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)
# non-sig p-value = 0.2472 rho = -0.1514865
cor.test(Whole_Blood_3$Zn, Whole_Blood_3$Cu,
         alternative = ("two.sided"),
         method = ("spearman"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)
# sig p-value = 0.006288 rho = 0.3505418


###### PCA ######

# going to remove outliers from data first

# a function for identifying outliers
outliers <- function(x) {
  Q1 <- quantile(x, probs = 0.25, na.rm = TRUE) # calculate first quantile
  Q3 <- quantile(x, probs = 0.75, na.rm = TRUE) # calculate third quantile
  IQR <- Q3 - Q1 # calculate IQR
  
  upper_limit = Q3 + 1.5*IQR # define upper limit
  lower_limit = Q1 - 1.5*IQR # define lower limit
  
  x > upper_limit | x < lower_limit # return true or false
}

remove_outliers <- function(Whole_Blood, cols = names(Whole_Blood)) {
  # for loop to traverse in columns vector
  for (col in cols) {
    # remove observation if it satisfied outlier function
    Whole_Blood_2 <- Whole_Blood[!outliers(Whole_Blood[[col]]), ]
  }
  # return dataframe
  Whole_Blood_2
}

# run function and create new dataframe
Whole_Blood_2 <- remove_outliers(Whole_Blood, c("Cu", "Zn", "As", "Hg", "Pb"))

# individual manual IQRs 
# Cu
hist(Whole_Blood$Cu)
Q <- quantile(Whole_Blood$Cu, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$Cu, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_2 <- subset(Whole_Blood, Whole_Blood$Cu > Lower & Whole_Blood$Cu < Upper)
hist(Whole_Blood_2$Cu)
shapiro.test(Whole_Blood_2$Cu)

#Zn
hist(Whole_Blood$Zn)
Q <- quantile(Whole_Blood$Zn, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$Zn, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_2 <- subset(Whole_Blood_2, Whole_Blood_2$Zn > Lower & Whole_Blood_2$Zn < Upper)
hist(Whole_Blood_2$Zn)
shapiro.test(Whole_Blood_2$Zn)

# As
hist(Whole_Blood$As)
Q <- quantile(Whole_Blood$As, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$As, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_2 <- subset(Whole_Blood_2, Whole_Blood_2$As > Lower & Whole_Blood_2$As < Upper)
hist(Whole_Blood_2$As)
shapiro.test(Whole_Blood_2$As) # normal p-value = 0.1319

# Hg
hist(Whole_Blood$Hg)
Q <- quantile(Whole_Blood$Hg, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$Hg, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_2 <- subset(Whole_Blood_2, Whole_Blood_2$Hg > Lower & Whole_Blood_2$Hg < Upper)
hist(Whole_Blood_2$Hg)
shapiro.test(Whole_Blood_2$Hg)

# Pb
hist(Whole_Blood$Pb)
Q <- quantile(Whole_Blood$Pb, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Whole_Blood$Pb, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Whole_Blood_2 <- subset(Whole_Blood_2, Whole_Blood_2$Pb > Lower & Whole_Blood_2$Pb < Upper)
hist(Whole_Blood_2$Pb)
shapiro.test(Whole_Blood_2$Pb)


# eliminate NAs and select the variables I'm interested in
Whole_Blood_3 <- Whole_Blood_2 %>%
  dplyr::select(Sex, Cu, Zn, As, Hg, Pb) %>%
  na.omit()

# normalize data
Whole_Blood_3[,-1] <- scale(Whole_Blood_3[,-1])

apply(Whole_Blood_3[,-1], 2, var)
Whole_Blood.pca <- prcomp(Whole_Blood_3[,-1], center = TRUE)
Whole_Blood.pca 
summary(Whole_Blood.pca)

fviz_eig(Whole_Blood.pca, addlabels=TRUE, hjust = -0.3) +
  geom_hline(yintercept = (1/ncol(Whole_Blood_3[,-1]))*100, linetype = "dashed")
# first two components most important, maybe 3

Whole_Blood.pca$rotation[,1:3]

cutoff <- sqrt(1/ncol(Whole_Blood_3[,-1]))
cutoff

foo <- data.frame(Whole_Blood.pca$rotation[,1:3]) %>%
  dplyr::mutate_all( function(x){replace(x, abs(x) < cutoff, NA)} )
row.names(foo) <- row.names(Whole_Blood.pca$rotation)
foo

ggbiplot(Whole_Blood.pca, groups = Whole_Blood_3$Sex)

 ###### NMDS ######
Whole_Blood_2 <- Whole_Blood %>%
  dplyr::select(Cu, As, Zn, Hg, Pb)
foo.dist <- vegdist(Whole_Blood_2, method = "bray", na.rm = TRUE)

Whole_Blood_2.mds <- metaMDS(comm = foo.dist, k = 2)
MDS_xy <- data.frame(Sex = Whole_Blood$Sex, Whole_Blood_2.mds$points)
ggplot(MDS_xy, aes(MDS1, MDS2, color = Sex)) + 
  geom_point(size =3) + 
  coord_equal() +
  theme_bw()

stressplot(Whole_Blood_2.mds)

vf <- envfit(Whole_Blood_2.mds, Whole_Blood_2, perm = 999, na.rm = TRUE)

spp.scrs <- data.frame(scores(vf, display = "vectors"))

ggplot(data = MDS_xy) +
  geom_point(mapping = aes(x = MDS1, y = MDS2, colour = Sex)) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1*.75, y = 0, yend = NMDS2*.75),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = spp.scrs, aes(x = NMDS1*.75, y = NMDS2*.75, label = rownames(spp.scrs)),
            size = 3) +
  scale_x_continuous(limits = c(-0.75,.5)) +
  coord_equal()

###### LDA ######
Whole_Blood_2 <- Whole_Blood %>%
  dplyr::select(Sex, Cu, As, Zn, Hg, Pb) %>%
  dplyr::filter(complete.cases(.))

Whole_Blood_2[,-1] <- scale(log10(Whole_Blood_2[,-1] + 1))

boxplot(Whole_Blood_2[,-1])

Whole_Blood.pca <- prcomp(~. -Sex, data = Whole_Blood_2, center = FALSE, scale. = FALSE)
Whole_Blood.pca
summary(Whole_Blood.pca)

fviz_eig(Whole_Blood.pca, addlabels = TRUE, ylim = c(0, 60)) +
  geom_hline(yintercept = (1/(ncol(Whole_Blood_2)-1))*100)

Whole_Blood.pca$rotation
cutoff <- sqrt(1/ncol(Whole_Blood_2[,-1]))
foo <- data.frame(Whole_Blood.pca$rotation[,1:3]) %>%
  dplyr::mutate_all( function(x){replace(x, abs(x) < cutoff, NA)} )
row.names(foo) <- row.names(Whole_Blood.pca$rotation)
foo

p1 <- ggbiplot(Whole_Blood.pca, groups = Whole_Blood_2$Sex)
p1

Whole_Blood_3 <- Whole_Blood_2[,-1] + abs(min(Whole_Blood_2[,-1]))


Whole_Blood.mds <- metaMDS(comm = Whole_Blood_3, distance = "euclid")
Whole_Blood.mds$stress

MDS_xy <- data.frame(Whole_Blood_2$Sex, Whole_Blood.mds$points)
MDS_xy$Sex <- Whole_Blood_2$Sex

stressplot(Whole_Blood.mds)

vf <- envfit(Whole_Blood.mds, Whole_Blood_3, perm = 999)

spp.scrs <- data.frame(scores(vf, display = "vectors"))

p2 <- ggplot(data = MDS_xy) +
  geom_point(mapping = aes(x = MDS1, y = MDS2, colour = Sex)) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1*.75, y = 0, yend = NMDS2*.75),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = spp.scrs, aes(x = NMDS1*.75, y = NMDS2*.75, label = rownames(spp.scrs)),
            size = 3) +
  scale_x_continuous(limits = c(-0.75,.5)) +
  coord_equal()
p2

ggarrange(p1,p2, nrow = 2, labels = c("PCA", "NMDS"))

my.lda <- lda(x = Whole_Blood_2[,-1], grouping  = Whole_Blood_2$Sex)
my.lda <- lda(Whole_Blood_2$Sex ~ as.matrix(Whole_Blood_2[,-1]))
my.lda

pred <- predict(my.lda)
names(pred)

struct.mat <-data.frame(LD1 = cor(Whole_Blood_2[,-1],pred$x[,1]),
                        LD2 = cor(Whole_Blood_2[,-1],pred$x[,2]))
struct.mat

d <- data.frame(type = Whole_Blood_2$Sex, LD1 = pred$x[,1], LD2 = pred$x[,2])

p3 <- ggplot(data = d, aes(x=LD1, y=LD2, col=type)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("LD1") + 
  ylab("LD2") +
  stat_density2d(bins = 15) +
  geom_point(alpha = .6) + 
  geom_text(data=struct.mat*9, aes(x=LD1, y=LD2, label=row.names(struct.mat)), size = 5, color="black") +
  geom_segment(data=struct.mat*9, aes(x=0, y=0, xend=LD1, yend=LD2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  coord_equal(1.1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3

###### PCA vs NMDS vs LDA ######
ggarrange(p1, p2, p3, ncol = 3, labels = c("PCA", "NMDS", "LDA"))
#### Inter-Elemental Relationships in Bone ####
# first to examine each element and remove outliers
# Cu
hist(Epiphysis$Cu)
q <- quantile(Epiphysis$Cu, probs = c(.25, .75), na.rm = TRUE)
IQR <- IQR(Epiphysis$Cu)
Lower <- q[1] - 1.5*IQR
Upper <- q[2] + 1.5*IQR
Epiphysis_2 <- subset(Epiphysis, Epiphysis$Cu > Lower & Epiphysis$Cu < Upper)
hist(Epiphysis_2$Cu)
shapiro.test(Epiphysis_2$Cu)
# need to normalize
hist((Epiphysis_2$Cu)^(1/3))
shapiro.test((Epiphysis_2$Cu)^(1/3)) # normal p-value = 0.6554
Epiphysis_3 <- Epiphysis_2 %>%
  mutate(Cu = (Cu^(1/3)))

# Zn
hist(Epiphysis$Zn)
q <- quantile(Epiphysis$Zn, probs = c(.25, .75), na.rm = TRUE)
IQR <- IQR(Epiphysis$Zn)
Lower <- q[1] - 1.5*IQR
Upper <- q[2] + 1.5*IQR
Epiphysis_2 <- subset(Epiphysis_2, Epiphysis_2$Zn > Lower & Epiphysis_2$Zn < Upper)
hist(Epiphysis_2$Zn)
shapiro.test(Epiphysis_2$Zn) # normal p-value = 0.4164

# As
hist(Epiphysis$As)
q <- quantile(Epiphysis$As, probs = c(.25, .75), na.rm = TRUE)
IQR <- IQR(Epiphysis$As)
Lower <- q[1] - 1.5*IQR
Upper <- q[2] + 1.5*IQR
Epiphysis_2 <- subset(Epiphysis_2, Epiphysis_2$As > Lower & Epiphysis_2$As < Upper)
hist(Epiphysis_2$As)
shapiro.test(Epiphysis_2$As) # normal p-value = 0.1346

# Hg
hist(Epiphysis$Hg)
q <- quantile(Epiphysis$Hg, probs = c(.25, .75), na.rm = TRUE)
IQR <- IQR(Epiphysis$Hg)
Lower <- q[1] - 1.5*IQR
Upper <- q[2] + 1.5*IQR
Epiphysis_2 <- subset(Epiphysis_2, Epiphysis_2$Hg > Lower & Epiphysis_2$Hg < Upper)
hist(Epiphysis_2$Hg)
shapiro.test(Epiphysis_2$Hg) # normal p-value = 0.6302

# Pb
hist(Epiphysis$Pb)
q <- quantile(Epiphysis$Pb, probs = c(.25, .75), na.rm = TRUE)
IQR <- IQR(Epiphysis$Pb)
Lower <- q[1] - 1.5*IQR
Upper <- q[2] + 1.5*IQR
Epiphysis_2 <- subset(Epiphysis_2, Epiphysis_2$Pb > Lower & Epiphysis_2$Pb < Upper)
hist(Epiphysis_2$Pb)
shapiro.test(Epiphysis_2$Pb) # normal p-value = 0.6997

# create a dataframe with only numeric values so chart.Correlation only plots
# whats necessary
Epiphysis_No_Outliers <- Epiphysis_2
Epiphysis_No_Outliers_2 <- Epiphysis_No_Outliers %>%
  dplyr::select(Cu, Zn, As, Hg, Pb)

# chart.Correlation for non-transformed no outliers epiphysis
chart.Correlation(Epiphysis_No_Outliers_2)
# shows correlation between Zn/Hg, Zn/Pb, and As/Pb

# create a dataframe with transformed data, no outliers
Epiphysis_No_Outliers_Transformed <- Epiphysis_2 %>%
  mutate(Cu = (Cu^(1/3)))
Epiphysis_No_Outliers_Transformed_2 <- Epiphysis_No_Outliers_Transformed %>%
  dplyr::select(Cu, Zn, As, Hg, Pb)
colnames(Epiphysis_No_Outliers_Transformed_2) = c('Cu^(1/3)', 'Zn', 'As', 'Hg', 'Pb')

# chart.Correlation
chart.Correlation(Epiphysis_No_Outliers_Transformed_2)
# shows correlation between Zn/Hg, Zn/Pb, and As/Pb
#### Mass-Wing length regression ####

# mass wing residuals
mod=lm(Mass ~ Wing,data=Correlated_Pb_Hg,na.action=na.exclude)
Correlated_Pb_Hg$condition_resid=resid(mod)

## mass/wing
Correlated_Pb_Hg$SMI=Correlated_Pb_Hg$Mass/Correlated_Pb_Hg$Wing

# Pb vs Condition residual
hist(Correlated_Pb_Hg$Pb) # potential outliers
Q <- quantile(Correlated_Pb_Hg$Pb, probs = c(.25, .75), na.rm = TRUE)
iqr <- IQR(Correlated_Pb_Hg$Pb, na.rm = TRUE)
Lower <- Q[1] - 1.5*iqr
Upper <- Q[2] + 1.5*iqr
Correlated_Pb_No_Outliers <- subset(Correlated_Pb_Hg, Correlated_Pb_Hg$Pb > Lower & Correlated_Pb_Hg$Pb < Upper)
hist(Correlated_Pb_No_Outliers$Pb)
shapiro.test(Correlated_Pb_No_Outliers$Pb) # non-normal p-value = 2.9188e-11
hist(log(Correlated_Pb_No_Outliers$Pb))
shapiro.test(log(Correlated_Pb_No_Outliers$Pb)) # non-normal p-value = 0.006016
hist(sqrt(Correlated_Pb_No_Outliers$Pb))
shapiro.test(sqrt(Correlated_Pb_No_Outliers$Pb)) # can't force normal so must use non-parametric tests

kruskal.test(data = Correlated_Pb_No_Outliers, condition_resid ~ Pb)
# non-sig p-value = 0.4512 so no sig difference in Pb by condition_resid

# Sex vs Condition residual
hist(Correlated_Pb_No_Outliers$condition_resid)
shapiro.test(Correlated_Pb_No_Outliers$condition_resid) # non-normal p-value = 5.489e-12
hist(log(Correlated_Pb_No_Outliers$condition_resid))
shapiro.test(log(Correlated_Pb_No_Outliers$condition_resid)) # normal p-value = 0.8193
Correlated_Pb_No_Outliers_2 <- Correlated_Pb_No_Outliers %>%
  mutate(condition_resid = log(condition_resid))

leveneTest(data = Correlated_Pb_No_Outliers_2, condition_resid ~ Sex)
# non-sig p-value = 0.4851 so can use ANOVA
my.aov_Resid_Sex <- aov(data = Correlated_Pb_No_Outliers_2, condition_resid ~ Sex)
summary(my.aov_Resid_Sex)
# sig p-value = 0.024 so sig difference in condition_resid by Sex

group_by(Correlated_Pb_No_Outliers_2, Sex) %>%
  summarise(
    count = n(),
    mean = mean(condition_resid, na.rm = TRUE),
    IQR = IQR(condition_resid, na.rm = TRUE)
  )
# log mean condition_resid F = 1.31, M = 0.693, U = 0.953
# mean condition_resid F = 3.706, M = 1.1, U = 2.593
# females are significantly fatter than males and unknowns
#### What else do hormones vary with? ####
###### Estradiol vs Metadata ######
# E2 vs Fat
hist(RBCs$E2)
shapiro.test(RBCs$E2)
RBCs_AK <- RBCs %>%
  filter(State == "AK")
hist(RBCs_AK$E2)
shapiro.test(RBCs_AK$E2) # normal p-value = 0.3928
leveneTest(data = RBCs_AK, E2 ~ Fat)
# non-sig p-value = 0.4663 so can use ANOVA
my.aov_E2_Fat <- aov(data = RBCs_AK, E2 ~ Fat)
summary(my.aov_E2_Fat)
# non-sig p-value = 0.53 so no sig difference in E2 by fat class

# E2 vs Age
# doesn't work because all birds are AHY

# E2 vs Sex
# can't do this bc E2 is all female birds

# E2 vs Fat
hist(RBCs$E2)
shapiro.test(RBCs$E2) # non-normal p-value = 0.0001
kruskal.test(data = RBCs, E2 ~ Fat)
# non-sig p-value = 0.1882 so no sig difference in E2 by fat class
#### Graphical Abstract ####

ggarrange(p1, p2, p3, ncol = 2, nrow = 2)
