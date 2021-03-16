#Loading the libraries. 

library(reshape2)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(markdown)

#Reading in transposed platereader data as csv. 

rawdata <- read.csv("../input/20210305_JEC21RNAimut_8DAY_30_TD_TRSP.csv")

#Tidying the data, changing time from s to hr, melt function to reshape the data to tidy format. 

rawdata_hrs <- mutate(rawdata, Time = Time/3600)  
reshaped <- melt(rawdata_hrs, id = c("Time", "Temp"), variable.name = "Well", value.name = "OD595")

#Loading the platemap to show where the strains lie in the wells. 

platemap <- read.csv("../input/20210305_JEC21RNAimut_8Day_30_SetupCSV.csv")

#Combining these two datasets by well. 

annotated <- inner_join(reshaped, platemap, by = "Well")

#Writing the annotated data set to a CSV file. 

#Replace this path with your path of choice. 

write.csv(annotated, "/Users/s2124188/Documents/Physical Rotations/Wallace & Bayne/Plate Reader Data/Run3_8_Day_Cultures.csv")

#Plotting all ODs un-normalised, by Biorep.

ggplot(data=annotated,aes(x=Time, y=OD595, color=Strain, group = Well)) +
  geom_line() +
  theme_bw() +
  facet_grid(BioRep ~ .) +
  labs(x = "Time (hrs)", y="Absorbance at 595 nm")

# Checking the stability of the blank well ODs, for use in normalisation. 

ggplot(data=filter(annotated,Strain==""),
       aes(x = Time, y = OD595, color = Strain, group = Well)) +
  geom_line() +
  theme_bw() +
  labs(x = "Time (hrs)", y = "Absorbance at 595nm")

# Checking OD of blank wells for removal of contaminated (B1)

annotated %>%
  group_by(Well) %>% 
  filter(Strain=="") %>% 
  summarise(maxOD = max(OD595)) %>% 
  arrange(desc(maxOD))

#Filtering annotated dataset for B1 contaminated blank well. 

annotated <- annotated %>% 
  filter(Well != "B1")

#Calculating the median OD and other staistics for blank wells for normalisation. 

blank_OD_summary <- annotated %>% 
  filter(Strain=="") %>% 
  group_by(Medium) %>% 
  summarise(OD_median = median(OD595),
            OD_mean = mean(OD595),
            OD_max = max(OD595),
            OD_min = min(OD595))

#Normalising the data and plotting the normalised ODs. 

normalisedOD <- annotated %>% 
  left_join(blank_OD_summary, by = "Medium") %>% 
  mutate(OD_corrected = OD595 - OD_median)

ggplot(data = filter(normalisedOD,Strain !=""),
       aes(x = Time, y = OD_corrected, color = Strain, group = Well)) +
  geom_line() +
  facet_grid(Medium ~ .) +
  theme_bw() +
  labs(x = "Time (hrs)", y = "Absorbance 595nm")       

#Calculating mean and other statistics across TRs for each BR. 

#With the whole dataset. 
normalisedOD_stats <- normalisedOD %>% 
  group_by(Time, BioRep, Strain, Condition) %>% 
  summarise(OD_mean=mean(OD595)) %>% 
  ungroup()



normalisedOD_stats_2 <- normalisedOD %>% 
  group_by(Well, Time, BioRep, Strain, Condition) %>% 
  summarise(OD_median=median(OD595),
            OD_mean=mean(OD595),
            OD_max=max(OD595),
            OD_min=min(OD595),
            OD_SD=sd(OD595))

#Convert BioRep to a 'character' as the ineger is confusing ggplot.

condition_names <- c('H' = "Hypoxia", 'N' = "Normoxia")


#Adjusted line graph plotting the mean across technical replicates (without the SD)

ggplot(data=filter(normalisedOD_stats, Strain!=""),
       aes(x = Time, y = OD_mean, color = Strain, group = ID)) +
  geom_line(size = 1) +
  facet_grid(Condition ~ ., labeller = as_labeller(condition_names)) +
  theme_bw() +
  labs(x = "Time (hrs)", y = "Absorbance at 595nm")

ggplot(data=normalisedOD_stats %>%
         filter(Strain != ""), 
       aes(x=Time, y=OD_mean, color=Strain, group = ID )) + 
  geom_line() + 
  facet_wrap(Condition ~ .) +
  theme_bw() +
  scale_y_continuous(trans = 'log10')+
  labs(x = "Time(Hrs)") +
  labs(y = "Absorbance 595nm") +
  labs(title = "Growth in YNB") +
  theme(text = element_text(size = 20))+
  scale_color_brewer(palette="Set1")

