#Loading the libraries. 

library(reshape2)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(markdown)

#Reading in transposed platereader data as csv. 

rawdata <- read.csv("../input/20210309_JEC21RNAi_Mut_YNB_YNBS_37_TD_TRSP.csv")

#Tidying the data, changing time from s to hr, melt function to reshape the data to tidy format. 

rawdata_hrs <- mutate(rawdata, Time = Time/3600)  
reshaped <- melt(rawdata_hrs, id = c("Time", "Temp"), variable.name = "Well", value.name = "OD595")

#Loading the platemap to show where the strains lie in the wells. 

platemap <- read.csv("../input/20210224_JEC21RNAimut_YNB+-Serum_37_SetupCSV.csv")

#Combining these two datasets by well. 

annotated <- inner_join(reshaped, platemap, by = "Well")

#Writing the annotated data set to a CSV file. 

#Adjust this filepath to desired path on your system. 

write.csv(annotated, "/Users/s2124188/Documents/Physical Rotations/Wallace & Bayne/Plate Reader Data/Run_2_YNP_SERUM/Output_Data/20210224_YNB+-Serum_JEC21RNAimut_37_EH_TRSP_Annotated")

#Plotting!

#Plot all ODs un-normalised by medium. 

ggplot(data=annotated,aes(x=Time, y=OD595, color=Strain, group = Well)) +
  geom_line() +
  theme_bw() +
  facet_grid(Medium ~ .) +
  labs(x = "Time (hrs)", y="Absorbance at 595 nm")

#Plot all ODs un-normlaised by strain. 

ggplot(data=annotated,aes(x=Time, y=OD595, color=Medium, group = Well)) +
  geom_line() +
  theme_bw() +
  facet_grid(Strain ~ .) +
  labs(x = "Time (hrs)", y="Absorbance at 595 nm")

#We can now check the OD of specific wells, looking at outliers with high OD across all wells.

annotated %>%
  group_by(Well) %>% 
  summarise(maxOD = max(OD595)) %>% 
  arrange(desc(maxOD))

#Also checking in specific strains or medium by the following, can change the Strain== or Medium== to check different ones. 

annotated %>%
  group_by(Well) %>% 
  filter(Strain=="") %>% 
  filter(Medium=="YNB") %>% 
  summarise(maxOD = max(OD595)) %>% 
  arrange(desc(maxOD))

# Checking the stability of the blank well ODs, for use in normalisation. 

ggplot(data=filter(annotated,Strain==""),
       aes(x = Time, y = OD595, color = Strain, group = Well)) +
  geom_line() +
  theme_bw() +
  facet_grid(Medium ~ .) +
  labs(x = "Time (hrs)", y = "Absorbance at 600nm")

#Checking a particular blank well if it looks weird. 

ggplot(data = filter(annotated,Well=="A6"),
       aes(x=Time, y=OD595, color=Strain, group = Well)) + 
  geom_line() + 
  theme_bw()+
  facet_grid(Medium ~ .) +
  labs(x="Time (hrs)", y="Absorbance at 595 nm")

#Calculating the median OD and other staistics for blank wells for normalisation. 

blank_OD_summary <- annotated %>% 
  filter(Strain=="") %>% 
  group_by(Medium) %>% 
  summarise(OD_median = median(OD595),
            OD_mean = mean(OD595),
            OD_max = max(OD595),
            OD_min = min(OD595))

#As above, but removing a weird blank for normalisation if something looks odd. 

blank_OD_summary <- annotated %>%
  filter(Strain=="") %>%
  filter(Well!= "A6")%>%
  filter(Well!= "A1") %>% 
  group_by(Medium) %>%
  summarise(OD_median=median(OD595),
            OD_mean=mean(OD595),
            OD_max=max(OD595),
            OD_min=min(OD595))

#Removing comtaminated wells from annotated dataset. 

annotated_2 <- annotated %>% 
  filter(Well!= "E7") %>% 
  filter(Well!= "F7")

#Normalising the data and plotting the normalised ODs. 

normalisedOD <- annotated_2 %>% 
  left_join(blank_OD_summary, by = "Medium") %>% 
  mutate(OD_corrected = OD595 - OD_median)

ggplot(data = filter(normalisedOD,Strain !=""),
       aes(x = Time, y = OD_corrected, color = Strain, group = Well)) +
  geom_line() +
  facet_grid(Medium ~ .) +
  theme_bw() +
  labs(x = "Time (hrs)", y = "Absorbance 600nm")

#Below are some more specific plots looking at different aspects of the dataset. 

#Growth in YNB

ggplot(data=normalisedOD %>%
         filter(Strain != "", Strain != "H99", Medium =="YNB"), 
       aes(x=Time, y=OD_corrected, color=Strain, group = Well)) + 
  geom_line() + 
  theme_bw() +
  labs(x = "Time(Hrs)") +
  labs(y = "Absorbance 595nm") +
  labs(title = "Growth in YNB") +
  theme(text = element_text(size = 20))+
  scale_color_brewer(palette="Set1")

#Growth in YNBS by replicate, log scale.

ggplot(data=normalisedOD %>%
         filter(Strain != "", Medium =="YNBS"), 
       aes(x=Time, y=OD_corrected, color=Strain, group = Well)) + 
  geom_line() + 
  facet_wrap(BioRep ~ .)+
  theme_bw() +
  scale_y_continuous(trans = 'log10')+
  labs(x = "Time(Hrs)") +
  labs(y = "Absorbance 595nm") +
  labs(title = "Growth in YNBS") +
  theme(text = element_text(size = 20))+
  scale_color_brewer(palette="Set1")

#Growth in YNB by replicate, log scale.

ggplot(data=normalisedOD %>%
         filter(Strain != "", Medium =="YNB"), 
       aes(x=Time, y=OD_corrected, color=Strain, group = Well)) + 
  geom_line() + 
  facet_wrap(BioRep ~ .)+
  theme_bw() +
  scale_y_continuous(trans = 'log10')+
  labs(x = "Time(Hrs)") +
  labs(y = "Absorbance 595nm") +
  labs(title = "Growth in YNB") +
  theme(text = element_text(size = 20))+
  scale_color_brewer(palette="Set1")

#All strains and media by Bioreplicate. 

ggplot(data=normalisedOD %>%
         filter(Strain != ""),
       aes(Time, OD_corrected, colour = Strain, group = Well)) + 
  geom_line() + 
  geom_hline(aes(yintercept=0.5), linetype="dotted", colour="black") +
  facet_wrap(BioRep ~ .) +
  theme_bw() +
  scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
  labs(x = "Time(Hrs)") +
  labs(y = "Absorbance 595nm") +
  theme(text = element_text(size = 10))+
  scale_color_brewer(palette="Set1")

sd()

