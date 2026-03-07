library(tidyverse)
library(vascr)

drugsinEGM1 <- vascr_import(
  "ECIS", "ECIS_260223_MFT_1_CG_drugscreenpaper_rep.abp", 
  "ECIS_260223_MFT_1_CG_drugscreenpaper_rep_RbA.csv", "rep1")




# Wellmaps ----------------------------------------------------------------



drugsinEGM1key <- tribble(
  ~SampleID, ~Row, ~Column, ~Sample, #triplicate treatments
  1, "A", "1 2 3",  "low edaravone",
  2, "B", "1 2 3",  "high edaravone",
  
  3, "A", "4 5 6",  "low dipyridamole",
  4, "B", "4 5 6",  "high dipyridamole",
  
  5, "A", "7 8 9",  "low imatinib",
  6, "B", "7 8 9",  "high imatinib",
  
  7, "A", "10 11 12",  "low ticagrelor",
  8, "B", "10 11 12",  "high ticagrelor",
  
  9, "C", "1 2 3", "low marimastat",
  10, "D", "1 2 3", "high marimastat",
  
  11, "C", "4 5 6", "low butylphthalide",
  12, "D", "4 5 6", "high butylphthalide",
  
  13, "C", "7 8 9", "low fasudil",
  14, "D", "7 8 9", "high fasudil",
  
  15, "C", "10 11 12", "low fingolimod",
  16, "D", "10 11 12", "high fingolimod",
  
  17, "E", "1 2 3", "low SCH79797",
  18, "F", "1 2 3", "high SCH79797",
  
  19, "E", "4 5 6", "low phenothiazine",
  20, "F", "4 5 6", "high phenothiazine",
  
  21, "E", "7 8 9", "high astaxanthin",
  22, "F", "7 8 9", "low astaxanthin",
  
  23, "G", "4 5 6", "low VPA",
  24, "H", "4 5 6", "high VPA",
  
  25, "G", "7 8 9", "low combo",
  26, "H", "7 8 9", "high combo",
  
  27, "G", "10 11 12", "low ibuprofen", # 
  28, "H", "10 11 12", "high ibuprofen", # JH's low conc
  
  29, "G", "1 2 3", "vehicle",
  30, "H", "1 2", "combo vehicle")
  
  

# Plots -------------------------------------------------------------------

drugsinEGM1labeled<- vascr_apply_map(drugsinEGM1, drugsinEGM1key)
  


drugsinEGM1plotdata <- drugsinEGM1labeled %>%
  vascr_subset(unit = "Rb") %>%
  vascr_zero_time(65) %>%
  vascr_resample_time(500) %>%
  vascr_normalise(-2, divide = TRUE) 

# edav
drugsinEGM1plotdata %>% 
  vascr_subset(time=c(-4,24),
  sampleid = c(29, 1,2)) %>% 
  vascr_summarise(level = "experiment") %>% vascr_plot_line() + ylim(0.3,1.1)

# dipyridamole
drugsinEGM1plotdata %>% 
  vascr_subset(time=c(-4,24),
               sampleid = c(29, 3,4)) %>% 
  vascr_summarise(level = "experiment") %>% vascr_plot_line() + ylim(0.3,1.1)

# imatinib
drugsinEGM1plotdata %>% 
  vascr_subset(time=c(-4,24),
               sampleid = c(29, 5,6)) %>% 
  vascr_summarise(level = "experiment") %>% vascr_plot_line() + ylim(0.3,1.1)

# ticagrelor
drugsinEGM1plotdata %>% 
  vascr_subset(time=c(-4,24),
               sampleid = c(29, 7,8)) %>% 
  vascr_summarise(level = "experiment") %>% vascr_plot_line() + ylim(0.3,1.1)

# butylphthalide and marimastat mixed up?
drugsinEGM1plotdata %>% 
  vascr_subset(time=c(-4,24),
               sampleid = c(29, 9,10)) %>% 
  vascr_summarise(level = "experiment") %>% vascr_plot_line() + ylim(0.3,1.1)

# marimastat
drugsinEGM1plotdata %>% 
  vascr_subset(time=c(-4,24),
               sampleid = c(29, 11,12)) %>% 
  vascr_summarise(level = "experiment") %>% vascr_plot_line() + ylim(0.3,1.1)


#fasudil
drugsinEGM1plotdata %>% 
  vascr_subset(time=c(-4,24),
               sampleid = c(29, 13,14)) %>% 
  vascr_summarise(level = "experiment") %>% vascr_plot_line() + ylim(0.3,1.1)

#fingolimod
drugsinEGM1plotdata %>% 
  vascr_subset(time=c(-4,24),
               sampleid = c(29, 15,16)) %>% 
  vascr_summarise(level = "experiment") %>% vascr_plot_line() + ylim(0.3,1.1)

# SCH79797
drugsinEGM1plotdata %>% 
  vascr_subset(time=c(-4,24),
               sampleid = c(29, 17,18)) %>% 
  vascr_summarise(level = "experiment") %>% vascr_plot_line() + ylim(0.3,1.1)

# phenothiazine
drugsinEGM1plotdata %>% 
  vascr_subset(time=c(-4,24),
               sampleid = c(29, 19,20)) %>% 
  vascr_summarise(level = "experiment") %>% vascr_plot_line() + ylim(0.3,1.1)

# astaxanthin
drugsinEGM1plotdata %>% 
  vascr_subset(time=c(-4,24),
               sampleid = c(29, 21,22)) %>% 
  vascr_summarise(level = "experiment") %>% vascr_plot_line() + ylim(0.3,1.1)

# VPA
drugsinEGM1plotdata %>% 
  vascr_subset(time=c(-4,24),
               sampleid = c(29, 23,24)) %>% 
  vascr_summarise(level = "experiment") %>% vascr_plot_line() + ylim(0.3,1.1)


# Combo
drugsinEGM1plotdata %>% 
  vascr_subset(time=c(-4,24),
               sampleid = c(29, 25,26,24,21)) %>% # vehicle has one bad well. But mean sits on top of other vehicle #30
  vascr_summarise(level = "experiment") %>% vascr_plot_line() + ylim(0.7,1.1)

# Ibuprofen
drugsinEGM1plotdata %>% 
  vascr_subset(time=c(-4,24),
               sampleid = c(29, 27,28)) %>% 
  vascr_summarise(level = "experiment") %>% vascr_plot_line() + ylim(0.3,1.1)

