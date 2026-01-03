library(tidyverse)
library(reticulate)
library(exiftoolr)

# This script transforms raw images into tsv files. This is done because raw
# data is too large to upload to the github repo. Original files will be made 
# available on demand.

# Original files and directories with the raw data were later compressed and 
#removed

reticulate::source_python("Scripts/analize.py")

# Gain vs Integration Time ------------------------------------------------

extract = function(path,gain){
  all = list.files(path,full.names = T) %>% map(get_spectrum) %>% bind_cols()
  all = t(all)
  all = as_tibble(all)
  colnames(all) = paste0("P",1:512)
  names = list.files(path)
  time = names %>% str_extract(pattern = regex("[:digit:]+(?=.dng)"))
  all$time = time
  all$gain = gain
  return(all)
}
all0.5 = extract(path = "Characterisation/Data/Gain vs Integration Time/0.5/",gain = 0.5)
all1 = extract(path = "Characterisation/Data/Gain vs Integration Time/1/",gain = 1)
all1.5 = extract(path = "Characterisation/Data/Gain vs Integration Time/1.5/",gain = 1.5)

all = bind_rows(all0.5,all1,all1.5)

write_tsv(all,file = "Characterisation/Data/Gain vs Integration Time/results.tsv")


# integration time --------------------------------------------------------

all1 = extract(path = "Characterisation/Data/Integration Time/signalVsTime/",gain = 1)
all2 = extract(path = "Characterisation/Data/Integration Time/SignalVsTime2/",gain = 1)
all3 = extract(path = "Characterisation/Data/Integration Time/SingalVsTime Long/",gain = 1)

all = bind_rows(all1,all2,all3)

write_tsv(all,file = "Characterisation/Data/Integration Time/results.tsv")

# integration time vs PWM -------------------------------------------------
extract2 = function(path){
  all = list.files(path,full.names = T) %>% map(get_spectrum) %>% unlist() %>% matrix(ncol = 512,byrow = T)
  all = as_tibble(all)
  colnames(all) = paste0("P",1:512)
  names = list.files(path)
  time = names %>% str_extract(pattern = regex("[:digit:]+(?=.dng)"))
  pwm = names %>% str_extract(pattern = regex("[:digit:].[:digit:]+(?=.measurement)"))
  all$time = time
  all$gain = pwm
  return(all)
}

all1 = extract2(path = "Characterisation/Data/integration time and PWM/TimeAndPWM2/")
all2 = extract2(path = "Characterisation/Data/integration time and PWM/TimeAndPWM3/")

all = bind_rows(all1,all2)

write_tsv(all,file = "Characterisation/Data/integration time and PWM/results.tsv")

# signal variation in time ------------------------------------------------

all = extract(path = "Characterisation/Data/rep/",gain = 1)
metadata = list.files("Characterisation/Data/rep/repetitionInTIme",full.names = T) %>% 
  map(exif_read) %>% bind_rows()
write_tsv(metadata,file = "Characterisation/Data/rep/metadata.tsv")
all = mutate(all,date = metadata$DateTimeOriginal)
write_tsv(all,file = "Characterisation/Data/rep/results.tsv")


