library(tidyverse)
library(reticulate)
library(figpatch)
library(patchwork)



image = imager::load.image("Characterisation/Data/measurement0.jpg")
RawImage = tibble(Red = image[,100,1],Green = image[,100,2],Blue = image[,100,3]) %>% 
  mutate(pixel = 1:length(image[,100,1] %>% unlist())) %>% 
  pivot_longer(cols = 1:3,names_to = "Channel",values_to = "Values") %>% 
  ggplot(aes(pixel,Values,colour = Channel)) +
  geom_line() + 
  scale_colour_manual(values = c("blue","green","red")) +
  theme_classic() +
  ylab("Value (AU)") +
  xlab("Pixel") +
  labs(tag = "c)") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.53),    
    legend.direction = "vertical",   
    legend.box = "vertical",
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

source_python(normalizePath("Design/Software/PEF0/src/ImageProcessor.py"))
Signal = ProcessThis$get_signal("Characterisation/Data/measurement0.dng")

FirstRowPlot = tibble(Signal = Signal[1,] %>% unlist()) %>% 
  mutate(Pixel = 1:2592) %>% 
  ggplot(aes(Pixel,Signal)) + 
  geom_line() +
  theme_classic() +
  xlab("Pixel") + 
  ylab("Signal (AU)") +
  labs(tag = "d)") +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

RowAveragedPlot = tibble(Signal = Signal[1:450,] %>% colMeans() %>% unlist()) %>% 
  mutate(Pixel = 1:2592) %>% 
  ggplot(aes(Pixel,Signal)) + 
  geom_line() +
  theme_classic() +
  xlab("Pixel") + 
  ylab("Signal (AU)") +
  labs(tag = "e)")
  
Summed = Signal %>% rowSums()
additive_CV = function(x){
  lengthMat = length(x)
  CVs = list()
  for(i in 2:lengthMat){
    selected = x[1:i]
    CVs[i-1] = sd(selected)/mean(selected)
  }
  return(CVs %>% unlist())
}

CVPlot = tibble("CV" = additive_CV(Summed),"Row" = 1:(nrow(Signal)-1)) %>%
  ggplot(aes(CV,Row)) + 
  geom_line() +
  theme_classic() +
  labs(tag = "a)") + 
  scale_y_reverse() +
  xlab("Coefficient of variance")


FirstSumm = Signal[1:450,1:2560] %>% colMeans()
FinalSpectrum = tibble(signal = FirstSumm %>% matrix(byrow = T,ncol = 5) %>% rowSums(),pixel = 1:512)

FinalSpectrumPlot = FinalSpectrum %>% 
  ggplot(aes(pixel,signal)) + 
  geom_line() +
  theme_classic() +
  xlab("Pixel") + 
  ylab("Signal (AU)") +
  labs(tag = "f)")

img = fig(path = "Characterisation/Data/measurement0.jpg",aspect.ratio = "free")
img = fig_tag(img,tag = "b)",pos = "topleft")

layout = "ABBBBBCCCC
          ABBBBBDDDD
          ABBBBBEEEE
          ABBBBBFFFF"

CVPlot + img + RawImage + FirstRowPlot + RowAveragedPlot + FinalSpectrumPlot + plot_layout(design = layout)
