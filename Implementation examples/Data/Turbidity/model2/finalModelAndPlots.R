library(jsonlite)
library(tidyverse)
library(pls)
library(caret)
library(patchwork)

# colour palettes ---------------------------------------------------------


palette1 = c("#264653",
             "#2A9D8F",
             "#E9C46A",
             "#F4A261",
             "#E76F51")

palette2 = c("#FCECC9",
             "#FCB0B3",
             "#F93943",
             "#7EB2DD",
             "#445E93")

palette3 = c("#3A3335",
             "#D81E5B",
             "#F0544F",
             "#FDF0D5",
             "#C6D8D3")

palette4 = c("#156064",
             "#00C49A",
             "#F8E16C",
             "#FFC2B4",
             "#FB8F67")


# user functions ----------------------------------------------------------


read_spectra = function(path){
  File = read_json(path)
  
  TurbiSpectra = map(File,function(x){x %>% .[["spectrum"]] %>% unlist()}) %>% bind_cols() %>% t() %>% as_tibble()
  ID = map(File,function(x){x %>% .[["tag"]] %>% strsplit("-") %>% unlist()}) %>% bind_cols() %>% t() %>% as_tibble()
  SiD = map(File,function(x){x %>% .[["ID"]] %>% unlist()}) %>% bind_cols() %>% t() %>% as_tibble()
  
  colnames(TurbiSpectra) = 1:512
  colnames(ID) = c("NTU","LED")
  TurbiSpectra = bind_cols(TurbiSpectra,ID) 
  TurbiSpectra$ID = SiD$V1  
  return(TurbiSpectra)
}
read_noise = function(path){
  File = read_json(path)
  
  TurbiSpectra = map(File,function(x){x %>% .[["noise"]] %>% unlist()}) %>% bind_cols() %>% t() %>% as_tibble()
  ID = map(File,function(x){x %>% .[["tag"]] %>% strsplit("-") %>% unlist()}) %>% bind_cols() %>% t() %>% as_tibble()
  SiD = map(File,function(x){x %>% .[["ID"]] %>% unlist()}) %>% bind_cols() %>% t() %>% as_tibble()
  
  colnames(TurbiSpectra) = 1:512
  colnames(ID) = c("NTU","LED")
  TurbiSpectra = bind_cols(TurbiSpectra,ID) 
  TurbiSpectra$ID = SiD$V1  
  return(TurbiSpectra)
}
### Noise
noise = read_noise("Implementation/Data/Turbidity/model2/data/CurvaDef.json")
noise = noise %>% filter(LED == "WA")
colnames(noise)[1:512] = paste0("p",1:512)
noise = noise %>% mutate_at(vars(starts_with("p")),function(x){-log(x)})
noise %>% select(starts_with("p")) %>% select(1:500) %>% prcomp() %>% biplot()
sdNoise = noise %>% select(starts_with("p")) %>%  .[1,] %>% sd()


# load signals ------------------------------------------------------------


calibration = read_spectra("Implementation/Data/Turbidity/model2/data/Curva.json")
test = read_spectra("Implementation/Data/Turbidity/model2/data/TEST.json")
calibration2 = read_spectra(path = "Implementation/Data/Turbidity/model2/data/CurvaDef.json")
blanco = read_spectra("Implementation/Data/Turbidity/model2/data/Blanco.json")

colnames(calibration)[1:512] = paste0("p",1:512)
colnames(calibration2)[1:512] = paste0("p",1:512)
colnames(test)[1:512] = paste0("p",1:512)
colnames(blanco)[1:512] = paste0("p",1:512)

blanco %>% filter(LED == "WA") %>% select(starts_with("p")) %>% prcomp(scale. = TRUE) %>% biplot()

AllData = bind_rows(calibration,calibration2,test,blanco)
AllData$NTU = AllData$NTU %>% as.numeric()

WA = AllData %>% filter(LED == "WA")

WA = WA %>% select(-LED,-ID)


# exploration -------------------------------------------------------------


logWa = WA %>% mutate_at(vars(starts_with("p")),function(x){-log(x)}) %>% filter(NTU >= 10 )
ggplot(logWa,aes(NTU,p150)) + geom_point()
notOut = ((!(logWa$p150 > -8 & logWa$NTU == 100)) & (!(logWa$p150 < -6 & logWa$NTU > 500)) & (logWa$NTU != 60)) %>% which() #Select Outliers
logWa = logWa[notOut,] # Remove outliers


# Transform to absorbance -------------------------------------------------


ggplot(logWa,aes(NTU,p150)) + geom_point() #check again


blanco %>% filter(LED == "WA") %>% select(starts_with("p")) %>% rowSums() %>% which.max()
invert = blanco %>% filter(LED == "WA") %>% 
  select(starts_with("p")) %>%
  .[9,] %>%
  unlist() %>%
  diag(ncol = 512) %>% 
  solve()
WAVals = WA %>% filter(NTU >= 10 ) %>% select(starts_with("p")) %>% as.matrix() %>% .[notOut,]
abs = -log10(WAVals %*% invert)


# explore data multidim ---------------------------------------------------

PCA = abs[,170:420] %>% prcomp()
Loadings = PCA$rotation %>% as_tibble() %>% mutate(pixel = 170:420)
Loadings_long = Loadings %>% select(c(1:6,49)) %>% 
  pivot_longer(starts_with("PC"),
               names_to = "Loading",
               values_to = "Value")
Loadings_long %>% ggplot(aes(pixel,Value)) + geom_line() + facet_wrap(~Loading)

# Absorbance spectra plot -------------------------------------------------

turbiAbs = abs[,170:420] %>% 
  as_tibble() %>% 
  rename_with(.cols = everything(),.fn = function(x){paste0("V",170:420)}) %>% 
  mutate(NTU = logWa$NTU,
         ID = 1:nrow(abs)) %>% 
  pivot_longer(starts_with("V"),
               names_to = "Pixel",
               values_to = "abs") %>% 
  mutate(Pixel = Pixel %>% 
           str_remove("V") %>% 
           as.numeric()
  ) %>% 
  ggplot(aes(Pixel,abs,group = ID,colour = NTU)) + 
  geom_line() +
  scale_colour_gradient(low = palette1[1],
                        high = palette1[3],
                        name = "Turbidity\n(NTU)",
                        breaks = c(100,300,500)) +
  theme_classic() +
  scale_x_continuous(breaks = c(200,300,400)) +
  ylab("Absorbance (AU)") +
  theme(legend.position = "none") 
# Score plot --------------------------------------------------------------

Scores = PCA$x %>% as_tibble() %>% 
  mutate(NTU =  logWa$NTU)
turbiScore = Scores %>% 
  ggplot(aes(PC1,PC2,colour = NTU)) +
  scale_colour_gradient(low = palette1[1],
                        high = palette1[3],
                        name = "Turbidity\n(NTU)",
                        breaks = c(100,300,500))+
  geom_point(size = 3) + 
  theme_classic() + 
  scale_x_continuous(breaks = c(-20,0,20)) +
  scale_y_continuous(breaks = c(-3,0,3)) +
  theme(legend.position = "none")


# Model selection ---------------------------------------------------------
set.seed(18)
train_ctrl = trainControl(method="repeatedcv", # type of resampling in this case Cross-Validated
                          number=10,
                          repeats = 100, # number of folds
                          search = "grid" # we are performing a "grid-search"
)



modelNTU = train(NTU ~ .,
                  data = abs[,170:420] %>% 
                   as_tibble() %>% 
                   mutate(NTU = logWa$NTU),
                  method = "pls",
                  tuneLength = 10,
                  trControl = train_ctrl
)
#modelNTU %>% plot()
#write_rds(modelNTU,file = "Implementation/Data/Turbidity/model2/output/FinalModel")
# Coefficient plot --------------------------------------------------------
coefData = tibble(C1 = modelNTU$finalModel$coefficients %>% .[1:length(170:420),,1],
                  C2 = modelNTU$finalModel$coefficients %>% .[1:length(170:420),,2],
                  C3 = modelNTU$finalModel$coefficients %>% .[1:length(170:420),,3],
                  pixel = 170:420
) %>% 
  pivot_longer(starts_with("C"),names_to = "coef",values_to = "value")

CoefPlot = ggplot(coefData,aes(pixel,value,colour = coef)) +
  scale_colour_manual(values = c(5,2,4),name = "Regression\nCoefficient")+
  geom_line(linewidth = 0.7) +
  theme_classic() +
  geom_hline(yintercept = 0,linetype = "dotted",colour = palette1[2]) + 
  scale_x_continuous(breaks = c(200,300,400)) +
  # scale_y_continuous(breaks = c(0,3,6)) +
  ylab("Regression Coefficient") +
  xlab("Pixel") +
  theme(legend.position = "bottom")

# Prediction Plot ---------------------------------------------------------

predVsReal = tibble(predicted = predict(modelNTU),
                    real = logWa$NTU)
predVsRealPlot = predVsReal %>% ggplot(aes(predicted,real,colour = real)) +
  scale_colour_gradient(low = palette1[1],
                        high = palette1[3],
                        name = "Turbidity\n(NTU)",
                        breaks = c(100,300,500))+
  geom_point(size = 2) + 
  geom_abline(slope = 1,intercept = 0,linetype = "dotted") +
  theme_classic() +
  xlab("Predicted Turbidity (NTU)") +
  ylab("Real Turbidity (NTU)") +
  scale_y_continuous(breaks = c(0,300,600)) +
  scale_x_continuous(breaks = c(0,300,600)) +
  theme(legend.position = "bottom")
ErrorDens = modelNTU$resample[201:300,] %>% 
  ggplot(aes(MAE)) + 
  geom_density(fill = palette1[4],alpha = 0.7) +
  theme_classic() +
  geom_vline(xintercept = modelNTU$results[3,4],linetype = "dashed") +
  geom_label(data = tibble(x = modelNTU$results[3,4] %>% unlist(),
                           y = 0.3,
                           lab = paste("MAE:",modelNTU$results[3,4] %>% unlist() %>% round(2)
                                       )
                           ),
             aes(x = x, y = y, label = lab))+
  ylab(label = "Density") +
  ylim(c(0,0.350))
Prediction = predVsRealPlot #+ inset_element(ErrorDens,left = 0.01,bottom = 0.5,right = 0.5,top =1 )


layout = "AAAB
          CCCD"
TurbidityPanel = turbiAbs + turbiScore + CoefPlot + Prediction + plot_layout(design = layout)
