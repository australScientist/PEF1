library(tidyverse)
library(patchwork)

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
# repetitivity ------------------------------------------------------------

repetitivity = read_tsv(file = "Characterisation/Data/rep/results.tsv")
repetitivity = repetitivity %>% mutate(date = date %>% ymd_hms())
repetitivityLong = pivot_longer(repetitivity,cols = starts_with("P"),names_to = "pixel",values_to = "value") %>% mutate(pixel = str_extract(pixel,"[:digit:]+") %>% as.numeric())

###
allPixels = ggplot() +
  geom_line(data = repetitivityLong,aes(pixel,value,group = date),alpha = 0.05) + 
  theme_classic() + 
  xlim(c(0,550)) + 
  geom_rect(
    data=tibble(
      x = 325,
      y = 28000,
      xend = 350,
      yend = 30000),
    aes(
      xmin = x, 
      ymin = y,
      xmax = xend, 
      ymax = yend),
    colour = "red3",
    alpha = 0.12
    ) + ylab("Signal (AU)")
p325_350 = repetitivityLong %>% 
  filter(pixel > 323 & pixel <351) %>% 
  ggplot(aes(pixel,value,group = date)) + geom_line(alpha = 0.025) + theme_classic() +
  # theme(axis.line.y = element_blank(),
  #       axis.title.y = element_blank(),axis.ticks.y = element_blank())
  scale_y_continuous(breaks = c(26000,29000)) +
  xlim(c(325,380)) + 
  scale_x_continuous(breaks = c(325,338,350)) + 
  geom_vline(xintercept = 338) + ylab("Signal (AU)")
histogram = repetitivityLong %>% 
  filter(pixel == 338) %>% 
  ggplot(aes(y = value)) + geom_histogram(bins = 20) + theme_classic() +
  scale_y_continuous(breaks = c(28667,29152,29543))  + ylab("Signal (AU)") #+ 
  # theme(axis.title = element_blank(),
  #       axis.line = element_blank(),
  #       axis.text = element_blank(),
  #       axis.ticks = element_blank())

allPixels + inset_element(p325_350,left = 0.65,bottom = 0.65,right = 1,top = 1) +
  inset_element(histogram,left = 0.65,bottom = 0,right = 1,top = 0.5)

allPixels + (p325_350/histogram) + plot_layout(widths = c(8,3))


allPixels = ggplot() +
  geom_line(data = repetitivityLong,aes(pixel,value,group = date),alpha = 0.05) + 
  theme_classic() + 
  xlim(c(0,550)) + 
  geom_segment(
    data=tibble(
      x = 338,
      y = 28000,
      xend = 338,
      yend = 30000),
    aes(
      x = x, 
      y = y,
      xend = xend, 
      yend = yend),
    colour = "red3",
    alpha = 1
  )

histogram = repetitivityLong %>% 
  filter(pixel == 338) %>% 
  ggplot(aes(y = value)) + geom_density() + geom_rug() + theme_classic() +
  scale_y_continuous(breaks = c(28667,29152,29543),position = "right") 
# theme(axis.title.y = element_blank(),
#       axis.line.y = element_blank(),
#       axis.text.y = element_blank(),
#       axis.ticks.y = element_blank())

allPixels + inset_element(histogram,left = 0.65,bottom = 0.65,right = 1,top = 1) 
allPixels = ggplot() +
  geom_line(data = repetitivityLong,aes(pixel,value,group = date),alpha = 0.05) + 
  theme_classic() + 
  xlim(c(0,550)) + 
  geom_segment(
    data=tibble(
      x = 160,
      y = 24000,
      xend = 160,
      yend = 28000),
    aes(
      x = x, 
      y = y,
      xend = xend, 
      yend = yend),
    colour = "red3",
    alpha = 1
  )
histogram = repetitivityLong %>% 
  filter(pixel == 160) %>% 
  ggplot(aes(y = value)) + geom_rug() + theme_classic() +
 theme(axis.title.y = element_blank())
#       axis.line.y = element_blank(),
#       axis.text.y = element_blank(),
#       axis.ticks.y = element_blank())

allPixels + inset_element(histogram,left = 0.02,bottom = 0.65,right = 0.4,top = 1) 

VanValen_cv = function(x){ # Van Valen CV is used because high colinearity makes det(cov) = 0
  covariance = cov(x)
  meanVector = colMeans(x)
  traceCov = diag(covariance) %>% sum()
  denominator = (t(meanVector) %*% meanVector) %>% unlist()
  numerator = traceCov
  CV = sqrt(numerator/denominator)
  return(CV)
}
CoefVar = VanValen_cv(repetitivity %>% select(starts_with("P"))) %>% .[1,1] %>% unlist()

PCArep = repetitivity %>% select(starts_with("P")) %>% prcomp() 
PCArep$x[,1] %>% hist()
sd(PCArep$x[,1])
plot(PCArep$x[,1] ~ repetitivity$date,type = "l")

repetitivityLong %>% 
  group_by(date) %>% 
  summarise(value = sum(value)) %>% 
  ggplot(aes(date,value)) +
  geom_line() + theme_minimal() +  ylim(c(0,6500000)) 

PCA = repetitivity[1:512] %>% 
  prcomp(scale. = T)

PCA$x %>% as_tibble() %>% 
  ggplot(aes(PC1,PC2)) +
  geom_point(alpha = 0.6) +
  # geom_path() +
  geom_density2d() + 
  theme_classic()
# signal vs integration time ----------------------------------------------

signalVsTime = read_tsv("Characterisation/Data/Integration Time/results.tsv") %>%
  filter(time != 0)
signalVsTime_long = pivot_longer(
  signalVsTime,
  cols = starts_with("P"),
  names_to = "pixel",
  values_to = "value") %>%  
  mutate(
    pixel = str_extract(pixel,
                        "[:digit:]+") %>% 
      as.numeric()
  )
# signalVsTime_long %>% 
#   ggplot(aes(pixel,value,colour = time,group = time)) + geom_line() + theme_minimal()

SignalTimePlot = signalVsTime %>% ggplot(aes(time/1000,P275)) + 
  geom_line() + 
  theme_minimal() +
  ylab("Intensity Pixel 275 (AU)") +
  xlab("Integration Time (ms)") +
  geom_segment(data = tibble(xstart = c(0,29,129),
                             xend = c(24,124,970),
                             y = 0,
                             regime = c("Linear","Overflow","Saturation")
                             ),
               aes(x = xstart,
                   xend = xend,
                   y = y,
                   yend = y,
                   colour = regime),
               linewidth = 3) + 
  scale_colour_manual(values = palette1[c(1,3,5)])
#signal should be less than 35000 time ~ 25000
signalVsTime_selected = signalVsTime %>% filter(time <= 25000) %>% mutate("rep" = c(rep(1,26),rep(2,26)))
signalVsTime_long = pivot_longer(
  signalVsTime_selected,
  cols = starts_with("P"),
  names_to = "pixel",
  values_to = "value") %>%  
  mutate(
    pixel = str_extract(pixel,
                        "[:digit:]+") %>% 
      as.numeric()
  )
spectraVsTime = signalVsTime_long %>% mutate(rep = paste("Camera",rep),time = time/1000) %>% 
  filter(time < 26) %>% 
  ggplot(aes(pixel,value,colour = time,group = time)) + 
  geom_line() + 
  theme_minimal() +
  facet_wrap(~rep) +
  ylab("Signal (AU)") +
  xlab("Pixel") +
  theme_classic() +
  scale_colour_gradient(name = "Integration\nTime (ms)",
                         low = palette2[5],
                         high =  palette2[1],breaks = c(6,14,22)) +
  scale_x_continuous(breaks = c(0,250,500)) 

pca_Time = signalVsTime_selected %>% select(starts_with("P")) %>% prcomp()

pca_Time_axis = pca_Time$x %>% 
  as_tibble() %>% 
  mutate(
    time = signalVsTime_selected$time,
    rep = signalVsTime_selected$rep %>% as_factor)

# ggplot(pca_Time_axis,aes(x= time, y = PC1,colour = rep)) + geom_point() + geom_smooth(method = "lm",se = F) + theme_minimal()

TimePCAReg = ggplot(pca_Time_axis %>% mutate(time = time/1000),aes(x= time, y = PC1)) + 
  geom_point(colour = palette1[1],size = 2) + 
  geom_line(aes(colour = rep),
            stat = "smooth",
            method = "lm",
            se = F,
            linewidth = 1.6,
            alpha = 0.5) + 
  theme_classic() +
  scale_colour_manual(values = palette1[c(2,4)],name = "Camera") + 
  xlab("Integration Time (ms)") 
layout = "AAAA
          BBCC
          BBCC"
SignalTimePlot + spectraVsTime + TimePCAReg + plot_layout(design = layout)

timeModel = lm(PC1 ~ time*rep,data = pca_Time_axis %>% mutate(time = time/1000))
timeModelSum = timeModel %>% summary()
timeRsq = timeModelSum$adj.r.squared %>% round(5)

# a mixed effect modesummary()# a mixed effect model to extract info about the variability in beta
lme4::lmer(PC1~time + (1+time|rep),pca_Time_axis)


# Gain vs Integration time ------------------------------------------------
gainTime = read_tsv("Characterisation/Data/Gain vs Integration Time/results.tsv")
signalVsTime_long = pivot_longer(
  gainTime,
  cols = starts_with("P"),
  names_to = "pixel",
  values_to = "value") %>%  
  mutate(
    pixel = str_extract(pixel,
                        "[:digit:]+") %>% 
      as.numeric()
  )
signalVsTime_long %>% ggplot(aes(pixel,value,colour = time,group = paste(time,gain))) + geom_line() + theme_minimal()
signalVsTime_long %>% filter(pixel == 50) %>% ggplot(aes(time,value,colour = gain %>% as_factor(),group = paste(gain))) + geom_line()

# time vs PWM -------------------------------------------------------------

pwmTime = read_tsv("Characterisation/Data/integration time and PWM/results.tsv")
signalVsTime_long = pivot_longer(
  pwmTime,
  cols = starts_with("P"),
  names_to = "pixel",
  values_to = "value") %>%  
  mutate(
    pixel = str_extract(pixel,
                        "[:digit:]+") %>% 
      as.numeric()
  )
signalVsTime_long %>% ggplot(aes(pixel,value,colour = time,group = paste(time,gain))) + geom_line() + theme_minimal()
signalVsTime_long %>% filter(pixel == 50) %>% ggplot(aes(time,value,colour = (gain),group = paste(gain,time))) + geom_point()


# Factorial ---------------------------------------------------------------


Results = jsonlite::read_json("Characterisation/Data/factorialExperiment/Lu1Ca1Ci1.json")
Spectra = Results %>% map(function(x){x[["spectrum"]]})
Tags = Results %>% map(function(x){x[["tag"]]})
Results[[1]][["spectrum"]]

Spectra = Spectra %>% map(function(x){x %>% unlist %>% return()}) %>% unlist() %>% matrix(byrow = TRUE,ncol = 512) %>% as_tibble()
colnames(Spectra) = 1:512
Spectra$Tag = Tags %>% map(function(x){x %>% unlist %>% return()}) %>% unlist()
Spectra = separate(Spectra,Tag,into = c("camera","light","percent"),sep = "-")
GatheredSpectra = gather(Spectra,"pixel","Value",1:512)

GatheredSpectra$pixel = as.numeric(GatheredSpectra$pixel)
GatheredSpectra$percent = as.numeric(GatheredSpectra$percent)
GatheredSpectra$camera = as.factor(GatheredSpectra$camera)
GatheredSpectra$light = as.factor(GatheredSpectra$light)
filter(GatheredSpectra,camera == 1, light == 1) %>% ggplot(aes(pixel,Value,colour = percent %>% as.factor())) + geom_line() #+ facet_wrap(~light)

GatheredSpectra %>% ggplot(aes(pixel,Value,colour = percent %>% as.factor())) + geom_line() + facet_grid(rows = vars(light),cols = vars(camera))

### mean spectrum
PercMean = Spectra %>% group_by(percent) %>% select(-camera,-light) %>%  summarise_all(mean)
gather(PercMean,pixel,value,2:513) %>% ggplot(aes(pixel %>% as.numeric(),value %>% as.numeric(),colour = as.factor(percent))) + geom_line()

filter(GatheredSpectra,camera == 1,light == 1,pixel == 300) %>% ggplot(aes(Value,percent)) + geom_line() 
nestedSpectra = Spectra %>% group_by(light,camera) %>% nest()
library(pls)
PLSR = function(x){
  model = pls::plsr(percent ~ ., data = x, ncomp = 1, validation = "LOO")
  return(model)
}
models = nestedSpectra$data %>% map(PLSR)
par(mfrow= c(4,5))
for(i in 1:20){
  plot(RMSEP(models[[i]]), legendpos = "topright")
}

par(mfrow= c(4,5))
for(i in 1:20){
  plot(models[[i]])
}
par(mfrow= c(1,1))
par(mfrow= c(1,1))
coeffs = models %>% map(coef) %>% map(as.vector) %>% bind_cols() %>% t() %>% as_tibble()
colnames(coeffs) = 1:512
coefsPCA = prcomp(coeffs)

colnames(Spectra)[514] = "slit"
vegan::adonis2(Spectra[,1:512]~Spectra$slit+Spectra$percent+Spectra$camera)

models %>% map(summary)




loadings = models %>% map(loadings) %>% map(as.vector) %>% bind_cols() %>% t() %>% as_tibble()
colnames(loadings) = 1:512
loadings$instrument = (1:10) %>% as.factor()
gather(loadings[-9,],pixel,value,1:512) %>% ggplot(aes(pixel %>% as.numeric(),value %>% as.numeric(),colour = instrument)) + geom_line()
averageLoadings = loadings[-9,-513] %>% colMeans()
project = function(spectra,referenceSpectra){
  spectra = spectra %>% unlist()
  referenceSpectra = referenceSpectra %>% unlist()
  normRef = sqrt(referenceSpectra %*% referenceSpectra) %>% as.vector()
  scalarProjection = ((spectra %*% referenceSpectra)/normRef) %>% as.vector()
  vectorProjection = scalarProjection*referenceSpectra/normRef
  return(vectorProjection)
}
project(loadings[2,-513],averageLoadings) %>% plot()
par(mfrow = c(1,3))
project(Spectra[3,1:512],loadings[2,-513]) %>% plot()
Spectra[3,1:512] %>% unlist() %>% plot()
loadings[2,-513] %>% unlist() %>% plot()

projected = apply(Spectra[1:512],1,project,loadings[2,-513]) %>% t() %>% as_tibble()
projected$light = Spectra$light
projected$camera = Spectra$camera
projected$percent = Spectra$percent
filter(GatheredSpectra,camera == 1,light == 1,pixel == 300) %>% ggplot(aes(Value,percent)) + geom_line() 
nestedSpectra = projected %>% group_by(light,camera) %>% nest()
library(pls)
PLSR = function(x){
  model = pls::plsr(percent ~ ., data = x, ncomp = 1, validation = "LOO")
  return(model)
}
models = nestedSpectra$data %>% map(PLSR)
par(mfrow= c(2,5))
for(i in 1:10){
  plot(models[[i]])
}

loadings2 = models %>% map(loadings) %>% map(as.vector) %>% bind_cols() %>% t() %>% as_tibble()
colnames(loadings2) = 1:512
loadings2$instrument = (1:10) %>% as.factor()
par(mfrow=c(1,1))
gather(loadings2[-9,],pixel,value,1:512) %>% ggplot(aes(pixel %>% as.numeric(),value %>% as.numeric(),colour = instrument)) + geom_line()

diag(c(rep(1,100),rep(1.5,100),rep(1,312)))
generate_data = function(along,original,center,width){
  original = original %>% as.matrix()
  cols = ncol(original)
  multiplier = 1:cols %>% map(function(x){exp(-((x-center)/width)^2)*along+1}) %>% as_vector()
  matrixMult = diag(multiplier)
  new = original %*% matrixMult
  new = new %>% as.vector()
  return(new)
}
generate_data(1,original = artifitialData[,1:512],center = 150,width = 20)

artifitial = map(seq(0,1,by=0.1),generate_data,original = artifitialData[,1:512],center = 150,width = 20) %>% bind_cols %>% t() %>% as_tibble
colnames(artifitial) = 1:512
artifitial$artInt = seq(0,1,by=0.1)
projectedArt = apply(artifitial[1:512],1,project,loadings1[2,-513]) %>% t() %>% as_tibble()
projectedArt$artInt = artifitial$artInt
par(mfrow=c(1,2))
gather(artifitial,pixel,value,1:512) %>% ggplot(aes(pixel %>% as.numeric(),value %>% as.numeric(),colour = artInt %>% as.factor())) + geom_line()
gather(projectedArt,pixel,value,1:512) %>% ggplot(aes(pixel %>% as.numeric(),value %>% as.numeric(),colour = artInt %>% as.factor())) + geom_line()

### PDS

PDS<-function(masterSpectra, slaveSpectra, MWsize, Ncomp, wavelength){
  
  require(pls)
  
  #Loop Initialization:
  i<-MWsize
  k<-i-1
  #Creation of an empty P matrix:
  P<-matrix(0,nrow=ncol(masterSpectra),ncol=ncol(masterSpectra)-(2*i)+2)
  InterceptReg<-c()
  
  while(i<=(ncol(masterSpectra)-k)){
    
    #PLS regression:
    fit<- plsr(masterSpectra[,i] ~ as.matrix(slaveSpectra[,(i-k):(i+k)]),
               ncomp=Ncomp, scale=F, method="oscorespls")
    
    #Extraction of the regression coefficients:
    coefReg<-as.numeric(coef(fit, ncomp=Ncomp, intercept = TRUE))
    InterceptReg<-c(InterceptReg,coefReg[1])
    coefReg<-coefReg[2:length(coefReg)]
    
    #Add coefficients to the transfer matrix:
    P[(i-k):(i+k),i-k]<-t(coefReg)
    
    rm(coefReg,fit)
    i<-i+1
    
    #Diplay progression:
    cat("\r",paste(round(i/ncol(masterSpectra)*100)," %",sep=""))}
  
  P<-data.frame(matrix(0,nrow=ncol(masterSpectra),ncol=k), P,
                matrix(0,nrow=ncol(masterSpectra),ncol=k))
  InterceptReg<-c(rep(0,k),InterceptReg,rep(0,k)) 
  
  Output<-list(P = P , Intercept = InterceptReg)
  
  return(Output)}
Master = Spectra %>% filter(camera == 1,light == 1) %>% .[,1:512] %>% as.matrix()
Slave = Spectra %>% filter(camera == 2,light == 4) %>% .[,1:512] %>% as.matrix()
ResPDS = PDS(masterSpectra = Master,slaveSpectra = Slave,wavelength = 1:512,MWsize = 10,Ncomp = 2)
ResPDS[["Intercept"]] %>% plot()
SlaveCor = Slave %*% as.matrix(ResPDS$P)
SlaveCor[4,] %>% as.vector() %>% plot()
SlaveCor = sweep(SlaveCor,2,as.numeric(t(ResPDS$Intercept)),"+")
SlaveCor[4,] %>% as.vector() %>% plot()

## Spectral Space Transformation:

SP1.1 = filter(Spectra,camera ==1, light == 1)
SP2.3 = filter(Spectra,camera == 2, light == 3)
Xcomb = bind_cols(SP2.3[1:512],SP1.1[1:512]) %>% as.matrix()
SVD = svd(Xcomb)
U = SVD$u
D = SVD$d
V = SVD$v
library(MASS)
(as.matrix(SP2.3[1,1:512])) %*% ginv(t(V[513:1024,1])) %*% (t(V[1:512,1]))
p1p2 = ginv(t(V[513:1024,1])) %*% (t(V[1:512,1]))
p2p2 = ginv(t(V[513:1024,1])) %*% (t(V[513:1024,1]))
xtest = as.matrix(SP1.1[21,1:512])
xtrans = xtest %*% p1p2 + xtest - xtest %*% p2p2
par(mfrow=c(3,1))
xtrans %>% as.vector() %>% plot()
xtest %>% as.vector() %>% plot()
as.matrix(SP2.3[21,1:512]) %>% as.vector() %>% plot()
par(mfrow=c(3,1))
xtest = generate_data(1,as.matrix(SP1.1[21,1:512]),center = 220,width = 5)
xtrans = xtest %*% p1p2 + xtest - xtest %*% p2p2
xtrans %>% as.vector() %>% plot()
xtest %>% as.vector() %>% plot()
as.matrix(SP2.3[21,1:512]) %>% as.vector() %>% plot()


minimal = filter(Spectra,camera == 1, light == 1) %>% select(-light,-camera)

model = pls::plsr(percent ~ ., data = minimal, ncomp = 10, validation = "LOO")
model %>% RMSEP() %>% plot
model = pls::plsr(percent ~ ., data = minimal, ncomp = 1, validation = "LOO")
model %>% plot()
model %>% plot(plottype="coef")


oliveoil %>% as.data.frame()
