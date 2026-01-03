library(jsonlite)
library(tidyverse)
library(pls)
library(caret)
library(gridExtra)

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

calibration = read_spectra("Implementation/Data/Turbidity/model2/data/Curva.json")
test = read_spectra("Implementation/Data/Turbidity/model2/data/TEST.json")
calibration2 = read_spectra(path = "Implementation/Data/Turbidity/model2/data/CurvaDef.json")
blanco = read_spectra("Implementation/Data/Turbidity/model2/data/Blanco.json")

colnames(calibration)[1:512] = paste0("p",1:512)
colnames(calibration2)[1:512] = paste0("p",1:512)
colnames(test)[1:512] = paste0("p",1:512)
colnames(blanco)[1:512] = paste0("p",1:512)


AllData = bind_rows(calibration,calibration2,test,blanco)
AllData$NTU = AllData$NTU %>% as.numeric()
AllData = AllData %>% select(-ID)

AllData = AllData %>% filter(NTU <501,NTU >19)
allTidy = AllData %>% group_by(LED) %>% nest()
allTidy$data = allTidy$data %>% map(function (x) {x[-c(1,3,29:38),]})
allTidy$Gathered = allTidy$data %>%map(function(x){gather(x,pixel,value,1:512) %>% mutate(pixel = pixel %>% str_extract("\\d+") %>% as.numeric())})
allTidy$Gathered = allTidy$Gathered %>% map(function(x){x %>% mutate(NTU = NTU %>% as_factor) })

allTidy$spectralPlots = map(allTidy$Gathered,function(x){ggplot(x,aes(pixel,value,colour = NTU)) + geom_point(size = 0.1) + ylab("Intensidad (UA)")})
allTidy$spectralPlots[[1]] = allTidy$spectralPlots[[1]] + ggtitle("Transmitancia")
allTidy$spectralPlots[[2]] = allTidy$spectralPlots[[2]] + ggtitle("Scattering 365nm")
allTidy$spectralPlots[[3]] = allTidy$spectralPlots[[3]] + ggtitle("Scattering 395nm")
allTidy$spectralPlots[[4]] = allTidy$spectralPlots[[4]] + ggtitle("Scattering 440nm")
allTidy$spectralPlots[[5]] = allTidy$spectralPlots[[5]] + ggtitle("Scattering 470nm")
allTidy$spectralPlots[[6]] = allTidy$spectralPlots[[6]] + ggtitle("Scattering 530nm")
allTidy$spectralPlots[[7]] = allTidy$spectralPlots[[7]] + ggtitle("Scattering 600nm")
allTidy$spectralPlots[[8]] = allTidy$spectralPlots[[8]] + ggtitle("Scattering 660nm")
grid.arrange(allTidy$spectralPlots[[1]],
             allTidy$spectralPlots[[2]],
             allTidy$spectralPlots[[3]],
             allTidy$spectralPlots[[4]],
             allTidy$spectralPlots[[5]],
             allTidy$spectralPlots[[6]],
             allTidy$spectralPlots[[7]],
             allTidy$spectralPlots[[8]],
             nrow = 4)

plot(allTidy$data[[1]]$NTU,log10(1/allTidy$data[[1]]$p150))
plot(allTidy$data[[3]]$NTU,allTidy$data[[3]]$p110)
plot(allTidy$data[[3]]$NTU,allTidy$data[[3]]$p110)
plot(allTidy$data[[4]]$NTU,allTidy$data[[4]]$p160)
plot(allTidy$data[[5]]$NTU,allTidy$data[[5]]$p180)
plot(allTidy$data[[6]]$NTU,allTidy$data[[6]]$p260)
plot(allTidy$data[[7]]$NTU,allTidy$data[[7]]$p360)
plot(allTidy$data[[8]]$NTU,allTidy$data[[8]]$p430)

allTidy$PCA = allTidy$data %>% map(function(x){x[1:512] %>% prcomp(scale = F)})
par(mfrow=c(1,1))

allTidy$PCA[[3]]$rotation[,1] %>% plot()
allTidy$PCA[[4]]$rotation[,1] %>% plot(main = "loadings PC1 (99,45%)")
allTidy$PCA[[5]]$rotation[,1] %>% plot()
allTidy$PCA[[6]]$rotation[,1] %>% plot()
allTidy$PCA[[7]]$rotation[,1] %>% plot()
allTidy$PCA[[8]]$rotation[,1] %>% plot()

allTidy$PCA[[3]] %>% biplot(var.axes = F) # 7 outlier posiblemente 2
allTidy$PCA[[4]] %>% biplot(var.axes = F) # posibles outliers 7,2,4,6
allTidy$PCA[[5]] %>% biplot(var.axes = F)
allTidy$PCA[[6]] %>% biplot(var.axes = F)
allTidy$PCA[[7]] %>% biplot(var.axes = F)
allTidy$PCA[[8]] %>% biplot(var.axes = F)

allTidy$data[[4]][c(-7,-2, -6 ,-4),1:512] %>% prcomp(scale = T) %>% biplot()

allTidy$PLSR1 = map(allTidy$data, function(x) {plsr(NTU ~ ., data = x,ncomp = 10, validation = "LOO")})
## PLot NCOMP
par(mfrow=c(2,1))
allTidy$PLSR1[[1]] %>% selectNcomp(method = "randomization",plot = T)
allTidy$PLSR1[[1]] %>% selectNcomp(method = "onesigma",plot = T)


allTidy$PLSR1[[3]] %>% selectNcomp(method = "randomization",plot = T)
allTidy$PLSR1[[3]] %>% selectNcomp(method = "onesigma",plot = T)
allTidy$PLSR1[[3]] %>% plot(plottype = "coef", ncomp=1:3)
plsr(NTU ~ .,data = allTidy$data[[3]] %>% select(NTU,90:150)) %>% plot()

allTidy$PLSR1[[4]] %>% selectNcomp(method = "randomization",plot = T)
allTidy$PLSR1[[4]] %>% selectNcomp(method = "onesigma",plot = T)
allTidy$PLSR1[[4]] %>% plot(plottype = "coef", ncomp=1:3)
plsr(NTU ~ .,data = allTidy$data[[4]] %>% select(NTU,0:100)) %>% plot() 
plsr(NTU ~ ., 
     data = allTidy$data[[4]] %>% 
       filter(NTU != 500,NTU != 480) %>% 
       .[-1,] %>% select(NTU,40:60)) %>% summary()


allTidy$PLSR1[[5]] %>% selectNcomp(method = "randomization",plot = T)
allTidy$PLSR1[[5]] %>% selectNcomp(method = "onesigma",plot = T)

allTidy$PLSR1[[6]] %>% selectNcomp(method = "randomization",plot = T)
allTidy$PLSR1[[6]] %>% selectNcomp(method = "onesigma",plot = T)
allTidy$PLSR1[[6]] %>% plot(plottype = "coef", ncomp=1:3)
plsr(NTU ~ .,data = allTidy$data[[6]] %>% select(NTU,100:200)) %>% plot()

allTidy$PLSR1[[7]] %>% selectNcomp(method = "randomization",plot = T)
allTidy$PLSR1[[7]] %>% selectNcomp(method = "onesigma",plot = T)
allTidy$PLSR1[[7]] %>% plot(plottype = "coef", ncomp=1:3)
plsr(NTU ~ .,data = allTidy$data[[7]] %>% select(NTU,0:100)) %>% plot()

allTidy$PLSR1[[8]] %>% selectNcomp(method = "randomization",plot = T)
allTidy$PLSR1[[8]] %>% selectNcomp(method = "onesigma",plot = T)
allTidy$PLSR1[[8]] %>% plot(plottype = "coef", ncomp=1:3)
plsr(NTU ~ .,data = allTidy$data[[4]] %>% select(NTU,400:500)) %>% plot()


allBound = allTidy$data %>% map(function(x){select(x,-NTU)}) %>% bind_cols()
allBound = allBound[513:4096]
allBound$NTU = allTidy$data[[1]]$NTU
plsr(NTU ~ .,data = allBound) %>% plot(plottype = "coef", ncomp=1:3)
plsr(NTU ~ .,data = allTidy$data[[4]] %>% select(NTU,0:100)) %>% plot()
allTidy$Max = c(160, 210,105,150,60,250,360,430)

TestLM = list()
TestLM[[1]] = lm(NTU ~ p150 , data = allTidy$data[[4]] %>% filter(NTU != 500,NTU != 480) %>% .[-1,])
TestLM[[2]] = lm(NTU ~ p150 + I(p150^2), data = allTidy$data[[4]] %>% filter(NTU != 500,NTU != 480) %>% .[-1,])
TestLM[[3]] = lm(NTU ~ p150 + I(p150^2) + I(p150^3), data = allTidy$data[[4]] %>% filter(NTU != 500,NTU != 480) %>% .[-1,])
TestLM[[4]] = lm(NTU ~ p150 + I(p150^2) + I(p150^3) + I(p150^4), data = allTidy$data[[4]] %>% filter(NTU != 500,NTU != 480) %>% .[-1,])

plot(predict(TestLM[[1]]),rstudent(TestLM[[1]]))
plot(predict(TestLM[[2]]),rstudent(TestLM[[2]]))
plot(predict(TestLM[[3]]),rstudent(TestLM[[3]]))
plot(predict(TestLM[[4]]),rstudent(TestLM[[4]]))
anova(TestLM[[3]],TestLM[[4]])
anova(TestLM[[4]])
summary(TestLM[[4]])

mod1summary = mod1 %>% summary()
mod1summary$sigma
mod1$residuals %>% shapiro.test()

Control = caret::trainControl(method = "repeatedcv",repeats = 200,number = 8)
mod2 = caret::train(NTU ~ p150 + I(p150^2) + I(p150^3) + I(p150^4),
                    data = allTidy$data[[4]] %>% filter(NTU != 500,
                                                        NTU != 480,
                                                        NTU != 420,
                                                        NTU != 60,
                                                        NTU != 110,
                                                        NTU != 200,
                                                        NTU != 350) %>% .[c(-1),],
                    method = "lm",trControl = Control)

mod2 %>% summary()
mod2$finalModel %>% anova()
mod2$results[,2]/mean(allTidy$data[[4]] %>% filter(NTU != 500,NTU != 480,NTU != 420,NTU != 60,NTU != 60) %>% .[-1,] %>% .$NTU)
plot(mod2$finalModel)

mod2$finalModel$residuals %>% hist()
mod2$finalModel$residuals %>% shapiro.test()

ggplot(data = allTidy$data[[4]] %>% filter(NTU != 500,NTU != 480) %>% .[-1,],
       aes(NTU,p150)) + geom_point() + geom_smooth(method = "lm",formula = y ~ x + I(x^2) + I(x^3))


### 530 nm

ggplot(data = allTidy$data[[6]] %>% filter(NTU != 500,NTU != 480) %>% .[-1,],
       aes(NTU,p250)) + geom_point() #+ geom_smooth(method = "lm",formula = y ~ x + I(x^2) + I(x^3))

TestLM2 = list()
TestLM2[[1]] = lm(NTU ~ p250 , data = allTidy$data[[6]] %>% filter(NTU != 500,NTU != 480) %>% .[-1,])
TestLM2[[2]] = lm(NTU ~ p250 + I(p250^2), data = allTidy$data[[6]] %>% filter(NTU != 500,NTU != 480) %>% .[-1,])
TestLM2[[3]] = lm(NTU ~ p250 + I(p250^2) + I(p250^3), data = allTidy$data[[6]] %>% filter(NTU != 500,NTU != 480) %>% .[-1,])
TestLM2[[4]] = lm(NTU ~ p250 + I(p250^2) + I(p250^3) + I(p250^4), data = allTidy$data[[6]] %>% filter(NTU != 500,NTU != 480) %>% .[-1,])

plot(predict(TestLM2[[1]]),rstudent(TestLM2[[1]]))
plot(predict(TestLM2[[2]]),rstudent(TestLM2[[2]]))
plot(predict(TestLM2[[3]]),rstudent(TestLM2[[3]]))
plot(predict(TestLM2[[4]]),rstudent(TestLM2[[4]]))
anova(TestLM[[3]],TestLM[[4]])
anova(TestLM[[2]])
summary(TestLM[[3]])

mod3 = caret::train(NTU ~ p250 + I(p250^2),
                    data = allTidy$data[[6]] %>% filter(NTU != 500,NTU != 480) %>% .[-1,],
                    method = "lm",trControl = Control)

mod3 %>% summary()
mod3$finalModel %>% anova()
mod3$results[,2]/mean(allTidy$data[[6]] %>% filter(NTU != 500,NTU != 480) %>% .[-1,] %>% .$NTU)
plot(mod3$finalModel)

mod2$finalModel$residuals %>% hist()
mod2$finalModel$residuals %>% shapiro.test()
##### FUNCTION PMAE

for (i in 1:1000){
  testIndex = sample(1:32,size = 8)
  train440 = clean440[-testIndex,]
  test440 = clean440[testIndex,]
  model440 = lm(NTU ~ p150 + I(p150^2) + I(p150^3) + I(p150^4),
                data = clean440)
  prediction = predict(model440,newdata = test440)
  #plot(test440$NTU,prediction)
  pmae[i] = abs(prediction-test440$NTU) %>% sum()/length(prediction)/mean(test440$NTU)
}
pma


### Dato 440

clean440 =allTidy$data[[4]][-c(1,3,29:38),]
pmae = list() %>% as.vector()
for (i in 1:100){
  testIndex = sample(1:32,size = 8)
  train440 = clean440[-testIndex,]
  test440 = clean440[testIndex,]
  model440 = lm(NTU ~ p150 + I(p150^2) + I(p150^3) + I(p150^4),
                data = clean440)
  prediction = predict(model440,newdata = test440)
  #plot(test440$NTU,prediction)
  pmae[i] = abs(prediction-test440$NTU) %>% sum()/length(prediction)/mean(test440$NTU)
}
pmae = pmae %>% unlist()
pmae %>% densityplot()
pmae[pmae<0.12] %>% mean()
ggplot(clean440,aes(NTU,p150)) + geom_point()
Control = caret::trainControl(method = "repeatedcv",repeats = 200,number = 8)
MODELS = list()

MODELS[[1]] = caret::train(NTU ~ p150 + I(p150^2) + I(p150^3) + I(p150^4),
                    data = clean440,
                    method = "lm",trControl = Control)
MODELS[[1]]$times


### Dato 530

clean530 =allTidy$data[[6]][-c(1,3,29:38),]
ggplot(clean530,aes(NTU,p260)) + geom_point()

model530 = lm(NTU ~ p260 + I(p260^2) + I(p260^3)+ I(p260^4),data = clean530)

Control = caret::trainControl(method = "repeatedcv",repeats = 200,number = 8)
MODELS = list()

MODELS[[2]] = caret::train(NTU ~ p260 + I(p260^2) + I(p260^3)+ I(p260^4),data = clean530,
                           method = "lm",trControl = Control)
## Dato 600
clean600 =allTidy$data[[7]][-c(1,3,29:38),]
ggplot(clean600,aes(NTU,p360)) + geom_point()
lm(clean600$NTU ~ clean600$p360 + I(clean600$p360^2) + I(clean600$p360^3))  %>% plot()
MODELS[[3]] = caret::train(NTU ~ p360 + I(p360^2),data = clean600,
                           method = "lm",trControl = Control)

## Dato 600
clean660 =allTidy$data[[8]][-c(1,3,29:38),]
ggplot(clean660,aes(NTU,p430)) + geom_point()

MODELS[[4]] = caret::train(NTU ~ p430,data = clean660,
                           method = "lm",trControl = Control)

WA = allTidy %>% filter(LED == "WA") %>% select(data) %>% unnest() %>% ungroup() %>% select(-LED)
