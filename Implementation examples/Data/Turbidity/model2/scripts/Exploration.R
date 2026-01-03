library(jsonlite)
library(tidyverse)
library(pls)
library(caret)

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

calibration = read_spectra("guido/Proyecto_ecofisiometro/Software/calibration/Models/Turbidity/model2/data/Curva.json")
test = read_spectra("guido/Proyecto_ecofisiometro/Software/calibration/Models/Turbidity/model2/data/TEST.json")
test$NTU = test$NTU %>% as.numeric()
colnames(test)[1:512] = paste0("p",1:512) 

calibration2 = calibration
colnames(calibration2)[1:512] = paste0("p",1:512) 
calibration2$NTU = calibration2$NTU %>% as.numeric() 
calibration2 = calibration2 %>% filter( NTU <510, NTU != 100)

plot(unlist(calibration2[30,1:512]))

calibration2 %>% filter(LED == "395") %>% ggplot(aes(NTU,p110)) + geom_point()
mod1 = plsr(NTU~.,data = calibration2 %>% filter(LED == 470,NTU != 100) %>% select(-ID,-LED),validation = "LOO")
plot(RMSEP(mod1),legend="topright")
mod1 = plsr(NTU~.,data = calibration2 %>% filter(LED == 395) %>% select(-ID,-LED),ncomp = 2,validation = "LOO")

plot(mod1,ncomp = 1,asp = 1, line = TRUE)
plot(mod1,plottype = "scores",asp = 1, line = TRUE)
plot(mod1,"loadings",comps = 1:2,legendpos = "topright")
plot(mod1,"coef")

mod2 = plsr(NTU~.,data = calibration2 %>% filter(LED == 395,NTU != 100) %>% select(90:175,NTU),ncomp = 2,validation = "LOO")
plot(mod2,ncomp = 1,asp = 1, line = TRUE)
plot(mod2,plottype = "scores",asp = 1, line = TRUE)
plot(mod2,"loadings",comps = 1:2,legendpos = "topright")
plot(mod2,"coef")
plot(mod2)
predict(mod2,test %>% filter(LED == 395,NTU >20 , NTU <510) %>% select(90:175,NTU))

modUni = lm(NTU ~ p110 + I(p110^2) ,data = calibration2 %>% filter(LED == 395))
modUni %>% plot()
modUni$residuals %>% shapiro.test()
modUni %>% summary()
calibration2 %>% filter(LED == "395") %>% ggplot(aes(NTU,p110)) + geom_point() + 
  geom_smooth(method = "lm",formula = y ~ x+ I(x^2))
predictionUni = predict(modUni,newdata = test %>% filter(LED == 395,NTU <500, NTU > 20))
plot(predictionUni,test %>% filter(LED == 395,NTU <500,NTU > 20) %>% .$NTU)

(predictionUni-(test %>% filter(LED == 395,NTU <500,NTU > 20) %>% .$NTU))/(test %>% filter(LED == 395,NTU <500,NTU > 20) %>% .$NTU)

(sqrt((predictionUni-(test %>% filter(LED == 395,NTU <500,NTU > 20) %>% .$NTU))^2 %>% sum())/length(predictionUni))/mean(test %>% filter(LED == 395,NTU <500,NTU > 20) %>% .$NTU)

tunning = lm(test %>% filter(LED == 395,NTU <500) %>% .$NTU ~ predictionUni)

predict(tunning,tibble(predictionUni))
(sqrt((predict(tunning,tibble(predictionUni))-(test %>% filter(LED == 395,NTU <500) %>% .$NTU))^2 %>% sum())/length(predict(tunning,tibble(predictionUni))))/mean(test %>% filter(LED == 395,NTU <500) %>% .$NTU)



blanco = read_spectra("guido/Proyecto_ecofisiometro/Software/calibration/Models/Turbidity/model2/data/Blanco.json")
colnames(blanco)[1:512] = paste0("p",1:512)

modUni2 = lm(p110~ NTU + I(NTU^2),data = calibration2 %>% filter(LED == 395))
modUni2 %>% plot()
modUni$coefficients
filter(blanco,LED == 395) %>% .$p110 %>% sd()/4
filter(blanco,LED == 395) %>% .$p110 %>% boxplot()
filter(blanco,LED == 395) %>% .$p110 %>% summary()
filter(blanco,LED == 395) %>% .$p110 %>% shapiro.test()



modUni3 = lm(NTU ~ p110 + I(p110^2),data = alles)
modUni3 %>% plot()
summary(modUni3)
alles %>% ggplot(aes(NTU,p110)) + geom_point() + geom_smooth(method = "lm",formula = y ~ x+I(x^2))

Control = caret::trainControl(method = "repeatedcv",repeats = 200,number = 3)
modUni4 = caret::train(NTU ~ p110 + I(p110^2),data = alles,method = "lm",trControl = Control)
modUni4 %>% summary()
modUni4$results[,2]/mean(alles$NTU)
############ NEW
CalDef = read_spectra("guido/Proyecto_ecofisiometro/Software/calibration/Models/Turbidity/model2/data/CurvaDef.json")

CalDef$NTU = CalDef$NTU %>% as.numeric()
Filter395 = CalDef %>% filter(LED == 395)
Filter395$NTU = Filter395$NTU %>% as.numeric()
colnames(Filter395)[1:512] = paste0("p",1:512)
ggplot(Filter395%>% filter(NTU != 480),aes(NTU,p110)) + geom_point() + geom_smooth(method = "lm", formula = y ~ x + I(x^2))
Filter395 = Filter395 %>% filter(NTU != 480)
ModDef1 = lm(NTU ~ p110 + I(p110^2),data = Filter395)
ModDef1 %>% plot()
ModDef1$residuals %>% shapiro.test()

alles = bind_rows(test %>% filter(LED ==395, NTU >20 ,NTU <501),
                  calibration2 %>% filter(LED ==395, NTU >20 ,NTU <510),
                  CalDef %>% Filter395(LED))

Control = caret::trainControl(method = "repeatedcv",repeats = 200,number = 8)
modUni4 = caret::train(NTU ~ p110 + I(p110^2),data = alles,method = "lm",trControl = Control)
modUni4
modUni4$results[,2]/mean(alles$NTU)
modUni4 %>% summary()
modUni4$finalModel %>% plot()

modUni4$finalModel %>% lmtest::gqtest(data = alles) 
modUni4$finalModel %>% lmtest::bpte



ggplot(alles,aes(NTU,p110)) + geom_point() + geom_smooth(method = "lm", formula = y ~ x + I(x^2))

glm(NTU ~ p110 + I(p110^2),family = "poisson",data = alles)
caret::train(NTU ~ p110 + I(p110^2),data = alles,method = "glm",trControl = Control,family = "poisson")


todo = bind_rows(calibration2,CalDef,test)
coso = tibble("signal395"= alles[90:130] %>% rowMeans(),"NTU" = alles$NTU)
models = list()
models[["average395"]] = lm(NTU ~ signal395 + I(signal395^2),coso)
