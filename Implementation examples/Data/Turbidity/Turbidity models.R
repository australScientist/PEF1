library(tidyverse)
library(jsonlite)
library(pls)


Turbidity = read_json("Implementation/Data/Turbidity/TURBIDITYModified.json")

TurbiSpectra = map(Turbidity,function(x){x %>% .[["spectrum"]] %>% unlist()}) %>% bind_cols() %>% t() %>% as_tibble()
ID = map(Turbidity,function(x){x %>% .[["tag"]] %>% strsplit("-") %>% unlist()}) %>% bind_cols() %>% t() %>% as_tibble()
SiD = map(Turbidity,function(x){x %>% .[["ID"]] %>% unlist()}) %>% bind_cols() %>% t() %>% as_tibble()

colnames(TurbiSpectra) = 1:512
colnames(ID) = c("NTU","LED")
TurbiSpectra = bind_cols(TurbiSpectra,ID) 
TurbiSpectra$ID = SiD$V1

MergedTurbiSpectra %>% filter(LED == "280") %>% ggplot(aes(Pixel,Value,colour = as.factor(ID))) + scale_y_log10() + geom_point()  + facet_wrap(~NTU)

table(TurbiSpectra$LED[-(164:169)],TurbiSpectra$NTU[-(164:169)])
CleanTurbiSpectra = TurbiSpectra[-(164:169),]
CleanTurbiSpectra = CleanTurbiSpectra %>% filter(LED != "280")
CleanTurbiSpectra[,513:515] %>% View
CleanTurbiSpectra %>% nrow()/8


colnames(WASpectra)[1:512] = paste0("p",1:512)

###
WASpectra = filter(CleanTurbiSpectra, LED == "WA")
WAFilter = !(WASpectra$ID %in% c(82,118,136,127,429,258))

filterSpectra = function(led,Coso = CleanTurbiSpectra,filtro = WAFilter){
  print(led)
  Sp = filter(Coso, LED == led)
  print(nrow(Sp))
  Sp = Sp[filtro,]
  print(nrow((Sp)))
  return(Sp)
}

CleanTurbiSpectra = map(unique(CleanTurbiSpectra$LED),filterSpectra) %>% bind_rows()
CleanTurbiSpectra$NTU = as.numeric(CleanTurbiSpectra$NTU)
###

MergedTurbiSpectra = gather(CleanTurbiSpectra,"Pixel","Value",1:512)
MergedTurbiSpectra$NTU = MergedTurbiSpectra$NTU %>% as.numeric()
MergedTurbiSpectra$Pixel = MergedTurbiSpectra$Pixel %>% as.numeric()
MergedTurbiSpectra$Value = MergedTurbiSpectra$Value %>% as.numeric()
MergedTurbiSpectra$LED = MergedTurbiSpectra$LED %>% as.factor()

MergedTurbiSpectra %>% ggplot(aes(Pixel,Value,colour = as.factor(NTU))) +
  geom_line() +
  scale_y_log10() + 
  facet_wrap(~ LED)

MergedTurbiSpectra %>% filter(LED == "660",NTU <601 &NTU > 1) %>% 
  ggplot(aes(Pixel,Value,colour = as.factor(NTU))) + 
  scale_y_log10() +
  geom_point()  

CleanTurbiSpectra %>% filter(LED == "660",NTU <601 & NTU > -1 ) %>% 
  ggplot(aes(NTU,`435`)) + geom_point() #+ scale_y_log10()
#### MODELOS

##simple Models
# LED 395 nm
LM395 = lm(NTU~ poly(`100`,4) ,
              data = CleanTurbiSpectra %>% filter(LED == "395",
                                                  NTU < 601) %>% .[-22,])
summary(LM395)
LM395 %>% plot()
model.matrix(LM395)
plot(LM395$fitted.values,CleanTurbiSpectra %>% filter(LED == "395",
                                                         NTU < 601) %>% .[513] %>% unlist())

LM440 = lm(NTU~`150` + I(`150`^2) + I(`150`^3) ,
           data = CleanTurbiSpectra %>% filter(LED == "440",
                                               NTU < 601) %>% .[-22,])

summary(LM440)
LM440 %>% plot()
model.matrix(LM440)
#plot(LM395$fitted.values,CleanTurbiSpectra %>% filter(LED == "440",

LM470 = lm(NTU~`180` + I(`180`^2) ,
           data = CleanTurbiSpectra %>% filter(LED == "470",
                                               NTU < 601) %>% .[-22,])

summary(LM470)
LM470 %>% plot()

LM530 = lm(NTU~`250` + I(`250`^2) ,
           data = CleanTurbiSpectra %>% filter(LED == "530",
                                               NTU < 601) %>% .[-22,])

summary(LM530)
LM530 %>% plot()


AIC(LM395,LM440,LM470,LM530)
## Bagging
(predict(LM395) + predict(LM440) + predict(LM470) + predict(LM530))/4


## Limit of Quantification 
LOQ = LM395 %>% 
CleanTurbiSpectra %>% 
  filter(LED == "395",
         NTU < 601) %>% 
  ggplot(aes(NTU,`100`)) + 
  geom_point() + 
  geom_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3)+ I(x^4))

PCA395 = CleanTurbiSpectra %>% 
  filter(LED == "395",NTU < 601) %>% 
  .[,1:512] %>% 
  prcomp(scale = TRUE)


PCA395 %>% summary()
PCA395$scale %>% View() #%>% ncol()
PCA395
Reconstruct = PCA395$rotation[,1:3] %*% t(PCA395$x[,1:3])
Reconstruct = Reconstruct[1:110,] %>% t() %>% as.tibble() 
Reconstruct$NTU = CleanTurbiSpectra %>% 
  filter(LED == "395",NTU < 601) %>% 
  .[,513] %>% unlist()
GatheredRec = Reconstruct %>% gather("pixel","value",1:110) %>% as.tibble()
GatheredRec$pixel = GatheredRec$pixel %>% as.numeric()
GatheredRec$value = GatheredRec$value %>% as.numeric()

GatheredRec %>% ggplot(aes(pixel,value,colour = NTU)) + geom_point()
plot(Reconstruct$`90`,Reconstruct$NTU)

SImplePCR = lm(NTU ~ `90` + I(`90`^2)+ I(`90`^3)+ I(`90`^4),Reconstruct)
summary(SImplePCR)


PLSR = plsr(NTU~ ., 
            data = CleanTurbiSpectra %>% 
              filter(LED == "395",NTU < 601) %>% .[c(50:200,513)],
            validation = "LOO")
PLSR %>% summary()
plot(RMSEP(PLSR))

PLSR = plsr(NTU~ ., 
            data = CleanTurbiSpectra %>% 
              filter(LED == "395",NTU < 601) %>% .[c(50:200,513)],
            ncomp = 5)
PLSR %>% summary()
PLSR %>% plot("loadings")
PLSR %>% plot("coef",ncomp = 1:4)

AIC(simpleLM,PLSR)

########  TEST
Test = read_json("Test samples 1.json")
TurbiTest = map(Test,function(x){x %>% .[["spectrum"]] %>% unlist()}) %>% bind_cols() %>% t() %>% as_tibble()
TurbID = map(Test,function(x){x %>% .[["tag"]] %>% strsplit("-") %>% unlist()}) %>% bind_cols() %>% t() %>% as_tibble()

colnames(TurbiTest) = 1:512
TurbiTest$ID = TurbID$V2 
TurbiTest$NTU = TurbID$V1

TurbiTest %>% filter(ID == 395)

Prediction = tibble(
  D395 = TurbiTest %>% filter(ID == 395) %>% .$NTU %>% as.numeric() - predict(LM395,TurbiTest %>% filter(ID == 395))*1.425967,
  D440 = TurbiTest %>% filter(ID == 440) %>% .$NTU %>% as.numeric() - predict(LM395,TurbiTest %>% filter(ID == 440)),
  D470 = TurbiTest %>% filter(ID == 470) %>% .$NTU %>% as.numeric() - predict(LM440,TurbiTest %>% filter(ID == 470)),
  D530 = TurbiTest %>% filter(ID == 530) %>% .$NTU %>% as.numeric() - predict(LM470,TurbiTest %>% filter(ID == 530))
)

lm()

differences = TurbiTest %>% filter(ID == 395) %>% .$NTU %>% as.numeric() - predict(LM530,TurbiTest %>% filter(ID == 530))
Tur

predict(PLSR,CleanTurbiSpectra %>% 
          filter(LED == "395",NTU < 601) %>% .[c(50:200,513)])

LM395_2 = lm(NTU~ poly(`100`,4) ,
           data = TurbiTest %>% filter(ID == "395"))
LM395$coefficients[2]/LM395_2$coefficients[2]

corrected = lm(mod1 ~ poly(mod2,2),
               tibble("mod1" = predict(LM395,
                                       tibble("100" = seq(1,12000,by =100))),
                      "mod2" = predict(LM395_2,
                                       tibble("100" = seq(1,12000,by =100))
                      )
                  )
               )


