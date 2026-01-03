library(tidyverse)
library(jsonlite)
library(prospectr)
library(pls)

Block1 = read_json(path = "Implementation/Data/Algae/new mixture/1block.json")
Block2 = read_json(path = "Implementation/Data/Algae/new mixture/2block.json")
Block3 = read_json(path = "Implementation/Data/Algae/new mixture/4block.json")

Test1 = read_json(path = "Implementation/Data/Algae/new mixture/1test.json")
Test2 = read_json(path = "Implementation/Data/Algae/new mixture/2test.json")
Test3 = read_json(path = "Implementation/Data/Algae/new mixture/3test.json")
Test4 = read_json(path = "Implementation/Data/Algae/new mixture/4test.json")


getSpectra = function(x){
  compData = x$tag %>% str_split("-") %>% .[[1]]
  Spectrum = x$spectrum %>% unlist() %>% t() %>% as_tibble()
  colnames(Spectrum) = 1:512 
  compData = tibble(
    "Scenedesmus" = compData[1] %>% as.numeric(),
    "Dolichospermum" = compData[2] %>% as.numeric(),
    "PBS" = compData[3] %>% as.numeric(),
    "LED" = compData[4]
    )
  out = bind_cols(compData,Spectrum)
}
## train sets
B1 = map(Block1,getSpectra)
B1 = B1 %>% bind_rows()
B2 = map(Block2,getSpectra)
B2 = B2 %>% bind_rows()
B3 = map(Block3,getSpectra)
B3 = B3 %>% bind_rows()
TrainData = bind_rows(B1,B2,B3)
rm(B1,B2,B3)

Train = list()
Train[["WA"]] = TrainData %>% filter(LED == "WA")
Train[["365"]] = TrainData %>% filter(LED == "365")
Train[["395"]] = TrainData %>% filter(LED == "395")
Train[["440"]] = TrainData %>% filter(LED == "440")
Train[["470"]] = TrainData %>% filter(LED == "470")
Train[["530"]] = TrainData %>% filter(LED == "530")
Train[["600"]] = TrainData %>% filter(LED == "600")
Train[["660"]] = TrainData %>% filter(LED == "660")

## Test sets
T1 = map(Test1,getSpectra)
T1 = T1 %>% bind_rows()
T2 = map(Test2,getSpectra)
T2 = T2 %>% bind_rows()
T3 = map(Test3,getSpectra)
T3 = T3 %>% bind_rows()
T4 = map(Test4,getSpectra)
T4 = T4 %>% bind_rows()
TestData = bind_rows(T1,T2,T3,T4)
rm(T1,T2,T3,T4)
Test = list()
Test[["WA"]] = TestData %>% filter(LED == "WA")
Test[["365"]] = TestData %>% filter(LED == "365")
Test[["395"]] = TestData %>% filter(LED == "395")
Test[["440"]] = TestData %>% filter(LED == "440")
Test[["470"]] = TestData %>% filter(LED == "470")
Test[["530"]] = TestData %>% filter(LED == "530")
Test[["600"]] = TestData %>% filter(LED == "600")
Test[["660"]] = TestData %>% filter(LED == "660")

## PRocessing and exploration

Train[["WAT"]] = bind_cols(Train[["WA"]][1:4],Train[["WA"]][5:516] %>% map(function(x){log(1/x)}) %>% bind_cols()) 

#El primer punto fue tomado con otra celda, se elimina.

plot(1:512,Train[["WAT"]][2,5:516] %>% unlist,type = "l",ylim = c(-10,2),col = "red")
map(3:56,function(x){lines(1:512,Train[["WAT"]][x,5:516] %>% unlist)})

Train[["WAT"]][5:400] %>% prcomp(center = T,scale. = F) %>% biplot()

# posibles outliers
Train[["WAT"]][,5:512] %>% prcomp(center = T,scale. = F) %>% biplot(arrow.len = 0,var.axes = F,ylabs = NULL)

plot(1:512,Train[["WAT"]][57,5:516] %>% unlist,type = "l",ylim = c(-10,2),col = "red")
map(c(2:56,58),function(x){lines(1:512,Train[["WAT"]][x,5:516] %>% unlist)})
#57 Salió mucho más baja la intensidad (log(1/x) más grande) y más ruidosa la última parte

plot(1:512,Train[["WA"]][2,5:516] %>% unlist,type = "l",ylim = c(0,40000),col = "red")
map(3:58,function(x){lines(1:512,Train[["WA"]][x,5:516] %>% unlist)})

plot(1:512,Train[["WA"]][38,5:516] %>% unlist,type = "l",ylim = c(0,40000),col = "red")
plot(1:512,Train[["WA"]][57,5:516] %>% unlist,type = "l",ylim = c(0,40000),col = "red")

plot(1:512,Train[["WA"]][2,5:516] %>% unlist,type = "l",ylim = c(0,40000),col = "black")
map(c(2:37,39:56,58),function(x){lines(1:512,Train[["WA"]][x,5:516] %>% unlist,col = if(Train[["WA"]][x,3] == 1){"green"} else {"black"})})

Train[["WA"]] %>% filter(PBS == 1) %>% .[,5:516] %>% rowSums() %>% which.min()
#El último espectro de los blancos que da raro. 
#Tiene la misma forma que todos los otros blancos, pero mucha menos intensidad.

plot(1:512,Train[["WA"]][2,5:516] %>% unlist,type = "l",ylim = c(0,40000),col = "black")
map(c(2:37,39:56),function(x){lines(1:512,Train[["WA"]][x,5:516] %>% unlist,col = if(Train[["WA"]][x,3] == 1){"green"} else {"black"})})

plot(1:512,Train[["WAT"]][2,5:516] %>% unlist,type = "l",ylim = c(-11,-3),col = "black")
map(c(2:37,39:56),function(x){lines(1:512,Train[["WAT"]][x,5:516] %>% unlist,col = if(Train[["WA"]][x,3] == 1){"green"} else {"black"})})
# remove outliers
Train[["WAclean"]] = Train[["WA"]][c(2:37,39:56),]
Train[["395clean"]] = Train[["395"]][c(2:37,39:56),]
Train[["440clean"]] = Train[["440"]][c(2:37,39:56),]
Train[["470clean"]] = Train[["470"]][c(2:37,39:56),]
Train[["530clean"]] = Train[["530"]][c(2:37,39:56),]
Train[["600clean"]] = Train[["600"]][c(2:37,39:56),]
Train[["660clean"]] = Train[["660"]][c(2:37,39:56),]

#absorbance respect to max intensity spectra
log10(Train[["WA"]][49,5:516]/Train[["WA"]][2,5:516]) %>% unlist() %>% plot(type = "l",ylim = c(-1,1))
map(c(3:37,39:56),
    function(x){log10(Train[["WA"]][49,5:516]/Train[["WA"]][x,5:516]) %>% 
        unlist() %>% 
        lines(type = "l",col = if(Train[["WA"]][x,3] == 1){"green"} else {"black"})})

#Absorbance respect to mean blank spectra
meanBlank = Train[["WA"]] %>%
  filter(PBS == 1) %>%
  .[5:516] %>% 
  colMeans() %>% unname()
log10(meanBlank / Train[["WA"]][2,5:516] %>% unlist) %>% plot(type = "l",ylim = c(-0.1,1))
map(c(3:37,39:56),
    function(x){log10(Train[["WA"]] %>%
                        filter(PBS == 1) %>%
                        .[5:516] %>%
                        colMeans() %>%
                        unname() / Train[["WA"]][x,5:516]) %>% 
        unlist() %>% 
        lines(type = "l",col = if(Train[["WA"]][x,3] == 1){"green"} else {"black"})})

## Generating white absorbance data
Train[["WAabs"]] = bind_cols(Train[["WAclean"]][1:4],
                           1:54 %>% 
                             map(
                               function(x,MOB){
                                 log10(MOB/Train[["WAclean"]][x,5:516])
                                 }
                               ,MOB = meanBlank) %>% bind_rows()) 

### Exploring the test set
plot(1:512,Test[["WA"]][2,5:516] %>% unlist,type = "l",ylim = c(0,40000),col = "red")
map(3:58,function(x){lines(1:512,Test[["WA"]][x,5:516] %>% unlist)})
Test[["WA"]][,5:516] %>% prcomp(center = T,scale. = F) %>% biplot

## No signs of outliers
Test[["WAabs"]] = bind_cols(Test[["WA"]][1:4],
                            1:nrow(Test[["WA"]])  %>% 
                               map(
                                 function(x,MOB){
                                   log10(MOB/Test[["WA"]][x,5:516])
                                 }
                                 ,MOB = meanBlank) %>% bind_rows()) 

plot(1:512,Test[["WAabs"]][2,5:516] %>% unlist,type = "l",ylim = c(0,1),col = "red")
map(3:58,function(x){lines(1:512,Test[["WAabs"]][x,5:516] %>% unlist,col = "red")})
map(3:58,function(x){lines(1:512,Train[["WAabs"]][x,5:516] %>% unlist)})

Train[["WAT"]][-57,5:512] %>% prcomp(center = T,scale. = F) %>% biplot(arrow.len = 0,var.axes = F,ylabs = NULL)
plot(1:512,Train[["WAT"]][38,5:516] %>% unlist,type = "l",ylim = c(-10,2),col = "red")
lines(1:512,Train[["WAT"]][39,5:516] %>% unlist)
lines(1:512,Train[["WAT"]][37,5:516] %>% unlist)
#38 es outlier, está anotado en el cuaderno. La señal salió garrapateada, se remueve

PCA = list()
PCA[["WA"]] = Train[["WAT"]][c(-58,-57,-38,-1),5:400] %>% scale(center = T,scale = F) %>% prcomp(center = T,scale. = F) 
PCA[["WA"]] %>% biplot(arrow.len = 0,var.axes = F,ylabs = NULL)
PCA[["WA"]] %>% plot
# scores first component
PCA[["WA"]]$rotation[,1] %>% unlist() %>% plot()
# scores second component 
PCA[["WA"]]$rotation[,2] %>% unlist() %>% plot()
# scores third component 
PCA[["WA"]]$rotation[,3] %>% unlist() %>% plot()
# scores fourth component 
PCA[["WA"]]$rotation[,4] %>% unlist() %>% plot()
# scores fifth component 
PCA[["WA"]]$rotation[,5] %>% unlist() %>% plot()
PCA[["395"]] = Train[["395clean"]][5:516] %>% scale(center = T,scale = F) %>% prcomp() %>% biplot() 
PCA[["440"]] = Train[["440clean"]][5:516] %>% scale(center = T,scale = F) %>% prcomp() %>% biplot() 
PCA[["470"]] = Train[["470clean"]][5:516] %>% scale(center = T,scale = F) %>% prcomp() %>% biplot() 

plsr()
baseline = baseline(Train[["WAT"]][c(-57,-38),5:516],1:512)
plot(1:512,baseline[1,] %>% unlist,type = "l",ylim = c(0,2),col = "red")
map(c(2:56),function(x){lines(1:512,baseline[x,] %>% unlist)})


compData = Block1[[49]]$tag %>% str_split("-") %>% .[[1]]
Spectrum = Block1[[1]]$spectrum %>% unlist() %>% as_tibble() %>% t()%>% as_tibble()


plsr(Dolichospermum~ ., data=Train[["WAabs"]][c(2,5:516)],ncomp = 8) %>% RMSEP() %>% plot()


derivative = savitzkyGolay(X = Train[["WAabs"]][5:516],m = 1,p = 3,w = 9) %>% as_tibble()
colnames(derivative) = paste("d",5:504)
union = bind_cols(Train[["WAabs"]][10:512],derivative)

manyThings = bind_cols(Train[["WAabs"]][5:516],
                       #Train[["395clean"]][5:180],
                       #Train[["440clean"]][5:230])
                       Train[["470clean"]][5:255])
                       #Train[["530clean"]][100:300],
                       #Train[["600clean"]][305:410],
                       #Train[["660clean"]][400:516])
par(mfrow = c(2,1))

baseline = baseline(Train[["WAabs"]][,5:516],1:512) %>% as_tibble()
#plot without baseline
plot(1:512,baseline[1,] %>% unlist,type = "l",ylim = c(0,0.4),col = "black")
map(c(2:37,39:56),
    function(x){
      lines(1:512,
            baseline[x,] %>%
              unlist,
            col = if(Train[["WAabs"]][x,3] == 1){
              "green"
              } else if (Train[["WAabs"]][x,2] == 1) {
                "orange"
              } else if (Train[["WAabs"]][x,1] == 1){
                "yellow"
              } else {
                "black"
              }
            )
      })
#plot with baseline
plot(1:512,Train[["WAabs"]][1,5:516] %>% unlist,type = "l",ylim = c(0,0.8),col = "black")
map(c(2:37,39:56),
    function(x){
      lines(1:512,
            Train[["WAabs"]][x,5:516] %>%
              unlist,
            col = if(Train[["WAabs"]][x,3] == 1){
              "green"
            } else if (Train[["WAabs"]][x,2] == 1) {
              "orange"
            } else if (Train[["WAabs"]][x,1] == 1){
              "yellow"
            } else {
              "black"
            }
      )
    })
PLSmodels = list()

PLSmodels[["absExplore"]] = plsr(Train[["WAabs"]]$Scenedesmus~ ., data=Train[["WAabs"]][5:430] ,validation = "LOO",ncomp = 20) 
PLSmodels[["DebasedExplore"]] = plsr(Train[["WAabs"]]$Scenedesmus~ ., data=baseline[,1:425] ,validation = "LOO",ncomp = 20) 

PLSmodels[["absExplore"]] %>% plot("validation")
PLSmodels[["DebasedExplore"]] %>% plot("validation")

PLSmodels[["absExplore"]] = plsr(Train[["WAabs"]]$Scenedesmus~ ., data=Train[["WAabs"]][5:516] ,validation = "LOO",ncomp = 20) 

plsr(Train[["WAabs"]]$Scenedesmus~ ., data=manyThings,ncomp = 5) %>% plot()

plsr(Train[["WAabs"]]$Scenedesmus~ ., data= ,validation = "LOO",ncomp = 20) %>% plot("validation")
plsr(Train[["WAabs"]]$Scenedesmus~ ., data=manyThings,ncomp = 5) %>% plot()


PLSmodels[["absDef"]] = plsr(Train[["WAabs"]]$Scenedesmus~ ., data=Train[["WAabs"]][5:430],ncomp = 6) 
tuning = tibble("real" = Test[["WAabs"]]$Scenedesmus[0:20],"predicted" = predict(PLSmodels[["absDef"]],newdata = Test[["WAabs"]][0:20,])[,,6])
plot(real ~ predicted,tuning)
mod2 = lm(real ~ predicted,data = tuning)

newPred = tibble("predicted" = predict(PLSmodels[["absDef"]],newdata = Test[["WAabs"]][21:56,])[,,6])
plot(Test[["WAabs"]]$Scenedesmus[21:56] ~ predict(mod2,newdata = newPred)) 

par(mfrow = c(2,1))
plot(predicted ~ real,tuning)
plot(predict(mod2,newdata = newPred) ~ Test[["WAabs"]]$Scenedesmus[21:56]) 
par(mfrow = c(1,1))

plsr(Train[["WAabs"]]$Scenedesmus~ ., data=manyThings,validation = "LOO",ncomp = 20) %>% plot("validation")

plsr(Train[["WAabs"]]$Scenedesmus~ ., data=manyThings,ncomp = 5) %>% plot()

## 395
plsr(Train[["WAabs"]]$Scenedesmus~ ., data=Train[["395clean"]][5:180],ncomp = 20,validation = "LOO") %>% plot("validation")
plsr(Train[["WAabs"]]$Scenedesmus~ ., data=Train[["395clean"]][5:180],ncomp = 3) %>% plot()
plsr(Train[["WAabs"]]$Dolichospermum~ ., data=Train[["395clean"]][5:180],ncomp = 20,validation = "LOO") %>% plot("validation")
plsr(Train[["WAabs"]]$Dolichospermum~ ., data=Train[["395clean"]][5:180],ncomp = 3) %>% plot()

## 440
plsr(Train[["WAabs"]]$Scenedesmus~ ., data=Train[["440clean"]][5:230],ncomp = 20,validation = "LOO") %>% plot("validation")
plsr(Train[["WAabs"]]$Scenedesmus~ ., data=Train[["440clean"]][5:230],ncomp = 4) %>% plot()
plsr(Train[["WAabs"]]$Dolichospermum~ ., data=Train[["440clean"]][5:230],ncomp = 20,validation = "LOO") %>% plot("validation")
plsr(Train[["WAabs"]]$Dolichospermum~ ., data=Train[["440clean"]][5:230],ncomp = 8) %>% plot()

## 470
plsr(Train[["WAabs"]]$Scenedesmus~ ., data=Train[["470clean"]][5:255],ncomp = 20,validation = "LOO") %>% plot("validation")
plsr(Train[["WAabs"]]$Scenedesmus~ ., data=Train[["470clean"]][5:255],ncomp = 3) %>% plot()
plsr(Train[["WAabs"]]$Dolichospermum~ ., data=Train[["470clean"]][5:255],ncomp = 20,validation = "LOO") %>% plot("validation")
plsr(Train[["WAabs"]]$Dolichospermum~ ., data=Train[["470clean"]][5:255],ncomp = 3) %>% plot()
## 530
plsr(Train[["WAabs"]]$Scenedesmus~ ., data=Train[["530clean"]][100:300],ncomp = 20,validation = "LOO") %>% plot("validation")
plsr(Train[["WAabs"]]$Scenedesmus~ ., data=Train[["530clean"]][100:300],ncomp = 4) %>% plot()
plsr(Train[["WAabs"]]$Dolichospermum~ ., data=Train[["530clean"]][100:300],ncomp = 20,validation = "LOO") %>% plot("validation")
plsr(Train[["WAabs"]]$Dolichospermum~ ., data=Train[["530clean"]][100:300],ncomp = 4) %>% plot()
## 600
plsr(Train[["WAabs"]]$Dolichospermum~ ., data=Train[["600clean"]][305:410],ncomp = 20,validation = "LOO") %>% plot("validation")
plsr(Train[["WAabs"]]$Dolichospermum~ ., data=Train[["600clean"]][305:410],ncomp = 3) %>% plot()
## 660
plsr(Train[["WAabs"]]$Scenedesmus~ ., data=Train[["660clean"]][400:516],ncomp = 20,validation = "LOO") %>% plot("validation")
plsr(Train[["WAabs"]]$Scenedesmus~ ., data=Train[["660clean"]][400:516],ncomp = 2) %>% plot()
plsr(Train[["WAabs"]]$Dolichospermum~ ., data=Train[["660clean"]][400:516],ncomp = 20,validation = "LOO") %>% plot("validation")
plsr(Train[["WAabs"]]$Dolichospermum~ ., data=Train[["660clean"]][400:516],ncomp = 2) %>% plot()


manyThings = bind_cols(Train[["WAabs"]][5:516],
                       #Train[["395clean"]][5:180])
                       Train[["440clean"]][5:230])
                       #Train[["470clean"]][5:255])
#Train[["530clean"]][100:300],
#Train[["600clean"]][305:410],
#Train[["660clean"]][400:516])

### union WAabs + 395

plsr(Train[["WAabs"]]$Scenedesmus~ ., data=manyThings[5:688],ncomp = 20,validation = "LOO") %>% plot("validation")
plsr(Train[["WAabs"]]$Scenedesmus~ ., data=manyThings[5:688],ncomp = 4) %>% plot()
plsr(Train[["WAabs"]]$Dolichospermum~ ., data=manyThings[5:688],ncomp = 20,validation = "LOO") %>% plot("validation")
plsr(Train[["WAabs"]]$Dolichospermum~ ., data=manyThings[5:688],ncomp = 8) %>% plot()

vector = c(
  Train[["WAabs"]][1,5:516] %>% unlist(),
  Train[["440clean"]][1,5:516] %>% unlist(),
  Train[["470clean"]][1,5:516] %>% unlist(),
  Train[["530clean"]][1,5:516] %>% unlist(),
  Train[["600clean"]][1,5:516] %>% unlist(),
  Train[["660clean"]][1,5:516] %>% unlist()
  ) %>% matrix(nrow = 1)

outerLong = outer(t(vector),vector) %>%
  as_tibble() %>% 
  mutate("row" = 1:3072) %>% 
  pivot_longer(starts_with("V"),names_to = "col") %>%
  mutate(col = col %>% str_remove("V") %>% as.numeric())

ggplot(outerLong,aes(row,col,fill = value)) + geom_raster()
