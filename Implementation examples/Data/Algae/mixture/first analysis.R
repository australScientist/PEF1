library(tidyverse)
library(jsonlite)
library(factoextra)
library(pls)

M1 = read_json("data/mezcla1.json")
M2 = read_json("data/mezcla2.json")
M3 = read_json("data/mezcla3.json")
M4 = read_json("data/mezcla4.json")
M5 = read_json("data/mezcla5.json")

extract_spectra = function(x){
  spectra = x %>% map(function(y){y$spectrum})
  spectra = spectra %>% unlist %>% matrix(byrow = T,ncol = 512) %>% as_tibble()
  colnames(spectra) = paste0("p",1:512)
  return(spectra)
}
extract_data = function(x){
  tags = x %>% map(function(y){y$tag %>% strsplit("-") %>% unlist})
  proportions = tags %>% map(function(y){y[1]}) %>% map(function(x){x %>% strsplit(",")})
  scenedesmus = proportions %>% map(function(x){x%>% unlist() %>% .[1]}) %>% unlist()
  dolichospermum = proportions %>% map(function(x){x%>% unlist() %>% .[2]}) %>% unlist()
  LED = tags %>% map(function(y){y%>% unlist() %>% .[2]}) %>% unlist()
  tags = tibble("Scenedesmus" = scenedesmus,"Dolichospermum" = dolichospermum,"LED" = LED)
  return(tags)
}

SpM1 = extract_spectra(M1)
SpM2 = extract_spectra(M2)
SpM3 = extract_spectra(M3)
SpM4 = extract_spectra(M4)
SpM5 = extract_spectra(M5)
Spectra = bind_rows(SpM1,SpM2,SpM3,SpM4,SpM5)
tag1 = extract_data(M1)
tag2 = extract_data(M2)
tag3 = extract_data(M3)
tag4 = extract_data(M4)
tag5 = extract_data(M5)
Tags = bind_rows(tag1,tag2,tag3,tag4,tag5)
Tags$Scenedesmus = ifelse(Tags$Scenedesmus == 900,90,Tags$Scenedesmus)
Tags$Scenedesmus = ifelse(Tags$Scenedesmus == 25,50,Tags$Scenedesmus)
Tags$Dolichospermum = ifelse(Tags$Dolichospermum == 1000,100,Tags$Dolichospermum)
Tags$Dolichospermum = ifelse(Tags$Dolichospermum == 25,50,Tags$Dolichospermum)
Tags$Scenedesmus = Tags$Scenedesmus %>% as.numeric()
Tags$Dolichospermum = Tags$Dolichospermum %>% as.numeric()

Tags$LED = Tags$LED %>% as_factor()
allData = bind_cols(Spectra,Tags)


gatheredData = allData %>% gather("pixel","value",1:512)
gatheredData$pixel = gatheredData$pixel %>% str_remove("p") %>%as.numeric()
gatheredData = gatheredData %>% mutate(value = ifelse(LED == "WA", log(1/value),value))
gatheredData %>% ggplot(aes(pixel,value,colour = Scenedesmus)) + geom_point(size = 0.25) + facet_wrap(~ LED,scales = "free")

WA = allData %>% filter(LED == "WA")
WA[1:512] = WA[1:512] %>% mutate_all(function(x){log(1/x)}) %>% apply(1,function(x){x/max(abs(x))}) %>% t()

PCAwa = prcomp(WA[WA$Scenedesmus + WA$Dolichospermum != 0,1:512],scale. = T)
PCAwa %>% fviz_pca_biplot(col.ind = WA$Scenedesmus[WA$Scenedesmus + WA$Dolichospermum != 0],geom.var = "",
                          gradient.cols = c("#207557","#bd71c7","#e3f542"))
PCAwa$rotation[,1] %>% plot()

PLSR  =plsr(Scenedesmus ~ .,data = WA %>% select(-Dolichospermum,-LED,200:450),ncomp = 5)
plot(PLSR %>% predict(type = "response") %>% .[,,5]  ~ WA$Scenedesmus)
### Dentro del mismo rango de concentraciones iniciales (pimeros 14 puntos)
PCAwa = prcomp(WA[1:27,1:512],scale. = T)
PCAwa %>% fviz_pca_biplot(col.ind = ScenedesmusConcentration,geom.var = "",
                          gradient.cols = c("#207557","#bd71c7","#e3f542"))
PCAwa$rotation[,1] %>% plot()

PLSR  =plsr(DolichospermumConcentration ~ .,data = WA[1:27,] %>% select(-Dolichospermum,-LED),ncomp = 5)

plot(PLSR %>% predict(type = "response") %>% .[,,5]  ~ DolichospermumConcentration)

### Dentro del mismo rango de concentraciones 2 iniciales (siguientes 13 puntos)
PCAwa = prcomp(WA[15:27,1:512],scale. = T)
PCAwa %>% fviz_pca_biplot(col.ind = WA$Scenedesmus[15:27],geom.var = "",
                          gradient.cols = c("#207557","#bd71c7","#e3f542"))
PCAwa$rotation[,1] %>% plot()

PLSR  =plsr(Scenedesmus ~ .,data = WA[15:27,] %>% select(-Dolichospermum,-LED),ncomp = 5)

plot(PLSR %>% predict(type = "response") %>% .[,,5]  ~ WA$Scenedesmus[15:27])

### Dentro del mismo rango de concentraciones 2 iniciales (siguientes 14 puntos)
PCAwa = prcomp(WA[28:42,1:512],scale. = T)
PCAwa %>% fviz_pca_biplot(col.ind = WA$Scenedesmus[28:42],geom.var = "",
                          gradient.cols = c("#207557","#bd71c7","#e3f542"))
PCAwa$rotation[,1] %>% plot()

PLSR  =plsr(Scenedesmus ~ .,data = WA[28:42,] %>% select(-Dolichospermum,-LED),ncomp = 5)

plot(PLSR %>% predict(type = "response") %>% .[,,5]  ~ WA$Scenedesmus[28:42])

### Dentro del mismo rango de concentraciones 2 iniciales (siguientes 14 puntos)
PCAwa = prcomp(WA[43:57,1:512],scale. = T)
PCAwa %>% fviz_pca_biplot(col.ind = WA$Scenedesmus[43:57],geom.var = "",
                          gradient.cols = c("#207557","#bd71c7","#e3f542"))
PCAwa$rotation[,1] %>% plot()

PLSR  =plsr(Scenedesmus ~ .,data = WA[43:57,] %>% select(-Dolichospermum,-LED),ncomp = 5)

plot(PLSR %>% predict(type = "response") %>% .[,,5]  ~ WA$Scenedesmus[43:57])

### Dentro del mismo rango de concentraciones 2 iniciales (siguientes 14 puntos)
PCAwa = prcomp(WA[58:1:512],scale. = T)
PCAwa %>% fviz_pca_biplot(col.ind = WA$Scenedesmus[43:57],geom.var = "",
                          gradient.cols = c("#207557","#bd71c7","#e3f542"))
PCAwa$rotation[,1] %>% plot()

PLSR  =plsr(Scenedesmus ~ .,data = WA[43:57,] %>% select(-Dolichospermum,-LED),ncomp = 5)

plot(PLSR %>% predict(type = "response") %>% .[,,5]  ~ WA$Scenedesmus[43:57])

# QUé pasa si usamos sólo las concentraciones preparadas con el mismo cultivo, el mismo día

### Dentro del mismo rango de concentraciones 2 iniciales (siguientes 14 puntos)
PCAwa = prcomp(WA[1:27,1:512],scale. = T)
PCAwa %>% fviz_pca_biplot(col.ind = WA$Scenedesmus[1:27],geom.var = "",
                          gradient.cols = c("#207557","#bd71c7","#e3f542"))
PCAwa$rotation[,1] %>% plot()

PLSR  =plsr(Scenedesmus ~ .,data = WA[1:27,] %>% select(-Dolichospermum,-LED),ncomp = 5)

plot(PLSR %>% predict(type = "response") %>% .[,,5]  ~ WA$Scenedesmus[1:27])

##### led 470

led440 = allData %>% filter(LED == "440")
led440[1:512] = led440[1:512] %>% apply(1,function(x){x/max(abs(x))}) %>% t()


PCA440 = prcomp(led440[led440$Scenedesmus + led440$Dolichospermum != 0,1:512],scale. = T)
PCA440 %>% fviz_pca_biplot(col.ind = led440$Scenedesmus[led440$Scenedesmus + led440$Dolichospermum != 0],geom.var = "",
                          gradient.cols = c("#207557","#bd71c7","#e3f542"))
PCAwa$rotation[,1] %>% plot()

PLSR  =plsr(Scenedesmus ~ .,data = led440 %>% select(-Dolichospermum,-LED,200:450),ncomp = 5)
plot(PLSR %>% predict(type = "response") %>% .[,,5]  ~ WA$Scenedesmus)

