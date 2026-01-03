library(tidyverse)
library(jsonlite)
library(pls)
library(mdatools)
library(prospectr)
library(factoextra)

get_data = function(path,rows = NULL){
  rawData = read_json(path)
  if (!is.null(rows)){
    rawData = rawData[rows]
  }
  Spectra = map(rawData,function(x){x %>% .[["spectrum"]] %>% unlist()}) %>% bind_cols() %>% t() %>% as_tibble()
  ID = map(rawData,function(x){x %>% .[["tag"]] %>% strsplit("-") %>% unlist()}) %>% 
   bind_cols() %>% t() %>% as_tibble()
  SiD = map(rawData,function(x){x %>% .[["ID"]] %>% unlist()}) %>% bind_cols() %>% t() %>% as_tibble()
  
  colnames(Spectra) = 1:512
  colnames(ID) = c("Tag","LED")
  ID$ID = SiD$V1
  Data = list(Spectra,ID)
  bind = possibly(bind_cols,list(Spectra,ID))
  Data = bind(Spectra,ID) 
  return(Data)
}

filterSpectra = function(led,Coso = CleanTurbiSpectra,filtro = WAFilter){
  print(led)
  Sp = filter(Coso, LED == led)
  print(nrow(Sp))
  Sp = Sp[filtro,]
  print(nrow((Sp)))
  return(Sp)
}
"Implementation/Data/Algae/new mixture/"
Algae = get_data("Implementation/Data/Algae/mixture/data/../Models/Algae/ALGAS JUEVES MODIFIED.json")
Monday = get_data("../Models/Algae/1results Modified.json",rows = 641:882)
Monday = Monday[Monday$Tag %in% (Monday$Tag %>% unique %>% .[1:17]),]


### EXPLORTE ALGAE

gatheredData = gather(Algae,"pixel","value",1:512)
gatheredData$pixel = gatheredData$pixel %>% as.numeric()
gatheredData$value = gatheredData$value %>% as.numeric()

gatheredData %>% filter(LED == "WA") %>% ggplot(aes(pixel,value, colour = as.factor(ID))) + geom_point(size = 0.5) +
  facet_grid(rows = vars(LED),cols = vars(Tag),scales = "free")
#Clean Algae
WAFilter = !(filter(Algae,LED == "WA") %>% .$ID %in% c(1,10,19,28,127,154))
CleanAlgae = map(unique(Algae$LED),filterSpectra,Coso = Algae) %>% bind_rows()
### EXPLORTE MONDAY

gatheredData = gather(Monday,"pixel","value",1:512)
gatheredData$pixel = gatheredData$pixel %>% as.numeric()
gatheredData$value = gatheredData$value %>% as.numeric()

gatheredData %>% filter(LED == "WA",ID != 802,ID != 649) %>% ggplot(aes(pixel,value, colour = as.factor(ID))) + geom_point(size = 0.5) +
  facet_grid(rows = vars(LED),cols = vars(Tag),scales = "free")

#Clean Monday
Monday = filter(Monday,Tag != "scenedesmus 1")
WAFilter = !(filter(Monday,LED == "WA") %>% .$ID %in% c(649,802))
CleanMonday = map(unique(Monday$LED),filterSpectra,Coso = Monday) %>% bind_rows()
CleanMonday$Tag = CleanMonday$Tag %>% str_split(" ") %>% map(function(x)x[1]) %>% unlist()

######### Joining Forces
Data = bind_rows(CleanAlgae,CleanMonday)
Data$Tag = ifelse(Data$Tag == "dolichospermun","Dolichospermum",Data$Tag)
Data = filter(Data,LED != 10)
## RAAW
gatheredData = gather(Data,"pixel","value",1:512)
gatheredData$pixel = gatheredData$pixel %>% as.numeric()
gatheredData$value = gatheredData$value %>% as.numeric()

gatheredData %>% filter(LED == "WA",) %>% ggplot(aes(pixel,value, colour = as.factor(ID))) + geom_point(size = 0.5) +
  facet_grid(rows = vars(LED),cols = vars(Tag),scales = "free")
## Scaling
ScaledSpectra = apply(Data[1:512],2,function(x){scale(x,center = T,scale = T)})  %>% as_tibble()
colnames(ScaledSpectra) = 1:512
ScaledSpectra$Tag = Data$Tag
ScaledSpectra$LED = Data$LED

gatheredData = gather(ScaledSpectra,"pixel","value",1:512)
gatheredData$pixel = gatheredData$pixel %>% as.numeric()
gatheredData$value = gatheredData$value %>% as.numeric()
gatheredData %>% ggplot(aes(pixel,value)) + geom_point(size = 0.5) +
  facet_grid(rows = vars(LED),cols = vars(Tag),scales = "free")

# Normalized
NormalizedSpectra = apply(Data[1:512],1,function(x){x/max(x)}) %>% t()  %>% as_tibble()
colnames(NormalizedSpectra) = 1:512
NormalizedSpectra$Tag = Data$Tag
NormalizedSpectra$LED = Data$LED

gatheredData = gather(NormalizedSpectra,"pixel","value",1:512)
gatheredData$pixel = gatheredData$pixel %>% as.numeric()
gatheredData$value = gatheredData$value %>% as.numeric()
gatheredData %>% filter(LED == "WA") %>% ggplot(aes(pixel,value)) + geom_point(size = 0.5) +
  facet_grid(rows = vars(LED),cols = vars(Tag),scales = "free")

#derivative 

Derivative = savitzkyGolay(Data[1:512],m = 1,p = 3,w = 15) %>% as_tibble()
colnames(NormalizedSpectra) = 1:498
Derivative$Tag = Data$Tag
Derivative$LED = Data$LED
Derivative$ID = Data$ID

gatheredData = gather(Derivative,"pixel","value",1:498)
gatheredData$pixel = gatheredData$pixel %>% as.numeric()
gatheredData$value = gatheredData$value %>% as.numeric()
gatheredData %>% filter(LED == "WA") %>% ggplot(aes(pixel,value,colour = as.factor(ID))) + geom_line(size = 0.5) +
  facet_grid(rows = vars(LED),cols = vars(Tag),scales = "free")

## DATA MERGED
DataMerged = bind_cols(filter(Data,LED == "WA") %>% .[1:512],
                       filter(Data,LED == "395") %>% .[1:512],
                       filter(Data,LED == "440") %>% .[1:512],
                       filter(Data,LED == "470") %>% .[1:512],
                       filter(Data,LED == "530") %>% .[1:512],
                       filter(Data,LED == "600") %>% .[1:512],
                       filter(Data,LED == "660") %>% .[1:512],
                       filter(Data,LED == "660") %>% .[513:515])
colnames(DataMerged) = c(1:3584,"tag","LED","ID")
gatheredData = gather(DataMerged,"pixel","value",1:3584)
gatheredData$pixel = gatheredData$pixel %>% as.numeric()
gatheredData$value = gatheredData$value %>% as.numeric()
gatheredData  %>% ggplot(aes(pixel,value,colour = tag)) + geom_point(size = 0.5) #+

OnlyData = DataMerged %>% .[1:3584] %>% as.matrix()
rownames(OnlyData) = DataMerged %>% .$tag
DataScaled = OnlyData %>% apply(2,scale)
rownames(DataScaled) = paste(DataMerged %>% .$tag,1:29)
clust = DataScaled %>% dist() %>% hclust()
dendrogram = fviz_dend(clust,k=2,rect = TRUE,rect_fill = TRUE,palette = "npg")
dendrogram$layers[[2]]$data$angle = 45
dendrogram
groups = cutree(clust,2)
fviz_cluster(list(data=DataScaled,cluster=groups),
             ellipse.type = "t",ellipse = T)

clust = kmeans(DataScaled,centers = 2,iter.max = 25)
fviz_cluster(object = clust,data = DataScaled,ellipse = T,ellipse.type = "norm")

predict(clust,DataMerged %>% filter(tag == "SanRoque") %>% .[1:3584] %>% as.matrix())

plot(clust)

PCA = prcomp(DataScaled)
PCA$rotation[,4] %>% unlist() %>% plot()
plot(PCA$x[,1:2],)
ggplot(PCA$x %>% as_tibble(),aes(PC1,PC2,colour=PCA$x %>% rownames())) + geom_point()

Classess = DataMerged %>% filter(tag != "SanRoque") %>% .$tag %>% as.factor()

PLSDA = plsda(x = DataScaled,
                c = Classess,center = F,scale = F,ncomp = 3,cv = 1)
plotYVariance(PLSDA)
plot(PLSDA)
mdatools::plotPredictions(PLSDA)
plotRegcoeffs(PLSDA)
plotXYScores(PLSDA)
PLSDA %>% summary() 

ResSanRoque = predict(PLSDA, DataMerged %>% filter(tag == "SanRoque") %>% .[1:3584] %>% as.matrix())
plotPredictions(ResSanRoque)
summary(PLSDA$coeffs,ncomp = 3)

Classess %>% summary()

## DATA MERGED
NormMerged = bind_cols(filter(NormalizedSpectra,LED == "WA") %>% .[1:512],
                       filter(NormalizedSpectra,LED == "395") %>% .[1:512],
                       filter(NormalizedSpectra,LED == "440") %>% .[1:512],
                       filter(NormalizedSpectra,LED == "470") %>% .[1:512],
                       filter(NormalizedSpectra,LED == "530") %>% .[1:512],
                       filter(NormalizedSpectra,LED == "600") %>% .[1:512],
                       filter(NormalizedSpectra,LED == "660") %>% .[1:512],
                       filter(NormalizedSpectra,LED == "660") %>% .[513:514])
colnames(NormMerged) = c(1:3584,"tag","LED")
gatheredData = gather(NormMerged,"pixel","value",1:3584)
gatheredData$pixel = gatheredData$pixel %>% as.numeric()
gatheredData$value = gatheredData$value %>% as.numeric()
gatheredData  %>% ggplot(aes(pixel,value,colour = tag)) + geom_point(size = 0.5)

OnlyData = NormMerged %>% filter(tag != "SanRoque") %>% .[1:3584] %>% as.matrix()
rownames(OnlyData) = NormMerged %>% filter(tag != "SanRoque") %>% .$tag
NormScaled = OnlyData %>% apply(2,scale)
rownames(NormScaled) = Classess
clust = NormScaled %>% dist() %>% hclust()

plot(clust)

PCA = prcomp(DataScaled)
PCA$rotation[,4] %>% unlist() %>% plot()
plot(PCA$x[,1:2])

Classess = DataMerged %>% filter(tag != "SanRoque") %>% .$tag %>% as.factor()

PLSDA = plsda(x = NormScaled,
              c = Classess,center = F,scale = F,ncomp = 3)
plotYVariance(PLSDA)
plot(PLSDA)
mdatools::plotPredictions(PLSDA)
plotRegcoeffs(PLSDA)
plotXYScores(PLSDA)
PLSDA %>% 
  
  Classess %>% summary()
