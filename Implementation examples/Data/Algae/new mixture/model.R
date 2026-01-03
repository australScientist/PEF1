library(tidyverse)
library(jsonlite)
library(prospectr)
library(pls)
library(caret)
library(patchwork)
library(baseline)
library(ggrepel)

sampler = function(var,data){
  x = data %>% select(starts_with(var)) %>% unlist() %>% unname()
  indices = (!(x %>% max() == x) )%>% which()
  indices = sample(indices,size = round(length(indices)*0.25))
  return(indices)
}

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

# Load data ---------------------------------------------------------------

B1 = read_json(path = "Implementation/Data/Algae/new mixture/1block.json")
B2 = read_json(path = "Implementation/Data/Algae/new mixture/2block.json")
B3 = read_json(path = "Implementation/Data/Algae/new mixture/4block.json")
B4 = read_json(path = "Implementation/Data/Algae/new mixture/1test.json")
B5 = read_json(path = "Implementation/Data/Algae/new mixture/2test.json")
B6 = read_json(path = "Implementation/Data/Algae/new mixture/3test.json")
B7 = read_json(path = "Implementation/Data/Algae/new mixture/4test.json")

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
B1 = map(B1,getSpectra)
B1 = B1 %>% bind_rows()
B2 = map(B2,getSpectra)
B2 = B2 %>% bind_rows()
B3 = map(B3,getSpectra)
B3 = B3 %>% bind_rows()
B4 = map(B4,getSpectra)
B4 = B4 %>% bind_rows()
B5 = map(B5,getSpectra)
B5 = B5 %>% bind_rows()
B6= map(B6,getSpectra)
B6 = B6 %>% bind_rows()
B7 = map(B7,getSpectra)
B7 = B7 %>% bind_rows()
Data = bind_rows(B1,B2,B3,B4,B5,B6,B7)
colnames(Data)[5:516] = paste0("p",1:512)

rm(B1,B2,B3,B4,B5,B6,B7)
DataSep = list()
DataSep[["WA"]] = Data %>% filter(LED == "WA")
DataSep[["365"]] = Data %>% filter(LED == "365")
DataSep[["395"]] = Data %>% filter(LED == "395")
DataSep[["440"]] = Data %>% filter(LED == "440")
DataSep[["470"]] = Data %>% filter(LED == "470")
DataSep[["530"]] = Data %>% filter(LED == "530")
DataSep[["600"]] = Data %>% filter(LED == "600")
DataSep[["660"]] = Data %>% filter(LED == "660")


DataSep[["WA"]] %>%
  filter(PBS != 1) %>% 
  select(starts_with("p",ignore.case = F)) %>%
  .[-c(36,50,90,80,96),] %>%  # multivariate outliers
  prcomp() %>% biplot()

Absorbance = DataSep[["WA"]] %>%
  filter(PBS != 1) %>% 
  select(starts_with("p",ignore.case = F)) %>%
  .[-c(36,50,90,80,96),]

meanBlank = DataSep[["WA"]] %>%
  filter(PBS == 1) %>% 
  .[-c(8,14),] %>% 
  select(starts_with("p",ignore.case = F)) %>% 
  colMeans()

blankMatrix = (1/meanBlank) %>% diag()

Absorbance = -log10(as.matrix(Absorbance) %*% blankMatrix) %>% as_tibble()
Absorbance %>% prcomp() %>% biplot()
Absorbance = cbind(DataSep[["WA"]] %>%
        filter(PBS != 1) %>% 
        select(!starts_with("p",ignore.case = F)) %>%
        .[-c(36,50,90,80,96),],Absorbance)

Blanks = DataSep[["WA"]] %>%
  filter(PBS == 1) %>% 
  .[-c(8,14),] %>% 
  select(starts_with("p",ignore.case = F))

BlanksAbs = -log10(as.matrix(Blanks) %*% blankMatrix) %>% as_tibble()
BlanksAbs = bind_cols(Blanks = DataSep[["WA"]] %>%
            filter(PBS == 1) %>% 
            .[-c(8,14),] %>% 
            select(!starts_with("p",ignore.case = F)),BlanksAbs)
Absorbance = bind_rows(Absorbance,BlanksAbs)

colours = rgb(Absorbance$PBS,Absorbance$Scenedesmus,Absorbance$Dolichospermum)

AbsPCA = Absorbance %>% select(starts_with("V",ignore.case = F)) %>% prcomp()

AbsPCA$x %>% as_tibble() %>% ggplot(aes(PC1,PC2)) + geom_point(colour = colours,size = 3)
AbsPCA$rotation %>% 
  as_tibble() %>% 
  mutate(pix = 1:512) %>%
  pivot_longer(starts_with("PC")) %>% 
  filter(name == "PC4") %>% 
  ggplot(aes(pix,value)) + geom_line()



AbsPCA$x %>% as_tibble() %>% ggplot(aes(PC1,PC2)) + geom_point(colour = colours)
AbsPCA$x %>% as_tibble() %>% ggplot(aes(PC1,PC2)) + geom_point(colour = colours,size = 3)

SC470 = DataSep[["470"]] %>%
  filter(PBS != 1) %>% 
  select(starts_with("p",ignore.case = F)) %>%
  .[-c(36,50,90,80,96),]

# Projection plot ---------------------------------------------------------


PLS2 = plsr(as.matrix(AbsNorm[c(1,2,3)]) ~ .,data = Absorbance %>% select(starts_with("V")),validation = "LOO",ncomp = 10)

visPLS2 = PLS2$scores[,2:3] %>% as_tibble() %>% mutate(colours = colours)
visPLS2_yScores = PLS2$Yloading %>% as.matrix() %>% .[,2:3]
genus = rownames(visPLS2_yScores)
visPLS2_yScores = visPLS2_yScores %>% as_tibble() %>% mutate("Genus" = genus)

projectionPlot = ggplot(visPLS2,aes(x = `Comp 2`,y = `Comp 3`)) + 
  geom_point(size = 3,colour = visPLS2$colours) +
  geom_segment(data = visPLS2_yScores,aes(x = 0,y = 0,xend = `Comp 2`,yend = `Comp 3`),
               arrow = arrow(type = "closed",length = unit(9,units = "pt"))) +
  # geom_label(data = visPLS2_yScores,aes(x = `Comp 2`*1.2,y = `Comp 3`*1.2,label = Genus)) + 
  theme_classic() + xlim(c(-1.2,1)) + xlab("PLS 2") + ylab("PLS3") +
  labs(tag = "b)")
#visPLS2_xScores = 


# Absorbance plot ---------------------------------------------------------

AbsorbanceLong = Absorbance %>% 
  mutate(ID = 1:nrow(Absorbance)) %>% 
  pivot_longer(cols = starts_with("V"),
               names_to = "Pixel",
               values_to = "Value") %>% 
  mutate(Pixel = Pixel %>% str_remove("V") %>% as.numeric(),
         colour = rgb(PBS,Dolichospermum,Scenedesmus))

  
AbsorbancePlot = AbsorbanceLong %>% ggplot(aes(Pixel,Value,colour = colour,group = ID)) +
  geom_line() +
  scale_colour_identity() +
  theme_classic() +
  ylab("Absorbance (AU)")

# colour Triangle ---------------------------------------------------------

n = 100
triangle_data = expand.grid(x = seq(0, 1, length.out = n), 
                             y = seq(0, 1, length.out = n)) %>% filter(x + y <= 1)

# Define barycentric coordinates for RGB blending
triangle_data = triangle_data %>%
  mutate(
    r = (1 - x - y) %>% round(4),  # Red proportion
    g = x,          # Green proportion
    b = y,          # Blue proportion
    color = rgb(r, g, b) 
  )

# 
Refference = ggplot(triangle_data, aes(x = x, y = y)) +
  geom_raster(aes(fill = color)) +
  scale_fill_identity() + 
  coord_fixed() +
  theme_void() +
  geom_label(data = tibble(x = c(0,0.3,1),
                          y = c(0,1,0),
                          label = c("PBS","Dolichospermum","Scenedesmus")),
             aes(x = x,y = y,label = label),size = 3) +
  ylim(c(-0.1,1.1)) +
  xlim(c(-0.25,1.35))

absPlot2 = AbsorbancePlot + inset_element(Refference, 0.4, 0.6, 0.8, 1)
absPlot2 = absPlot2 

# Train -------------------------------------------------------------------
set.seed(18)
scenedesmus = (26+18+34+22+33+38+38+28)*10000/(2*4)*5
dolichospermum = 2e6*(0.23662+0.24116)/2+98619
Absorbance = mutate(Absorbance,Scenedesmus = Scenedesmus * scenedesmus,
       Dolichospermum = Dolichospermum * dolichospermum)
AbsNorm = mutate(Absorbance,
                 Scenedesmus = ifelse(Absorbance$Scenedesmus != 0,
                                      Absorbance$Scenedesmus/(Absorbance$Scenedesmus + Absorbance$Dolichospermum),
                                  0),
                 Dolichospermum = ifelse(Absorbance$Dolichospermum != 0 ,
                                         Absorbance$Dolichospermum/(Absorbance$Scenedesmus + Absorbance$Dolichospermum),
                                         0)
                 )

indices = createDataPartition(Absorbance$Scenedesmus, 
                              p = .8, 
                              list = FALSE, 
                              times = 1) %>% .[,1] %>% unlist()
AbsTrain = Absorbance[indices,]
AbsTest = Absorbance[-indices,]

AbsNormTrain = AbsNorm[indices,]
AbsNormTest = AbsNorm[-indices,]

train_ctrl = trainControl(method="repeatedcv", # type of resampling in this case Cross-Validated
                          number=10,
                          repeats = 10, # number of folds
                          search = "grid" # we are performing a "grid-search"
)



modelScen = train(Scenedesmus ~ .,
      data = AbsTrain %>% select(starts_with("V"),Scenedesmus) %>% as_tibble(),
      method = "pls",
      tuneLength = 15,
      trControl = train_ctrl
)
plot(predict(modelScen,AbsTest)~AbsTest$Scenedesmus) 

modelDolico = train(Dolichospermum ~ .,
              data = AbsTrain %>% select(starts_with("V"),Dolichospermum) %>% as_tibble(),
              method = "pls",
              tuneLength = 15,
              trControl = train_ctrl
)
plot(predict(modelScen,AbsTest)~AbsTest$Scenedesmus) 


modelNormScen = train(Scenedesmus ~ .,
                  data = AbsNormTrain %>% select(starts_with("V"),Scenedesmus) %>% as_tibble(),
                  method = "pls",
                  tuneLength = 15,
                  trControl = train_ctrl
)


modelNormDolico = train(Dolichospermum ~ .,
                    data = AbsNormTrain %>% select(starts_with("V"),Dolichospermum) %>% as_tibble(),
                    method = "pls",
                    tuneLength = 15,
                    trControl = train_ctrl
)


MAPEtestScen = tibble("Predicted" = predict(modelNormScen,AbsNormTest),
       "Real" = AbsNormTest$Scenedesmus,
       "Set" = "Test") %>% mutate(MAE = abs(Predicted - Real)*100)%>% 
  select(MAE) %>% unlist() %>% mean()%>% 
  format(digits = 2)

MAPEtestDol = tibble("Predicted" = predict(modelNormDolico,AbsNormTest),
                      "Real" = AbsNormTest$Dolichospermum,
                      "Set" = "Test") %>% mutate(MAE = abs(Predicted - Real)*100)%>% 
  select(MAE) %>% unlist() %>% mean()%>% 
  format(digits = 2)

MAEtestScen = tibble("Predicted" = predict(modelScen,AbsTest),
                      "Real" = AbsTest$Scenedesmus,
                      "Set" = "Test") %>% mutate(MAE = abs(Predicted - Real))%>% 
  select(MAE) %>% unlist() %>% mean()%>% 
  format(digits = 2)

MAEtestDol = tibble("Predicted" = predict(modelDolico,AbsTest),
                     "Real" = AbsTest$Dolichospermum,
                     "Set" = "Test") %>% mutate(MAE = abs(Predicted - Real))%>% 
  select(MAE) %>% unlist() %>% mean()%>% 
  format(digits = 2)



plot(predict(modelDolico,AbsTest)~AbsTest$Dolichospermum) 
mockTest = AbsTest
colnames(mockTest)[1:2] = c("Dolichospermum","Scenedesmus")

plot(predict(model2,mockTest)~mockTest$Dolichospermum) 
plot(predict(model,mockTest)~mockTest$Scenedesmus) 

RegCoefDolicho = modelDolico$finalModel$coefficients %>% 
  .[1:512,,1:4]%>% 
  as_tibble() %>% 
  rename_with(.cols = everything(),
              .fn = function(x){paste0("C",1:4)}) %>% 
  mutate(pixel = 1:512,"Genus" = "Dolichospermum") %>% 
  pivot_longer(starts_with("C"),names_to = "Component",values_to = "Value")

CoefDolicho = RegCoefDolicho %>% ggplot(aes(pixel,Value,colour = Component)) +
  geom_line(linewidth = 1.3)  + 
  scale_colour_manual(values = c("steelblue1",
                                 "navyblue",
                                 "blue",
                                 "royalblue2"
                                 )) +
  theme_classic() +
  labs(tag = "c)") +
  theme(
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  ylab("Coefficient") +
  geom_label(data = tibble(x = 0,
                           y = -5e4,
                           lab = " Dolichospermum "),
             aes(x = x,
                 y = y,
                 label = lab),
             fill = "black",
             colour = "white",
             label.r = unit(0,"mm"),
             hjust = 0,
             fontface = "bold") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.15),    
    legend.direction = "horizontal",   
    legend.box = "horizontal"         
  )

RegCoefScen = modelScen$finalModel$coefficients %>% 
  .[1:512,,1:4]%>% 
  as_tibble() %>% 
  rename_with(.cols = everything(),
              .fn = function(x){paste0("C",1:4)}) %>% 
  mutate(pixel = 1:512,"Genus" = "Scenedesmus") %>% 
  pivot_longer(starts_with("C"),names_to = "Component",values_to = "Value")

CoefScen = RegCoefScen %>% ggplot(aes(pixel,Value,colour = Component)) +
  geom_line(linewidth = 1.3) +
  theme_classic()+
  scale_colour_manual(values = c("palegreen1",
                                 "green4",
                                 "olivedrab2",
                                 "green"))+
  labs(tag = "d)")  +
  ylab("Coefficient") +
  geom_label(data = tibble(x = 0,
                           y = 2e5,
                           lab = " Scenedesmus "),
             aes(x = x,
                 y = y,
                 label = lab),
             fill = "black",
             colour = "white",
             label.r = unit(0,"mm"),
             hjust = 0,
             fontface = "bold") +
  scale_x_continuous(breaks = c(100,250,400)) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.15),    
    legend.direction = "horizontal",   
    legend.box = "horizontal"         
  )

CoefScen




regressionPlot = bind_rows(RegCoefScen,RegCoefDolicho) %>% 
  ggplot(aes(pixel,value,colour = Genus)) + 
  geom_line(linewidth = 1.3) +
  theme_classic() +
  scale_colour_manual(values = c("green2","blue2")) +
  geom_hline(yintercept = 0,linetype = "dashed")+
  ylab("Regression Coefficient") +
  xlab("Pixel") + theme(legend.position = "top")+
  labs(tag = "c)")

absPlot2 + projectionPlot + regressionPlot+ plot_layout(design = "AAB
                                                        CCB")
predRealScen = bind_rows(tibble("Predicted" = predict(modelScen),
                      "Real" = AbsTrain$Scenedesmus,
                      "Set" = "Train"),
                      tibble("Predicted" = predict(modelScen,AbsTest),
                             "Real" = AbsTest$Scenedesmus,
                             "Set" = "Test"))
RVsPscen = ggplot(predRealScen) +
  geom_point(aes(Real,Predicted),
             colour = "green2",
             alpha = ifelse(predRealScen$Set == "Train",0.3,1),
             size = ifelse(predRealScen$Set == "Train",1.5,2.5)) +
  geom_abline(slope = 1,intercept = 0,linetype = "dotted") +
  theme_classic()+
  labs(tag = "e)")

predRealDol = bind_rows(tibble("Predicted" = predict(modelDolico),
                                "Real" = AbsTrain$Dolichospermum,
                                "Set" = "Train"),
                         tibble("Predicted" = predict(modelDolico,AbsTest),
                                "Real" = AbsTest$Dolichospermum,
                                "Set" = "Test"))
RVsPdol = ggplot(predRealDol) +
  geom_point(aes(Real,Predicted),
             colour = "blue",
             alpha = ifelse(predRealScen$Set == "Train",0.3,1),
             size = ifelse(predRealScen$Set == "Train",1.5,2.5)) +
  geom_abline(slope = 1,intercept = 0,linetype = "dotted") +
  theme_classic()+
  labs(tag = "f)")

layout = paste0(str_dup("A",9),str_dup("B",3),"\n",
                str_dup("A",9),str_dup("B",3),"\n",
                str_dup("A",9),str_dup("B",3),"\n",
                str_dup("A",9),str_dup("B",3),"\n",
                str_dup("A",9),str_dup("B",3),"\n",
                str_dup("A",9),str_dup("B",3),"\n",
                str_dup("D",9),str_dup("C",3),"\n",
                str_dup("D",9),str_dup("C",3),"\n",
                str_dup("D",9),str_dup("C",3),"\n",
                str_dup("F",9),str_dup("E",3),"\n",
                str_dup("F",9),str_dup("E",3),"\n",
                str_dup("F",9),str_dup("E",3),"\n")

absPlot2 + projectionPlot + RVsPdol  + CoefDolicho + RVsPscen +CoefScen +plot_layout(design = layout)


CoefDolicho/CoefScen