library(tidyverse)
library(reticulate)

reticulate::source_python("analize.py")


spectra = read_csv("Data/Spectra.csv",col_names = 1:512)


spectra[2:6,300] %>% unlist() %>% plot(type="l")
spectra[3,300] %>% unlist() %>% plot(type="l")
spectra[4,300] %>% unlist() %>% plot(type="l")
spectra[5,300] %>% unlist() %>% plot(type="l")

spectra = read_csv("Data/Spectra0-50000C1000.csv",col_names = 1:512)

spectra[2:50,10] %>% unlist() %>% plot(type="l")
spectra[3,] %>% unlist() %>% plot(type="l")
spectra[4,] %>% unlist() %>% plot(type="l")
spectra[40,] %>% unlist() %>% plot(type="l")

spectra = read_csv("Data/SpectraPWMat20000.csv",col_names = 1:512)
colnames(spectra) = 1:512
list.files("Data/PWM and time/20000/")
spectra = spectra[spectra %>% rowSums() %>% sort() %>% match(spectra %>% rowSums()),]
spectra$pwm = seq(0.1,0.95,0.05)
ggplot(spectra %>% gather("pix","value",1:512),aes(pix %>% as.numeric(),value,colour= pwm %>% as.factor())) + geom_line()

PCR = pls::pcr(pwm ~ . ,data=spectra,scale= TRUE,center=TRUE,ncomp=1) 
PLSR = pls::plsr(pwm ~ . ,data=spectra,scale= TRUE,center=TRUE,ncomp=1)
PLSR %>% plot()
PLSR %>% summary()

spectra[1:18,10] %>% unlist() %>% sort() %>% summary() 
PCA = prcomp((spectra),center = TRUE,scale. = TRUE)
summary(PCA)
PCA$x[,1] %>% plot()
spectra[3,] %>% unlist() %>% plot(type="l")
spectra[4,] %>% unlist() %>% plot(type="l")
spectra[40,] %>% unlist() %>% plot(type="l")

spectra = read_csv("Data/SpectraPWMat40000.csv",col_names = 1:512)

list.files("Data/PWM and time/40000/")
spectra[1:18,10] %>% unlist() %>% sort() %>% plot()
spectra[2:50,10] %>% unlist() %>% plot(type="l")
spectra[3,] %>% unlist() %>% plot(type="l")
spectra[4,] %>% unlist() %>% plot(type="l")
spectra[40,] %>% unlist() %>% plot(type="l")

# PWM an time to find time behaviour
all = paste0("Data/TimeAndPWM2/",list.files("Data/TimeAndPWM2/")) %>% map(py$get_spectrum) %>% bind_cols()
all = t(all)
all = as_tibble(all)
colnames(all) = paste0("P",1:512)
names = list.files("Data/TimeAndPWM2/")
PWM = names %>% str_extract(pattern = regex("...(?=m)"))
expTime = names %>% str_extract(pattern = regex("(?<=t).*(?=\\.)"))
all$PWM = PWM
all$expTime = expTime %>% as.numeric()
ggplot(all ,aes(x = expTime,y = P344,colour = PWM)) + geom_line() + geom_point()
ggplot(all %>% filter(expTime <46000),aes(x = expTime,y = P344,colour = PWM)) + geom_point()
mod = lm(expTime ~ P344*as.numeric(PWM) ,data = all %>% filter(expTime <46000))

models = all %>% filter(expTime <41000) %>% group_by(PWM) %>% nest() %>% .$data %>% map(function(x){lm(P344 ~ expTime,data = x)})
coefs = models %>% map(coef) %>% bind_rows()
coefs$PWM = c(0,0.2,0.4,0.6,0.8,1)
ggplot(coefs,aes(expTime,PWM)) + geom_point()

# New PWM an time to find time behaviour
all = paste0("Data/TimeAndPWM3/",list.files("Data/TimeAndPWM3/")) %>% map(py$get_spectrum) %>% bind_cols()
all = t(all)
all = as_tibble(all)
colnames(all) = paste0("P",1:512)
names = list.files("Data/TimeAndPWM3/")
PWM = names %>% str_extract(pattern = regex(".*(?=measurement)"))
expTime = names %>% str_extract(pattern = regex("(?<=t).*(?=\\.)"))
all$PWM = PWM %>% as.numeric() %>% round(digits = 4)
all$expTime = expTime %>% as.numeric() %>% round(digits = 3)
ggplot(all %>% filter(PWM<0.3) ,aes(x = expTime,y = P344,colour = PWM)) + geom_point()
ggplot(all %>% filter(expTime <46000),aes(x = expTime,y = P344,colour = PWM %>% as_factor())) + geom_point() + geom_line()
mod = lm(expTime ~ P344*as.numeric(PWM) ,data = all %>% filter(expTime <46000))

models = all %>% filter(P344 <34000) %>% group_by(PWM) %>% nest() %>% .$data %>% map(function(x){lm(P344 ~ expTime,data = x)})
coefs = models %>% map(coef) %>% bind_rows()
coefs$PWM = all %>% filter(P344 <34000) %>% group_by(PWM) %>% nest() %>% .$PWM
ggplot(coefs,aes(PWM,expTime)) + geom_point()
ggplot(all %>% filter(expTime == 21000) ,aes(x = PWM,y = P344)) + geom_point() + geom_line()


# How does the least intense light behave
filter(all,PWM == 0.015) %>% ggplot(aes(expTime,P344)) + geom_point()
mod = lm(P344 ~ expTime, data = filter(all,PWM == 0.015))
mod %>% summary()
#is it enough to use 4 points?
modRestricted = mod = lm(P344 ~ expTime, data = filter(all,PWM == 0.015 & expTime <= 21000))
modRestricted %>% summary()
coef(mod)/coef(modRestricted)
## 
prediction = predict(modRestricted,newdata = filter(all,PWM == 0.015)) %>% unname() %>% unlist()
real = filter(all,PWM == 0.015) %>% select(P344) %>% unlist()
plot(real,prediction)
lm(real~prediction) %>% summary()


#NEW DATA
library(rjson)

#Results = jsonlite::read_json("Data/results.json")
Results = jsonlite::read_json("Data/factorialExperiment/Lu1Ca1Ci1.json")
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
