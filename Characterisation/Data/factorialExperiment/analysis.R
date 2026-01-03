library(tidyverse)
library(jsonlite)

files = list.files("Characterisation/Data/factorialExperiment/Results/",full.names = T) %>% str_subset(".json") 
# Replacing infinity so that can be  read
for(i in files){
  read_file(i) %>% str_replace_all("Infinity","0") %>% write_file(file = i %>% str_replace("Raw",replacement = "Processed"))
}

files = list.files("Calibration/Ecophysiometer/Raw",full.names = T) %>% str_subset(".json") 
# Replacing infinity so that can be  read
for(i in files){
  read_file(i) %>% str_extract_all("\"tag\": \"B") %>% View()
}
files = list.files("Characterisation/Data/factorialExperiment/Results/processed/",full.names = T) %>% str_subset(".json") 
for(i in files){
  read_file(i) %>% str_replace_all("NaN","0") %>% write_file(file = i %>% str_replace("Raw",replacement = "Processed"))
}

#checking that all files can be read
files = list.files("Characterisation/Data/factorialExperiment/Results/processed/",full.names = T) %>% str_subset(".json") 

for(i in files){
  print(i)
  read_json(i)
}

# Processing files --------------------------------------------------------



read_spectra = function(file){
  require(jsonlite)
  require(tidyverse)
  
  json = read_json(path = file)
  spectra = map(json,function(x){x$spectrum %>% unlist()}) %>% 
    unlist() %>%
    matrix(byrow = T,ncol = 512) %>%
    as_tibble() %>% 
    rename_all(function(x){x %>% str_replace("V","p")})
  Tags = map(json,function(x){x$tag %>% unlist()}) %>% 
    unlist()
  Led = map(json,function(x){x$tag %>% unlist()}) %>% 
    unlist() %>% str_extract("(?<=-)[:alnum:]+")
  Block = file %>% str_extract("[:alnum:]+(?=.json)")
  output = mutate(spectra,Tags,Led,Block) %>% select(Led,Tags,Block,everything()) %>% separate(Tags,into = c("Board","Case","Cell","Sensor","Slit","PWM"),sep = "-")
  return(output)
}
factorial = map(files,.f = read_spectra) %>% bind_rows()



factorialLong = factorial %>% 
  pivot_longer(cols = starts_with("p",ignore.case = F),names_to = "pixel",values_to = "signal") %>% 
  mutate(pixel = pixel %>% str_remove("p") %>% as.numeric())
factorialLong %>% ggplot(aes(pixel,signal,group = paste(Board,Case,Cell,Sensor,Slit,Led))) + geom_line() + facet_wrap(PWM)
ggplot(factorial,aes(PWM %>% as.numeric(),p50,group = paste(Board,Case,Cell,Sensor,Slit))) + geom_line(alpha = 0.2)

factorialLong = factorial %>% filter(p300<35000,p250 <35000,p240 <35000) %>% 
  pivot_longer(cols = starts_with("p",ignore.case = F),names_to = "pixel",values_to = "signal") %>% 
  mutate(pixel = pixel %>% str_remove("p") %>% as.numeric())
factorialLong %>% ggplot(aes(pixel,signal,group = paste(Board,Case,Cell,Sensor,Slit,PWM))) + geom_line(alpha = 0.05)

tagged = factorial %>% unite(col = "TAG",sep = "-",c("Board","Case","Cell","Sensor","Slit"))
factorialLong = tagged %>% 
  pivot_longer(cols = starts_with("p",ignore.case = F),names_to = "pixel",values_to = "signal") %>% 
  mutate(pixel = pixel %>% str_remove("p") %>% as.numeric())
factorialLong %>% ggplot(aes(pixel,signal,group = paste(PWM))) + geom_line(alpha = 0.05) + facet_wrap(~ TAG)
nested = factorialLong %>% group_by(TAG) %>% nest()

ggplot(nested$data[[190]],aes(pixel,signal,group = PWM)) + geom_line()

factorial %>% select(starts_with("p",ignore.case = F)) %>% summarise_all(.funs = function(x){sd(x)/mean(x)}) %>% unlist() %>% plot
factorial = factorial %>% mutate(PWM= PWM %>% as.numeric())
ggplot(factorial,aes(PWM %>% as.numeric(),p50,group = paste(Board,Case,Cell,Sensor,Slit))) + geom_line(alpha = 0.2)

# easy to find outliers

outliers = bind_rows(
  factorial %>% filter((PWM == 0.1 & p50> 2500)),
  factorial %>% filter((PWM == 0 & p50 >500)),
  factorial %>% filter((PWM == 0.01 & p50 >500)),
  factorial %>% filter((PWM == 0.2 & p50 >3000)),
  factorial %>% filter((PWM == 1 & p50 >17000)),
  factorial %>% filter((PWM == 0.5 & p50 >8000)),
  factorial %>% filter((PWM == 0.9 & p50 >14900)),
  factorial %>% filter((PWM == 0.4 & p50 >6600)),
  factorial %>% filter((PWM == 0.3 & p50 >4850)),
  factorial %>% filter((PWM == 0.6 & p50 >10000)),
  factorial %>% filter((PWM == 0.7 & p50 >11500)),
  factorial %>% filter((PWM == 0.7 & p50  < 50))
)

anti_join(factorial %>% 
            mutate(
              tag = paste(
                Led,
                Board,
                Case,
                Cell,
                Sensor,
                Slit,
                sep = "-")
              ),
          outliers %>%
            mutate(
              tag = paste(
                Led,
                Board,
                Case,
                Cell,
                Sensor,
                Slit,
                sep = "-"
                )
              ),
          by = join_by(tag == tag)) %>%
            ggplot(aes(PWM %>% as.numeric(),p50,group = paste(Board,Case,Cell,Sensor,Slit))) + geom_line(alpha = 0.2) + 
  geom_hline(yintercept = 11500) +
  geom_vline(xintercept = 0.7)

clean = anti_join(factorial %>% 
                    mutate(
                      tag = paste(
                        Led,
                        Board,
                        Case,
                        Cell,
                        Sensor,
                        Slit,
                        sep = "-")
                    ),
                  outliers %>%
                    mutate(
                      tag = paste(
                        Led,
                        Board,
                        Case,
                        Cell,
                        Sensor,
                        Slit,
                        sep = "-"
                      )
                    ),
                  by = join_by(tag == tag)) %>% 
  mutate(
    Board = ifelse(Board == "11","1",Board),
    Board = ifelse(Board == "22","2",Board)) %>% 
  filter(Board != "4")

ggplot(clean) + 
  geom_line(aes(PWM,p50,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) + theme_classic()
ggplot(clean %>% filter(PWM > 0.01)) + 
  geom_line(aes(PWM,p50,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) + theme_classic()

ggplot(clean) + 
  geom_line(aes(PWM,p50,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,p50,colour = Led),method = "lm",se = F) + theme_classic()
ggplot(clean) + 
  geom_line(aes(PWM,p50,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,p50,colour = Board),method = "lm",se = F) + theme_classic()
ggplot(clean) + 
  geom_line(aes(PWM,p50,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,p50,colour = Case),method = "lm",se = F) + theme_classic()
ggplot(clean) + 
  geom_line(aes(PWM,p50,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,p50,colour = Cell),method = "lm",se = F) + theme_classic()
ggplot(clean) + 
  geom_line(aes(PWM,p50,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,p50,colour = Sensor),method = "lm",se = F) + theme_classic()
ggplot(clean) + 
  geom_line(aes(PWM,p50,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,p50,colour = Slit),method = "lm",se = F) + theme_classic()
  
ggplot(clean) + 
  geom_line(aes(PWM,p50, colour = Led,group = paste(Board,Case,Cell,Sensor,Slit)),alpha = 0.3) +
  geom_smooth(aes(PWM,p50,colour = Led),method = "lm",se = F) + theme_classic()
ggplot(clean) + 
  geom_line(aes(PWM,p50,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,p50,colour = Board),method = "lm",se = F) + theme_classic()
ggplot(clean) + 
  geom_line(aes(PWM,p50,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,p50,colour = Case),method = "lm",se = F) + theme_classic()
ggplot(clean) + 
  geom_line(aes(PWM,p50,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,p50,colour = Cell),method = "lm",se = F) + theme_classic()
ggplot(clean) + 
  geom_line(aes(PWM,p50,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,p50,colour = Sensor),method = "lm",se = F) + theme_classic()
ggplot(clean) + 
  geom_line(aes(PWM,p50,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,p50,colour = Slit),method = "lm",se = F) + theme_classic()

factorialLong = clean %>% 
  pivot_longer(cols = starts_with("p",ignore.case = F),names_to = "pixel",values_to = "signal") %>% 
  mutate(pixel = pixel %>% str_remove("p") %>% as.numeric())
factorialLong %>% ggplot(aes(pixel,signal,group = paste(Board,Case,Cell,Sensor,Slit,Led))) + geom_line(alpha = 0.05) + facet_wrap(~ PWM)
factorialLong %>% filter(PWM == 0.2)%>% ggplot(aes(pixel,signal,group = paste(Board,Case,Cell,Sensor,Slit,Led))) + geom_line(alpha = 0.05) 
clean %>% select(starts_with("p")) %>% prcomp() %>% summary()

PC = clean %>% select(starts_with("p")) %>% select(1:75) %>% prcomp(,scale. = T) %>% .$x %>% .[,1] %>% unlist()
reduced = clean %>% select(1:7) %>% mutate("PC" = PC)
ggplot(reduced) + geom_line(aes(PWM,PC,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) + geom_vline(xintercept = 0.3) + geom_hline(yintercept = 0)
ggplot(reduced) + geom_line(aes(PWM,PC,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04)

ggplot(reduced) + 
  geom_line(aes(PWM,PC,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,PC,colour = Board),method = "lm",se = F) + theme_classic()
ggplot(reduced) + 
  geom_line(aes(PWM,PC,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,PC,colour = Case),method = "lm",se = F) + theme_classic()
ggplot(reduced) + 
  geom_line(aes(PWM,PC,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,PC,colour = Cell),method = "lm",se = F) + theme_classic()
ggplot(reduced) + 
  geom_line(aes(PWM,PC,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,PC,colour = Sensor),method = "lm",se = F) + theme_classic()
ggplot(reduced) + 
  geom_line(aes(PWM,PC,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,PC,colour = Slit),method = "lm",se = F) + theme_classic()
ggplot(reduced) + 
  geom_line(aes(PWM,PC,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) +
  geom_smooth(aes(PWM,PC,colour = Led),method = "lm",se = F) + theme_classic()




reduced %>% filter((PWM < 0.3 & PC > 0)) %>% .$Led
## Multivarite
PERMANOVA1 = vegan::adonis2(formula = clean %>% select(starts_with("p")) %>% dist()~ ., data = clean[,1:7],parallel = 31,permutations = 100000)
saveRDS(PERMANOVA1,file = "PERMANOVA1")
PERMANOVA2 = vegan::adonis2(formula = clean %>% select(starts_with("p")) %>% vegan::decostand("hellinger") %>% dist()~ ., data = clean[,1:7],parallel = 31,permutations = 100000)
saveRDS(PERMANOVA1,file = "PERMANOVA2")

library(lme4)
vainilla = lm(p50 ~ PWM+ Led + Board + Case + Cell + Sensor + PWM*Led + PWM*Board + PWM*Case + PWM*Cell + PWM*Sensor,claen)
mixed = lmer(p50 ~ PWM + (1+PWM|Led)+ (1+PWM|Board)+ (1+PWM|Case)+ (1+PWM|Cell)+ (1+PWM|Sensor)+ (1+PWM|Slit),data = claen)
mixed = lmer(p50 ~ PWM + (1 |Led)+ (1 |Board)+ (1 + PWM |Case)+ (1 |Cell)+ (1 |Sensor)+ (1 |Slit),data = claen)

mixed = lmer(p50 ~ PWM + (PWM|Led)+ (PWM|Board)+ (PWM|Case)+ (PWM|Cell)+ (PWM|Sensor)+ (PWM|Slit),data = clean)
mixed = lmer(p50 ~ PWM + (PWM|Led)+ (PWM|Board)+ (PWM|Case)+ (PWM|Cell)+ (PWM|Sensor)+ (PWM|Slit),data = clean)

mixedRest = lmer(p50 ~ (1|Led)+ (1|Board)+ (1|Case)+ (1|Cell)+ (1|Sensor)+ (1|Slit),data = clean %>% filter(PWM == 0.2))
mixedRest = lmer(p50 ~ (1|Led)+ (1|Board)+ (1|Case)+ (1|Cell)+ (1|Sensor)+ (1|Slit) + (1|PWM),data = clean %>% filter(PWM < 0.4))

mixedRest = lmer(p50 ~ (1|Led)+ (1|Board)+ (1|Case)+ (1|Cell)+ (1|Sensor)+ (1|Slit) + (1|PWM),data = clean)
GmixedRest = glmer(p50 ~ (1|Led)+ (1|Board)+ (1|Case)+ (1|Cell)+ (1|Sensor)+ (1|Slit) + (1|PWM),data = clean %>% filter((PWM >0.01 & PWM <0.6)),family = Gamma(link = "log"))
GmixedRest = glmer(MValue ~ (1|Led)+ (1|Board)+ (1|Case)+ (1|Cell)+ (1|Sensor)+ (1|Slit) + (1|PWM),data = clean %>% mutate(MValue = rowMeans(clean %>% select(starts_with("p")) %>% select(1:75))) %>% filter(PWM >0.01),family = Gamma(link = "log"))
GmixedRest %>% plot()

GmixedRest %>% plot()

mixedRest = lmer(PC ~ (1|Led) + (1|Board) + (1|Case)+ (1|Cell)+ (1|Sensor)+ (1|Slit) + (1|PWM),data = reduced)
mixedRest = lmer(PC ~ (1|Led) + (1|Case)+ (1|Cell)+ (1|Sensor)+ (1|Slit) + (1|PWM),data = reduced)

mixedRest = glmer(
  PC ~ (1|Led) + (1|Case)+ (1|Cell)+ (1|Sensor)+ (1|Slit) + (1|PWM),
  data = reduced %>% filter(PWM >0.01) %>% mutate(PC = PC - min(PC) +1)%>% filter(PWM >0.01),
  family = Gamma(link = "log"))

mixed = brm(p50 ~ PWM + (PWM|Led)+ (PWM|Board)+ (PWM|Case)+ (PWM|Cell)+ (PWM|Sensor)+ (PWM|Slit),data = clean)
mixed = brm(p50 ~ PWM + (PWM|Led)+ (PWM|Board)+ (PWM|Case)+ (PWM|Cell)+ (PWM|Sensor)+ (PWM|Slit),data = clean)

factorial %>% filter(!(PWM == 0 & p50 >500))
factorial %>% 
  filter(!(Led == 1 & Board == 4 & Case == 1 & Cell == 4 & Sensor == 2 & Slit == 1)) %>% 
  filter(!(Led == 4 & Board == 1 & Case == 4 & Cell == 1 & Sensor == 1 & Slit == 1)) %>% 
  filter(!(Led == 3 & Board == 2 & Case == 3 & Cell == 2 & Sensor == 3 & Slit == 1))
  ggplot(aes(PWM %>% as.numeric(),p50,group = paste(Board,Case,Cell,Sensor,Slit))) + geom_line(alpha = 0.2)


factorial %>% ggplot(aes(x = p50+PWM)) + geom_histogram() + scale_x_log10()

factorial %>% filter(PWM + p50> 1000) %>% View()
factorial %>% filter(PWM + p50> 100) %>% select(starts_with("p")) %>% .[1,] %>% unlist() %>% plot
factorial 


# Load necessary packages
library(brms)
library(rstan)

# Hypothetical data
df <- data.frame(
  light_intensity = rep(seq(0, 100, 10), times = 60),
  signal_output = rnorm(600, mean = 50, sd = 10),
  part1 = rep(letters[1:6], each = 100),
  part2 = rep(LETTERS[1:6], each = 100),
  part3 = rep(1:6, each = 100),
  part4 = rep(c("A", "B", "C", "D", "E", "F"), each = 100),
  part5 = rep(c("x", "y", "z", "u", "v", "w"), each = 100),
  part6 = rep(c(10, 20, 30, 40, 50, 60), each = 100)
)

# Specify the model formula
model_formula <- bf(p50 ~ PWM + (PWM|Led)+ (PWM|Board)+ (PWM|Case)+ (PWM|Cell)+ (PWM|Sensor)+ (PWM|Slit))

# Fit the hierarchical Bayesian model
fit <- brm(
  formula = model_formula,
  data = clean,
  family = gaussian(),
  prior = c(
    prior(normal(0, 10), class = "b"),
    prior(cauchy(0, 2), class = "sd")
  ),
  iter = 2000, chains = 4, cores = 4
)

# Summarize the model
summary(fit)

# Plot the posterior distributions
plot(fit)
saveRDS(object = fit,"Characterisation/Data/factorialExperiment/Results/bayesian")
