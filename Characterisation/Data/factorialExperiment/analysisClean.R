library(tidyverse)
library(lme4)
library(patchwork)
library(treemapify)
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

# Load data ---------------------------------------------------------------


files = list.files("Characterisation/Data/factorialExperiment/Results/processed/",full.names = T) %>% str_subset(".json") 

factorial = purrr::map(files,.f = read_spectra) %>% bind_rows()


factorialLong = factorial %>% 
  pivot_longer(cols = starts_with("p",ignore.case = F),names_to = "pixel",values_to = "signal") %>% 
  mutate(pixel = pixel %>% str_remove("p") %>% as.numeric())

factorialLong %>% ggplot(aes(pixel,signal,group = paste(Board,Case,Cell,Sensor,Slit,Led))) +
  geom_line(alpha = 0.01) + 
  facet_wrap(~PWM) +
  theme_classic()

factorialLong %>% filter(pixel < 75) %>% ggplot(aes(pixel,signal,group = paste(Board,Case,Cell,Sensor,Slit,Led))) +
  geom_line(alpha = 0.01) + 
  facet_wrap(~PWM) +
  theme_classic()

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

# Variable Selection ------------------------------------------------------

selected = clean %>%
  filter(PWM >0.01) %>%
  .[1:7] %>% 
  mutate("Sum" = clean %>% filter(PWM >0.01) %>% select(starts_with("p",ignore.case = FALSE)) %>% .[1:75] %>% rowSums())


ggplot(selected) + 
  geom_line(aes(PWM,Sum,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04) + 
  theme_classic()
#Modifying
selected = selected %>%
  mutate("SumScaled" = Sum/max(Sum),
         PWM = PWM %>% as.numeric())


ggplot(selected) + 
  geom_line(aes(
    PWM,
    SumScaled,
    group = paste(
      Board,
      Case,
      Cell,
      Sensor,
      Slit
      ),
    colour = Led),
    alpha = 0.04) +
  theme_classic()  + theme(legend.position = "none",axis.title = element_blank()) +
  
  ggplot(selected %>% filter(PWM == 1)) +
  geom_density(aes(y = SumScaled,fill = Led),alpha = 0.4) + theme_void() + plot_layout(design = paste0(str_dup("A",9),"B"))

# PWMLED = ggplot(selected) + 
#   geom_line(aes(
#     PWM,
#     SumScaled,
#     group = paste(
#       Board,Case,Cell,Sensor,Slit,Led
#       )
#     ),alpha = 0.04,linewidth = 0.3) +
#   geom_line(
#     aes(
#       PWM,
#       SumScaled,
#       colour = Led
#       ),stat = "smooth",
#     method = "lm",
#     se = F,
#     formula = y ~ x, linewidth = 1.6
#     ) + theme_classic() + theme(legend.position = "none",axis.title = element_blank())  +
#   scale_colour_manual(values = palette4[2:5]) + labs(tag = "d)")
#   
# PWMLED_DENS =  ggplot(selected %>% filter(PWM == 1)) +
#   geom_density(aes(y = SumScaled,fill = Led),alpha = 0.4) + theme_void()   +
#   scale_colour_manual(values = palette4[2:5])+
#   scale_colour_manual(values = palette4[2:5])


# Model -------------------------------------------------------------------
library(rethinking)
# Generative Exploration
id = 1:3
a = vector()
a[id] = dnorm(0,3)

size = 1e5
bSe = rlnorm(n = size,meanlog = log(1/sqrt(1+1)),sdlog = sqrt(log(1+1/1)))
bSl=rlnorm(n = size,meanlog = log(1/sqrt(1+1)),sdlog = sqrt(log(1+1/1)))
bCe=rlnorm(n = size,meanlog = log(1/sqrt(1+1)),sdlog = sqrt(log(1+1/1)))
bCa= rlnorm(n = size,meanlog = log(1/sqrt(1+1)),sdlog = sqrt(log(1+1/1)))
bLight= rlnorm(n = size,log(0.5/sqrt(2)),log(1.25)) #rbeta(n = size,shape1 = 2,shape2 = 2) 

bSe = rbeta(n = size,shape1 = 2,shape2 = 2) 
bSl=rbeta(n = size,shape1 = 2,shape2 = 2) 
bCe=rbeta(n = size,shape1 = 2,shape2 = 2) 
bCa= rbeta(n = size,shape1 = 2,shape2 = 2) 
bLight= rbeta(n = size,shape1 = 2,shape2 = 2) 


mu = (bSe*bSl*bCe*bCa*bLight %*% t(seq(0.1,1,by=0.1))) %>%
  as_tibble() %>% mutate(ID = 1:length(bSe)) %>% 
  pivot_longer(starts_with("V")) %>% 
  mutate(name = rep(seq(0.1,1,by=0.1),size)) 

sigma = runif(length(mu$value),min = 0,max = 0.0001)#rexp(length(mu$value),10)

y = map_vec(mu$value,rnorm, n = 1,sd = 0.02)

mu$signal = y

ggplot(mu,aes(name,value,group = ID)) + geom_line()


# Model Try

BayesianModel = alist(
  SumScaled ~ dnorm(mu,sigma),
  mu <- bBoard[Board]*bCase[Case]*bCell[Cell]*bSensor[Sensor]*bSlit[Slit]*bLight*PWM,
  bBoard[Board] ~ dlnorm(meanlog = log(1.0/sqrt(1+1)),sdlog = sqrt(log(1+1.0/1))),
  bCase[Case] ~ dlnorm(meanlog = log(1.0/sqrt(1+1)),sdlog = sqrt(log(1+1.0/1))),
  bCell[Cell] ~ dlnorm(meanlog = log(1.0/sqrt(1+1)),sdlog = sqrt(log(1+1.0/1))),
  bSensor[Sensor] ~ dlnorm(meanlog = log(1.0/sqrt(1+1)),sdlog = sqrt(log(1+1.0/1))),
  bSlit[Slit] ~ dlnorm(meanlog = log(1.0/sqrt(1+1)),sdlog = sqrt(log(1+1.0/1))),
  bLight ~ dlnorm(meanlog = log(1.0/sqrt(1+1)),sdlog = sqrt(log(1+1.0/1))),
  sigma ~ dlnorm(0,1)
)
# BayesianModel =  alist(
#   
#   Sum ~ dnorm(mu,sigma),
#   mu <- bBoard[Board]*bCase[Case]*bCell[Cell]*bSensor[Sensor]*bSlit[Slit]*bLight*PWM,
#   bBoard[Board] ~ dlnorm(0, 1),
#   bCase[Case] ~ dlnorm(0, 1),
#   bCell[Cell] ~ dlnorm(0, 1),
#   bSensor[Sensor] ~ dlnorm(0, 1),
#   bSlit[Slit] ~ dlnorm(0, 1),
#   bLight ~ dlnorm(0, 1),
#   sigma ~ dexp(1)  # Precision parameter
# )

# modelQuap =  rethinking::quap(BayesianModel,
#                               data = selected %>% mutate_all(as.numeric),
#                               control = list(iter = 1e10))


modelTrain = rethinking::ulam(BayesianModel,
                              data = selected %>% mutate_all(as.numeric),
                              cores = 32,
                              chains = 64,
                               start = list(
                                 bBoard = rep(1, nlevels(selected$Board %>% as.factor())),
                                 bCase = rep(1, nlevels(selected$Case %>% as.factor())),
                                 bCell = rep(1, nlevels(selected$Cell %>% as.factor())),
                                 bSensor = rep(1, nlevels(selected$Sensor %>% as.factor())),
                                 bSlit = rep(1, nlevels(selected$Slit %>% as.factor())),
                                 bLight = 1,
                                 sigma = 1
                               ),iter = 10000,max_treedepth = 15,
                               control = list(adapt_delta = 0.99))

saveRDS(modelTrain,file = "Characterisation/Data/factorialExperiment/BayesianModel")
library(rethinking)
post <- extract.samples(modelTrain)

summary(post)



# Plots data --------------------------------------------------------------

PWMBOARD = ggplot(selected) + 
  geom_line(aes(PWM,SumScaled,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04,linewidth = 0.3) +
  geom_line(aes(PWM,SumScaled,colour = Board),stat = "smooth",method = "lm",se = F, linewidth = 1.6) +
  theme_classic() + 
  # theme(legend.position = "none",
  #       axis.title = element_blank())  +
  scale_colour_manual(values = palette4[2:5]) + labs(tag = "a)") + 
  scale_x_continuous(breaks = c(0,0.5,1)) + 
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme(legend.position = "inside",
        legend.box = "inside", 
        legend.position.inside = c(0.05,0.9),
        legend.direction = "horizontal",
        legend.title.position = "left",
        legend.title.align = 0.5,
        legend.spacing.y = unit(-1, "cm"),
        legend.box.background = element_rect(
          fill = NA,colour = NA),
        legend.justification = "left"
  ) +
  xlab("Light intensity")+
  ylab("Signal (AU)")


PWMBOARD_DENS = ggplot(selected %>% filter(PWM == 1)) +
  geom_density(aes(y = SumScaled,fill = Board),alpha = 0.4) + theme_void() +
  scale_fill_manual(values = palette4[2:5])+
  theme(legend.position = "none",
        axis.title = element_blank())


PWMCASE = ggplot(selected) + 
  geom_line(aes(PWM,SumScaled,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04,linewidth = 0.3) +
  geom_line(aes(PWM,SumScaled,colour = Case),stat = "smooth", linewidth = 1.6,method = "lm",se = F) + 
  theme_classic() + 
  # theme(legend.position = "none",axis.title = element_blank())  +
  scale_colour_manual(values = palette4[2:5]) + labs(tag = "b)") +
  scale_x_continuous(breaks = c(0,0.5,1)) + 
  scale_y_continuous(breaks = c(0,0.5,1))  +
  theme(legend.position = "inside",
        legend.box = "inside", 
        legend.position.inside = c(0.05,0.9),
        legend.direction = "horizontal",
        legend.title.position = "left",
        legend.title.align = 0.5,
        legend.spacing.y = unit(-1, "cm"),
        legend.box.background = element_rect(
          fill = NA,colour = NA),
        legend.justification = "left"
  ) +
  xlab("Light intensity")+
  ylab("Signal (AU)")

PWMCASE_DENS =  ggplot(selected %>% filter(PWM == 1)) +
  geom_density(aes(y = SumScaled,fill = Case),alpha = 0.4,linewidth = 0.3) + theme_void() +
  scale_fill_manual(values = palette4[2:5]) +
  theme(legend.position = "none",
        axis.title = element_blank())

PWMCELL = ggplot(selected) + 
  geom_line(aes(PWM,SumScaled,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04,linewidth = 0.3) +
  geom_line(aes(PWM,SumScaled,colour = Cell),stat = "smooth",method = "lm",se = F, linewidth = 1.6) + 
  theme_classic() +
  # theme(legend.position = "none",axis.title.x = element_blank())  +
  scale_colour_manual(values = palette4[2:5]) + ylab("Scaled Signal (AU)")  +
  labs(tag = "c)") +
  scale_x_continuous(breaks = c(0,0.5,1)) + scale_y_continuous(breaks = c(0,0.5,1)) +
  theme(legend.position = "inside",
        legend.box = "inside", 
        legend.position.inside = c(0.05,0.9),
        legend.direction = "horizontal",
        legend.title.position = "left",
        legend.title.align = 0.5,
        legend.spacing.y = unit(-1, "cm"),
        legend.box.background = element_rect(
          fill = NA,colour = NA),
        legend.justification = "left"
  ) +
  xlab("Light intensity")+
  ylab("Signal (AU)")

PWMCELL_DENS = ggplot(selected %>% filter(PWM == 1)) +
  geom_density(aes(y = SumScaled,fill = Cell),alpha = 0.4) + theme_void() + plot_layout(design = paste0(str_dup("A",9),"B"))  +
  scale_fill_manual(values = palette4[2:5]) +
  theme(legend.position = "none",
        axis.title = element_blank())

PWMSENSOR = ggplot(selected) + 
  geom_line(aes(PWM,SumScaled,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04,linewidth = 0.3) +
  geom_line(aes(PWM,SumScaled,colour = Sensor),stat = "smooth",method = "lm",se = F, linewidth = 1.6) +
  theme_classic() + 
  # theme(legend.position = "none",axis.title.y = element_blank())  +
  scale_colour_manual(values = palette4[2:5]) + labs(tag = "e)") +
  scale_x_continuous(breaks = c(0,0.5,1)) + 
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme(legend.position = "inside",
        legend.box = "inside", 
        legend.position.inside = c(0.05,0.9),
        legend.direction = "horizontal",
        legend.title.position = "left",
        legend.title.align = 0.5,
        legend.spacing.y = unit(-1, "cm"),
        legend.box.background = element_rect(
          fill = NA,colour = NA),
        legend.justification = "left"
  ) +
  xlab("Light intensity")+
  ylab("Signal (AU)")

PWMSENSOR_DENS = ggplot(selected %>% filter(PWM == 1)) +
  geom_density(aes(y = SumScaled,fill = Sensor),alpha = 0.4) + theme_void() + plot_layout(design = paste0(str_dup("A",9),"B"))  +
  scale_fill_manual(values = palette4[2:5]) +
  theme(legend.position = "none",
        axis.title = element_blank())

PWMSLIT = ggplot(selected) + 
  geom_line(aes(PWM,SumScaled,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04,linewidth = 0.3) +
  geom_smooth(aes(PWM,SumScaled,colour = Slit),stat = "smooth",method = "lm",se = F, linewidth = 1.6) + 
  theme_classic()  + 
  # theme(legend.position = "none",axis.title.y = element_blank())  +
  scale_colour_manual(values = palette4[2:5]) + labs(tag = "d)") + 
  scale_x_continuous(breaks = c(0,0.5,1)) + 
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme(legend.position = "inside",
        legend.box = "inside", 
        legend.position.inside = c(0.05,0.9),
        legend.direction = "horizontal",
        legend.title.position = "left",
        legend.title.align = 0.5,
        legend.spacing.y = unit(-1, "cm"),
        legend.box.background = element_rect(
          fill = NA,colour = NA),
        legend.justification = "left"
  ) +
  xlab("Light intensity")+
  ylab("Signal (AU)")

PWMSLIT_DENS = ggplot(selected %>% filter(PWM == 1)) +
  geom_density(aes(y = SumScaled,fill = Slit),alpha = 0.4) + theme_void() + plot_layout(design = paste0(str_dup("A",9),"B"))  +
  scale_fill_manual(values = palette4[2:5]) +
  theme(legend.position = "none",
        axis.title = element_blank())


# plots effects model -----------------------------------------------------


# Convert samples to a data frame
df <- as.data.frame(post$bCase)  # Replace `bBoard` with your parameter of interest

# Reshape for ggplot (if needed for multiple columns)
colnames(df) <- paste0("Slit", 1:ncol(df))  # Name each column for boards
df_long <- tidyr::pivot_longer(df, cols = everything(), names_to = "Slit", values_to = "Posterior")

# Plot densities
ggplot(df_long, aes(x = Posterior, fill = Slit)) +
  geom_density(alpha = 0.5) +
  labs(
    # title = "Posterior Densities for bBoard",
       x = "Parameter Value",
       y = "Density") +
  theme_minimal()

PWMBOARD_DENS = ggplot(df_long, aes(y = Posterior, fill = Slit)) +
  geom_density(alpha = 0.4) + theme_void() +
  scale_colour_manual(values = palette4[2:5])
# PWMSLIT + PWMBOARD_DENS

# Combine posterior samples into one long data frame
params <- list(
  bBoard = as.data.frame(post$bBoard),
  bCase = as.data.frame(post$bCase),
  bCell = as.data.frame(post$bCell),
  bLight = as.data.frame(post$bSensor),
  bLight = as.data.frame(post$bLight),
  bSlit = as.data.frame(post$bSlit),
  bSensor = as.data.frame(post$bSensor)
)
df_all <- purrr::map_dfr(
  params, ~ tidyr::pivot_longer(
    .x, cols = everything(),
    names_to = "Component",
    values_to = "Posterior"),
  .id = "Parameter")

# Plot
effectDens = ggplot(df_all, aes(x = Posterior, fill = Parameter)) +
  geom_density(alpha = 0.4) +
  # facet_wrap(~ Parameter, scales = "free") +
  labs(x = "Parameter Value",
       y = "Density") +
  theme_classic() + 
  xlim(c(0,6)) + 
  # scale_x_continuous(breaks = c(1,4,7))+
  scale_fill_manual(values = c(palette4, palette2[5],palette3[1])) + 
  theme(legend.position.inside = c(0.6,0.9),
        legend.position = "inside",
        legend.box = "inside",
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.title.align = 0.5,
        legend.spacing.y = unit(-1, "cm"),
        legend.box.background = element_rect(
          fill = NA,colour = NA),
        legend.justification = "centre"
  ) + 
  labs(tag = "f)")

# Variance exploration
# Compute variances for each component
variances <- list(
  bBoard = apply(post$bBoard, 2, var),  # Variance of each board
  bCase = apply(post$bCase, 2, var),    # Variance of each case
  bCell = apply(post$bCell, 2, var),    # Variance of each cell
  bSensor = apply(post$bSensor, 2, var),# Variance of each sensor
  bSlit = apply(post$bSlit, 2, var),      # Variance of each slit
  bLight = apply(post$bLight, 2, var)# Variance of PWM (light intensity)
)



# Compute the mean variance across subcomponents (if hierarchical)
mean_variances <- sapply(variances, mean)
# Add residual variance
residual_variance <- var(post$sigma)

# Update total variance
total_variance <- sum(mean_variances) + residual_variance

# Recompute relative variances including residual
relative_variances <- c(mean_variances, Residual = residual_variance) / total_variance
# 
# # Total variance across all components
# total_variance <- sum(mean_variances)
# 
# # Relative variances as proportions
# relative_variances <- mean_variances / total_variance
# 
# Convert to percentages (optional)
relative_variances_percent <- relative_variances * 100
# Create a data frame for plotting
df_var <- data.frame(
  Component = names(relative_variances),
  RelativeVariance = relative_variances_percent
)

BarVariance = ggplot(df_var%>% filter(Component != "Residual"), 
       aes(x = reorder(Component, -RelativeVariance), 
           y = RelativeVariance,
           fill = Component)) +
  geom_bar(stat = "identity") +
  labs(x = "Component",
       y = "Variance Contribution (%)") +
  theme_classic()  +
  scale_fill_manual(values = c(palette4, palette2[5],palette3[1])) + 
  geom_text(aes(x = reorder(Component, -RelativeVariance),
                 y = -1,label = RelativeVariance %>% round(2) %>% paste("%"))) +
  theme(legend.position = "none")  + 
  ylim(c(-1,55)) + 
  labs(tag = "g)")
        
variances2 <- list(
  bBoard = apply(post$bBoard, 2, var),  # Variance of each board
  bCase = apply(post$bCase, 2, var),    # Variance of each case
  bCell = apply(post$bCell, 2, var),    # Variance of each cell
  bSensor = apply(post$bSensor, 2, var),# Variance of each sensor
  # bLight = apply(post$bLight, 2, var),# Variance of PWM (light intensity)
  bSlit = apply(post$bSlit, 2, var)      # Variance of each slit
)



# Compute the mean variance across subcomponents (if hierarchical)
mean_variances2 <- sapply(variances2, mean)
# Add residual variance

# Update total variance
total_variance2 <- sum(mean_variances2) + residual_variance

# Recompute relative variances including residual
relative_variances2 <- c(mean_variances2, Residual = residual_variance) / total_variance2
# 
# # Total variance across all components
# total_variance <- sum(mean_variances)
# 
# # Relative variances as proportions
# relative_variances <- mean_variances / total_variance
# 
# Convert to percentages (optional)
relative_variances_percent2 <- relative_variances2 * 100
# Create a data frame for plotting
df_var2 <- data.frame(
  Component = names(relative_variances2),
  RelativeVariance = relative_variances_percent2
)

pieVariance2 = ggplot(
  df_var2 %>% 
    filter(Component != "Residual"),
  aes(x = "", 
      y = RelativeVariance,
      fill = Component)
  ) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "right")  +
  scale_fill_manual(values = c(palette4[-4], palette2[5],palette3[1])) + 
  theme_void() +
  theme(
    legend.position = "none"
  )

triangle = ggplot(tibble(x = c(0,2,1),
                         y = c(0,0,0.3))) + 
  # ylim(c(0,1)) +
  geom_polygon(aes(x = x,y = y),fill = "lightgrey") + 
  theme_void() + 
  scale_fill_discrete()+ 
  theme(legend.position = "none")
triangle = ggplot() + 
  geom_segment(aes(x = 0,y = 0.9915,xend = 10,yend = 0.9915 ),linewidth = 3,colour = "gray")+ 
  # geom_segment(aes(x = 0,y = 1,xend = 0,yend = -0.2 ),linewidth = 3)+ 
  # geom_segment(aes(x = 10,y = 1,xend = 10,yend = -0.2 ),linewidth = 3)+ 
  geom_segment(aes(x = 5,y = 1,xend = 5,yend = 3),
               arrow = arrow(length = unit(10,units = "point")),
               linewidth = 1.3) + 
  theme_void()

plotVars = BarVariance + 
  pieVariance2 %>% inset_element(left = 0.15,bottom =0.3 ,right = 1 ,top = 1) +
  triangle %>% inset_element(left = 0.15,bottom =0.25 ,right = 1 ,top = 0.35)
plotVars

layout = paste0(str_dup("A",7),"B",str_dup("C",7),"D","\n",
                str_dup("E",7),"F",str_dup("G",7),"H","\n",
                str_dup("I",7),"J",str_dup("K",8),"\n",
                str_dup("L",8),str_dup("K",8),"\n",
                str_dup("L",8),str_dup("K",8),"\n"
                )

PWMBOARD + 
  PWMBOARD_DENS + 
  PWMCASE + 
  PWMCASE_DENS + 
  PWMCELL +
  PWMCELL_DENS +
  PWMSLIT +
  PWMSLIT_DENS +
  PWMSENSOR +
  PWMSENSOR_DENS + 
  plotVars +
  effectDens +
  plot_layout(design = layout)

ggplot(df_var2) +
  geom_treemap_text(aes(
    area = RelativeVariance,
    fill = Component,
    label = relative_variances)
    ) +
  theme_minimal()  +
  scale_fill_manual(values = c(palette4, palette2[5],palette3[1]))



# test --------------------------------------------------------------------
(post$bBoard[,1] * post$bCase[,1] * post$bCell[,1] *post$bSensor[,1] *post$bSlit[,1] *post$bLight[,1]) %>% hist()

# comparing effect sizes --------------------------------------------------

# Compute posterior means for each beta
effect_sizes <- list(
  bBoard = apply(post$bBoard, 2, mean),
  bCase = apply(post$bCase, 2, mean),
  bCell = apply(post$bCell, 2, mean),
  bSensor = apply(post$bSensor, 2, mean),
  bSlit = apply(post$bSlit, 2, mean),
  bLight = mean(post$bLight)  # Single value
)

# Calculate average effect size per component (if hierarchical)
mean_effect_sizes <- sapply(effect_sizes, mean)

# Normalize effect sizes to sum to 1
relative_effect_sizes <- mean_effect_sizes / sum(mean_effect_sizes)

# Convert to percentages (optional)
relative_effect_sizes_percent <- relative_effect_sizes * 100
# Data frame for plotting
df_effect <- data.frame(
  Component = names(relative_effect_sizes),
  RelativeEffect = relative_effect_sizes_percent
)

# Plot with ggplot2
library(ggplot2)
ggplot(df_effect, aes(x = reorder(Component, -RelativeEffect), y = RelativeEffect)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Relative Effect Sizes of Components",
       x = "Component",
       y = "Effect Size Contribution (%)") +
  theme_minimal()
# Combine posterior samples into a data frame for plotting
params <- list(
  bBoard = as.data.frame(post$bBoard),
  bCase = as.data.frame(post$bCase),
  bCell = as.data.frame(post$bCell),
  bSensor = as.data.frame(post$bSensor),
  bSlit = as.data.frame(post$bSlit),
  bLight = as.data.frame(post$bLight)  # Convert scalar to single-column DF
)
df_beta <- purrr::map_dfr(params, ~ tidyr::pivot_longer(.x, cols = everything(), names_to = "Subcomponent", values_to = "Beta"), .id = "Component")

# Plot density
ggplot(df_beta, aes(x = Beta, fill = Component)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Component, scales = "free") +
  labs(title = "Posterior Densities of Effect Sizes",
       x = "Beta Value",
       y = "Density") +
  theme_minimal()



# Experimental Plot -------------------------------------------------------

# Define range for PWM
PWM_range <- seq(0, 1, length.out = 100)

# Create new data frame for predictions (one component at a time)
new_data <- expand.grid(
  PWM = PWM_range,
  Board = unique(selected$Board),
  Case = unique(selected$Case),
  Cell = unique(selected$Cell),
  Sensor = unique(selected$Sensor),
  Slit = unique(selected$Slit),
  Led = unique(selected$Led)
)

# Posterior predictions
posterior_predictions <- rethinking::link(
  BayesianModel,
  data = new_data
)

# Compute means and credible intervals
pred_means <- apply(posterior_predictions, 2, mean)
pred_ci <- apply(posterior_predictions, 2, quantile, probs = c(0.05, 0.95))

new_data$Mean <- pred_means
new_data$Lower <- pred_ci[1, ]
new_data$Upper <- pred_ci[2, ]
PWMBOARD <- ggplot(selected) +
  geom_line(aes(PWM, SumScaled, group = paste(Board, Case, Cell, Sensor, Slit, Led)),
            alpha = 0.04, linewidth = 0.3) +
  geom_line(data = subset(new_data, Board == "1"),  # Filter specific Board
            aes(x = PWM, y = Mean, colour = Board), linewidth = 1.6) +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_blank()) +
  scale_colour_manual(values = palette4[2:5]) +
  labs(tag = "a)") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1))
# Combine all beta parameters into a single data frame
beta_df <- data.frame(
  Parameter = c(rep("bBoard", ncol(post$bBoard)),
                rep("bCase", ncol(post$bCase)),
                rep("bCell", ncol(post$bCell)),
                rep("bSensor", ncol(post$bSensor)),
                rep("bSlit", ncol(post$bSlit)),
                "bLight"),
  Value = c(as.vector(post$bBoard),
            as.vector(post$bCase),
            as.vector(post$bCell),
            as.vector(post$bSensor),
            as.vector(post$bSlit),
            post$bLight)
)

# Plot
beta_density <- ggplot(beta_df, aes(x = Value, fill = Parameter)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Density of Parameters", x = "Parameter Value", y = "Density")
# Calculate variances of each component
variances <- sapply(list(
  bBoard = apply(post$bBoard, 2, var),
  bCase = apply(post$bCase, 2, var),
  bCell = apply(post$bCell, 2, var),
  bSensor = apply(post$bSensor, 2, var),
  bSlit = apply(post$bSlit, 2, var)
), mean)  # Average across subcomponents

# Normalize variances
relative_variances <- variances / sum(variances)

# Data frame for plotting
variance_df <- data.frame(
  Component = names(relative_variances),
  Variance = relative_variances
)

# Plot
variance_plot <- ggplot(variance_df, aes(x = reorder(Component, -Variance), y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Relative Variances", x = "Component", y = "Relative Variance") +
  theme_minimal()
# Calculate mean effect size
effect_sizes <- sapply(list(
  bBoard = apply(post$bBoard, 2, mean),
  bCase = apply(post$bCase, 2, mean),
  bCell = apply(post$bCell, 2, mean),
  bSensor = apply(post$bSensor, 2, mean),
  bSlit = apply(post$bSlit, 2, mean),
  bLight = mean(post$bLight)
), mean)

# Normalize effect sizes
relative_effect_sizes <- effect_sizes / sum(effect_sizes)

# Data frame for plotting
effect_df <- data.frame(
  Component = names(relative_effect_sizes),
  EffectSize = relative_effect_sizes
)

# Plot
effect_size_plot <- ggplot(effect_df, aes(x = reorder(Component, -EffectSize), y = EffectSize)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  labs(title = "Relative Effect Sizes", x = "Component", y = "Relative Effect Size") +
  theme_minimal()
library(patchwork)

final_plot <- PWMBOARD + PWMCASE + PWMCELL + PWMSLIT + PWMSENSOR +
  beta_density + variance_plot + effect_size_plot +
  plot_layout(ncol = 4)
print(final_plot)


# Other models ------------------------------------------------------------


GmixedRest = glmer(
  SumScaled ~ (1|Led)+ (1|Board)+ (1|Case)+ (1|Cell)+ (1|Sensor)+ (1|Slit) + (1|PWM),
  data = selected ,family = Gamma(link = "log"))

GmixedSum = glmer(
  SumScaled ~ (1|Led) + (1|Case)+ (1|Cell)+ (1|Sensor)+ (1|Slit) + (1|PWM),
  data = selected ,family = Gamma(link = "inverse"))

GmixedSum = glmer(
  SumScaled ~ (1|Led) + (1|Case)+ (1|Cell)+ (1|Sensor)+ (1|Slit)+ (1|PWM),
  data = selected ,family = Gamma(link = "log"))

GmixedSum %>% plot

GmixedSum %>% summary()

GmixedSum = glmer(
  SumScaled ~ PWM + (PWM|Case)+ (PWM|Cell)+ (PWM|Sensor)+ (PWM|Slit) ,
  data = selected,
  family = gaussian())
GmixedSum = glmer(
  SumScaled ~ PWM + (PWM|Case)+ (PWM|Cell)+ (PWM|Sensor)+ (PWM|Slit) ,
  data = selected,
  family = gaussian())


GmixedSum %>% plot
GmixedSum %>% summary()
GmixedSum %>% resid() %>% qqnorm()
GmixedSum %>% resid() %>% qqline()
GmixedSum %>% resid() %>% shapiro.test()
GmixedSum %>% resid() %>% hist()

GmixedSum = glmer(
  SumScaled ~ logPWM+ (1|Led) + (1|Case)+ (1|Cell)+ (1|Sensor)+ (1|Slit),
  data = selected %>% mutate(logPWM = log(PWM) %>% scale()) ,
  family = gaussian(link = "log"))

GmixedSum = glmmTMB::glmmTMB(SumScaled ~ PWM + (1|Case)+ (1|Cell)+ (1|Sensor)+ (1|Slit), 
                             family = gaussian(), 
                             dispformula = ~ PWM, 
                             data = selected)# %>% plot

GmixedSum %>% resid() %>% shapiro.test()

plot(GmixedSum %>% resid()~GmixedSum %>% predict())


GmixedSum %>% VarCorr() %>% as_tibble() %>% mutate("ICC" = vcov/sum(vcov))
selected[-(GmixedSum %>% residuals(type = "pearson") < -2),]

GmixedSum %>% residuals(type = "pearson") %>% qqnorm()

outliers = selected[(GmixedSum %>% residuals(type = "pearson") < -2),]
outliers = inner_join(selected%>%
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
            ),outliers %>%
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
          by = join_by(tag == tag))



ggplot(selected) + 
  geom_line(aes(PWM,SumScaled,group = paste(Board,Case,Cell,Sensor,Slit,Led)),alpha = 0.04)+ 
  geom_line(data = outliers,
            aes(PWM.x,SumScaled.x,group = paste(Case.x,Cell.x,Slit.x,Led.x))) + 
  theme_classic()
outliers = selected[(GmixedSum %>% residuals(type = "pearson") < -2),]
selected2 = anti_join(selected%>%
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
                         ),outliers %>%
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
                       by = join_by(tag == tag))
GmixedSum2 = lmer(
  log(SumScaled) ~ (1|Led) + (1|Case)+ (1|Cell)+ (1|Sensor)+ (1|Slit) + (1|PWM),
  data = selected2 ,REML = T)
GmixedSum2 %>% sjPlot::plot_model(type = "re")
GmixedSum2 %>% sjPlot::plot_model(type = "diag")

GmixedSum %>% summary()
GmixedSum %>% VarCorr() %>% as_tibble() %>% mutate("ICC" = vcov/sum(vcov),"Percentage of Variation" = ICC * 100)


# STAN model --------------------------------------------------------------

library(rstan)

# Define the Stan model
stan_model_code <- "
data {
  int<lower=1> N;                   // Number of observations
  int<lower=1> J_board;             // Number of boards
  int<lower=1> J_case;              // Number of cases
  int<lower=1> J_cell;              // Number of cells
  int<lower=1> J_sensor;            // Number of sensors
  int<lower=1> J_slit;              // Number of slits

  int<lower=1, upper=J_board> board[N];   // Board index for each observation
  int<lower=1, upper=J_case> case[N];     // Case index for each observation
  int<lower=1, upper=J_cell> cell[N];     // Cell index for each observation
  int<lower=1, upper=J_sensor> sensor[N]; // Sensor index for each observation
  int<lower=1, upper=J_slit> slit[N];     // Slit index for each observation
  real<lower=0> pwm[N];                   // PWM input (scaled between 0.01 and 1)
  real<lower=0> y[N];                     // Scaled response (SumScaled)
}

parameters {
  vector<lower=0>[J_board] b_board;   // Board effects
  vector<lower=0>[J_case] b_case;     // Case effects
  vector<lower=0>[J_cell] b_cell;     // Cell effects
  vector<lower=0>[J_sensor] b_sensor; // Sensor effects
  vector<lower=0>[J_slit] b_slit;     // Slit effects
  real<lower=0> b_light;              // Light effect
  real<lower=0> sigma;                // Standard deviation
}

model {
  // Priors
  b_board ~ lognormal(log(0.5 / sqrt(2)), log(1.25));
  b_case ~ lognormal(log(0.5 / sqrt(2)), log(1.25));
  b_cell ~ lognormal(log(0.5 / sqrt(2)), log(1.25));
  b_sensor ~ lognormal(log(0.5 / sqrt(2)), log(1.25));
  b_slit ~ lognormal(log(0.5 / sqrt(2)), log(1.25));
  b_light ~ lognormal(log(0.5 / sqrt(2)), log(1.25));
  sigma ~ exponential(1);

  // Likelihood
  for (n in 1:N) {
    y[n] ~ normal(
      b_board[board[n]] *
      b_case[case[n]] *
      b_cell[cell[n]] *
      b_sensor[sensor[n]] *
      b_slit[slit[n]] *
      b_light *
      pwm[n],
      sigma
    );
  }
}
"
# Compile the model
stan_model <- stan_model(model_code = stan_model_code)

# Prepare the data for Stan
stan_data <- list(
  N = nrow(selected),                      # Number of observations
  J_board = length(unique(selected$Board)),# Number of boards
  J_case = length(unique(selected$Case)),  # Number of cases
  J_cell = length(unique(selected$Cell)),  # Number of cells
  J_sensor = length(unique(selected$Sensor)),# Number of sensors
  J_slit = length(unique(selected$Slit)),  # Number of slits
  board = as.numeric(selected$Board),      # Board indices
  case = as.numeric(selected$Case),        # Case indices
  cell = as.numeric(selected$Cell),        # Cell indices
  sensor = as.numeric(selected$Sensor),    # Sensor indices
  slit = as.numeric(selected$Slit),        # Slit indices
  pwm = selected$PWM,                      # PWM inputs
  y = selected$SumScaled                   # Scaled response
)

# Fit the model
fit <- sampling(
  object = stan_model,
  data = stan_data,
  iter = 2000,
  chains = 4,
  cores = 4,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

