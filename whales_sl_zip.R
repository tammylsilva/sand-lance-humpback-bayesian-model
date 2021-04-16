
#Bayesian model
#Predict humpback whale counts based on sand lance counts from Seaboss data
#Tammy Silva
#10/18/2018

#load packages
library(rjags)
library(coda) 
library(MCMCvis) 
library(tidyverse) 
library(tidybayes)


###Load and format data

#read in and format visual observations data
vis_obs<-read.csv("./data/seaboss/vis_obs_clean_july2020.csv", header=TRUE)

#pull out just sl, great shearwater and humpback whale observations
vis_obs<-dplyr::select(vis_obs, local_date, local_time, site, sl, contains("GRSH"), contains("HUWH"))

#add up columns for total shearwater and whale sightings
vis_obs$b_total <- rowSums(dplyr::select(vis_obs, contains("GRSH")))
vis_obs$w_total <- rowSums(dplyr::select(vis_obs, contains("HUWH")))

#select relevant data columns
data<-select(vis_obs, local_date, site, sl, w_total)

#remove all rows with no data 'nd' for SL count
data<-data[!data$sl == "nd", ]

#for some reason, sl is read in as a factor (bc of 'nd') so need to convert to character and the integer
data$sl <- as.integer(as.character(data$sl))

#format local_date as dates
data$local_date<-as.Date(data$local_date, format="%m/%d/%Y")

#create a column with just the year
data$year<-as.factor(format(data$local_date,format="%Y"))

#create a column with just the month
data$month<-as.factor(format(data$local_date,format="%m"))

#create a new column for season
data<-mutate(data, season = ifelse(month %in% 08:12, "Fall", ifelse(month=="07", "Summer",ifelse(month=="08","Summer", ifelse(month=="09", "Fall", "Spring"))))) #i don't know why 08:12 doesn't work

#separate the site into the block and the station number
data<-separate(data, site, into = c('block', 'stat'), sep = 1, remove=FALSE,convert = TRUE)

#keep only standard sites: N1-N16, C1-C14, S1-S14 (no N2?)
data<- data %>% filter(block != "R") %>% droplevels # get rid of R sites
data<-data %>% filter(stat < 17) %>% droplevels
data<-data %>% filter(site !="C15") %>% droplevels
data<-data %>% filter(site !="C16") %>% droplevels
data<-data %>% filter(site !="S15") %>% droplevels
data<-data %>% filter(site !="S16") %>% droplevels

#remove november 2018 since not complete cruise
data <- filter(data, local_date!='2018-11-09')
data <- filter(data, local_date!='2018-11-30')

#assign a cruise number
data <- data %>% 
  mutate(cruise_num = group_indices_(data, .dots=c("season", "year")))
data$cruise_num<-as.factor(data$cruise)

#remove duplicate vis obs on sites in fall 2015
data<- data[-c(77:78),]

#remove the 44 sl data point
#data<- data %>% filter(sl<44)

#remove cruises where sl=0 or whales=0 - fall13, fall14, fall16, fall17, fall17, summer18, summer19
data <- filter(data, cruise_num!='1') %>% droplevels()
data <- filter(data, cruise_num!='6') %>% droplevels()
data <- filter(data, cruise_num!='8') %>% droplevels()
data <- filter(data, cruise_num!='9') %>% droplevels()
data <- filter(data, cruise_num!='5') %>% droplevels()
data <- filter(data, cruise_num!='13') %>% droplevels()
data <- filter(data, cruise_num!='12') %>% droplevels()
data <- filter(data, cruise_num!='11') %>% droplevels()

#rename some cols
colnames(data)[1] <- c("date") 
colnames(data)[6] <- c("whales")

###summarize some things for data summary table in manuscript
# data %>% group_by(cruise_num) %>% 
#   summarise(sl=sum(sl), whales=sum(whales),n_distinct(site))
# 
# data$block <- as.factor(data$block)
# 
# fall2014 <- filter(data, cruise_num==2)
# spring2015 <- filter(data, cruise_num==7)
# fall2015 <- filter(data, cruise_num==3)
# fall2016 <- filter(data, cruise_num==4)
# spring2018 <- filter(data, cruise_num==10)
# 
# data %>% group_by(cruise_num, block) %>% summarise(slP = sum(sl>0), wP = sum(whales>0))

#format data for model
data$sl<-as.integer(data$sl)
data$whales<-as.integer(data$whales)
data$site_num<-as.numeric(as.factor(data$site)) #give site a number (separate column) for getting a coefficient
data$season<-as.numeric(as.factor(data$season))#-1  #comment out the -1 if put summer back in (more than 2 factors)
data$region<-as.numeric(as.factor(data$block))
data$year_num <- as.numeric(data$year)


data <- as_tibble(data)


###Create data list for model

#number of observations
n<-length(data$whales)


win_data <- list(n=n, whales=data$whales, sl=data$sl,site=data$site_num, block=data$region, season=data$season,
                 year=data$year_num)
win_data

#################################################################################

###Run the model
###ZIP GLMM - Zero inflated poisson mixed effects model
modelstring="model{
#likelihood
for (i in 1:n) {

#sand lance model

sl[i] ~ dpois(mu_sl[i])
mu_sl[i] <- max(lambda_sl[i]*z_sl[i], 0.00001) #only adjusts if z=0
log(lambda_sl[i]) <- block_sl[block[i]] + year_sl[year[i]] + site_sl[site[i]]
z_sl[i] ~ dbern(beta_sl[season[i]])
#z_sl[i] ~ dbern(beta_sl)

#whale model

whales[i] ~ dpois(mu_w[i])
mu_w[i] <- max(lambda_w[i]*z_w[i], 0.00001) #only adjusts if z=0
log(lambda_w[i]) <-beta*log(lambda_sl[i])  + year_w[year[i]]  + site_w[site[i]]
#log(lambda_w[i]) <-beta*log(lambda_sl[i]) + site_w[site[i]]
#z_w[i] ~ dbern(beta_w[cruise[i]])
z_w[i] ~ dbern(beta_w[season[i]])
#z_w[i] ~ dbern(beta_w)

}

#priors

#sand lance model

#probability of 0 inflation for each season for sl (prob of success - count>0)
beta_sl[1] ~ dunif(0,1)
beta_sl[2] ~ dunif(0,1)
#beta_sl ~ dunif(0,1)

#variance for random effects
alpha1 ~ dgamma(0.01,0.01) #Gelman - uniforms on standard deviation, or inverse gamma on precision, gamma on variance
alpha2 ~ dgamma(0.01,0.02) #smaller variance help autocorrelation?

#site level effect sl
for (i in 1:43){
site_sl[i] ~ dnorm(0, alpha1)} #instead of 0.01, include parameter that determines precision
                           #constrain site, site specific deviations


#year effects sl
for (i in 1:4){
year_sl[i] ~ dnorm(0, alpha2)} 
#year_sl[i] ~ dnorm(0, 0.01)} #trying year as fixed effect to reduce autocorrelation   
#instead of 0.01, include parameter that determines precision,#constrain year, deviations from block mean
      
#block sl
for (i in 1:3){
block_sl[i] ~ dnorm(0, 0.01)} #these effects independent, individual block means (Harrison et al 2018)

#whale model

beta ~ dnorm(0, 0.01) #parameter for sand lance, norm prior with mean 0 and prec=0.01 = 1/var, var=100, sd=10

#probability of 0 inflation for each season for w
beta_w[1] ~ dunif(0,1)
beta_w[2] ~ dunif(0,1)
#beta_w ~ dunif(0,1)

#variance for random effects
alpha1w ~ dgamma(0.01,0.02) #Gelman - uniforms on standard deviation, or inverse gamma on precision, gamma on variance
alpha2w ~ dgamma(0.01,0.02) #also see Zuur et al 2012

#site level effect w
for (i in 1:43){
site_w[i] ~ dnorm(0, alpha1w)} #instead of 0.01, include parameter that determines precision
                               #constrain site, site specific deviations

#year effects w
for (i in 1:4){
year_w[i] ~ dnorm(0, alpha2w)}
#year_w[i] ~ dnorm(0, 0.01)} 
 #instead of 0.01, include parameter that determines precision, #constrain year, deviations from block mean
    
#block w
#for (i in 1:3){
#block_w[i] ~ dnorm(0, 0.01)} #these effects independent, individual block means

}"

###Zero inflated negative binomial model - ZINB GLMM

#modelstring="model{
#likelihood
#for (i in 1:n) {

#sand lance model

#negative binomial / ZINB - r = overdispersion parameter, p(mu here) = success parameter
#sl[i] ~ dnegbin(mu_sl[i],r)
#mu_sl[i] <- r/(r+lambda_sl[i]) #for negative binomial
#mu_sl[i] <- r/(r + mu_eff[i])
#mu_eff[i] <- max(lambda_sl[i] * z_sl[i], 0.000001) #only adjusts if z=0
#log(lambda_sl[i]) <- block_sl[block[i]] + site_sl[site[i]] + year_sl[year[i]] 
#z_sl[i] ~ dbern(beta_sl[i])
#z_sl[i] ~ dbern(beta_sl[season[i]])

#whale model

#whales[i] ~ dpois(mu_w[i])
#mu_w[i] <- max(lambda_w[i]*z_w[i], 0.00001) #only adjusts if z=0
#log(lambda_w[i]) <-beta*log(lambda_sl[i])  + year_w[year[i]]  + site_w[site[i]]
#log(lambda_w[i]) <-beta*log(lambda_sl[i]) + site_w[site[i]]
#z_w[i] ~ dbern(beta_w[cruise[i]])
#z_w[i] ~ dbern(beta_w[season[i]])
#z_w[i] ~ dbern(beta_w)


#negative binomial / ZINB
#whales[i] ~ dnegbin(mu_w[i],rw)
#mu_w[i] <- rw/(rw+lambda_w[i]) #for negative binomial
#mu_w[i] <- rw/(rw + mu_eff_w[i])
#mu_eff_w[i] <- max(lambda_w[i]*z_w[i], 0.000001) #only adjusts if z=0
#log(lambda_w[i]) <- beta*log(lambda_sl[i]) + site_w[site[i]]  + year_w[year[i]] + theta_w*season[i]
#z_w[i] ~ dbern(beta_w[i])
#z_w[i] ~ dbern(beta_w[season[i]])

#}

#priors

#sand lance model

#o inf negbin
#r~dunif(0,50)

#probability of 0 inflation for each season for sl
#beta_sl[1] ~ dunif(0,1)
#beta_sl[2] ~ dunif(0,1)
#beta_sl ~ dunif(0,1)

#variance for random effects
#alpha1 ~ dgamma(0.01,0.01) #Gelman - uniforms on standard deviation, or inverse gamma on precision, gamma on variance
#alpha2 ~ dgamma(0.01,0.02) #smaller variance help autocorrelation?

#site level effect sl
#for (i in 1:43){
#site_sl[i] ~ dnorm(0, alpha1)} #instead of 0.01, include parameter that determines precision
                           #constrain site, site specific deviations

##year effects sl
#for (i in 1:4){
#year_sl[i] ~ dnorm(0, alpha2)} 
#year_sl[i] ~ dnorm(0, 0.01)} #trying year as fixed effect to reduce autocorrelation   
#instead of 0.01, include parameter that determines precision,#constrain year, deviations from block mean
      
#block sl
#for (i in 1:3){
#block_sl[i] ~ dnorm(0, 0.01)} #these effects independent, individual block means (Harrison et al 2018)

#whale model

#o inf negbin
#rw~dunif(0,50)

#beta ~ dnorm(0, 0.01) #parameter for sand lance, norm prior with mean 0 and prec=0.01 = 1/var, var=100, sd=10

#probability of 0 inflation for each season for w
#beta_w[1] ~ dunif(0,1)
#beta_w[2] ~ dunif(0,1)
#beta_w ~ dunif(0,1)

#variance for random effects
#alpha1w ~ dgamma(0.01,0.02) #Gelman - uniforms on standard deviation, or inverse gamma on precision, gamma on variance
#alpha2w ~ dgamma(0.01,0.02) #also see Zuur et al 2012

#site level effect w
#for (i in 1:43){
#site_w[i] ~ dnorm(0, alpha1w)} #instead of 0.01, include parameter that determines precision
                               #constrain site, site specific deviations

#year effects w
#for (i in 1:4){
#year_w[i] ~ dnorm(0, alpha2w)}
#year_w[i] ~ dnorm(0, 0.01)} 
 #instead of 0.01, include parameter that determines precision, #constrain year, deviations from block mean
    
#block w
#for (i in 1:3){
#block_w[i] ~ dnorm(0, 0.01)} #these effects independent, individual block means

#}"

#adaptations, burn in and iterations
n.adapt=50000 #number of iterations JAGS uses to choose sampler and assure optimum mixing of MCMC chains
n.update=50000 #burn-in
n.iter=1000000 #number of iterations to be stored in the chain as samples from the posterior distribution

#create a model object
model <- jags.model(textConnection(modelstring), data=win_data, n.adapt=n.adapt, n.chains=4)#if use dummy var & make list prior

#burn-in the chains 
update(model, n.iter=n.update)

#extract random samples from the posterior distribution of the parameters of a jags model
#samps<-coda.samples(model,variable.names=c("beta","block_sl","year_sl","year_w","site_sl","site_w"), 
                    #n.iter=n.iter,thin=1000)

samps2<-coda.samples(model,variable.names=c("beta_sl", "beta_w","alpha1","alpha2","alpha1w","alpha2w"),
                     n.iter=n.iter,thin=1000)
#samps3<-coda.samples(model,variable.names=c("lambda_w","lambda_sl"),n.iter=n.iter,thin=1000)
samps4<-coda.samples(model,variable.names=c("z_sl","z_w"),n.iter=n.iter,thin=1000)
samps5<-coda.samples(model,variable.names=c("mu_w","mu_sl"),n.iter=n.iter,thin=1000)

samps1 <- coda.samples(model,variable.names=c("beta","block_sl","year_sl","year_w","site_sl","site_w","lambda_sl",
                                              "lambda_w"), n.iter=n.iter,thin=1000)

#save model output
#save.image(file = "model6_FINAL0828.RData")

##################################################################

###Model checking

#load in model output
load("model6_FINAL_Oct2.RData")   

#trace plots - save manually as pdfs
MCMCtrace(samps5)

#save model summaries as dataframes - includes Rhat and eff samples
summ_samps1 <- as.data.frame(MCMCsummary(samps1, round=2,n.eff=TRUE,func=median,func_name="median"))
summ_samps2 <- as.data.frame(MCMCsummary(samps2, round=2,n.eff=TRUE,func=median,func_name="median"))
#summ_samps3 <- as.data.frame(MCMCsummary(samps3, round=2,n.eff=TRUE))
summ_samps4 <- as.data.frame(MCMCsummary(samps4, round=2,n.eff=TRUE, func=median,func_name="median"))
summ_samps5 <- as.data.frame(MCMCsummary(samps5, round=2,n.eff=TRUE, func=median,func_name="median"))

#how many lambda_sl Rhats are greater than 1.2?
summ_samps1 %>%
  rownames_to_column() %>%
  filter(Rhat>1.2) %>%
#count them
summarise(sum(Rhat>1.2))

#how many z Rhats are greater than 1.2?
summ_samps4 %>%
  rownames_to_column() %>%
  filter(Rhat>1.2) %>%
  #count them
  summarise(sum(Rhat>1.2))

#how many mu Rhats are greater than 1.2?
summ_samps5 %>%
  rownames_to_column() %>%
  filter(Rhat>1.2) %>%
  #count them
  summarise(sum(Rhat>1.2))

#count the NaN Rhat in z and mu
summ_samps4[is.na(summ_samps4$Rhat),] #find which rows
sum(is.na(summ_samps4$Rhat)) #count them
summ_samps5[is.na(summ_samps5$Rhat),] #find which rows
sum(is.na(summ_samps5$Rhat)) #count them

#create tables of model parameter results for manuscript
#library(gt) #does not export as csv or word! 
#library(reportr)

#change row names so they match parameter names
summ_samps1_short <- summ_samps1 %>% 
  rownames_to_column(var="parameter") %>%
  select(-Rhat,-n.eff,-"50%") 
  
summ_samps1_short<- summ_samps1_short %>%
  filter(grepl("beta|block|year|site",parameter)) %>%
  mutate(parameter= recode(parameter, "block_sl[1]" = "block_C","block_sl[2]" = "block_N",
                  "block_sl[3]" = "block_S", "year_sl[1]"= "year2014_sl","year_sl[2]"= "year2015_sl",
                  "year_sl[3]"= "year2016_sl", "year_sl[4]"= "year2018_sl", "year_w[1]"= "year2014_w",
                  "year_w[2]"= "year2015_w", "year_w[3]"= "year2016_w", "year_w[4]"= "year2018_w",
                  
                  "site_sl[1]"="site_sl_C1","site_sl[2]"="site_sl_C10","site_sl[3]"="site_sl_C11",
                  "site_sl[4]"="site_sl_C12","site_sl[5]"="site_sl_C13","site_sl[6]"="site_sl_C14",
                  "site_sl[7]"="site_sl_C2","site_sl[8]"="site_sl_C3","site_sl[9]"="site_sl_C4",
                  "site_sl[10]"="site_sl_C5","site_sl[11]"="site_sl_C6","site_sl[12]"="site_sl_C7",
                  "site_sl[13]"="site_sl_C8", "site_sl[14]"="site_sl_C9","site_sl[15]"="site_sl_N1",
                  "site_sl[16]"="site_sl_N10","site_sl[17]"="site_sl_N11",
                  "site_sl[18]"="site_sl_N12","site_sl[19]"="site_sl_N13","site_sl[20]"="site_sl_N14",
                  "site_sl[21]"="site_sl_N15","site_sl[22]"="site_sl_N16","site_sl[23]"="site_sl_N3",
                  "site_sl[24]"="site_sl_N4","site_sl[25]"="site_sl_N5","site_sl[26]"="site_sl_N6",
                  "site_sl[27]"="site_sl_N7", "site_sl[28]"="site_sl_N8","site_sl[29]"="site_sl_N9",
                  "site_sl[30]"="site_sl_S1","site_sl[31]"="site_sl_S10","site_sl[32]"="site_sl_S11",
                  "site_sl[33]"="site_sl_S12","site_sl[34]"="site_sl_S13","site_sl[35]"="site_sl_S14",
                  "site_sl[36]"="site_sl_S2","site_sl[37]"="site_sl_S3","site_sl[38]"="site_sl_S4",
                  "site_sl[39]"="site_sl_S5","site_sl[40]"="site_sl_S6","site_sl[41]"="site_sl_S7",
                  "site_sl[42]"="site_sl_S8", "site_sl[43]"="site_sl_S9",
                  
                  "site_w[1]"="site_w_C1","site_w[2]"="site_w_C10","site_w[3]"="site_w_C11",
                  "site_w[4]"="site_w_C12","site_w[5]"="site_w_C13","site_w[6]"="site_w_C14",
                  "site_w[7]"="site_w_C2","site_w[8]"="site_w_C3","site_w[9]"="site_w_C4",
                  "site_w[10]"="site_w_C5","site_w[11]"="site_w_C6","site_w[12]"="site_w_C7",
                  "site_w[13]"="site_w_C8", "site_w[14]"="site_w_C9","site_w[15]"="site_w_N1",
                  "site_w[16]"="site_w_N10","site_w[17]"="site_w_N11",
                  "site_w[18]"="site_w_N12","site_w[19]"="site_w_N13","site_w[20]"="site_w_N14",
                  "site_w[21]"="site_w_N15","site_w[22]"="site_w_N16","site_w[23]"="site_w_N3",
                  "site_w[24]"="site_w_N4","site_w[25]"="site_w_N5","site_w[26]"="site_w_N6",
                  "site_w[27]"="site_w_N7", "site_w[28]"="site_w_N8","site_w[29]"="site_w_N9",
                  "site_w[30]"="site_w_S1","site_w[31]"="site_w_S10","site_w[32]"="site_w_S11",
                  "site_w[33]"="site_w_S12","site_w[34]"="site_w_S13","site_w[35]"="site_w_S14",
                  "site_w[36]"="site_w_S2","site_w[37]"="site_w_S3","site_w[38]"="site_w_S4",
                  "site_w[39]"="site_w_S5","site_w[40]"="site_w_S6","site_w[41]"="site_w_S7",
                  "site_w[42]"="site_w_S8", "site_w[43]"="site_w_S9")) %>%
          relocate(median,.before=mean)

summ_samps2_short <- summ_samps2 %>% 
  rownames_to_column(var="parameter") %>%
  select(-Rhat,-n.eff,-"50%") %>%
  mutate(parameter= recode(parameter, "beta_sl[1]" = "beta_sl_fall",
            "beta_sl[2]" ="beta_sl_spring","beta_w[1]"= "beta_w_fall",
            "beta_w[2]"= "beta_w_spring")) %>%
  relocate(median,.before=mean)

#put samps2 and samps1 tables together
results_summary <- rbind(summ_samps1_short, summ_samps2_short)

#save as a csv
write_csv(results_summary,"results_summary_update_10_5.csv")
#export this table

#makes a nice table in gt package but cannot export
#summ_samps2_short %>% gt %>% #just like ggplot, gt serves as first function call to make tables in gt
  #tab_spanner(label="95% credible interval",columns=vars("2.5%","97.5%")) %>%
    #cols_move_to_end(columns=vars(mean,sd,"2.5%","97.5%"))

#autocorrelation plots - save as pdfs
pdf("./results/bayesian_model/autocorr_samps1.pdf")
autocorr.plot(samps1)
dev.off()

pdf("./results/bayesian_model/autocorr_samps2.pdf")
autocorr.plot(samps2)
dev.off()

pdf("./results/bayesian_model/autocorr_samps3.pdf")
autocorr.plot(samps3)
dev.off()

pdf("./results/bayesian_model/autocorr_samps4.pdf")
autocorr.plot(samps4)
dev.off()

pdf("./results/bayesian_model/autocorr_samps5.pdf")
autocorr.plot(samps5)
dev.off()

#Posterior predictive checks for ZIP GLMM (code by Gavin Fay)

#load data, grab the posteriors for the means for each data point, 
#rejig the mcmc output  

library(janitor)

load("./results/bayesian_model/model6.RData")
ppchecks <- do.call(rbind, samps5) %>%    #there's got to be a tidyverse way to do this....
  as_tibble() %>% 
  janitor::clean_names() %>% #names columns based on the data point name
  rowid_to_column() %>%  #takes the row id and makes it it's own column
  pivot_longer(values_to = "value",
               names_to = "label",
               cols = -rowid) %>% 
  separate(label, into = c("param", "species", "idat")) %>% 
  select(rowid, species, idat, value) %>% 
  mutate(newdat = rpois(nrow(.),value)) %>% 
  I()
ppchecks

#compute the posterior predictive distributions for some summary statistics by species, 
#mean, variance, & proportion of zeroes

ppcheck_summaries <- ppchecks %>% 
  group_by(species, rowid) %>% 
  mutate(is_zero = ifelse(newdat == 0, 1, 0)) %>% 
  summarise(mu = mean(newdat, na.rm = TRUE),
            var = var(newdat, na.rm = TRUE), 
            n = length(newdat),
            nzero = sum(is_zero)) %>% 
  mutate(propzero = nzero/n) %>% 
  ungroup() %>% 
  select(species, rowid, mu, var, propzero) %>% 
  pivot_longer(names_to = "stat",
               values_to = "value",
               cols = c("mu", "var", "propzero")) %>%
  I()
ppcheck_summaries

#Create summaries of the data by species, mean, variance, & proportion of zeroes

data_summaries <- data %>% 
  select(sl, whales) %>% 
  rowid_to_column() %>% 
  pivot_longer(names_to = "species",
               values_to = "count",
               cols = -rowid) %>% 
  mutate(species = case_when(
    species == "whales" ~ "w",
    TRUE ~ species)) %>% 
  group_by(species) %>% 
  mutate(is_zero = ifelse(count == 0, 1, 0)) %>% 
  summarise(mu = mean(count, na.rm = TRUE),
            var = var(count, na.rm = TRUE), 
            n = length(count),
            nzero = sum(is_zero)) %>% 
  mutate(propzero = nzero/n) %>% 
  ungroup() %>% 
  select(species, mu, var, propzero) %>% 
  pivot_longer(names_to = "stat",
               values_to = "value",
               cols = c("mu", "var", "propzero")) %>%
  I()
data_summaries  

#Calculate the Bayesian p-values

bpvals <- ppcheck_summaries %>% 
  group_by(species, stat) %>% 
  nest() %>% #takes the other cols(rowid, value) and makes a col for each eithin each species-stat group (so 6 - sl mu, var and prop0, whales mu var prop0)
  mutate(cdf = map(data, ~ecdf(.x$value))) %>% #this data means the data from bpval$data (the values) , ecdf finds fraction of observations less or equal to value
  left_join(data_summaries) %>%  #line above says take each value (mu, var, etc) from each spp and make a histo of the values
  mutate(pval = map2_dbl(value, cdf, ~.y(.x))) %>% #proportion of draws less than data value
  I()                                              #think of this like f(x) - the way ecdf plots look
bpvals

#plot the distributions & values for the summary statistics and save plot

tiff("./figures/Bayesian_manuscript/post_pred_checks.tiff", width = 10, height = 8, units = 'in', res = 300)
ppcheck_summaries %>% 
  ggplot() +
  aes(x = value, fill = stat) +
  geom_histogram(col="white") +
  geom_vline(data = data_summaries, aes(xintercept = value), col = "black",size=1) +
  facet_grid(species~stat, scales = "free") +
  theme_bw(base_size=14)
# NULL
dev.off()

#summarize posterior predictive checks for block, season and year

#add block, season and year in here for ppchecks
#add a idat column to data
data_for_ppchecks<- data %>% rowid_to_column(var="idat")
data_for_ppchecks$idat <- as.character(data_for_ppchecks$idat)
data_for_ppchecks <- data_for_ppchecks %>% select(idat, block, year, season,site)

#merge with ppchecks based on idat

ppchecks <- full_join(ppchecks, data_for_ppchecks, by="idat")

ppcheck_summaries_block <- ppchecks %>% 
  group_by(species, rowid, block) %>% 
  mutate(is_zero = ifelse(newdat == 0, 1, 0)) %>% 
  summarise(mu = mean(newdat, na.rm = TRUE),
            var = var(newdat, na.rm = TRUE), 
            n = length(newdat),
            nzero = sum(is_zero)) %>% 
  mutate(propzero = nzero/n) %>% 
  ungroup() %>% 
  select(species, rowid,block, mu, var, propzero) %>% 
  pivot_longer(names_to = "stat",
               values_to = "value",
               cols = c("mu", "var", "propzero")) %>%
  I()
ppcheck_summaries_block

ppcheck_summaries_year <- ppchecks %>% 
  group_by(species, rowid, year) %>% 
  mutate(is_zero = ifelse(newdat == 0, 1, 0)) %>% 
  summarise(mu = mean(newdat, na.rm = TRUE),
            var = var(newdat, na.rm = TRUE), 
            n = length(newdat),
            nzero = sum(is_zero)) %>% 
  mutate(propzero = nzero/n) %>% 
  ungroup() %>% 
  select(species, rowid,year, mu, var, propzero) %>% 
  pivot_longer(names_to = "stat",
               values_to = "value",
               cols = c("mu", "var", "propzero")) %>%
  I()
ppcheck_summaries_year

ppcheck_summaries_season <- ppchecks %>% 
  group_by(species, rowid, season) %>% 
  mutate(is_zero = ifelse(newdat == 0, 1, 0)) %>% 
  summarise(mu = mean(newdat, na.rm = TRUE),
            var = var(newdat, na.rm = TRUE), 
            n = length(newdat),
            nzero = sum(is_zero)) %>% 
  mutate(propzero = nzero/n) %>% 
  ungroup() %>% 
  select(species, rowid,season, mu, var, propzero) %>% 
  pivot_longer(names_to = "stat",
               values_to = "value",
               cols = c("mu", "var", "propzero")) %>%
  I()
ppcheck_summaries_season

#data summary by block
data_summaries_block <- data %>% 
  select(sl, whales,block) %>% 
  rowid_to_column() %>% 
  pivot_longer(names_to = "species",
               values_to = "count",
               cols = c(-rowid,-block)) %>% 
  mutate(species = case_when(
    species == "whales" ~ "w",
    TRUE ~ species)) %>% 
  group_by(species,block) %>% 
  mutate(is_zero = ifelse(count == 0, 1, 0)) %>% 
  summarise(mu = mean(count, na.rm = TRUE),
            var = var(count, na.rm = TRUE), 
            n = length(count),
            nzero = sum(is_zero)) %>% 
  mutate(propzero = nzero/n) %>% 
  ungroup() %>% 
  select(species,block, mu, var, propzero) %>% 
  pivot_longer(names_to = "stat",
               values_to = "value",
               cols = c("mu", "var", "propzero")) %>%
  I()
data_summaries_block

#data summary by year
data_summaries_year <- data %>% 
  select(sl, whales,year) %>% 
  rowid_to_column() %>% 
  pivot_longer(names_to = "species",
               values_to = "count",
               cols = c(-rowid,-year)) %>% 
  mutate(species = case_when(
    species == "whales" ~ "w",
    TRUE ~ species)) %>% 
  group_by(species,year) %>% 
  mutate(is_zero = ifelse(count == 0, 1, 0)) %>% 
  summarise(mu = mean(count, na.rm = TRUE),
            var = var(count, na.rm = TRUE), 
            n = length(count),
            nzero = sum(is_zero)) %>% 
  mutate(propzero = nzero/n) %>% 
  ungroup() %>% 
  select(species,year, mu, var, propzero) %>% 
  pivot_longer(names_to = "stat",
               values_to = "value",
               cols = c("mu", "var", "propzero")) %>%
  I()
data_summaries_year

#data summary by season

data_summaries_season <- data %>% 
  select(sl, whales,season) %>% 
  rowid_to_column() %>% 
  pivot_longer(names_to = "species",
               values_to = "count",
               cols = c(-rowid,-season)) %>% 
  mutate(species = case_when(
    species == "whales" ~ "w",
    TRUE ~ species)) %>% 
  group_by(species,season) %>% 
  mutate(is_zero = ifelse(count == 0, 1, 0)) %>% 
  summarise(mu = mean(count, na.rm = TRUE),
            var = var(count, na.rm = TRUE), 
            n = length(count),
            nzero = sum(is_zero)) %>% 
  mutate(propzero = nzero/n) %>% 
  ungroup() %>% 
  select(species,season, mu, var, propzero) %>% 
  pivot_longer(names_to = "stat",
               values_to = "value",
               cols = c("mu", "var", "propzero")) %>%
  I()
data_summaries_season

#bpvals

#by block
bpvals_block <- ppcheck_summaries_block %>% 
  group_by(species,block, stat) %>% 
  nest() %>% 
  mutate(cdf = map(data, ~ecdf(.x$value))) %>% #ecdf finds fraction of observations less or equal to value
  left_join(data_summaries_block) %>% 
  mutate(pval = map2_dbl(value, cdf, ~.y(.x))) %>% 
  I()
bpvals_block

#by year
bpvals_year <- ppcheck_summaries_year %>% 
  group_by(species,year, stat) %>% 
  nest() %>% 
  mutate(cdf = map(data, ~ecdf(.x$value))) %>% #ecdf finds fraction of observations less or equal to value
  left_join(data_summaries_year) %>% 
  mutate(pval = map2_dbl(value, cdf, ~.y(.x))) %>% 
  I()
bpvals_year

#by season
bpvals_season <- ppcheck_summaries_season %>% 
  group_by(species,season, stat) %>% 
  nest() %>% 
  mutate(cdf = map(data, ~ecdf(.x$value))) %>% #ecdf finds fraction of observations less or equal to value
  left_join(data_summaries_season) %>% 
  mutate(pval = map2_dbl(value, cdf, ~.y(.x))) %>% 
  I()
bpvals_season

#plot and save post pred check figures for block, season and year
tiff("./figures/Bayesian_manuscript/post_pred_checks_block.tiff", width = 10, height = 8, units = 'in', res = 300)
ppcheck_summaries_block %>% 
  ggplot() +
  aes(x = value, fill = stat) +
  geom_histogram(col="white") +
  geom_vline(data = data_summaries_block, aes(xintercept = value), col = "black",size=1) +
  facet_grid(species+block~stat, scales = "free") +
  theme_bw(base_size=14)
NULL
dev.off()

tiff("./figures/Bayesian_manuscript/post_pred_checks_year.tiff", width = 10, height = 8, units = 'in', res = 300)
ppcheck_summaries_year %>% 
  ggplot() +
  aes(x = value, fill = stat) +
  geom_histogram(col="white") +
  geom_vline(data = data_summaries_year, aes(xintercept = value), col = "black",size=1) +
  facet_grid(species+year~stat, scales = "free") +
  theme_bw(base_size=14)
NULL
dev.off()

tiff("./figures/Bayesian_manuscript/post_pred_checks_season.tiff", width = 10, height = 8, units = 'in', res = 300)
ppcheck_summaries_season %>% 
  ggplot() +
  aes(x = value, fill = stat) +
  geom_histogram(col="white") +
  geom_vline(data = data_summaries_season, aes(xintercept = value), col = "black",size=1) +
  facet_grid(species+season~stat, scales = "free") +
  theme_bw(base_size=14)
NULL
dev.off()