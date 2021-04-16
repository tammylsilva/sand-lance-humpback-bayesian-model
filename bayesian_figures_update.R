
#This script makes figures from sand lance / humpback Bayesian model output
#8/20/2020
#Tammy Silva

library(tidyverse)
library(tidybayes)
library(viridis)
library(raster)
library(ggsn) #scale bars and north arrows

#load in model results
#ZIP GLMM
load("./results/bayesian_model/model6_FINAL_Oct2.RData") 

###site and environmental data for maps
#read in site data
sites<-read.csv(file="./data/seaboss/sites.csv", header=TRUE)

#remove rows with NA (no positions for stations)
sites<-na.omit(sites)

#get rid of LSB and R sites
sites <- filter(sites, !grepl('LSB|R', site))

#split site into 2 columns-letters vs numbers
sites <- separate(sites, site, into = c("text", "num"), sep = "(?<=[A-Za-z])(?=[0-9])")
sites$num <- as.numeric(sites$num)

#only keep the standard 44 sites and remove all others
sites <- filter(sites, num < 17)

#put columns back together to make sites normal again
sites <- unite(sites, site, c("text", "num"),sep = "")

#sites to factor
sites$site <- as.factor(sites$site)

sites<-separate(sites, site, into = c('block', 'stat'), sep = 1, remove=FALSE,convert = TRUE)

#remove extra S and C sites
sites <- filter(sites, site !="S15")
sites <- filter(sites, site !="S16")
sites <- filter(sites, site !="C15")
sites <- filter(sites, site !="C16")

#make sites lon minus
sites$lon<- -sites$lon

#read in bathymetry data
str_name='./data/environmental/mikesbathym.tif'
depth_raster<-raster(str_name)

#crop the depth_raster to the extent of sbnms
depth_crop<-crop(depth_raster, sbnms)
depth_mask<-mask(depth_crop, sbnms)

#for ggplot
depth_mask_fort <- as.data.frame(depth_mask,xy=TRUE)

#set colors for depth
library(RColorBrewer)
cols<-colorRampPalette(brewer.pal(9, "Blues"))(100)

#######################################################################################

###exploratory data figures

###visualize the entire raw filtered data set (cleaned up but all cruises) (n=439)

#plot histogram of whales to see zero inflation (total or color by a factor)

tiff("./figures/Bayesian_manuscript/whale_hist_alldata_season.tiff", width = 6, height = 5, units = 'in', res = 300)
ggplot(full_filt_data,aes(whales, fill=season)) + 
  geom_histogram(binwidth=1) + 
  scale_fill_viridis(discrete=T) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,400)) + 
  labs(x="Whales", y="Count") + 
  theme_bw(base_size=18)
dev.off()

#plot histogram of sand lance to see zero inflation (total or color by a factor)

tiff("./figures/Bayesian_manuscript/sl_hist_alldata_season.tiff", width = 6, height = 5, units = 'in', res = 300)
ggplot(full_filt_data,aes(sl,fill=season)) + 
  geom_histogram(binwidth=1) + 
  scale_fill_viridis(discrete=T) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,400)) + 
  labs(x="Sand lance", y="Count") + 
  theme_bw(base_size=18)
dev.off()

#plot total number of whales by block, year and season

#fall whales counts

tiff("./figures/Bayesian_manuscript/fall_whale_counts.tiff", width = 6, height = 5, units = 'in', res = 300)
full_filt_data%>%
  filter(season=='Fall') %>%
 ggplot(aes(x=year, y=whales, fill=block)) + geom_bar(stat="identity") +  
  scale_fill_viridis_d(name="Block",drop=FALSE)+ scale_x_discrete(drop=FALSE) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,60)) + 
  labs(x="Year", y="Whales") + 
  theme_bw(base_size=18) + geom_text(x=2, y=55,label="Fall",size=7)
dev.off()


#spring whales counts

tiff("./figures/Bayesian_manuscript/spring_whale_counts.tiff", width = 6, height = 5, units = 'in', res = 300)
full_filt_data%>%
  filter(season=='Spring') %>%
ggplot(aes(x=year, y=whales, fill=block)) + geom_bar(stat="identity") +  
  scale_fill_viridis_d(name="Block",drop=FALSE)+ scale_x_discrete(drop=FALSE) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,60)) + 
  labs(x="Year", y="Whales") + 
  theme_bw(base_size=18) + geom_text(x=2, y=55,label="Spring",size=7)
dev.off()

#plot total number of sand lance by block, year and season

#fall sand lance counts

tiff("./figures/Bayesian_manuscript/fall_sl_counts.tiff", width = 6, height = 5, units = 'in', res = 300)
full_filt_data%>%
  filter(season=='Fall') %>%
  ggplot(aes(x=year, y=sl, fill=block)) + geom_bar(stat="identity") +  
  scale_fill_viridis_d(name="Block",drop=FALSE)+ scale_x_discrete(drop=FALSE) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,130)) + 
  labs(x="Year", y="Sand lance") + 
  theme_bw(base_size=18) + geom_text(x=2, y=120,label="Fall",size=7)
dev.off()

#spring whales counts

tiff("./figures/Bayesian_manuscript/spring_sl_counts.tiff", width = 6, height = 5, units = 'in', res = 300)
full_filt_data%>%
  filter(season=='Spring') %>%
  ggplot(aes(x=year, y=sl, fill=block)) + geom_bar(stat="identity") +  
  scale_fill_viridis_d(name="Block",drop=FALSE)+ scale_x_discrete(drop=FALSE) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,130)) + 
  labs(x="Year", y="Sand lance") + 
  theme_bw(base_size=18) + geom_text(x=2, y=120,label="Spring",size=7)
dev.off()


###visualize raw data (n=164) used in the model

#plot histogram of whales to see zero inflation (total or color by a factor)

tiff("./figures/Bayesian_manuscript/whale_hist.tiff", width = 6, height = 5, units = 'in', res = 300)
ggplot(data,aes(whales))+#, fill=block)) + 
  geom_histogram(binwidth=1) + 
  #scale_fill_viridis(discrete=T) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,130)) + 
  labs(x="Whales", y="Count") + 
  theme_bw(base_size=18)
dev.off()

#plot histogram of sand lance to see zero inflation (total or color by a factor)

tiff("./figures/Bayesian_manuscript/sl_hist.tiff", width = 6, height = 5, units = 'in', res = 300)
ggplot(data,aes(sl, fill=block)) + 
  geom_histogram(binwidth=1) + 
 scale_fill_viridis(discrete=T) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,130)) + 
  #scale_x_continuous(expand=c(0,0), limits=c(0,30)) +
  labs(x="Sand lance", y="Count") + 
  theme_bw(base_size=18)
dev.off()

#zoom in of higher counts
ggplot(data,aes(sl, fill=block)) + 
  geom_histogram(binwidth=1) + 
  scale_fill_viridis(discrete=T) +
  scale_y_continuous(expand = c(0, 0),breaks=c(0,1), limits=c(0,1)) + 
  scale_x_continuous(limits=c(15,45)) +
  labs(x="Sand lance", y="Count") + 
  theme_bw(base_size=18)
ggsave("./figures/Bayesian_manuscript/sl_highcounts_zoom.tiff",width=4, height=4, units="in")


#plot counts at each site - change for sl or whales
counts_sites <- data %>% group_by(site) %>% summarise(sl=sum(sl), whales=sum(whales)) %>% 
  ggplot(aes(x=factor(site,levels=c("S14","S13","S12","S11","S10","S9","S8","S7","S6","S5","S4","S3","S2","S1",
                                    "C14","C13","C12","C11","C10","C9","C8","C7","C6","C5","C4","C3","C2","C1",
                                    "N16","N15","N14","N13","N12","N11","N10","N9","N8","N7","N6","N5","N4","N3","N1")),
                      y=whales)) + geom_bar(stat="identity") + scale_fill_viridis(discrete=T) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,60))  + 
  labs(x="Sand lance", y="Count") + 
  theme_bw(base_size=18) + coord_flip()
ggsave("./figures/Bayesian_manuscript/w_counts_by_site.tiff",width=8, height=10, units="in")

#####################################################################

#summarize model parameters figures

#model checking and plotting with tidybayes
#https://mjskay.github.io/tidybayes/articles/tidybayes.html for tutorial
#https://mjskay.github.io/tidybayes/ for quick summary of plot types

#extract model output

#this creates a tibble with values for each MCMC draw for each chain for each variable 
#variables that are indexed (block, site, year - use the indexing from the model)
samps_data<-samps %>% recover_types(data) %>% spread_draws(beta, block_sl[block],year_sl[year],year_w[year],
                                                           site_sl[site], site_w[site]) %>%
  median_qi() #this generates medians and upper and lower % ranges for every variable
  #can specify certain ones if you want
  
  
#stat plot showing 50% to 95% intervals  (but can change using .width argument)
#rather then summarize posterior before calling ggplot (with median call), we can
#use stat_pointinterval to perform the summary within ggplot
  
###Plots of model parameter estimates

#beta

samps1 %>% 
  spread_draws(beta) %>%
  ggplot(aes(x=beta)) +
  stat_pointinterval(.width=c(.95,.5)) +
  labs(x="Parameter estimate" , y="beta") +
  theme_bw(base_size=18) + 
  scale_y_continuous(breaks=NULL)
ggsave("./figures/Bayesian_manuscript/post_param_estimates_beta.tiff")

#block
  
  samps1 %>% recover_types(data) %>%
  spread_draws(block_sl[block]) %>%
  ggplot(aes(y=factor(block,level=c("S","C","N")),x=block_sl)) +
  stat_pointinterval(.width=c(.95,.5)) +
  labs(x="Parameter estimate" , y="Block") +
  theme_bw(base_size=18) #+ 
  #geom_vline(xintercept=0,linetype="dashed")
ggsave("./figures/Bayesian_manuscript/post_param_estimates_block.tiff")

#year_sl

samps1 %>% recover_types(data) %>%
  spread_draws(year_sl[year]) %>%
  ggplot(aes(y=factor(year,level=c("2018","2016","2015","2014")),x=year_sl)) +
  stat_pointinterval(.width=c(.95,.5)) +
  labs(x="Parameter estimate" , y="Year") +
  theme_bw(base_size=18) + 
  geom_vline(xintercept=0,linetype="dashed")
ggsave("./figures/Bayesian_manuscript/post_param_estimates_year_sl.tiff")

#year_w

samps1 %>% recover_types(data) %>%
  spread_draws(year_w[year]) %>%
  ggplot(aes(y=factor(year,level=c("2018","2016","2015","2014")),x=year_w)) +
  stat_pointinterval(.width=c(.95,.5)) +
  labs(x="Parameter estimate" , y="Year") +
  theme_bw(base_size=18) + 
  geom_vline(xintercept=0,linetype="dashed")
ggsave("./figures/Bayesian_manuscript/post_param_estimates_year_w.tiff")

#site_sl

samps1 %>% recover_types(data) %>%
  spread_draws(site_sl[site]) %>%
  ggplot(aes(y=factor(site,level=c("S14","S13","S12","S11","S10","S9","S8","S7","S6","S5","S4","S3","S2","S1",
                                   "C14","C13","C12","C11","C10","C9","C8","C7","C6","C5","C4","C3","C2","C1",
                                   "N16","N15","N14","N13","N12","N11","N10","N9","N8","N7","N6","N5","N4","N3","N1")),
             x=site_sl)) +
  stat_pointinterval(.width=c(.95,.5)) +
  labs(x="Parameter estimate" , y="Site") +
  theme_bw(base_size=18) + 
  geom_vline(xintercept=0,linetype="dashed")
ggsave("./figures/Bayesian_manuscript/post_param_estimates_site_sl.tiff",width=8, height=10, units="in")

#site_w

samps1 %>% recover_types(data) %>%
  spread_draws(site_w[site]) %>%
  ggplot(aes(y=factor(site,level=c("S14","S13","S12","S11","S10","S9","S8","S7","S6","S5","S4","S3","S2","S1",
                                   "C14","C13","C12","C11","C10","C9","C8","C7","C6","C5","C4","C3","C2","C1",
                                   "N16","N15","N14","N13","N12","N11","N10","N9","N8","N7","N6","N5","N4","N3","N1")),
             x=site_w)) +
  stat_pointinterval(.width=c(.95,.5)) +
  labs(x="Parameter estimate" , y="Site") +
  theme_bw(base_size=18) + 
  geom_vline(xintercept=0,linetype="dashed")
ggsave("./figures/Bayesian_manuscript/post_param_estimates_site_w.tiff",width=8, height=10, units="in")


###plot parameter estimates in response space

#block
samps1 %>% recover_types(data) %>%
  spread_draws(block_sl[block]) %>%
  ggplot(aes(y= factor(block,level=c("S","C","N")),x=exp(block_sl))) +
  stat_pointinterval(.width=c(.5,.95)) +
  labs(x="Predicted Number of Sand lance" , y="Block") +
  theme_bw(base_size=18) + 
  #geom_vline(xintercept=0,linetype="dashed")
  ggsave("./figures/Bayesian_manuscript/response_space_estimates_block.tiff")

#year_sl
#response space = exp(block + year)
#viridis color palette
#purple #440154FF
#yellow #FDE725FF
#green #21908CFF

samps1 %>% recover_types(data) %>%
  spread_draws(year_sl[year],block_sl[block]) %>%
  mutate(block_year = exp(block_sl + year_sl)) %>%
  ggplot(aes(y=factor(year,level=c("2018","2016","2015","2014")),x=block_year,group=block,color=factor(block))) +
  scale_color_viridis(discrete=TRUE,name="Block") +
  stat_pointinterval(.width=c(.95,.5),position="dodge") +
  labs(x="Predicted number of sand lance" , y="Year") +
  theme_bw(base_size=18) + 
  geom_vline(xintercept=exp(-2.60),linetype="dashed",color="#21908CFF",size=1) + #block N median (still exp right?)
  geom_vline(xintercept=exp(-0.31),linetype="dashed",color="#440154FF",size=1) + #block C median 
  geom_vline(xintercept=exp(1.32),linetype="dashed",color="#FDE725FF",size=1)  #block S median 
ggsave("./figures/Bayesian_manuscript/response_space_estimates_year_sl.tiff",width=5, height=5, units='in')

#year_w
#join betas, lambda_sl, year
#so year estimate of whales = exp(beta*log(lambda_sl) + year) - correct?
#to plot vertical lines with mean number of whales per block - exp(block*beta+year_w)

samps1 %>% recover_types(data) %>%
  spread_draws(year_sl[year],year_w[year],beta,lambda_sl[lambda_sl],block_sl[block]) %>%
  #mutate(whales_year = exp((beta*log(lambda_sl)) + year_w)) %>% #is this right?
  mutate(whales_year = exp(beta*(block_sl + year_sl) + year_w)) %>%  
  #ggplot(aes(y=factor(year,level=c("2018","2016","2015","2014")),x=whales_year)) +
  ggplot(aes(y=factor(year,level=c("2018","2016","2015","2014")),x=whales_year,group=block,color=factor(block))) +
  scale_color_viridis(discrete=TRUE,name="Block") +
  stat_pointinterval(.width=c(.95,.5),position="dodge") +
  labs(x="Predicted number of whales" , y="Year") +
  theme_bw(base_size=18) +
  geom_vline(xintercept=exp(-2.60 * 0.35),linetype="dashed",color="#21908CFF",size=1) + #block N median #whales, 0.35-median beta (still exp right?)
  geom_vline(xintercept=exp(-0.31* 0.35),linetype="dashed",color="#440154FF",size=1) + #block C median 
  geom_vline(xintercept=exp(1.32* 0.35),linetype="dashed",color="#FDE725FF",size=1)  #block S median 
ggsave("./figures/Bayesian_manuscript/response_space_estimates_year_w_block_whales_update.tiff",width=5, height=5, units='in')


 #site_sl exp(block +site)

samps1 %>% recover_types(data) %>%
  spread_draws(block_sl[block],site_sl[site]) %>%
  filter(site!='N15') %>% #to remove N15 b'c it's huge
  #filter(site=='N15') %>% #just N15
  mutate(block_site = exp(block_sl + site_sl)) %>%
  ggplot(aes(y=factor(site,level=c("S14","S13","S12","S11","S10","S9","S8","S7","S6","S5","S4","S3","S2","S1",
                                   "C14","C13","C12","C11","C10","C9","C8","C7","C6","C5","C4","C3","C2","C1",
                                   "N16","N15","N14","N13","N12","N11","N10","N9","N8","N7","N6","N5","N4","N3","N1"))
             ,x=block_site)) +
  stat_pointinterval(.width=c(.95,.5)) +
  labs(x="Predicted number of sand lance" , y="Site") +
  theme_bw(base_size=18) #+ 
  #geom_vline(xintercept=0,linetype="dashed")
ggsave("./figures/Bayesian_manuscript/response_space_estimates_site_sl_noN15.tiff",width=8, height=10, units="in")

#site_w

samps1 %>% recover_types(data) %>%
  spread_draws(site_w[site],beta,lambda_sl[lambda_sl]) %>%
  mutate(whales_site = exp((beta*log(lambda_sl)) + site_w)) %>%
  ggplot(aes(y=factor(site,level=c("S14","S13","S12","S11","S10","S9","S8","S7","S6","S5","S4","S3","S2","S1",
                                   "C14","C13","C12","C11","C10","C9","C8","C7","C6","C5","C4","C3","C2","C1",
                                   "N16","N15","N14","N13","N12","N11","N10","N9","N8","N7","N6","N5","N4","N3","N1"))
             ,x=whales_site)) +
  stat_pointinterval(.width=c(.95,.5)) +
  labs(x="Predicted number of whales" , y="Site") +
  theme_bw(base_size=18) #+ 
  #geom_vline(xintercept=0,linetype="dashed")
ggsave("./figures/Bayesian_manuscript/response_space_estimates_site_w.tiff",width=8, height=10, units="in")

#########################################################################

###calculate derived quantities for additional model inference

###get probability that fall occurrence is greater than spring occurrence
#create of df of mcmc draws of params
post_season_data <- as.data.frame(rbind(samps2[[1]], samps2[[2]],samps2[[3]], samps2[[4]]))

#get just season params
post_season_data <- post_season_data[,5:8]

#rename cols
colnames(post_season_data) <- c("sl_fall", "sl_spring", "w_fall", "w_spring")

#find numbers of draws where fall param > spring param
mean(post_season_data$sl_fall > post_season_data$sl_spring)
mean(post_season_data$w_fall > post_season_data$w_spring)
###let's see probability that S is greater than c or N by calculating a bayesian p value

#create a df of mcmc draws of params stored in samps
post_data <- as.data.frame(rbind(samps1[[1]], samps1[[2]],samps1[[3]], samps1[[4]]))

#pull out the block draws
block_data <-post_data %>% select(contains("block_sl")) 
colnames(block_data)<-c("C","N","S")

#find numbers of draws where params greater than another block
mean(block_data$C < block_data$S)
mean(block_data$N < block_data$C)
mean(block_data$N < block_data$S)

###look at site probabilities of sites having greater than average sand lance than the respective block

site_sl_data <- post_data %>% select(contains("site_sl")) 

head(site_sl_data)

site_sl_data_tib <-site_sl_data %>%  
  as_tibble() %>% 
  janitor::clean_names() %>% #names columns based on the data point name
  rowid_to_column() %>%  #takes the row id and makes it it's own column
  pivot_longer(values_to = "value",
               names_to = "label",
               cols = -rowid) %>% 
  separate(label, into = c("param","species", "site_num")) 
  
  prob_greater_than_avg <- site_sl_data_tib %>%
  group_by(site_num) %>%
  summarise(p = 1-ecdf(value)(0))
  
  prob_greater_than_avg$site_num <- as.numeric(prob_greater_than_avg$site_num)
  
  #get site names/num pairs from data to pair with lat lon and probabilities for plotting
  
  site_name_num <- data %>% dplyr::select(site, site_num) 
  site_name_num <- unique(site_name_num)
  site_name_num <- inner_join(site_name_num,sites,by="site")
  
  prob_greater_than_avg <- inner_join(prob_greater_than_avg,site_name_num, by="site_num")
  prob_greater_than_avg
  

###make a map and color by this probability
  #color gradient continuous because I'm mapping prob>average (not divergent colors bc I'm not plotting above & below here)
  
  #use mean param values for a divergent color scale -  red to blue with white as 0
  #color gradient, blue to red with 0.5 as white (where param would equal 0 and be equal to block mean)
  # mid <- 0.5
  # sp + scale_color_gradient(low = "blue", high = "red")#for a color gradient
  # sp + scale_color_gradient2(midpoint = mid, low = "blue", mid = "white",
  #                            high = "red", space = "Lab" ) #for a color gradient with a different midpoint
  
  g <- ggplot(data = states) +
    geom_polygon(aes(x=long,y=lat))+
    geom_path(data = sbnms_fort, 
              aes(x = long, y = lat, group = group),
              color = 'black', size = .4) +
    coord_map(xlim = c(-70.62,-70),ylim = c(42.11,42.8)) + #just SBNMS
    labs(x="Longitude",y="Latitude") +
    geom_contour(data=bathy_fort,
                 aes(x=x,y=y,z=z),
                 breaks=c(-40,-50),
                 color="grey") +
    geom_point(data=prob_greater_than_avg,aes(x=lon,y=lat,color=p,size=0.1)) +
    #to color all sites by prob >block mean
    # geom_text(data=s1,aes(x=lon,y=lat,label=site),size=3,vjust=-0.4) +
    # geom_text(data=s8,aes(x=lon,y=lat,label=site),size=3,vjust=1.4) +
    # geom_text(data=n1,aes(x=lon,y=lat,label=site),size=3,vjust=-0.4) +
    # geom_text(data=n8,aes(x=lon,y=lat,label=site),size=3,vjust=1.4) +
    # geom_text(data=c1,aes(x=lon,y=lat,label=site),size=3,hjust=-0.4) +
    # geom_text(data=c8,aes(x=lon,y=lat,label=site),size=3,hjust=1.4) +
    theme_bw(base_size=18) +
    theme(legend.title = element_blank(),legend.background = element_rect(fill=NA),legend.position=c(0.85,0.8)) +
    guides(colour = guide_colourbar(label.theme=element_text(size=18), barheight=10,
                                     ticks.linewidth = 2.5)) +
    guides(size=FALSE) 
  
  g <- g + annotate("text", x=-70.15, y=42.8, label="Probability of greater than \n block average sand lance", size=7)
   
    #guides(fill=guide_colourbar(legend.position="bottom",title.position="bottom"))
  ggsave("./figures/Bayesian_manuscript/prob_sites_greater_than_avg_sl_update.tiff",width=8,height=10, units="in")
    
    
###find probability of whale aggregations
  #does SW corner sites have aggregations of whales? SW corner: S1-S4, S11-S14

#already tooks mus(mean number of whales for each data point) and generated new counts for posterior predictive checks
  #use these new data, group by site and find proportion that has whales > 5 (or whatever aggregation)

ppchecks #with data for ppchecks included so I have sites in there

#say an aggregation =5 whales
prob_aggregation <- ppchecks %>%
  filter(species=='w') %>%
  group_by(site) %>%
  summarise(p_agg = round(mean(newdat>5),2))

prob_aggregation

#join this with sites so I can plot these probs on a map
prob_aggregation <- inner_join(prob_aggregation,sites, by="site")

#plot a map

wplot<- ggplot(data = states) +
  geom_polygon(aes(x=long,y=lat))+
  geom_path(data = sbnms_fort, 
            aes(x = long, y = lat, group = group),
            color = 'black', size = .4) +
  coord_map(xlim = c(-70.62,-70),ylim = c(42.11,42.8)) + #just SBNMS
  
  labs(x="Longitude",y="Latitude") +
  geom_contour(data=bathy_fort,
               aes(x=x,y=y,z=z),
               breaks=c(-40,-50),
               color="grey") +
  geom_point(data=prob_aggregation,aes(x=lon,y=lat,color=p_agg,size=0.1)) +#to color all sites by prob >block mean
  # geom_text(data=s1,aes(x=lon,y=lat,label=site),size=3,vjust=-0.4) +
  # geom_text(data=s8,aes(x=lon,y=lat,label=site),size=3,vjust=1.4) +
  # geom_text(data=n1,aes(x=lon,y=lat,label=site),size=3,vjust=-0.4) +
  # geom_text(data=n8,aes(x=lon,y=lat,label=site),size=3,vjust=1.4) +
  # geom_text(data=c1,aes(x=lon,y=lat,label=site),size=3,hjust=-0.4) +
  # geom_text(data=c8,aes(x=lon,y=lat,label=site),size=3,hjust=1.4) +
  theme_bw(base_size=18) +
  theme(legend.title = element_blank(),legend.background = element_rect(fill=NA),legend.position=c(0.85,0.77)) +
  guides(colour = guide_colourbar(label.theme=element_text(size=18), barheight=10,
                                  ticks.linewidth = 2.5)) +
  guides(size=FALSE) 

wplot <- wplot + annotate("text", x=-70.15, y=42.8, label="Probability of whale \n aggregation (>5)", size=7)

  
ggsave("./figures/Bayesian_manuscript/prob_sites_whales_aggreg5_update_2021_01_07.tiff",width=8,height=8, units="in")

#####################################################################

#other figures for manuscript

###Figure 1 - map of study area
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgdal)
library(marmap)
library(sp)
library(maps)
library(mapdata)
library(mapproj)
library(viridis)
library(scales) #to show viridis color numbers

#find color numbers for 3 viridis colors
#show_col(viridis_pal()(3))
#purple #440154FF
#yellow #FDE725FF
#green #21908CFF

#us east coast map for inset
#usa <- map_data("usa") # we already did this, but we can do it again
#ggplot() + geom_polygon(data = usa, aes(x=long, y = lat, group = group)) + 
  #coord_map(xlim = c(-82,-67),ylim = c(25,47.4)) +
  #()
#ggsave("./figures/Bayesian_manuscript/USEastCoast.tiff",width=8, height=8, units="in")

#MA data
states <- map_data("state") %>%
  filter(region=="massachusetts")

#read in sbnms boudaries
sbnms<-readOGR("./data/environmental/sbnms.shp")

#for ggplot
sbnms_fort <- fortify(sbnms)

#get bathy data from marmap
bathy = getNOAA.bathy(lon1 = -71, lon2 = -69, lat1 = 42, lat2 = 42.8, 
                  resolution = 1)

#fortify for ggplot
bathy_fort <- fortify(bathy)

#make separate tibbles of block for labeling points
s1 <- sites %>% filter(block=="S" & stat <8)
s8 <- sites %>% filter(block=="S" & stat >7)
c1 <- sites %>% filter(block=="C" & stat <8)
c8 <- sites %>% filter(block=="C" & stat >7)
n1 <- sites %>% filter(block=="N" & stat <9)
n8 <- sites %>% filter(block=="N" & stat >8)

#make the zoomed in map of sites
 site_map <- ggplot(data = states) +
   geom_polygon(aes(x=long,y=lat))+
  #geom_path(data = sbnms_fort, 
            #aes(x = long, y = lat, group = group),
            #color = 'black', size = .2) +
   #coord_map(xlim = c(-71,-70),ylim = c(41.8,42.8)) + #just MA, comment out site labels
   coord_map(xlim = c(-70.5,-70),ylim = c(42.13,42.45)) +#just entire Bank
   
   labs(x="Longitude",y="Latitude") +
   #coord_cartesian(xlim = c(-70.4,-70.23),ylim = c(42.15,42.2)) + #SW corner
    #geom_raster(data=depth_mask_fort, aes(x=x,y=y,fill=mikesbathym)) +
   #scale_fill_gradientn(colors=cols) +
  geom_contour(data=bathy_fort,
                aes(x=x,y=y,z=z),
                breaks=c(-40,-50),
                color="grey") +
  geom_point(data=sites,aes(x=lon,y=lat,color=block, size=0.1)) + #for all points colored by block 0.1 for zoom in
   #geom_point(data=sw_w_p_aggreg,aes(x=lon,y=lat,color=paggregation,size=0.1)) +
   
   scale_color_viridis(discrete=TRUE,name="Block") +
   geom_text(data=s1,aes(x=lon,y=lat,label=site),size=3,vjust=-0.4) +
   geom_text(data=s8,aes(x=lon,y=lat,label=site),size=3,vjust=1.4) +
   geom_text(data=n1,aes(x=lon,y=lat,label=site),size=3,vjust=-0.4) +
   geom_text(data=n8,aes(x=lon,y=lat,label=site),size=3,vjust=1.4) +
   geom_text(data=c1,aes(x=lon,y=lat,label=site),size=3,hjust=-0.4) +
   geom_text(data=c8,aes(x=lon,y=lat,label=site),size=3,hjust=1.4) +
   theme_bw(base_size=18) +
   guides(size=FALSE) +
   scalebar(x.min=-70.5,x.max=-70.4, y.min=42.14,y.max=42.22, st.bottom=TRUE,
            height=0.07,st.dist=0.1, dist=4, dist_unit="km", transform=TRUE, model="WGS84")

 #north2(site_map, x=0.3,y=0.48, scale=0.1,symbol=3) #add a north arrow
 #cant use ggsave with north arrow bc can't add north arrow to gg object
 tiff(file="./figures/Bayesian_manuscript/sites_update.tiff",width=8,height=8, units="in", res=300)
 north2(site_map, x=0.2,y=0.31, scale=0.07,symbol=3) #add a north arrow
 dev.off()
    #add a scale bar
 #ggsave("./figures/Bayesian_manuscript/sites_update.tiff",width=8, height=8, units="in")
 
 ######################################################################




#####################################################################################







