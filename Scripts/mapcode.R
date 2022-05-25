### Making maps of the turbinaria data
### Nyssa Silbiger ####
### 5/25/2022 ####

## load library#####
library(here)
library(tidyverse)
library(patchwork)
library(ggmap)
library(viridis)
library(maptools)
library(kriging)
library(ggnewscale)
library(wql)
library(glue)
library(gridExtra)


## Read in data
# Turb data
AugData<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Nutrients/Turb_NC.csv") %>%
  mutate(Season = "Dry")

MarData<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/March2022/Nutrients/Turb_NC.csv")%>%
  mutate(Season = "Wet")

# Lat Long data
Locations<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/Sandwich_Locations_Final.csv")%>%
  select(Location, CowTagID, lat, lon, Plate_Seep)


# Bring all turb data together

AllTurb<-bind_rows(AugData,MarData)%>%
  left_join(Locations)


# Base Maps
#API<-names(read_table(here("Data","API.txt")))
#register_google(key = API) ### use your own API

VarariBaseMap<-get_map(AllTurb %>% filter(Location == "Varari", Plate_Seep == "Plate") %>% select(lon,lat), maptype = 'satellite', zoom = 18)

# Cabral
CabralBaseMap<-get_map(AllTurb %>% filter(Location == "Cabral", Plate_Seep == "Plate") %>% select(lon,lat), maptype = 'satellite', zoom = 18)

### Make a spatial kriging file with polygon layers ####
### Bring in the polygons for the sites
#Varari
V_kml <- getKMLcoordinates(kmlfile=here("Data","Polygons","Varari_Polygon.kml"), ignoreAltitude=T)
#Cabral
C_kml <- getKMLcoordinates(kmlfile=here("Data","Polygons","Cabral2.kml"), ignoreAltitude=T)


# make a function to do all the krigings
Krig_function <-function(dat_in = data, Lat = "lat", Lon = "lon", Param = "Values", poly ) {
  
  dat <- dat_in[,c(Lon, Lat, Param)]
  names(dat) <- c('Lon', 'Lat', 'Param') 
  # VData <- AllChemData %>%
  #   filter(Location == location,
  #          Tide == tide,
  #          Day_Night == day_night,
  #          Date == date,
  #          Plate_Seep == "Plate") %>%
  #   drop_na({{parameter}})
  
  dat<-dat%>%
    drop_na()
  
  x <- dat$Lon
  y <- dat$Lat
  z <-dat$Param
  
  krig1 <- kriging(x, y, z, pixels=500,polygons=poly, lags = 3) ###pixels controls how fine or course you want the prediction data frame to be
  krig2 <- krig1$map
  return(krig2)
}

# And do it "safely"
Krig_function_safe<-safely(Krig_function) # skip the NAs without breaking the code

# plot map function
V_krig_map<-function(datakrig=preds){
  
  ggmap(VarariBaseMap)+
    geom_point(data=datakrig, aes(x=x, y=y, colour=pred), size=4, alpha=0.5) + 
    # geom_point(data = VData, aes(x=lon, y=lat))+
    scale_color_viridis_c(" ", option = "plasma")+
    coord_sf() +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +  
    theme(panel.grid.major = element_line(color = 'white', linetype = "dashed",size = 0.5),
          plot.background=element_rect(fill='white'))+
    ggtitle(glue("Varari: {.y}")) 
  #   ggtitle(paste("Varari",DN, TD))
}


## August map

# nest by all parameters, tides, day/Night, Date, etc to make it easy to plot all types of maps
# Varari
Varari_kriging<-AllTurb %>%
  select(Location, CowTagID,Season,lat, lon, Plate_Seep, del15N, N_percent, C_N)%>% # this is temporary until we get the temperature data entered
  filter(Plate_Seep == "Plate", # only plot the plates because the seep samples skew the maps
         Location == "Varari") %>%
  droplevels()%>%
  pivot_longer(cols = del15N:C_N, names_to = "Parameters", values_to = "Values") %>%
  select(lat, lon,Season, Parameters, Values) %>% # select the values that are important for the kriging
  group_nest(Season, Parameters) %>% # the parameters to group by
  # left_join(min_max)%>% # add in the mins and max values for the plots
  mutate(preds = map(data, ~Krig_function_safe(dat_in = .x, poly = V_kml)), # run the function for every nested group
         #   preds = map(preds, head, -1), # remove the error column
         #  preds = map(preds, flatten_df), # flatten back to a tibble 
         # mutate(preds = unlist(preds))
         longname = paste(Parameters, Season),
         plots = map2(preds, longname, ~ggmap(VarariBaseMap)+
                        geom_point(data=.x$result, aes(x=x, y=y, colour=pred), size=4, alpha=0.5) + 
                        geom_point(data = Locations %>% filter(Location == 'Varari', Plate_Seep=="Plate"), aes(x=lon, y=lat))+
                        scale_color_viridis_c(" ", option = "plasma")+
                        coord_sf() +
                        theme(axis.line=element_blank(),
                              axis.text.x=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks=element_blank(),
                              axis.title.x=element_blank(),
                              axis.title.y=element_blank()) +  
                        theme(panel.grid.major = element_line(color = 'white', linetype = "dashed",size = 0.5),
                              plot.background=element_rect(fill='white'))+
                        ggtitle(glue("Varari: {.y}"))))


for(i in 1:length(Varari_kriging$plots)){
  try({
    ggsave(plot = Varari_kriging$plots[[i]], file = here("Output",paste0(Varari_kriging$longname[i],"_Varari",".png")))}, silent = TRUE)
}                     


### Cabral #####
Cabral_kriging<-AllTurb %>%
  select(Location, CowTagID,Season,lat, lon, Plate_Seep, del15N, N_percent, C_N)%>% # this is temporary until we get the temperature data entered
  filter(Plate_Seep == "Plate", # only plot the plates because the seep samples skew the maps
         Location == "Cabral") %>%
  droplevels()%>%
  pivot_longer(cols = del15N:C_N, names_to = "Parameters", values_to = "Values") %>%
  select(lat, lon,Season, Parameters, Values) %>% # select the values that are important for the kriging
  group_nest(Season, Parameters) %>% # the parameters to group by
  mutate(preds = map(data, ~Krig_function_safe(dat_in = .x, poly = C_kml)), # run the function for every nested group
        longname = paste(Parameters, Season),
         plots = map2(preds, longname, ~ggmap(CabralBaseMap)+
                        geom_point(data=.x$result, aes(x=x, y=y, colour=pred), size=4, alpha=0.5) + 
                        geom_point(data = Locations %>% filter(Location == 'Cabral', Plate_Seep=="Plate"), aes(x=lon, y=lat))+
                        scale_color_viridis_c(" ", option = "plasma")+
                        coord_sf() +
                        theme(axis.line=element_blank(),
                              axis.text.x=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks=element_blank(),
                              axis.title.x=element_blank(),
                              axis.title.y=element_blank()) +  
                        theme(panel.grid.major = element_line(color = 'white', linetype = "dashed",size = 0.5),
                              plot.background=element_rect(fill='white'))+
                        ggtitle(glue("Cabral: {.y}"))))

for(i in 1:length(Cabral_kriging$plots)){
  try({
    ggsave(plot = Cabral_kriging$plots[[i]], file = here("Output",paste0(Cabral_kriging$longname[i],"_Cabral",".png")))}, silent = TRUE)
}                     


#### Comparisons of turb data from August to March

TurbLong<-AllTurb %>%
  select(Location, CowTagID,Season,lat, lon, Plate_Seep, del15N, N_percent, C_N)%>% # this is temporary until we get the temperature data entered
  pivot_longer(cols = del15N:C_N, names_to = "Parameters", values_to = "Values") 

shapes <-c("Plates" = 21, "Seeps" = 18, "Reef Crest" = 15)

ggplot(data = TurbLong %>% filter(Plate_Seep == "Plate"), aes(x = Location, y = Values, color = Season))+
  geom_jitter(aes(x = Location, y = Values, color = Season), width = 0.2)+
  geom_point(data = TurbLong %>% filter(Plate_Seep == "Seep"), aes(x = Location, y = Values), size = 5, shape = 18)+
  geom_point(data = TurbLong %>% filter(Plate_Seep == "Offshore"), aes(x = Location, y = Values), size = 5, shape = 15)+
  scale_shape_manual(name = " ",values = shapes)+
  labs(shape = "type")+
  facet_wrap(~Parameters, scales = "free")

   
