# A script to create PCAs for 24 hour biogeochem data for Cabral and Varari
# Edited on 12/21/2022
# Created by Nyssa Silbiger 

####################################

# load libraries #################
library(tidyverse)
library(ggfortify)
library(lubridate)
library(ggforce)
library(viridis)
library(patchwork)
library(here)
library(wesanderson)
library(broom)
library(ggtext)
library(glue)
library(htmlTable)

# load the 24 hour chemistry data #####################
#Data<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv")
Data<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC2.csv") %>% mutate(Season = "Dry")

Data_march<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/March2022/Allbiogeochemdata_QC_march.csv")%>% mutate(Season = "Wet") %>%
  mutate(DateTime = mdy_hms(paste(Date,as.character(Time))))

## Load the turb nutrient data
turbdata<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Nutrients/Turb_NC.csv") %>% mutate(Season = "Dry")

turbdata_march<- read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/March2022/Nutrients/Turb_NC.csv") %>% mutate(Season = "Wet")

turb_all<-bind_rows(turbdata, turbdata_march)

# Load the Seep Data
VarariSeep<-read_csv(here("Data","Varari","AllVarariSeepData.csv")) %>%
  rename_with(.cols = TempInSitu:PAR_calc,function(x){paste0(x,"_seep")}) %>% # rename columns to say seep at the end
  rename(DateTime = date, Location = Site, Season = Season_seep)

CabralSeep<-read_csv(here("Data","Cabral","AllCabralSeepData.csv")) %>%
  rename_with(.cols = TempInSitu:PAR_calc,function(x){paste0(x,"_seep")}) %>% # rename columns to say seep at the end
  rename(DateTime = date, Location = Site, Season = Season_seep)

### Danielles community composition data
CommData<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/June2022/BenthicData/Full_Metadata_DMB.csv") %>%
  select(Location,CowTagID, Rubble, Sand,LiveCoral,DeadCoral, dist_to_seep_m, adj_CT_depth_cm, meanRugosity) %>%
  drop_na()

# bind them together
SeepAll<-bind_rows(VarariSeep, CabralSeep)

# Bind with the chem data
Data<-Data %>%
  bind_rows(Data_march)%>%
  left_join(SeepAll) %>%
  left_join(CommData) %>%
  mutate(NewDay_Night = case_when( 
    Tide == "Low" & Day_Night == "Day" ~ "Night", # change day and night delineations... dusk and dawn are now night and day 
    Tide == "Low" & Day_Night == "Night" ~ "Day",
    Tide == "High" & Day_Night == "Day" ~ "Day",
    Tide == "High" & Day_Night == "Night" ~ "Night"
  ),
  Day_Night = NewDay_Night) %>% # replace it
  select(!NewDay_Night)


# # Load water level data 
# # Bring in the depth data
# WLPath<-here("Data", "Varari", "WL")
# files <- dir(path = WLPath,pattern = ".csv", full.names = TRUE)
# 
# WL_Varari<-files %>%
#   set_names()%>% # set's the id of each list to the file name
#   map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
#   select(Date = date, Depth)%>%
#   mutate(Location = "Varari")
# 
# WLPath<-here("Data", "Cabral", "WL")
# files <- dir(path = WLPath,pattern = ".csv", full.names = TRUE)

# #Cabral
# WL_Cabral<-files %>%
#   set_names()%>% # set's the id of each list to the file name
#   map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
#   select(Date = date, Depth)%>%
#   mutate(Location = "Cabral")

# WL_all <- bind_rows(WL_Varari, WL_Cabral) %>%
#   rename(DateTime = Date)

# # Bind with the chem data
# Data<-Data %>%
#   left_join(WL_all) 

# drop dataframes I dont need
rm(CabralSeep,VarariSeep, SeepAll)

### Varari #####
## There seems to be a contaminated nutrient sample for V2 Low tide on the 8/8/2021.  Remove this point
remove<-Data %>% filter(CowTagID=="V2", Tide == "Low",Day_Night=="Night", Date =="8/8/2021") 
removea<-Data %>% filter(CowTagID == "V17", Tide == "High", Day_Night == "Night", Season == "Wet")# Ammonia is an outlier

## also filter out all the data from the first low tide where the water level was super high
remove2<-Data %>% filter(Tide =="Low", Day_Night=="Night", Date == "8/6/2021")
remove_varari<-bind_rows(remove, removea, remove2) # bring the "bad data" into one tibble
 
# fdom or silicate data off for these
remove3<-Data %>% filter(CowTagID== "C4", Tide =="Low", Day_Night=="Day", Date == "8/9/2021")
remove4<-Data %>% filter(CowTagID== "C2", Tide =="High", Day_Night=="Night", Date == "3/31/2022")
remove5<-Data %>% filter(CowTagID== "C17", Tide =="Low", Day_Night=="Day", Date == "3/30/2022")

remove_cabral<- bind_rows(remove3, remove4,remove5) # bring the "bad data" into one tibble




## Make histograms of all the data 
# Varari
Data %>% filter(Plate_Seep == "Plate", Location == "Varari") %>%
  anti_join(remove_varari) %>%
  select(Location, Day_Night, Tide, Salinity:Season) %>%
  pivot_longer(Salinity:Lignin_Like)%>%
  ggplot(aes(x = value, fill =  paste(Day_Night, Tide)))+
  geom_density(alpha = 0.2)+
  facet_wrap(~name*Season, scales = "free")

# Cabral
Data %>% filter(Plate_Seep == "Plate", Location == "Cabral") %>%
  anti_join(remove_cabral) %>%
  select(Location, Day_Night, Tide, Salinity:Season) %>%
  pivot_longer(Salinity:Lignin_Like)%>%
  ggplot(aes(x = value, fill =  paste(Day_Night, Tide)))+
  geom_density(alpha = 0.2)+
  facet_wrap(~name*Season, scales = "free")


## create log transformed data for everything except for salinity, temperature, TA, and pH
Datalog<-Data %>%
  mutate_at(vars(Phosphate_umolL:Lignin_Like), .funs = function(x){log(x+0.001)})

remove_vararilog<-remove_varari %>% mutate_at(vars(Phosphate_umolL:Lignin_Like), .funs = function(x){log(x+0.001)})
remove_cabrallog<-remove_cabral %>% mutate_at(vars(Phosphate_umolL:Lignin_Like), .funs = function(x){log(x+0.001)})

# view it again-- a bit more normally distributes 
# Varari
Datalog %>% filter(Plate_Seep == "Plate", Location == "Varari") %>%
  anti_join(remove_vararilog) %>%
  select(Location, Day_Night, Tide, Salinity:Season) %>%
  pivot_longer(Salinity:Lignin_Like)%>%
  ggplot(aes(x = value, fill =  paste(Day_Night, Tide)))+
  geom_density(alpha = 0.2)+
  facet_wrap(~name*Season, scales = "free")

# Cabral
Datalog %>% filter(Plate_Seep == "Plate", Location == "Cabral") %>%
  anti_join(remove_cabrallog) %>%
  select(Location, Day_Night, Tide, Salinity:Season) %>%
  pivot_longer(Salinity:Lignin_Like)%>%
  ggplot(aes(x = value, fill =  paste(Day_Night, Tide)))+
  geom_density(alpha = 0.2)+
  facet_wrap(~name*Season, scales = "free")

# extract the params for the PCA
V_pca_Data_wet<-Datalog %>%
  filter(Season == "Wet") %>%
  anti_join(remove_vararilog)%>%
  #anti_join(remove2)%>%
 # filter(Location == "Varari", Tide %in% c("High","Low")) %>%
  # filter(DateTime %in% c(ymd_hms("2021-08-05 11:57:00"),ymd_hms("2021-08-05 00:00:00"),ymd_hms("2021-08-08 18:30:00"), ymd_hms("2021-08-08 07:30:00"), ymd_hms("2021-08-04 23:51:00"))) %>%
   filter(Location == "Varari", Plate_Seep=="Plate") %>%
  #select(Salinity,pH,Phosphate_umolL, Silicate_umolL, NN_umolL, Ammonia_umolL ) %>%
  select(Salinity,pH,Phosphate_umolL:NN_umolL,Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like, MarineHumic_Like, TA)%>%
  drop_na()

V_pca_Data<-Datalog %>%
  filter(Season == "Dry") %>%
  anti_join(remove_vararilog)%>%
 # anti_join(remove2)%>%
  # filter(Location == "Varari", Tide %in% c("High","Low")) %>%
  # filter(DateTime %in% c(ymd_hms("2021-08-05 11:57:00"),ymd_hms("2021-08-05 00:00:00"),ymd_hms("2021-08-08 18:30:00"), ymd_hms("2021-08-08 07:30:00"), ymd_hms("2021-08-04 23:51:00"))) %>%
  filter(Location == "Varari", Plate_Seep=="Plate") %>%
  #select(Salinity,pH,Phosphate_umolL, Silicate_umolL, NN_umolL, Ammonia_umolL ) %>%
  select(Salinity,pH,Phosphate_umolL:NN_umolL,Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like, MarineHumic_Like,TA)%>%
  drop_na()

# Run the PCA
pca_V <- prcomp(V_pca_Data, scale. = TRUE, center = TRUE)
pca_V_wet <- prcomp(V_pca_Data_wet, scale. = TRUE, center = TRUE)

#pca_V <- prcomp(V_pca_Data, scale. = TRUE, center = TRUE)

# calculate percent explained by each PC
perc.explained<-round(100*pca_V$sdev/sum(pca_V$sdev),1)
perc.explained_wet<-round(100*pca_V_wet$sdev/sum(pca_V_wet$sdev),1)

# Extract the scores and loadings
PC_scores <-as_tibble(pca_V$x[,1:2])
PC_scores_wet <-as_tibble(pca_V_wet$x[,1:2])

PC_loadings<-as_tibble(pca_V$rotation) %>%
  bind_cols(labels = rownames(pca_V$rotation))%>%
  mutate(groupings = case_when( # add groupings
    labels %in% c("Ammonia_umolL","NN_umolL","Phosphate_umolL","Silicate_umolL")~ "Nutrient Chemistry",
    labels == "Salinity" ~ "Salinity",
    labels %in% c("pH","TA") ~ "Carbonate Chemistry",
    labels %in% c("HIX","Tryptophan_Like","Tyrosine_Like","VisibleHumidic_Like", "MarineHumic_Like")~"fDOM"
  ),
  nicenames = case_when(labels == "TempInSitu_seep" ~ "Temperature",
                        labels == "pH" ~ "pH<sub>T</sub>",
                     #   labels == "Lignin_Like" ~"Lignin Like",
                    #    labels == "M_C" ~ "M:C",
                        labels == "Tyrosine_Like" ~ "Tyrosine Like",
                        labels == "Tryptophan_Like" ~ "Tryptophan Like",
                        labels == "HIX"~"HIX",
                        labels == "MarineHumic_Like" ~ "Marine Humic Like",
                        labels == "VisibleHumidic_Like" ~ "Visible Humic Like",
                        labels == "Ammonia_umolL" ~ "Ammonium",
                        labels == "TA" ~ "Total Alkalinity",
                        labels == "Phosphate_umolL" ~ "Phosphate",
                        labels == "NN_umolL" ~ "Nitrate+Nitrite",
                        labels == "Silicate_umolL" ~ "Silicate",
                        labels == "Salinity" ~"Salinity"))

PC_loadings_wet<-as_tibble(pca_V_wet$rotation) %>%
  bind_cols(labels = rownames(pca_V_wet$rotation))%>%
  mutate(groupings = case_when( # add groupings
    labels %in% c("Ammonia_umolL","NN_umolL","Phosphate_umolL","Silicate_umolL")~ "Nutrient Chemistry",
    labels == "Salinity" ~ "Salinity",
    labels %in% c("pH","TA") ~ "Carbonate Chemistry",
    labels %in% c("HIX","Tryptophan_Like","Tyrosine_Like","VisibleHumidic_Like", "MarineHumic_Like")~"fDOM"
  ),
  nicenames = case_when(labels == "TempInSitu_seep" ~ "Temperature",
                        labels == "pH" ~ "pH<sub>T</sub>",
                        #   labels == "Lignin_Like" ~"Lignin Like",
                        #    labels == "M_C" ~ "M:C",
                        labels == "Tyrosine_Like" ~ "Tyrosine Like",
                        labels == "Tryptophan_Like" ~ "Tryptophan Like",
                        labels == "HIX"~"HIX",
                        labels == "MarineHumic_Like" ~ "Marine Humic Like",
                        labels == "VisibleHumidic_Like" ~ "Visible Humic Like",
                        labels == "Ammonia_umolL" ~ "Ammonium",
                        labels == "TA" ~ "Total Alkalinity",
                        labels == "Phosphate_umolL" ~ "Phosphate",
                        labels == "NN_umolL" ~ "Nitrate+Nitrite",
                        labels == "Silicate_umolL" ~ "Silicate",
                        labels == "Salinity" ~"Salinity"))

# Put it with all the original data
V_pca_Data_all_wet<-Data %>%
  anti_join(remove_varari) %>%
  filter(Season == "Wet") %>%
  filter(Location == "Varari", Plate_Seep=="Plate") %>%
  drop_na(Salinity,pH,Phosphate_umolL:Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like)%>%
  bind_cols(PC_scores_wet)

V_pca_Data_all<-Data %>%
  filter(Season == "Dry") %>%
  anti_join(remove_varari)%>%
  #anti_join(remove2)%>%
  filter(Location == "Varari", Plate_Seep=="Plate") %>%
  drop_na(Salinity,pH,Phosphate_umolL:Lignin_Like )%>%
  bind_cols(PC_scores)

# scores plot
p1<-V_pca_Data_all %>%
  ggplot(aes(x = PC1, y = PC2, color = Tide, shape = Day_Night))+
  geom_point(size = 2) +
 # geom_point(data = V_pca_Data_all %>% filter(Plate_Seep=="Seep"), aes(x = PC1, y = PC2,shape = Day_Night ), color = "black")+
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  scale_shape_manual(values = c(22,16))+
  scale_colour_hue(l = 45)+
  scale_fill_hue(l = 45)+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  ggforce::geom_mark_ellipse(
    aes(fill = Tide, label = paste(Day_Night, Tide), color =Tide), 
    alpha = .15, show.legend = FALSE,  label.buffer = unit(1, "mm"))+
  labs(title = "Dry",
       x = paste0("PC1 ","(",perc.explained[1],"%)"),
       y = paste0("PC2 ","(",perc.explained[2],"%)"))+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18))

p1_wet<-V_pca_Data_all_wet %>%
  ggplot(aes(x = -PC1, y = -PC2, color = Tide, shape = Day_Night))+
  geom_point(size = 2) +
  # geom_point(data = V_pca_Data_all %>% filter(Plate_Seep=="Seep"), aes(x = PC1, y = PC2,shape = Day_Night ), color = "black")+
  coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
  scale_shape_manual(values = c(22,16))+
  scale_colour_hue(l = 45)+
  scale_fill_hue(l = 45)+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  ggforce::geom_mark_ellipse(
    aes(fill = Tide, label = paste(Day_Night, Tide), color =Tide), 
    alpha = .15, show.legend = FALSE,  label.buffer = unit(1, "mm"))+
  labs(title = "Wet",
       x = paste0("PC1 ","(",perc.explained_wet[1],"%)"),
       y = paste0("PC2 ","(",perc.explained_wet[2],"%)"))+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18))


# loadings plot 
p2<-PC_loadings %>%
  ggplot(aes(x=PC1, y=PC2, label=nicenames, color = groupings))+
  geom_richtext(aes(x = PC1*10+0.1, y = PC2*10+.1 ), show.legend = FALSE, size = 5, fill=NA, label.colour = NA) +
  geom_segment(data = PC_loadings, aes(x=0,y=0,xend=PC1*10,yend=PC2*10),size = 1.2,
      arrow=arrow(length=unit(0.1,"cm")))+
   # annotate("text", x = PC_loadings$PC1*10+0.1, y = PC_loadings$PC2*10+.1,
  #          label = PC_loadings$labels)+
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  labs(color ="")+
  scale_color_manual(values = wes_palette("Darjeeling1"))+
   theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.75, 0.85),
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

p2_wet<-PC_loadings_wet %>%
  ggplot(aes(x=-PC1, y=-PC2, label=nicenames, color = groupings))+
  geom_richtext(aes(x = -PC1*10+0.1, y = -PC2*10+.1 ), show.legend = FALSE, size = 5, fill=NA, label.colour = NA) +
  geom_segment(data = PC_loadings_wet, aes(x=0,y=0,xend=-PC1*10,yend=-PC2*10),size = 1.2,
               arrow=arrow(length=unit(0.1,"cm")))+
  # annotate("text", x = PC_loadings$PC1*10+0.1, y = PC_loadings$PC2*10+.1,
  #          label = PC_loadings$labels)+
  coord_cartesian(xlim = c(-8, 8), ylim = c(-6, 6)) +
  labs(color ="",
       x = "PC1",
       y = "PC2")+
  scale_color_manual(values = wes_palette("Darjeeling1"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
       # legend.position = c(0.75, 0.85),
        legend.position = "none",
      #  legend.text = element_markdown(size = 16),
      #  legend.key.size = unit(1, 'cm'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

VarariPCA<-(p1+p1_wet)/(p2+p2_wet)+ 
  patchwork::plot_annotation(#"Varari Plates", 
                             theme = theme(plot.title = element_text(size = rel(1.5), 
                                                                     face = "bold", hjust = 0.5, 
                                                                     margin = margin(t = 10, b = 20,
                                                                                     unit = "pt"))))

ggsave(plot = VarariPCA, filename = here("Output","VarariPCA.png"), width = 16, height = 16)

### Site level pca with ranges
V_pca_Data_site<-Datalog %>%
  anti_join(remove_vararilog)%>%
 # anti_join(remove2)%>%
  filter(Location == "Varari",
         Plate_Seep=="Plate") %>%
  group_by(CowTagID, Season)%>% # calculate the range by cowtag
  summarise_at(vars(Salinity,pH,TA,Phosphate_umolL:Lignin_Like), .funs = function(x) {max(x, na.rm = TRUE)-min(x, na.rm = TRUE)}) %>%
 # summarise_at(vars(Salinity,pH,Phosphate_umolL:Lignin_Like), .funs = function(x) {var(x, na.rm = TRUE)/mean(x, na.rm = TRUE)}) %>%
  #summarise_at(vars(Salinity,pH,Phosphate_umolL:Lignin_Like), .funs = function(x) {mean(x, na.rm = TRUE)}) %>%
  ungroup()%>%
  left_join(bind_rows(turbdata, turbdata_march))%>%
  select(CowTagID, Season, Salinity,pH,Phosphate_umolL, Silicate_umolL, NN_umolL, Ammonia_umolL, del15N, N_percent,VisibleHumidic_Like, MarineHumic_Like,Tryptophan_Like, Tyrosine_Like,TA ) %>%
 # select(Salinity,pH,Phosphate_umolL:Lignin_Like )%>%
  drop_na()

# Run the PCA
pca_V_site <- prcomp(V_pca_Data_site[,-c(1:2)], scale. = TRUE, center = TRUE)
# Extract the scores and loadings
PC_scores <-as_tibble(pca_V_site$x[,1:2])
PC_loadings<-as_tibble(pca_V_site$rotation) %>%
  bind_cols(labels = rownames(pca_V_site$rotation))

V_summaryplot<-PC_scores %>%
  bind_cols(V_pca_Data_site)%>%
  left_join(CommData)%>%
  ggplot(aes(x = PC1, y = PC2))+
 # geom_point(color = "red") +
  geom_point(aes(x = PC1, y = PC2, color = Season, size = -dist_to_seep_m), alpha = 0.2)+
  geom_text(aes(x = PC1, y = PC2,label = CowTagID, color = Season), show.legend = FALSE)+
  geom_segment(data = PC_loadings, aes(x=0,y=0,xend=PC1*10,yend=PC2*10),
               arrow=arrow(length=unit(0.1,"cm")), color = "grey")+
  annotate("text", x = PC_loadings$PC1*10+0.1, y = PC_loadings$PC2*10+.1,
           label = PC_loadings$labels)+
  labs(title = "Ranges for each plate",
       size = "Distance to Seep")+
  
  theme_bw()+
  theme(#legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())



### some plots of the ranges
V_pca_Data_site %>%
left_join(CommData)%>%
ggplot( aes(x = Silicate_umolL, y = NN_umolL, label = CowTagID, color = -dist_to_seep_m))+
  geom_point()+
  geom_label()+
  geom_smooth(method = "lm")+
  labs(x = "Silicate range",
       y = "NN range")+
  theme_bw()+
  facet_wrap(~Season)

V_pca_Data_site %>%
  left_join(CommData)%>%
ggplot(aes(y = N_percent, x = NN_umolL, label = CowTagID, color = -dist_to_seep_m))+
  geom_point()+
  geom_label()+
  geom_smooth(method = "lm")+
  labs(y = "%Tissue N from turbinaria",
       x = "NN diel range")+
  theme_bw()+
  facet_wrap(~Season)

V_pca_Data_site %>%
  left_join(CommData)%>%
ggplot( aes(y = Salinity, x = Silicate_umolL, color = -dist_to_seep_m))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(y = "Salinity",
       x = "Silicate diel range")+
  theme_bw()+
  facet_wrap(~Season)

V_pca_Data_site %>%
  left_join(CommData)%>%
ggplot(aes(x = Ammonia_umolL, y = pH))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(x = "Ammonia diel",
       y = "pH diel range")+
  theme_bw()+
  facet_wrap(~Season)

#### Cabral #####

#Extract the cabral data
C_pca_Data<-Datalog %>%
  anti_join(remove_cabrallog)%>%
  filter(Location == "Cabral", Plate_Seep=="Plate",
         Season == "Dry") %>%
  select(Salinity,pH,TA,Phosphate_umolL:NN_umolL,Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like, MarineHumic_Like)%>%
  #select(Salinity,pH,Phosphate_umolL:Lignin_Like )%>%
  #select(Salinity,pH,Phosphate_umolL, Silicate_umolL, NN_umolL, Ammonia_umolL ) %>%
  drop_na()

C_pca_Data_wet<-Datalog %>%
  anti_join(remove_cabrallog)%>%
 # anti_join(remove5)%>%
  filter(Location == "Cabral", Plate_Seep=="Plate",
         Season == "Wet") %>%
  select(Salinity,pH,TA,Phosphate_umolL:NN_umolL,Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like, MarineHumic_Like)%>%
  #select(Salinity,pH,Phosphate_umolL:Lignin_Like )%>%
  #select(Salinity,pH,Phosphate_umolL, Silicate_umolL, NN_umolL, Ammonia_umolL ) %>%
  drop_na()
# Run the PCA

pca_C <- prcomp(C_pca_Data, scale. = TRUE, center = TRUE)
pca_C_wet <- prcomp(C_pca_Data_wet, scale. = TRUE, center = TRUE)

# calculate percent explained by each PC
perc.explainedC<-round(100*pca_C$sdev/sum(pca_C$sdev),1)
perc.explainedC_wet<-round(100*pca_C_wet$sdev/sum(pca_C_wet$sdev),1)


# Extract the scores and loadings
PC_scoresC <-as_tibble(pca_C$x[,1:2])
PC_scoresC_wet <-as_tibble(pca_C_wet$x[,1:2])

PC_loadingsC<-as_tibble(pca_C$rotation) %>%
  bind_cols(labels = rownames(pca_C$rotation))%>%
  mutate(groupings = case_when( # add groupings
    labels %in% c("Ammonia_umolL","NN_umolL","Phosphate_umolL","Silicate_umolL")~ "Nutrient Chemistry",
    labels == "Salinity" ~ "Salinity",
    labels %in% c("TA","pH") ~ "Carbonate Chemistry",
    labels %in% c("HIX","Lignin_Like","M_C","MarineHumic_Like","Tryptophan_Like","Tyrosine_Like","VisibleHumidic_Like","MarineHumic_Like")~"fDOM"
  ),
  nicenames = case_when(labels == "TempInSitu_seep" ~ "Temperature",
                        labels == "pH" ~ "pH<sub>T</sub>",
                        #   labels == "Lignin_Like" ~"Lignin Like",
                        #    labels == "M_C" ~ "M:C",
                        labels == "Tyrosine_Like" ~ "Tyrosine Like",
                        labels == "Tryptophan_Like" ~ "Tryptophan Like",
                        labels == "HIX"~"HIX",
                        labels == "MarineHumic_Like" ~ "Marine Humic Like",
                        labels == "VisibleHumidic_Like" ~ "Visible Humic Like",
                        labels == "Ammonia_umolL" ~ "Ammonium",
                        labels == "TA" ~ "Total Alkalinity",
                        labels == "Phosphate_umolL" ~ "Phosphate",
                        labels == "NN_umolL" ~ "Nitrate+Nitrite",
                        labels == "Silicate_umolL" ~ "Silicate",
                        labels == "Salinity" ~"Salinity"))

PC_loadingsC_wet<-as_tibble(pca_C_wet$rotation) %>%
  bind_cols(labels = rownames(pca_C_wet$rotation))%>%
  mutate(groupings = case_when( # add groupings
    labels %in% c("Ammonia_umolL","NN_umolL","Phosphate_umolL","Silicate_umolL")~ "Nutrient Chemistry",
    labels == "Salinity" ~ "Salinity",
    labels  %in% c("TA","pH") ~ "Carbonate Chemistry",
    labels %in% c("HIX","Lignin_Like","M_C","MarineHumic_Like","Tryptophan_Like","Tyrosine_Like","VisibleHumidic_Like","MarineHumic_Like")~"fDOM"
  ),
  nicenames = case_when(labels == "TempInSitu_seep" ~ "Temperature",
                        labels == "pH" ~ "pH<sub>T</sub>",
                        #   labels == "Lignin_Like" ~"Lignin Like",
                        #    labels == "M_C" ~ "M:C",
                        labels == "Tyrosine_Like" ~ "Tyrosine Like",
                        labels == "Tryptophan_Like" ~ "Tryptophan Like",
                        labels == "HIX"~"HIX",
                        labels == "MarineHumic_Like" ~ "Marine Humic Like",
                        labels == "VisibleHumidic_Like" ~ "Visible Humic Like",
                        labels == "Ammonia_umolL" ~ "Ammonium",
                        labels == "TA" ~ "Total Alkalinity",
                        labels == "Phosphate_umolL" ~ "Phosphate",
                        labels == "NN_umolL" ~ "Nitrate+Nitrite",
                        labels == "Silicate_umolL" ~ "Silicate",
                        labels == "Salinity" ~"Salinity"))

# Put it with all the original data
C_pca_Data_all<-Data %>%
  anti_join(remove_cabral)%>%
  select(!Jamie_Plate_ID)%>% # Jamie's plates are all NA here
  filter(Location == "Cabral", Plate_Seep=="Plate",
         Season == "Dry") %>%
  drop_na(Salinity,pH,Phosphate_umolL:NN_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like) %>%
  bind_cols(PC_scoresC)  

C_pca_Data_all_wet<-Data %>%
  anti_join(remove_cabral)%>%
  select(!Jamie_Plate_ID)%>% # Jamie's plates are all NA here
  filter(Location == "Cabral", Plate_Seep=="Plate",
         Season == "Wet") %>%
  drop_na(Salinity,pH,Phosphate_umolL:NN_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like) %>%
  bind_cols(PC_scoresC_wet)  

# scores plot
p1c<-C_pca_Data_all %>%
  ggplot(aes(x = PC1, y = PC2, color = Tide, shape = Day_Night))+
  geom_point() +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  scale_shape_manual(values = c(22,16))+
  scale_colour_hue(l = 45)+
  scale_fill_hue(l = 45)+
  labs(title = "Dry",
       x = paste0("PC1 ","(",perc.explainedC[1],"%)"),
       y = paste0("PC2 ","(",perc.explainedC[2],"%)"))+
  ggforce::geom_mark_ellipse(
    aes(fill = Tide, label = paste(Day_Night, Tide), color = Tide), 
    alpha = .15, show.legend = FALSE,  label.buffer = unit(1, "mm"))+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18))

p1c_wet<-C_pca_Data_all_wet %>%
  ggplot(aes(x = PC1, y = -PC2, color = Tide, shape = Day_Night))+
  geom_point() +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  scale_shape_manual(values = c(22,16))+
  scale_colour_hue(l = 45)+
  scale_fill_hue(l = 45)+
  labs(title = "Wet",
       x = paste0("PC1 ","(",perc.explainedC_wet[1],"%)"),
       y = paste0("PC2 ","(",perc.explainedC_wet[2],"%)"))+
  ggforce::geom_mark_ellipse(
    aes(fill = Tide, label = paste(Day_Night, Tide), color = Tide), 
    alpha = .15, show.legend = FALSE,  label.buffer = unit(1, "mm"))+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18))

# scores plot with depth and light as continuous instead of discrete... missing depth data from lowtide at night :(
p1c_DL<-C_pca_Data_all %>%
  ggplot(aes(x = PC1, y = PC2, color = PAR_calc_seep+1, size = -Depth_seep))+
  geom_point() +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  scale_color_gradient(low = "black", high = "yellow", trans = "log")+
  theme_bw()+
  theme(#legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

# loadings plot 
p2c<-PC_loadingsC %>%
  ggplot(aes(x=PC1, y=PC2, label=nicenames, color = groupings))+
  geom_segment(aes(x=0,y=0,xend=PC1*10,yend=PC2*10),
               arrow=arrow(length=unit(0.1,"cm")))+
#  geom_text(aes(x = PC1*10+0.1, y = PC2*10+.1 ), show.legend = FALSE) +
  geom_richtext(aes(x = PC1*10+0.1, y = PC2*10+.1 ), show.legend = FALSE, size = 5, fill=NA, label.colour = NA)+
  scale_color_manual(values = wes_palette("Darjeeling1"))+
  #  annotate("text", x = PC_loadingsC$PC1*10+0.1, y = PC_loadingsC$PC2*10+.1,
  #        label = PC_loadingsC$labels)+
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  labs(color = "")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.75, 0.20),
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

p2c_wet<-PC_loadingsC_wet %>%
  ggplot(aes(x=PC1, y=-PC2, label=nicenames, color = groupings))+
  geom_segment(aes(x=0,y=0,xend=PC1*10,yend=-PC2*10),
               arrow=arrow(length=unit(0.1,"cm")))+
#  geom_text(aes(x = PC1*10+0.1, y = -PC2*10+.1 ), show.legend = FALSE) +
  geom_richtext(aes(x = PC1*10+0.1, y = -PC2*10+.1 ), show.legend = FALSE, size = 5, fill=NA, label.colour = NA)+
  scale_color_manual(values = wes_palette("Darjeeling1"))+
  #  annotate("text", x = PC_loadingsC$PC1*10+0.1, y = PC_loadingsC$PC2*10+.1,
  #        label = PC_loadingsC$labels)+
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  labs(color = "")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        #legend.position = c(0.75, 0.20)
        legend.position = "none"
        )

CabralPCA<-(p1c+p1c_wet)/(p2c+p2c_wet)+ 
  patchwork::plot_annotation(#"Cabral Plates", 
                             theme = theme(plot.title = element_text(size = rel(1.5), face = "bold", 
                                                                     hjust = 0.5, 
                                                                     margin = margin(t = 10, b = 20, 
                                                                                     unit = "pt"))))

ggsave(plot = CabralPCA, filename = here("Output","CabralPCA.png"), width = 16, height = 16)


#### Add Cabral Site level here ####
C_pca_Data_site<-Datalog %>%
  anti_join(remove_cabrallog)%>%
  # anti_join(remove2)%>%
  filter(Location == "Cabral",
         Plate_Seep=="Plate") %>%
  group_by(CowTagID, Season)%>% # calculate the range by cowtag
  summarise_at(vars(Salinity,pH,TA,Phosphate_umolL:Lignin_Like), .funs = function(x) {max(x, na.rm = TRUE)-min(x, na.rm = TRUE)}) %>%
  # summarise_at(vars(Salinity,pH,Phosphate_umolL:Lignin_Like), .funs = function(x) {var(x, na.rm = TRUE)/mean(x, na.rm = TRUE)}) %>%
  #summarise_at(vars(Salinity,pH,Phosphate_umolL:Lignin_Like), .funs = function(x) {mean(x, na.rm = TRUE)}) %>%
  ungroup()%>%
  left_join(bind_rows(turbdata, turbdata_march))%>%
  select(CowTagID, Season, Salinity,TA,pH,Phosphate_umolL, Silicate_umolL, NN_umolL, Ammonia_umolL, del15N, N_percent,VisibleHumidic_Like, MarineHumic_Like,Tryptophan_Like, Tyrosine_Like ) %>%
  # select(Salinity,pH,Phosphate_umolL:Lignin_Like )%>%
  drop_na()

# Run the PCA
pca_C_site <- prcomp(C_pca_Data_site[,-c(1:2)], scale. = TRUE, center = TRUE)
# Extract the scores and loadings
PC_scores <-as_tibble(pca_C_site$x[,1:2])
PC_loadings<-as_tibble(pca_C_site$rotation) %>%
  bind_cols(labels = rownames(pca_C_site$rotation))

C_summaryplot<-PC_scores %>%
  bind_cols(C_pca_Data_site)%>%
  left_join(CommData)%>%
  ggplot(aes(x = PC1, y = PC2))+
  # geom_point(color = "red") +
  geom_point(aes(x = PC1, y = PC2,label = CowTagID, color = Season, size = -dist_to_seep_m), alpha = 0.2)+
  geom_text(aes(x = PC1, y = PC2,label = CowTagID, color = Season), show.legend = FALSE)+
  geom_segment(data = PC_loadings, aes(x=0,y=0,xend=PC1*10,yend=PC2*10),
               arrow=arrow(length=unit(0.1,"cm")), color = "grey")+
  annotate("text", x = PC_loadings$PC1*10+0.1, y = PC_loadings$PC2*10+.1,
           label = PC_loadings$labels)+
  labs(title = "Ranges for each plate",
       size = "Distance to Seep")+
  
  theme_bw()+
  theme(#legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())

V_summaryplot+C_summaryplot + plot_annotation(tag_levels = "A")
ggsave(here("Output","SitePCASummary.png"), width = 12, height = 5)

### plot some of the ranges
### some plots of the ranges
C_pca_Data_site %>%
  left_join(CommData)%>%
  ggplot( aes(x = Silicate_umolL, y = NN_umolL, label = CowTagID, color = -dist_to_seep_m))+
  geom_point()+
  geom_label()+
  geom_smooth(method = "lm")+
  labs(x = "Silicate range",
       y = "NN range")+
  theme_bw()+
  facet_wrap(~Season)

C_pca_Data_site %>%
  left_join(CommData)%>%
  ggplot(aes(y = N_percent, x = NN_umolL, label = CowTagID, color = -dist_to_seep_m))+
  geom_point()+
  geom_label()+
  geom_smooth(method = "lm")+
  labs(y = "%Tissue N from turbinaria",
       x = "NN diel range")+
  theme_bw()+
  facet_wrap(~Season)

C_pca_Data_site %>%
  left_join(CommData)%>%
  ggplot( aes(y = Salinity, x = Silicate_umolL, color = -dist_to_seep_m))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(y = "Salinity",
       x = "Silicate diel range")+
  theme_bw()+
  facet_wrap(~Season)

C_pca_Data_site %>%
  left_join(CommData)%>%
  ggplot(aes(x = Ammonia_umolL, y = pH))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(x = "Ammonia diel",
       y = "pH diel range")+
  theme_bw()+
  facet_wrap(~Season)



#### Seep PCA #####
#Varari
V_pca_Seep<-Datalog %>%
  filter(Plate_Seep == "Seep" |Plate_Seep == "Spring")%>% # anchor by the spring because water was sampled at seep differently between seasons and sites
  filter(Location == "Varari", 
        # Plate_Seep=="Seep",
         Season == "Dry" ) %>%
  select(Salinity,pH,Phosphate_umolL:NN_umolL,Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like, MarineHumic_Like, TA) %>%
  #select(Salinity,pH,Phosphate_umolL:Lignin_Like ) %>%
  drop_na(Salinity,pH,Phosphate_umolL:NN_umolL,Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like, MarineHumic_Like, TA)

V_pca_Seep_wet<-Datalog %>%
  filter(Plate_Seep == "Seep" |Plate_Seep == "Spring")%>% # anchor by the spring because water was sampled at seep differently between 
  filter(Location == "Varari", 
         #Plate_Seep=="Seep",
         Season == "Wet") %>%
  select(Salinity,pH,Phosphate_umolL:NN_umolL,Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like, MarineHumic_Like, TA) %>%
  #select(Salinity,pH,Phosphate_umolL:Lignin_Like ) %>%
  drop_na(Salinity,pH,Phosphate_umolL:NN_umolL,Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like, MarineHumic_Like, TA)

# Run the PCA
pca_V_Seep <- prcomp(V_pca_Seep, scale. = TRUE, center = TRUE)
pca_V_Seep_wet <- prcomp(V_pca_Seep_wet, scale. = TRUE, center = TRUE)

# Extract the scores and loadings
PC_scores_Seep <-as_tibble(pca_V_Seep$x[,1:2])
PC_loadings_Seep<-as_tibble(pca_V_Seep$rotation) %>%
  bind_cols(labels = rownames(pca_V_Seep$rotation))

PC_scores_Seep_wet <-as_tibble(pca_V_Seep_wet$x[,1:2])
PC_loadings_Seep_wet<-as_tibble(pca_V_Seep_wet$rotation) %>%
  bind_cols(labels = rownames(pca_V_Seep_wet$rotation))

# calculate percent explained by each PC
perc.explainedSeep<-round(100*pca_V_Seep$sdev/sum(pca_V_Seep$sdev),1)
perc.explainedSeep_wet<-round(100*pca_V_Seep_wet$sdev/sum(pca_V_Seep_wet$sdev),1)

PC_loadings_Seep<-as_tibble(pca_V_Seep$rotation) %>%
  bind_cols(labels = rownames(pca_V_Seep$rotation))%>%
  mutate(groupings = case_when( # add groupings
    labels %in% c("Ammonia_umolL","NN_umolL","Phosphate_umolL","Silicate_umolL")~ "Nutrient Chemistry",
    labels == "Salinity" ~ "Salinity",
    labels %in% c("pH","TA") ~ "Carbonate Chemistry",
    labels %in% c("HIX","Lignin_Like","M_C","MarineHumic_Like","Tryptophan_Like","Tyrosine_Like","VisibleHumidic_Like")~"fDOM"
  ),
  nicenames = case_when(labels == "TempInSitu_seep" ~ "Temperature",
                        labels == "pH" ~ "pH<sub>T</sub>",
                        labels == "Lignin_Like" ~"Lignin Like",
                        labels == "M_C" ~ "M:C",
                        labels == "Tyrosine_Like" ~ "Tyrosine Like",
                        labels == "Tryptophan_Like" ~ "Tryptophan Like",
                        labels == "HIX"~"HIX",
                        labels == "MarineHumic_Like" ~ "Marine Humic Like",
                        labels == "VisibleHumidic_Like" ~ "Visible Humic Like",
                        labels == "Ammonia_umolL" ~ "Ammonium",
                        labels == "TA" ~ "Total Alkalinity",
                        labels == "Phosphate_umolL" ~ "Phosphate",
                        labels == "NN_umolL" ~ "Nitrate+Nitrite",
                        labels == "Silicate_umolL" ~ "Silicate",
                        labels == "Salinity" ~"Salinity"))

PC_loadings_Seep_wet<-as_tibble(pca_V_Seep_wet$rotation) %>%
  bind_cols(labels = rownames(pca_V_Seep_wet$rotation))%>%
  mutate(groupings = case_when( # add groupings
    labels %in% c("Ammonia_umolL","NN_umolL","Phosphate_umolL","Silicate_umolL")~ "Nutrient Chemistry",
    labels == "Salinity" ~ "Salinity",
    labels %in% c("pH","TA") ~ "Carbonate Chemistry",
    labels %in% c("HIX","Lignin_Like","M_C","MarineHumic_Like","Tryptophan_Like","Tyrosine_Like","VisibleHumidic_Like")~"fDOM"
  ),
  nicenames = case_when(labels == "TempInSitu_seep" ~ "Temperature",
                        labels == "pH" ~ "pH<sub>T</sub>",
                        labels == "Lignin_Like" ~"Lignin Like",
                        labels == "M_C" ~ "M:C",
                        labels == "Tyrosine_Like" ~ "Tyrosine Like",
                        labels == "Tryptophan_Like" ~ "Tryptophan Like",
                        labels == "HIX"~"HIX",
                        labels == "MarineHumic_Like" ~ "Marine Humic Like",
                        labels == "VisibleHumidic_Like" ~ "Visible Humic Like",
                        labels == "Ammonia_umolL" ~ "Ammonium",
                        labels == "TA" ~ "Total Alkalinity",
                        labels == "Phosphate_umolL" ~ "Phosphate",
                        labels == "NN_umolL" ~ "Nitrate+Nitrite",
                        labels == "Silicate_umolL" ~ "Silicate",
                        labels == "Salinity" ~"Salinity"))

# Put it with all the original data
V_pca_Data_all_Seep<-Data %>%
  filter(Plate_Seep == "Seep" |Plate_Seep == "Spring")%>% # anchor by the spring because water was sampled at seep differently between 
  filter(Location == "Varari", 
       #  Plate_Seep=="Seep", 
         Season == "Dry") %>%
  drop_na(Salinity,pH,Phosphate_umolL:NN_umolL,Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like) %>%
  bind_cols(PC_scores_Seep)

V_pca_Data_all_Seep_wet<-Data %>%
  filter(Plate_Seep == "Seep" |Plate_Seep == "Spring")%>% # anchor by the spring because water was sampled at seep differently between 
  filter(Location == "Varari", 
         #Plate_Seep=="Seep", 
         Season == "Wet") %>%
  drop_na(Salinity,pH,Phosphate_umolL:NN_umolL,Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like) %>%
  bind_cols(PC_scores_Seep_wet)

# scores plot
p1seep<-PC_loadings_Seep%>%
#  drop_na(Depth_seep)%>%
  ggplot(aes(x = PC1, y = PC2, label=nicenames, color = groupings))+
   #geom_point(data = V_pca_Data_all_Seep, inherit.aes = FALSE, aes(x = PC1, y = PC2, shape = Day_Night)) +
    coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  #scale_shape_manual(values = c(22,16))+
#  scale_color_gradient(low = "black", high = "yellow")+
#  scale_colour_viridis(limits = c(-1.1,-0.4))+
  geom_segment(data = PC_loadings_Seep, aes(x=0,y=0,xend=PC1*3,yend=PC2*3),
               arrow=arrow(length=unit(0.1,"cm")))+
   # annotate("text", x = PC_loadings_Seep$PC1*3+0.1, y = PC_loadings_Seep$PC2*3+.1,
  #          label = PC_loadings_Seep$labels)+
  geom_richtext(aes(x = PC1*3+0.1, y = PC2*3+0.1), show.legend = FALSE, size = 5, fill=NA, label.colour = NA)+
  scale_size(limits = c(-1.1,-0.4))+
  scale_color_manual(values = wes_palette("Darjeeling1"))+
  
#  guides(color=guide_legend(), size = guide_legend())+
  labs(#color = "PAR", size = "Water Depth", shape = "Day or Night",
       #title = "Varari Seep",
       color = "",
       x = paste0("PC1 ","(",perc.explainedSeep[1],"%)"),
       y = paste0("PC2 ","(",perc.explainedSeep[2],"%)"),
       title = "Varari (Dry)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_markdown(size = 16),
       # legend.key.size = unit(1, 'cm'),
        #legend.position = c(0.75, 0.85),
        legend.position = "none")

p1seep_wet<-PC_loadings_Seep_wet%>%
   ggplot(aes(x = PC1, y = PC2, label=nicenames, color = groupings))+
   coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) + 
 #   geom_point(data=V_pca_Data_all_Seep_wet, aes(x = -PC1, y = -PC2), inherit.aes = FALSE)+
   geom_segment(data = PC_loadings_Seep_wet, aes(x=0,y=0,xend=PC1*3,yend=PC2*3),
               arrow=arrow(length=unit(0.1,"cm")))+
   geom_richtext(aes(x = PC1*3+0.1, y = PC2*3+0.1), show.legend = FALSE, size = 5, fill=NA, label.colour = NA)+
  scale_size(limits = c(-1.1,-0.4))+
  scale_color_manual(values = wes_palette("Darjeeling1"))+
  labs(#color = "PAR", size = "Water Depth", shape = "Day or Night",
    #title = "Varari Seep",
    color = "",
    x = paste0("PC1 ","(",perc.explainedSeep_wet[1],"%)"),
    y = paste0("PC2 ","(",perc.explainedSeep_wet[2],"%)"),
    title = "Varari (Wet)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        legend.position = c(0.25, 0.25))

## Cabral Seep
C_pca_Seep<-Datalog %>%
  filter(Plate_Seep == "Seep" |Plate_Seep == "Spring")%>% # anchor by the spring because water was sampled at seep differently between 
  filter(Location == "Cabral", 
       #  Plate_Seep=="Seep",
         Season == "Dry") %>%
  select(Salinity,pH,Phosphate_umolL:NN_umolL,Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like, MarineHumic_Like, TA) %>%
  drop_na(Salinity,pH,Phosphate_umolL:NN_umolL,Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like,MarineHumic_Like, TA)

C_pca_Seep_wet<-Datalog %>%
  filter(Plate_Seep == "Seep" |Plate_Seep == "Spring")%>% # anchor by the spring because water was sampled at seep differently between 
  filter(Location == "Cabral", 
    #     Plate_Seep=="Seep", 
         Season == "Wet") %>%
  select(Salinity,pH,Phosphate_umolL:NN_umolL,Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like, MarineHumic_Like, TA) %>%
  drop_na(Salinity,pH,Phosphate_umolL:NN_umolL,Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like, MarineHumic_Like, TA)

# Run the PCA
pca_C_Seep <- prcomp(C_pca_Seep, scale. = TRUE, center = TRUE)
pca_C_Seep_wet <- prcomp(C_pca_Seep_wet, scale. = TRUE, center = TRUE)

# Extract the scores and loadings
PC_scores_SeepC <-as_tibble(pca_C_Seep$x[,1:2])
PC_loadings_SeepC<-as_tibble(pca_C_Seep$rotation) %>%
  bind_cols(labels = rownames(pca_C_Seep$rotation))

PC_scores_SeepC_wet <-as_tibble(pca_C_Seep_wet$x[,1:2])
PC_loadings_SeepC_wet<-as_tibble(pca_C_Seep_wet$rotation) %>%
  bind_cols(labels = rownames(pca_C_Seep_wet$rotation))

# calculate percent explained by each PC
perc.explainedSeepC<-round(100*pca_C_Seep$sdev/sum(pca_C_Seep$sdev),1)
perc.explainedSeepC_wet<-round(100*pca_C_Seep_wet$sdev/sum(pca_C_Seep_wet$sdev),1)


PC_loadings_SeepC<-as_tibble(pca_C_Seep$rotation) %>%
  bind_cols(labels = rownames(pca_C_Seep$rotation))%>%
  mutate(groupings = case_when( # add groupings
    labels %in% c("Ammonia_umolL","NN_umolL","Phosphate_umolL","Silicate_umolL")~ "Nutrient Chemistry",
    labels == "Salinity" ~ "Salinity",
    labels %in% c("pH","TA") ~ "Carbonate Chemistry",
    labels %in% c("HIX","Lignin_Like","M_C","MarineHumic_Like","Tryptophan_Like","Tyrosine_Like","VisibleHumidic_Like")~"fDOM"
  ),
  nicenames = case_when(labels == "TempInSitu_seep" ~ "Temperature",
                        labels == "pH" ~ "pH<sub>T</sub>",
                        labels == "Lignin_Like" ~"Lignin Like",
                        labels == "M_C" ~ "M:C",
                        labels == "Tyrosine_Like" ~ "Tyrosine Like",
                        labels == "Tryptophan_Like" ~ "Tryptophan Like",
                        labels == "HIX"~"HIX",
                        labels == "MarineHumic_Like" ~ "Marine Humic Like",
                        labels == "VisibleHumidic_Like" ~ "Visible Humic Like",
                        labels == "Ammonia_umolL" ~ "Ammonium",
                        labels == "TA" ~ "Total Alkalinity",
                        labels == "Phosphate_umolL" ~ "Phosphate",
                        labels == "NN_umolL" ~ "Nitrate+Nitrite",
                        labels == "Silicate_umolL" ~ "Silicate",
                        labels == "Salinity" ~"Salinity"))

PC_loadings_SeepC_wet<-as_tibble(pca_C_Seep_wet$rotation) %>%
  bind_cols(labels = rownames(pca_C_Seep_wet$rotation))%>%
  mutate(groupings = case_when( # add groupings
    labels %in% c("Ammonia_umolL","NN_umolL","Phosphate_umolL","Silicate_umolL")~ "Nutrient Chemistry",
    labels == "Salinity" ~ "Salinity",
    labels %in% c("pH","TA") ~ "Carbonate Chemistry",
    labels %in% c("HIX","Lignin_Like","M_C","MarineHumic_Like","Tryptophan_Like","Tyrosine_Like","VisibleHumidic_Like")~"fDOM"
  ),
  nicenames = case_when(labels == "TempInSitu_seep" ~ "Temperature",
                        labels == "pH" ~ "pH<sub>T</sub>",
                        labels == "Lignin_Like" ~"Lignin Like",
                        labels == "M_C" ~ "M:C",
                        labels == "Tyrosine_Like" ~ "Tyrosine Like",
                        labels == "Tryptophan_Like" ~ "Tryptophan Like",
                        labels == "HIX"~"HIX",
                        labels == "MarineHumic_Like" ~ "Marine Humic Like",
                        labels == "VisibleHumidic_Like" ~ "Visible Humic Like",
                        labels == "Ammonia_umolL" ~ "Ammonium",
                        labels == "TA" ~ "Total Alkalinity",
                        labels == "Phosphate_umolL" ~ "Phosphate",
                        labels == "NN_umolL" ~ "Nitrate+Nitrite",
                        labels == "Silicate_umolL" ~ "Silicate",
                        labels == "Salinity" ~"Salinity"))

# Put it with all the original data
C_pca_Data_all_Seep<-Data %>%
  filter(Plate_Seep == "Seep" |Plate_Seep == "Spring")%>% # anchor by the spring because water was sampled at seep differently between 
  filter(Location == "Cabral", 
  #       Plate_Seep=="Seep",
         Season == "Dry") %>%
  drop_na(Salinity,pH,Phosphate_umolL:Lignin_Like,TA) %>%
  bind_cols(PC_scores_SeepC)

C_pca_Data_all_Seep_wet<-Data %>%
  filter(Plate_Seep == "Seep" |Plate_Seep == "Spring")%>% # anchor by the spring because water was sampled at seep differently between 
  filter(Location == "Cabral", 
        # Plate_Seep=="Seep", 
         Season == "Wet") %>%
  drop_na(Salinity,pH,Phosphate_umolL:Lignin_Like,TA) %>%
  bind_cols(PC_scores_SeepC_wet)

# scores plot
p2seep<-PC_loadings_SeepC %>%
  ggplot(aes(x = -PC1, y = -PC2, label=nicenames, color = groupings))+
    coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  geom_segment(data = PC_loadings_SeepC, aes(x=0,y=0,xend=-PC1*3,yend=-PC2*3),
               arrow=arrow(length=unit(0.1,"cm")))+
  geom_richtext(aes(x = -PC1*3+0.1, y = -PC2*3+0.1), show.legend = FALSE, size = 5, fill=NA, label.colour = NA)+
  scale_size(limits = c(-.55,-0.05)) +
  guides(color=guide_legend(), size = guide_legend())+
  scale_color_manual(values = wes_palette("Darjeeling1"))+
  labs(#color = "PAR", size = "Water Depth"
    color = "",
    x = paste0("PC1 ","(",perc.explainedSeepC[1],"%)"),
    y = paste0("PC2 ","(",perc.explainedSeepC[2],"%)"),
       title = "Cabral Seep (Dry)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        axis.title = element_text(size = 18, hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 16),
        legend.position = "none")

p2seep_wet<-PC_loadings_SeepC_wet %>%
  ggplot(aes(x = PC1, y = PC2, label=nicenames, color = groupings))+
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  geom_segment(data = PC_loadings_SeepC_wet, aes(x=0,y=0,xend=PC1*3,yend=PC2*3),
               arrow=arrow(length=unit(0.1,"cm")))+
  geom_richtext(aes(x = PC1*3+0.1, y = PC2*3+0.1), show.legend = FALSE, size = 5, fill=NA, label.colour = NA)+
  scale_size(limits = c(-.55,-0.05)) +
  guides(color=guide_legend(), size = guide_legend())+
  scale_color_manual(values = wes_palette("Darjeeling1"))+
  labs(#color = "PAR", size = "Water Depth"
    color = "",
    x = paste0("PC1 ","(",perc.explainedSeepC_wet[1],"%)"),
    y = paste0("PC2 ","(",perc.explainedSeepC_wet[2],"%)"),
    title = "Cabral Seep (Wet)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        axis.title = element_text(size = 18, hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 16),
        legend.position = "none")

# save the plots
SeepPCA<-(p1seep+p1seep_wet)/(p2seep+p2seep_wet)+ 
  patchwork::plot_annotation(#"Seep", 
                             theme = theme(plot.title = element_text(size = rel(1.5), face = "bold", hjust = 0.5, 
                                                                     margin = margin(t = 10, b = 20, unit = "pt"))))

ggsave(plot = SeepPCA, filename = here("Output","SeepPCA.png"), width = 16, height = 16)



### plots of each parameter versus depth
Data %>%
  filter(Plate_Seep=="Seep", Salinity>10) %>% # I think I need to remove the bottle sample...
  pivot_longer(cols = Salinity:Ammonia_umolL) %>% 
  ggplot(aes(x = Depth_seep, y = value, color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_y_continuous(trans = "log10")+
  facet_wrap(~Location*name, scales = "free")

## Calculate the ranges for each parameter
ranges<-Data %>%
  filter(Plate_Seep=="Seep"| Plate_Seep == "Spring") %>% # just do the seep data
  select(Location,Season, Salinity, TA,pH,Phosphate_umolL:Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like, HIX, TempInSitu_seep, MarineHumic_Like, Lignin_Like, M_C)%>%
 # select(Location, Salinity, TA: Lignin_Like, TempInSitu_seep) %>%
  pivot_longer(cols = c(Salinity, TA: M_C)) %>%
  drop_na()%>%
  group_by(Location,name, Season)%>%
  summarise(min = round(min(value),2),
            max = round(max(value),2)) %>%
    mutate(unit = case_when(name == "TempInSitu_seep" ~ " &deg;C", # add units
                            name %in% c("Ammonia_umolL","NN_umolL","Silicate_umolL","Phosphate_umolL")~" &mu;mol L<sup>-1</sup>",
                            name %in% c("M_C","HIX")~ " Ratio",
                            name %in% c("Salinity")~ " psu",
                            name %in% c("pH")~ " pH total scale",
                            name == "TA"~ " &mu;mol kg<sup>-1</sup>",
                            name %in% c("Lignin_Like","Tyrosine_Like","Tryptophan_Like","MarineHumic_Like","VisibleHumidic_Like")~" Raman Units"),
      range = paste0("[",min," - ",max,unit,"]"),
      nicenames = case_when(name == "TempInSitu_seep" ~ "Temperature",
                            name == "pH" ~ "pH<sub>T</sub>",
                            name == "Lignin_Like" ~"Lignin Like",
                         #   name == "M_C" ~ "M:C",
                            name == "Tyrosine_Like" ~ "Tyrosine Like",
                            name == "Tryptophan_Like" ~ "Tryptophan Like",
                            name == "HIX"~"HIX",
                            name == "MarineHumic_Like" ~ "Marine Humic Like",
                            name == "VisibleHumidic_Like" ~ "Visible Humic Like",
                            name == "Ammonia_umolL" ~ "Ammonium",
                            name == "TA" ~ "Total Alkalinity",
                            name == "Phosphate_umolL" ~ "Phosphate",
                            name == "NN_umolL" ~ "Nitrate+Nitrite",
                            name == "Silicate_umolL" ~ "Silicate",
                            name == "Salinity" ~"Salinity"
                            
        
      ))
  


# Make a correlation plot of everything versus depth
# function for correlation coef
cor_fun <- function(data) cor.test(data$value, data$Depth_seep, method = "pearson")%>% tidy()

cortest<-Datalog %>%
  filter(Plate_Seep=="Seep") %>% # just do the seep data
  select(Location, Season, Salinity, TA: Lignin_Like, TempInSitu_seep, Depth_seep)%>%
  pivot_longer(cols = c(Salinity, TA: Lignin_Like, TempInSitu_seep)) %>%
  drop_na()%>%
  group_by(Location, name, Season)%>%
  nest() %>%
  mutate(model = map(data, cor_fun)) %>%
  select(Location, name,Season, model) %>%
  unnest(model)%>% # calculate correlation coefficient
  mutate(sig = ifelse(p.value<0.1, 1,0 ))# add a 1 if significant correlation (to 0.1 pvalue)

# make a plot of it
Vcortest<-cortest %>%
  filter(Location == "Varari") %>%
  ggplot(aes(x = fct_reorder(name, abs(estimate)), y = Season, color = estimate, size = abs(estimate)))+
  geom_point()+
  geom_point(aes(shape = factor(sig)), size = 2, color = "white")+
  scale_color_gradient2(low = "#005AB5", high = "#DC3220",mid = "black", midpoint = 0, limits = c(-1,1)  )+
  scale_size(range = c(0.1,10))+
  scale_shape_manual(values = c(NA,8))+
  coord_flip()+
  theme_bw()+
  guides(size = "none",
         shape = "none")+
  labs(x = "",
       y = "",
       color = "Correlation coefficient with water depth",
       title = "Varari")

Ccortest<-cortest %>%
  filter(Location == "Cabral") %>%
  ggplot(aes(x = fct_reorder(name, abs(estimate)), y = Season, color = estimate, size = abs(estimate)))+
  geom_point()+
  geom_point(aes(shape = factor(sig)), size = 2, color = "white")+
  scale_color_gradient2(low = "#005AB5", high = "#DC3220",mid = "black", midpoint = 0, limits = c(-1,1) )+
  scale_size(range = c(0.1,10))+
  scale_shape_manual(values = c(NA,8))+
  coord_flip()+
  theme_bw()+
  guides(size = "none",
         shape = "none")+
  labs(x = "",
       y = "",
       color = "Correlation coefficient with water depth",
       title = "Cabral")

Vcortest+Ccortest+plot_layout(guides = "collect")+ plot_annotation(tag_levels = "A")&theme(legend.position = 'bottom')

ggsave(here("Output","CorrelationPlot_seep.pdf"), height = 8, width = 10)

## same thing but everything versus salinity

cor_fun <- function(data) cor.test(data$value, data$Salinity, method = "pearson")%>% tidy()

cortest<-Datalog %>%
  filter(Plate_Seep=="Seep" | Plate_Seep == "Spring" ) %>% # just do the seep data
  select(Location,Season, Salinity, TA: Lignin_Like, TempInSitu_seep)%>%
  pivot_longer(cols = c(TA: Lignin_Like, TempInSitu_seep)) %>%
  drop_na()%>%
  group_by(Location, name, Season)%>%
  nest() %>%
  mutate(model = map(data, cor_fun)) %>%
  select(Location, name,Season, model) %>%
  unnest(model)%>% # calculate correlation coefficient
  mutate(sig = ifelse(p.value<0.05, 1,0 ))# add a 1 if significant correlation (to 0.1 pvalue)

Vcortest_salinity<-cortest %>%
  filter(!name %in% c("M_C","HIX"),
         Location == "Varari")%>%
  ggplot(aes(x = fct_reorder(name, abs(estimate)), y = Season, color = estimate, size = abs(estimate)))+
  geom_point()+
  geom_point(aes(shape = factor(sig)), size = 2, color = "white")+
  scale_color_gradient2(low = "#005AB5", high = "#DC3220",mid = "black", midpoint = 0,limits = c(-1,1) )+
  scale_size(range = c(0.1,10))+
  scale_shape_manual(values = c(NA,8))+
  coord_flip()+
  theme_bw()+
  guides(size = "none",
         shape = "none")+
  labs(x = "",
       y = "",
       color = "Correlation coefficient with Salinity",
       title = "Varari")

Ccortest_salinity<-cortest %>%
  filter(!name %in% c("M_C","HIX"),
         Location == "Cabral")%>%
  ggplot(aes(x = fct_reorder(name, abs(estimate)), y = Season, color = estimate, size = abs(estimate)))+
  geom_point()+
  geom_point(aes(shape = factor(sig)), size = 2, color = "white")+
  scale_color_gradient2(low = "#005AB5", high = "#DC3220",mid = "black", midpoint = 0,limits = c(-1,1) )+
  scale_size(range = c(0.1,10))+
  scale_shape_manual(values = c(NA,8))+
  coord_flip()+
  theme_bw()+
  guides(size = "none",
         shape = "none")+
  labs(x = "",
       y = "",
       color = "Correlation coefficient with Salinity",
       title = "Cabral")

Vcortest_salinity+Ccortest_salinity+plot_layout(guides = "collect")+ plot_annotation(tag_levels = "A")&theme(legend.position = 'bottom')

ggsave(here("Output","CorrelationPlot_seepSalinity.pdf"), height = 8, width = 10)


# just both with range

cortest %>%
  left_join(ranges) %>% # join with the ranges
 # filter(Location == "Varari") %>%
  filter(!name %in% c("HIX","M_C") )%>%
  ggplot(aes(x = fct_reorder(nicenames, estimate),#y = 1,
             y = c(rep(1,24),rep(1.8,24)),
             color = estimate, size = abs(estimate)))+
  geom_point()+
  geom_point(aes(shape = factor(sig)), size = 2, color = "white")+
  geom_richtext(aes(y = c(rep(1.4,24),rep(2.2,24)),
                    label = range), size = 5, color = "black", fill = NA, label.colour = NA)+
#  ylim (0.8,3)+
  scale_color_gradient2(low = "#005AB5", high = "#DC3220",mid = "black", midpoint = 0, limits = c(-1,1))+
  scale_size(range = c(0.1,10))+
  scale_shape_manual(values = c(NA,8))+
  coord_flip()+
  theme_bw()+
  guides(size = "none",
         shape = "none",
         colour = guide_colourbar(title.position="top", title.hjust = 0.5))+
  labs(x = "",
       y = "",
       color = "Correlation with Salinity"
     #  title = "Correlation with Salinity"
       )+
  scale_y_continuous(breaks = c(1,1.8), labels = c("Dry","Wet"), limits = c(0.8,2.5))+
  theme(#axis.ticks.x = element_blank(),
      #  axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_markdown(size = 14),
        axis.text.x = element_markdown(size = 14),
        legend.title = element_text(size=14),
        legend.text = element_text(size = 12),
      #  panel.grid.major.x = element_blank(),
      #  panel.grid.minor.x = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm")
        )+
  facet_wrap(~Location, ncol = 2)

ggsave(here("Output","CorrelationPlot_seepSalinityV.png"), height = 8, width = 16)

### Pure Spring Water
mean_spring<-Data %>%
  filter(CowTagID %in% c("VSPRING","Varari_Well","CSPRING_BEACH2","CSPRING_ROAD")) %>% # just do the seep data
  select(Location, Salinity,CowTagID, TA,pH,Phosphate_umolL:Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like, HIX, TempInSitu_seep)%>%
  # select(Location, Salinity, TA: Lignin_Like, TempInSitu_seep) %>%
  pivot_longer(cols = c(Salinity, TA: TempInSitu_seep)) %>%
  drop_na()%>%
  group_by(CowTagID, name)%>%
  summarise(mean = round(min(value),2)) %>%
  mutate(unit = case_when(name == "TempInSitu_seep" ~ " &deg;C", # add units
                          name %in% c("Ammonia_umolL","NN_umolL","Silicate_umolL","Phosphate_umolL")~" &mu;mol L<sup>-1</sup>",
                          name %in% c("M_C","HIX")~ " Ratio",
                          name %in% c("Salinity")~ " psu",
                          name %in% c("pH")~ " pH total scale",
                          name == "TA"~ " &mu;mol kg<sup>-1</sup>",
                          name %in% c("Lignin_Like","Tyrosine_Like","Tryptophan_Like","MarineHumic_Like","VisibleHumidic_Like")~" Raman Units"),
         nicenames = case_when(name == "TempInSitu_seep" ~ "Temperature",
                               name == "pH" ~ "pH<sub>T</sub>",
                               #    name == "Lignin_Like" ~"Lignin Like",
                               #   name == "M_C" ~ "M:C",
                               name == "Tyrosine_Like" ~ "Tyrosine Like",
                               name == "Tryptophan_Like" ~ "Tryptophan Like",
                               name == "HIX"~"HIX",
                               #  name == "MarineHumic_Like" ~ "Marine Humic Like",
                               name == "VisibleHumidic_Like" ~ "Visible Humic Like",
                               name == "Ammonia_umolL" ~ "Ammonium",
                               name == "TA" ~ "Total Alkalinity",
                               name == "Phosphate_umolL" ~ "Phosphate",
                               name == "NN_umolL" ~ "Nitrate+Nitrite",
                               name == "Silicate_umolL" ~ "Silicate",
                               name == "Salinity" ~"Salinity"
                               
                               
         )) %>%
  pivot_wider(values_from = mean, names_from  = CowTagID)



mx<-mean_spring %>%
  column_to_rownames(var = "nicenames") %>%
  select(unit,CSPRING_BEACH2,CSPRING_ROAD,VSPRING,Varari_Well ) %>%
  rename("Cabral Beach"=CSPRING_BEACH2, "Cabral Road"=CSPRING_ROAD,
         "Varari Beach"=VSPRING, "Varari Well"=Varari_Well) %>%
    htmlTable()

mx
   
### Some plots of Silicate vs pH and NN

Data %>%
  filter(Plate_Seep == "Plate")%>%
  anti_join(remove2) %>%
  anti_join(remove3) %>%
  anti_join(remove4) %>%
  anti_join(remove5) %>%
  ggplot(aes(x = log(Silicate_umolL), y = pH, color = Day_Night, shape = Tide))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free") 
ggsave(here("Output","SivspH.png"), width = 5, height = 5)

Data %>%
  filter(Plate_Seep == "Plate")%>%
  anti_join(remove2) %>%
  anti_join(remove3) %>%
  anti_join(remove4) %>%
  anti_join(remove5) %>%
  ggplot(aes(x = log(Silicate_umolL), y = log(NN_umolL), color = Day_Night, shape = Tide))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")
ggsave(here("Output","SivsNN.png"), width = 5, height = 5)

### Only keep certain dataframes for eco metab script
rm(list= ls()[!(ls() %in% c("Data", "Datalog", "remove_varari","remove_cabral","remove_vararilog","remove_cabrallog","turb_all"))])
