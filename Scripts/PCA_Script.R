# A script to create PCAs for 24 hour biogeochem data for Cabral and Varari
# Edited on 2/15/2022
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

# load the 24 hour chemistry data #####################
#Data<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv")
Data<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC2.csv")

## Load the turb nutrient data
turbdata<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Nutrients/Turb_NC.csv")

# Load the Seep Data
VarariSeep<-read_csv(here("Data","Varari","AllVarariSeepData.csv")) %>%
  rename_with(.cols = TempInSitu:PAR_calc,function(x){paste0(x,"_seep")}) %>% # rename columns to say seep at the end
  rename(DateTime = date, Location = Site)

CabralSeep<-read_csv(here("Data","Cabral","AllCabralSeepData.csv")) %>%
  rename_with(.cols = TempInSitu:PAR_calc,function(x){paste0(x,"_seep")}) %>% # rename columns to say seep at the end
  rename(DateTime = date, Location = Site)

# bind them together
SeepAll<-bind_rows(VarariSeep, CabralSeep)

# Bind with the chem data
Data<-Data %>%
  left_join(SeepAll) 

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
remove<-Data %>% filter(CowTagID=="V2", Tide == "Low",Day_Night=="Day", Date =="8/8/2021")

remove3<-Data %>% filter(CowTagID== "C4", Tide =="Low", Day_Night=="Night", Date == "8/9/2021")

## also filter out all the data from the first low tide where the water level was super high

remove2<-Data %>% filter(Tide =="Low", Day_Night=="Day", Date == "8/6/2021")

# extract the params for the PCA
V_pca_Data<-Data %>%
  anti_join(remove)%>%
  anti_join(remove2)%>%
 # filter(Location == "Varari", Tide %in% c("High","Low")) %>%
  filter(Location == "Varari", Plate_Seep=="Plate") %>%
  #select(Salinity,pH,Phosphate_umolL, Silicate_umolL, NN_umolL, Ammonia_umolL ) %>%
  select(Salinity,pH,Phosphate_umolL:Lignin_Like)%>%
  drop_na()

# Run the PCA
pca_V <- prcomp(V_pca_Data, scale. = TRUE, center = TRUE)

# Extract the scores and loadings
PC_scores <-as_tibble(pca_V$x[,1:2])
PC_loadings<-as_tibble(pca_V$rotation) %>%
  bind_cols(labels = rownames(pca_V$rotation))%>%
  mutate(groupings = case_when( # add groupings
    labels %in% c("Ammonia_umolL","NN_umolL","Phosphate_umolL","Silicate_umolL")~ "Nutrient Chemistry",
    labels == "Salinity" ~ "Salinity",
    labels == "pH" ~ "Carbonate Chemistry",
    labels %in% c("HIX","Lignin_Like","M_C","MarineHumic_Like","Tryptophan_Like","Tyrosine_Like","VisibleHumidic_Like")~"fDOM"
  ))

# Put it with all the original data
V_pca_Data_all<-Data %>%
  anti_join(remove)%>%
  anti_join(remove2)%>%
  filter(Location == "Varari", Plate_Seep=="Plate") %>%
#  filter(Location == "Varari", Tide %in% c("High","Low")) %>%
# drop_na(Salinity,pH,Phosphate_umolL, Silicate_umolL, NN_umolL, Ammonia_umolL) %>%
  drop_na(Salinity,pH,Phosphate_umolL:Lignin_Like )%>%
  #drop_na(Salinity,pH,Phosphate_umolL:Ammonia_umolL)%>%
  bind_cols(PC_scores)

# scores plot
p1<-V_pca_Data_all %>%
  ggplot(aes(x = PC1, y = PC2, color = Tide, shape = Day_Night))+
  geom_point(size = 2) +
 # geom_point(data = V_pca_Data_all %>% filter(Plate_Seep=="Seep"), aes(x = PC1, y = PC2,shape = Day_Night ), color = "black")+
  coord_cartesian(xlim = c(-6, 4), ylim = c(-6, 4)) +
  scale_shape_manual(values = c(22,16))+
  scale_colour_hue(l = 45)+
  scale_fill_hue(l = 45)+
  ggforce::geom_mark_ellipse(
    aes(fill = Tide, label = paste(Day_Night, Tide), color = Tide), 
    alpha = .15, show.legend = FALSE,  label.buffer = unit(1, "mm"))+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# scores plot with depth and light as continuous instead of discrete... missing depth data from lowtide at night :(
p1_DL<-V_pca_Data_all %>%
  ggplot(aes(x = PC1, y = PC2, color = PAR_calc_seep+1, size = -Depth_seep))+
  geom_point() +
  coord_cartesian(xlim = c(-4, 7), ylim = c(-6, 4)) +
  scale_color_gradient(low = "black", high = "yellow", trans = "log")+
  theme_bw()+
  theme(#legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        

# loadings plot 
p2<-PC_loadings %>%
  ggplot(aes(x=PC1, y=PC1, label=labels, color = groupings))+
    geom_segment(data = PC_loadings, aes(x=0,y=0,xend=PC1*10,yend=PC2*10),
      arrow=arrow(length=unit(0.1,"cm")))+
  geom_text(aes(x = PC1*10+0.1, y = PC2*10+.1 ), show.legend = FALSE) +
  # annotate("text", x = PC_loadings$PC1*10+0.1, y = PC_loadings$PC2*10+.1,
  #          label = PC_loadings$labels)+
  coord_cartesian(xlim = c(-6, 4), ylim = c(-6, 4)) +
  labs(color ="")+
  scale_color_manual(values = wes_palette("Darjeeling1"))+
   theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = c(0.75, 0.85))

VarariPCA<-p1+p2+ 
  patchwork::plot_annotation("Varari Plates", 
                             theme = theme(plot.title = element_text(size = rel(1.5), face = "bold", hjust = 0.5, 
                                                                     margin = margin(t = 10, b = 20, unit = "pt"))))

ggsave(plot = VarariPCA, filename = here("Output","VarariPCA.png"), width = 14, height = 8)

### Site level pca with variances
V_pca_Data_site<-Data %>%
  anti_join(remove)%>%
  anti_join(remove2)%>%
  filter(Location == "Varari", Plate_Seep=="Plate") %>%
  group_by(CowTagID)%>% # calculate the range by cowtag
  summarise_at(vars(Salinity,pH,Phosphate_umolL:Lignin_Like), .funs = function(x) {max(x, na.rm = TRUE)-min(x, na.rm = TRUE)}) %>%
  ungroup()%>%
  left_join(turbdata)%>%
  select(CowTagID, Salinity,pH,Phosphate_umolL, Silicate_umolL, NN_umolL, Ammonia_umolL, del15N, N_percent ) %>%
 # select(Salinity,pH,Phosphate_umolL:Lignin_Like )%>%
  drop_na()

# Run the PCA
pca_V_site <- prcomp(V_pca_Data_site[,-1], scale. = TRUE, center = TRUE)
# Extract the scores and loadings
PC_scores <-as_tibble(pca_V_site$x[,1:2])
PC_loadings<-as_tibble(pca_V_site$rotation) %>%
  bind_cols(labels = rownames(pca_V_site$rotation))

PC_scores %>%
  bind_cols(V_pca_Data_site)%>%
  ggplot(aes(x = PC1, y = PC2))+
 # geom_point(color = "red") +
  geom_text(aes(x = PC1, y = PC2,label = CowTagID),color = "red")+
  geom_segment(data = PC_loadings, aes(x=0,y=0,xend=PC1*10,yend=PC2*10),
               arrow=arrow(length=unit(0.1,"cm")), color = "grey")+
  annotate("text", x = PC_loadings$PC1*10+0.1, y = PC_loadings$PC2*10+.1,
           label = PC_loadings$labels)+
  
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


### some plots of the ranges
ggplot(V_pca_Data_site, aes(x = Silicate_umolL, y = NN_umolL, label = CowTagID, color = del15N))+
  geom_point()+
  geom_label()+
  geom_smooth(method = "lm")+
  labs(x = "Silicate range",
       y = "NN range")+
  theme_bw()


ggplot(V_pca_Data_site, aes(y = N_percent, x = NN_umolL))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(y = "%Tissue N from turbinaria",
       x = "NN diel range")+
  theme_bw()


ggplot(V_pca_Data_site, aes(y = Salinity, x = Silicate_umolL))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(y = "Salinity",
       x = "Silicate diel range")+
  theme_bw()


ggplot(V_pca_Data_site, aes(x = Ammonia_umolL, y = pH))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(x = "Ammonia diel",
       y = "pH diel range")+
  theme_bw()



#### Cabral #####

#Extract the cabral data
C_pca_Data<-Data %>%
  anti_join(remove3)%>%
  filter(Location == "Cabral", Plate_Seep=="Plate") %>%
  select(Salinity,pH,Phosphate_umolL:Lignin_Like )%>%
  #select(Salinity,pH,Phosphate_umolL, Silicate_umolL, NN_umolL, Ammonia_umolL ) %>%
  drop_na()

# Run the PCA
pca_C <- prcomp(C_pca_Data, scale. = TRUE, center = TRUE)

# Extract the scores and loadings
PC_scoresC <-as_tibble(pca_C$x[,1:2])
PC_loadingsC<-as_tibble(pca_C$rotation) %>%
  bind_cols(labels = rownames(pca_C$rotation))%>%
  mutate(groupings = case_when( # add groupings
    labels %in% c("Ammonia_umolL","NN_umolL","Phosphate_umolL","Silicate_umolL")~ "Nutrient Chemistry",
    labels == "Salinity" ~ "Salinity",
    labels == "pH" ~ "Carbonate Chemistry",
    labels %in% c("HIX","Lignin_Like","M_C","MarineHumic_Like","Tryptophan_Like","Tyrosine_Like","VisibleHumidic_Like")~"fDOM"
  ))


# Put it with all the original data
C_pca_Data_all<-Data %>%
  anti_join(remove3)%>%
  select(!Jamie_Plate_ID)%>% # Jamie's plates are all NA here
  filter(Location == "Cabral", Plate_Seep=="Plate") %>%
  drop_na(Salinity,pH,Phosphate_umolL:Lignin_Like) %>%
  bind_cols(PC_scoresC)  

# scores plot
p1c<-C_pca_Data_all %>%
  ggplot(aes(x = PC1, y = PC2, color = Tide, shape = Day_Night))+
  geom_point() +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  scale_shape_manual(values = c(22,16))+
  scale_colour_hue(l = 45)+
  scale_fill_hue(l = 45)+
  ggforce::geom_mark_ellipse(
    aes(fill = Tide, label = paste(Day_Night, Tide), color = Tide), 
    alpha = .15, show.legend = FALSE,  label.buffer = unit(1, "mm"))+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# scores plot with depth and light as continuous instead of discrete... missing depth data from lowtide at night :(
p1c_DL<-C_pca_Data_all %>%
  ggplot(aes(x = PC1, y = PC2, color = PAR_calc_seep+1, size = -Depth_seep))+
  geom_point() +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  scale_color_gradient(low = "black", high = "yellow", trans = "log")+
  theme_bw()+
  theme(#legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# loadings plot 
p2c<-PC_loadingsC %>%
  ggplot(aes(x=PC1, y=PC1, label=labels, color = groupings))+
  geom_segment(aes(x=0,y=0,xend=PC1*10,yend=PC2*10),
               arrow=arrow(length=unit(0.1,"cm")))+
  geom_text(aes(x = PC1*10+0.1, y = PC2*10+.1 ), show.legend = FALSE) +
  scale_color_manual(values = wes_palette("Darjeeling1"))+
  #  annotate("text", x = PC_loadingsC$PC1*10+0.1, y = PC_loadingsC$PC2*10+.1,
  #        label = PC_loadingsC$labels)+
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  labs(color = "")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = c(0.75, 0.20))

CabralPCA<-p1c+p2c+ 
  patchwork::plot_annotation("Cabral Plates", 
                             theme = theme(plot.title = element_text(size = rel(1.5), face = "bold", hjust = 0.5, 
                                                                     margin = margin(t = 10, b = 20, unit = "pt"))))

ggsave(plot = CabralPCA, filename = here("Output","CabralPCA.pdf"), width = 14, height = 8)

### Make a PCA of the seeps colored by water depth ####
#Varari
V_pca_Seep<-Data %>%
  filter(Location == "Varari", Plate_Seep=="Seep") %>%
  select(Salinity,pH,Phosphate_umolL:Lignin_Like ) %>%
  drop_na()

# Run the PCA
pca_V_Seep <- prcomp(V_pca_Seep, scale. = TRUE, center = TRUE)

# Extract the scores and loadings
PC_scores_Seep <-as_tibble(pca_V_Seep$x[,1:2])
PC_loadings_Seep<-as_tibble(pca_V_Seep$rotation) %>%
  bind_cols(labels = rownames(pca_V_Seep$rotation))

# Put it with all the original data
V_pca_Data_all_Seep<-Data %>%
  filter(Location == "Varari", Plate_Seep=="Seep") %>%
  drop_na(Salinity,pH,Phosphate_umolL:Lignin_Like ) %>%
  bind_cols(PC_scores_Seep)

# scores plot
p1seep<-V_pca_Data_all_Seep %>%
  drop_na(Depth_seep)%>%
  ggplot(aes(x = PC1, y = PC2, color = log(PAR_calc_seep+1)))+
   geom_point(aes(size = -Depth_seep)) +
  #  coord_cartesian(xlim = c(-4, 7), ylim = c(-4, 4)) +
  #scale_shape_manual(values = c(22,16))+
  scale_color_gradient(low = "black", high = "yellow")+
#  scale_colour_viridis(limits = c(-1.1,-0.4))+
   geom_segment(data = PC_loadings_Seep, aes(x=0,y=0,xend=PC1*3,yend=PC2*3),
               arrow=arrow(length=unit(0.1,"cm")), color = "grey")+
  annotate("text", x = PC_loadings_Seep$PC1*3+0.1, y = PC_loadings_Seep$PC2*3+.1,
           label = PC_loadings_Seep$labels)+
  scale_size(limits = c(-1.1,-0.4))+
  guides(color=guide_legend(), size = guide_legend())+
  labs(color = "PAR", size = "Water Depth", shape = "Day or Night",
       title = "Varari Seep")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## Cabral Seep
C_pca_Seep<-Data %>%
  filter(Location == "Cabral", Plate_Seep=="Seep") %>%
  select(Salinity,pH,Phosphate_umolL:Lignin_Like ) %>%
  drop_na()

# Run the PCA
pca_C_Seep <- prcomp(C_pca_Seep, scale. = TRUE, center = TRUE)

# Extract the scores and loadings
PC_scores_SeepC <-as_tibble(pca_C_Seep$x[,1:2])
PC_loadings_SeepC<-as_tibble(pca_C_Seep$rotation) %>%
  bind_cols(labels = rownames(pca_C_Seep$rotation))

# Put it with all the original data
C_pca_Data_all_Seep<-Data %>%
  filter(Location == "Cabral", Plate_Seep=="Seep") %>%
  drop_na(Salinity,pH,Phosphate_umolL:Lignin_Like) %>%
  bind_cols(PC_scores_SeepC)

# scores plot
p2seep<-C_pca_Data_all_Seep %>%
  drop_na(Depth_seep)%>%
  ggplot(aes(x = PC1, y = PC2, color = log(PAR_calc_seep+1)))+
  geom_point(aes(size = -Depth_seep)) +
  #  coord_cartesian(xlim = c(-4, 7), ylim = c(-4, 4)) +
  scale_color_gradient(low = "black", high = "yellow")+
  #scale_shape_manual(values = c(22,16))+
#  scale_colour_viridis(limits = c(-.55,-.05))+
#  scale_fill_viridis(trans = "reverse")+
  geom_segment(data = PC_loadings_SeepC, aes(x=0,y=0,xend=PC1*10,yend=PC2*10),
               arrow=arrow(length=unit(0.1,"cm")), color = "grey")+
  annotate("text", x = PC_loadings_SeepC$PC1*10+0.1, y = PC_loadings_SeepC$PC2*10+.1,
           label = PC_loadings_SeepC$labels)+
  scale_size(limits = c(-.55,-0.05)) +
  guides(color=guide_legend(), size = guide_legend())+
  labs(color = "PAR", size = "Water Depth", 
       title = "Cabral Seep")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# save the plots
SeepPCA<-p1seep+p2seep+ 
  patchwork::plot_annotation("Seep", 
                             theme = theme(plot.title = element_text(size = rel(1.5), face = "bold", hjust = 0.5, 
                                                                     margin = margin(t = 10, b = 20, unit = "pt"))))

ggsave(plot = SeepPCA, filename = here("Output","SeepPCA.pdf"), width = 12, height = 6)



### plots of each parameter versus depth
Data %>%
  filter(Plate_Seep=="Seep", Salinity>10) %>% # I think I need to remove the bottle sample...
  pivot_longer(cols = Salinity:Ammonia_umolL) %>% 
  ggplot(aes(x = Depth_seep, y = value))+
  geom_point(aes(color = Day_Night))+
  geom_smooth(method = "lm")+
  scale_y_continuous(trans = "log10")+
  facet_wrap(~Location*name, scales = "free")

## Calculate the ranges for each parameter
ranges<-Data %>%
  filter(Plate_Seep=="Seep") %>% # just do the seep data
  select(Location, Salinity, TA: Lignin_Like, TempInSitu_seep) %>%
  pivot_longer(cols = c(Salinity, TA: Lignin_Like, TempInSitu_seep)) %>%
  drop_na()%>%
  group_by(Location, name)%>%
  summarise(min = round(min(value),2),
            max = round(max(value),2)) %>%
  mutate(range = paste0("[",min," - ",max,"]"))
  

# Make a correlation plot of everything versus depth
# function for correlation coef
cor_fun <- function(data) cor.test(data$value, data$Depth_seep, method = "pearson")%>% tidy()

cortest<-Data %>%
  filter(Plate_Seep=="Seep") %>% # just do the seep data
  select(Location, Salinity, TA: Lignin_Like, TempInSitu_seep, Depth_seep)%>%
  pivot_longer(cols = c(Salinity, TA: Lignin_Like, TempInSitu_seep)) %>%
  drop_na()%>%
  group_by(Location, name)%>%
  nest() %>%
  mutate(model = map(data, cor_fun)) %>%
  select(Location, name, model) %>%
  unnest(model)%>% # calculate correlation coefficient
  mutate(sig = ifelse(p.value<0.1, 1,0 ))# add a 1 if significant correlation (to 0.1 pvalue)

# make a plot of it
cortest %>%
  ggplot(aes(x = fct_reorder(name, abs(estimate)), y = Location, color = estimate, size = abs(estimate)))+
  geom_point()+
  geom_point(aes(shape = factor(sig)), size = 2, color = "white")+
  scale_color_gradient2(low = "#005AB5", high = "#DC3220",mid = "black" )+
  scale_size(range = c(0.1,10))+
  scale_shape_manual(values = c(NA,8))+
  coord_flip()+
  theme_bw()+
  guides(size = "none",
         shape = "none")+
  labs(x = "",
       y = "",
       color = "Correlation coefficient",
       title = "Correlation with water depth")


ggsave(here("Output","CorrelationPlot_seep.pdf"), height = 8, width = 5)


# just Varari

cortest %>%
  filter(Location == "Varari") %>%
  ggplot(aes(x = fct_reorder(name, estimate, .desc = TRUE),y = 1, color = estimate, size = abs(estimate)))+
  geom_point()+
  geom_point(aes(shape = factor(sig)), size = 2, color = "white")+
  scale_color_gradient2(low = "#005AB5", high = "#DC3220",mid = "black" )+
  scale_size(range = c(0.1,10))+
  scale_shape_manual(values = c(NA,8))+
  coord_flip()+
  theme_bw()+
  guides(size = "none",
         shape = "none")+
  labs(x = "",
       y = "",
       color = "Correlation coefficient",
       title = "Correlation with water depth")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

## same thing but everything versus salinity

cor_fun <- function(data) cor.test(data$value, data$Salinity, method = "pearson")%>% tidy()

cortest<-Data %>%
  filter(Plate_Seep=="Seep") %>% # just do the seep data
  select(Location, Salinity, TA: Lignin_Like, TempInSitu_seep)%>%
  pivot_longer(cols = c(TA: Lignin_Like, TempInSitu_seep)) %>%
  drop_na()%>%
  group_by(Location, name)%>%
  nest() %>%
  mutate(model = map(data, cor_fun)) %>%
  select(Location, name, model) %>%
  unnest(model)%>% # calculate correlation coefficient
  mutate(sig = ifelse(p.value<0.05, 1,0 ))# add a 1 if significant correlation (to 0.1 pvalue)

cortest %>%
  ggplot(aes(x = fct_reorder(name, abs(estimate)), y = Location, color = estimate, size = abs(estimate)))+
  geom_point()+
  geom_point(aes(shape = factor(sig)), size = 2, color = "white")+
  scale_color_gradient2(low = "#005AB5", high = "#DC3220",mid = "black" )+
  scale_size(range = c(0.1,10))+
  scale_shape_manual(values = c(NA,8))+
  coord_flip()+
  theme_bw()+
  guides(size = "none",
         shape = "none")+
  labs(x = "",
       y = "",
       color = "Correlation coefficient",
       title = "Correlation with Salinity")

ggsave(here("Output","CorrelationPlot_seepSalinity.pdf"), height = 8, width = 5)


# just Varari

cortest %>%
  left_join(ranges) %>% # join with the ranges
  filter(Location == "Varari") %>%
  ggplot(aes(x = fct_reorder(name, estimate),y = 1, color = estimate, size = abs(estimate)))+
  geom_point()+
  geom_point(aes(shape = factor(sig)), size = 2, color = "white")+
  geom_text(aes(y = 1.15, label = range), size = 4, color = "black")+
  ylim (0.8,1.4)+
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
       color = "Correlation coefficient",
       title = "Correlation with Salinity")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
      #  panel.grid.major.x = element_blank(),
      #  panel.grid.minor.x = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm")
        )

ggsave(here("Output","CorrelationPlot_seepSalinityV.png"), height = 8, width = 7)

