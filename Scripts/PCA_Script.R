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

# load the 24 hour chemistry data #####################
Data<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv")

# Load water level data 
# Bring in the depth data
WLPath<-here("Data", "Varari", "WL")
files <- dir(path = WLPath,pattern = ".csv", full.names = TRUE)

WL_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(Date = date, Depth)%>%
  mutate(Location = "Varari")

WLPath<-here("Data", "Cabral", "WL")
files <- dir(path = WLPath,pattern = ".csv", full.names = TRUE)

#Cabral
WL_Cabral<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(Date = date, Depth)%>%
  mutate(Location = "Cabral")

WL_all <- bind_rows(WL_Varari, WL_Cabral) %>%
  rename(DateTime = Date)

# Bind with the chem data
Data<-Data %>%
  left_join(WL_all) 

### Varari #####
## There seems to be a contaminated nutrient sample for V2 Low tide on the 8/8/2021.  Remove this point
remove<-Data %>% filter(CowTagID=="V2", Tide ==" Low",Day_Night=="Day", Date == ymd("2021-08-08"))

remove3<-Data %>% filter(CowTagID== "C4", Tide =="Low", Day_Night=="Night", Date == ymd("2021-08-09"))

## also filter out all the data from the first low tide where the water level was super high

remove2<-Data %>% filter(Tide =="Low", Day_Night=="Day", Date == ymd("2021-08-06"))

# extract the params for the PCA
V_pca_Data<-Data %>%
  anti_join(remove)%>%
  anti_join(remove2)%>%
  filter(Location == "Varari", Plate_Seep=="Plate") %>%
  #select(Salinity,pH,Phosphate_umolL, Silicate_umolL, NN_umolL, Ammonia_umolL ) %>%
  select(Salinity,pH,Phosphate_umolL:Lignin_Like )%>%
  drop_na()

# Run the PCA
pca_V <- prcomp(V_pca_Data, scale. = TRUE, center = TRUE)

# Extract the scores and loadings
PC_scores <-as_tibble(pca_V$x[,1:2])
PC_loadings<-as_tibble(pca_V$rotation) %>%
  bind_cols(labels = rownames(pca_V$rotation))

# Put it with all the original data
V_pca_Data_all<-Data %>%
  anti_join(remove)%>%
  anti_join(remove2)%>%
  filter(Location == "Varari", Plate_Seep=="Plate") %>%
 # drop_na(Salinity,pH,Phosphate_umolL, Silicate_umolL, NN_umolL, Ammonia_umolL) %>%
  drop_na(Salinity,pH,Phosphate_umolL:Lignin_Like )%>%
  bind_cols(PC_scores)

# scores plot
p1<-V_pca_Data_all %>%
  ggplot(aes(x = PC1, y = PC2, color = Tide, shape = Day_Night))+
  geom_point() +
  coord_cartesian(xlim = c(-4, 7), ylim = c(-6, 4)) +
  scale_shape_manual(values = c(22,16))+
  scale_colour_hue(l = 45)+
  scale_fill_hue(l = 45)+
  ggforce::geom_mark_ellipse(
    aes(fill = Tide, label = paste(Day_Night, Tide), color = Tide), 
    alpha = .15, show.legend = FALSE,  label.buffer = unit(1, "mm"))+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        
# loadings plot 
p2<-PC_loadings %>%
  ggplot(aes(x=PC1, y=PC1, label=labels))+
    geom_segment(aes(x=0,y=0,xend=PC1*10,yend=PC2*10),
      arrow=arrow(length=unit(0.1,"cm")), color = "grey")+
  annotate("text", x = PC_loadings$PC1*10+0.1, y = PC_loadings$PC2*10+.1,
           label = PC_loadings$labels)+
  coord_cartesian(xlim = c(-4, 8), ylim = c(-6, 4)) +
   theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

VarariPCA<-p1+p2+ 
  patchwork::plot_annotation("Varari Plates", 
                             theme = theme(plot.title = element_text(size = rel(1.5), face = "bold", hjust = 0.5, 
                                                                     margin = margin(t = 10, b = 20, unit = "pt"))))

ggsave(plot = VarariPCA, filename = here("Output","VarariPCA.pdf"), width = 12, height = 6)
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
  bind_cols(labels = rownames(pca_C$rotation))

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

# loadings plot 
p2c<-PC_loadingsC %>%
  ggplot(aes(x=PC1, y=PC1, label=labels))+
  geom_segment(aes(x=0,y=0,xend=PC1*10,yend=PC2*10),
               arrow=arrow(length=unit(0.1,"cm")), color = "grey")+
  annotate("text", x = PC_loadingsC$PC1*10+0.1, y = PC_loadingsC$PC2*10+.1,
           label = PC_loadingsC$labels)+
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

CabralPCA<-p1c+p2c+ 
  patchwork::plot_annotation("Cabral Plates", 
                             theme = theme(plot.title = element_text(size = rel(1.5), face = "bold", hjust = 0.5, 
                                                                     margin = margin(t = 10, b = 20, unit = "pt"))))

ggsave(plot = CabralPCA, filename = here("Output","CabralPCA.pdf"), width = 12, height = 6)

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
  drop_na(Depth)%>%
  ggplot(aes(x = PC1, y = PC2, color = -Depth))+
   geom_point(aes(size = -Depth, shape = Day_Night)) +
  #  coord_cartesian(xlim = c(-4, 7), ylim = c(-4, 4)) +
  #scale_shape_manual(values = c(22,16))+
  scale_colour_viridis(limits = c(-1.1,-0.4))+
   geom_segment(data = PC_loadings_Seep, aes(x=0,y=0,xend=PC1*3,yend=PC2*3),
               arrow=arrow(length=unit(0.1,"cm")), color = "grey")+
  annotate("text", x = PC_loadings_Seep$PC1*3+0.1, y = PC_loadings_Seep$PC2*3+.1,
           label = PC_loadings_Seep$labels)+
  scale_size(limits = c(-1.1,-0.4))+
  guides(color=guide_legend(), size = guide_legend())+
  labs(color = "Water Depth", size = "Water Depth", shape = "Day or Night",
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
  drop_na(Depth)%>%
  ggplot(aes(x = PC1, y = PC2, color = -Depth))+
  geom_point(aes(size = -Depth, shape = Day_Night)) +
  #  coord_cartesian(xlim = c(-4, 7), ylim = c(-4, 4)) +
  #scale_shape_manual(values = c(22,16))+
  scale_colour_viridis(limits = c(-.55,-.05))+
#  scale_fill_viridis(trans = "reverse")+
  geom_segment(data = PC_loadings_SeepC, aes(x=0,y=0,xend=PC1*10,yend=PC2*10),
               arrow=arrow(length=unit(0.1,"cm")), color = "grey")+
  annotate("text", x = PC_loadings_SeepC$PC1*10+0.1, y = PC_loadings_SeepC$PC2*10+.1,
           label = PC_loadings_SeepC$labels)+
  scale_size(limits = c(-.55,-0.05)) +
  guides(color=guide_legend(), size = guide_legend())+
  labs(color = "Water Depth", size = "Water Depth", shape = "Day or Night",
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
  ggplot(aes(x = Depth, y = value))+
  geom_point(aes(color = Day_Night))+
  geom_smooth(method = "lm")+
  scale_y_continuous(trans = "log10")+
  facet_wrap(~Location*name, scales = "free")

# there is something wrong with how Depth is joining