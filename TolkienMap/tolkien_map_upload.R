library(sf)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(ggimage)
library(geomtextpath)
library(grid)
library(extrafont)
library(mapview)
library(tidyverse)

# Import map

# Link to geojson file
# https://github.com/robhop/fylker-og-kommuner/blob/main/Kommuner-L.geojson

kommuner_sf <- read_sf("Kommuner-L.geojson") # Sampling map shoreline

kommuner_sf <- st_transform(kommuner_sf, 4326)

# Information to subset data

# Districs ---------------------------------------------------------------------
# Districts Nordland
# https://no.wikipedia.org/wiki/Nordland#Kommuner
nordlandInfo <- read.csv("nordlandInfo.csv", sep = "\t")

# Subset to desired municipals
temp <- nordlandInfo[nordlandInfo$Distrikt %in% c("Salten", "Lofoten", "Ofoten"),]$Nr
Salten <- kommuner_sf[kommuner_sf$kommunenummer %in% temp,]


# Fylker -----------------------------------------------------------------------
# Add Fylke to dataset
kom_meta <- read.csv("kommuner.csv", sep = "\t")
norge_sf <- merge(x = kommuner_sf, y = kom_meta[,-c(2,3,5)], by.x = "id", by.y = "Nr", all.x = TRUE)

# Oslo is broken for some reason, adding manually
norge_sf[1,5] <- "Oslo"
norge_sf[1,6] <- 45412

sf_use_s2(FALSE) #change translation method

fylker_sf <- norge_sf %>% 
  group_by(Fylke) %>% 
  st_make_valid() %>%
  summarize() 

# pre-load ---------------------------------------------------------------------

# set location and run bounding box (bbox) below before running pre-load
norge_fill <- norge_sf %>% 
  st_make_valid() %>% # union does not work unless the shapes are considered valid first
  st_union()
norge_crop <- st_crop(norge_fill, bbox)


# Add extra line outside of coast
coast_line_1 <- st_buffer(x = norge_fill, dist = (bbox[3]-bbox[1]) / 100, joinStyle = "MITRE", mitreLimit = 1)
coast_line_2 <- st_buffer(x = norge_fill, dist = (bbox[3]-bbox[1]) / 200, joinStyle = "MITRE", mitreLimit = 1)

coast_crop_1 <- st_crop(coast_line_1, bbox)
coast_crop_2 <- st_crop(coast_line_2, bbox)

# Select area to map -----------------------------------------------------------

loc = "Vestvågøy"
loc = "Bodø"

df_sf <- norge_sf[norge_sf$kommunenavn == loc | norge_sf$Fylke == loc,] 

# SET BOUNDING BOX run before pre-load
bbox <-  df_sf %>% 
  st_as_sf(coords = c("X2","X1"), crs = 4326) %>% 
  st_bbox()

#combine geometries
sf_use_s2(FALSE) # Can't subset without this function turned off

# remove borders
df_sf <- df_sf %>% 
  st_make_valid() %>% # union does not work unless the shapes are considered valid first
  st_union()


# smooth map -------------------------------------------------------------------
sf_use_s2(TRUE) #change translation method

df_sf <- st_crop(df_sf, bbox)

norge_crop_round <- norge_crop %>% st_buffer(dist = 0.0001, endCapStyle = "ROUND")
coast_crop_2_round <- coast_crop_2 %>% st_buffer(dist = 0.0002, endCapStyle = "ROUND")
coast_crop_1_round <- coast_crop_1 %>% st_buffer(dist = 0.0002, endCapStyle = "ROUND")

sf_use_s2(FALSE) #change translation method


# Base map ---------------------------------------------------------------------

p <- ggplot() + 
  geom_sf(data = norge_fill, linewidth = 1, fill = "white", color = "black") + 
  geom_sf(data = coast_line_2, linewidth = 1, fill = "transparent", color = "gray20") +
  geom_sf(data = coast_line_1, linewidth = 1, fill = "transparent", color = "gray20") +
  geom_sf(data = df_sf, linewidth = 1.5, fill = "white", color = "black") +
  xlim(bbox[c(1,3)]) +
  ylim(bbox[c(2,4)]) +
  theme_void() +
  theme(legend.position = "bottom", 
        legend.key = element_rect(fill = NA),
        panel.background = element_rect(fill = "lightgray"),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major = element_blank(),
        axis.title=element_blank(),
        text=element_text(size=14, family="AnironC"),
        panel.grid.minor = element_blank()) +
  labs(shape = "Vegkategori:")
p



# N100 kartdata: ---------------------------------------------------------------

# Veier ####
n100_vei <- read_sf("Basisdata_18_Nordland_25833_N100Samferdsel_GML.gml",
                    "veglenke")

n100_vei <- st_transform(n100_vei, 4326)

n100_vei = st_intersection(n100_vei, df_sf)


# Elver ####
n100_elv <- read_sf("Basisdata_18_Nordland_25833_N100Arealdekke_GML.gml",
                    "Elv")

n100_elv <- st_transform(n100_elv, 4326)

n100_elv = st_intersection(n100_elv, df_sf)

# Mountains: terrengpunkt ####

#load data
mnts <- read_sf("Basisdata_18_Nordland_25833_N100Hoyde_GML.gml",
                "terrengpunkt")

# Transform coordinate translation
mnts <- st_transform(mnts, 4326)

# Mask of data to mapping area
mnts = st_intersection(mnts, df_sf)

mnts <- mnts %>% 
  dplyr::mutate(lon = st_coordinates(.)[,1],
                lat = st_coordinates(.)[,2])

mnts_1 <- c("tolkien_mnts/1_1.png", "tolkien_mnts/1_2.png", "tolkien_mnts/1_3.png", "tolkien_mnts/1_4.png", "tolkien_mnts/1_5.png")
mnts_2 <- c("tolkien_mnts/2_1.png", "tolkien_mnts/2_2.png", "tolkien_mnts/2_3.png", "tolkien_mnts/2_4.png", "tolkien_mnts/2_5.png", "tolkien_mnts/2_6.png")   
mnts_3 <- c("tolkien_mnts/3_1.png", "tolkien_mnts/3_2.png", "tolkien_mnts/3_3.png", "tolkien_mnts/3_4.png", "tolkien_mnts/3_5.png", "tolkien_mnts/3_6.png") 
mnts_4 <- c("tolkien_mnts/4_1.png", "tolkien_mnts/4_2.png", "tolkien_mnts/4_3.png", "tolkien_mnts/4_4.png", "tolkien_mnts/4_5.png")
mnts_5 <- c("tolkien_mnts/5_1.png", "tolkien_mnts/5_2.png", "tolkien_mnts/5_3.png", "tolkien_mnts/5_4.png", "tolkien_mnts/5_5.png", "tolkien_mnts/5_6.png")
mnts_6 <- c("tolkien_mnts/6_1")
mnts_7 <- c("tolkien_mnts/7_1.png", "tolkien_mnts/7_2.png")
mnts_8 <- c("tolkien_mnts/8_1.png", "tolkien_mnts/8_2.png", "tolkien_mnts/8_3.png")


# Randomly assign mountain icon to point
mnts <- mnts %>%
  mutate(IMAGE = ifelse(høyde > 100 , sample(mnts_8), 
                        ifelse(høyde < 100, sample(mnts_1), "nada")))

mnts <- mnts[order(-mnts$lat),]


# Stedsnavn ####
n100_stedsnavn <- read_sf("Basisdata_18_Nordland_25833_N100Stedsnavn_GML.gml",
                          "StedsnavnTekst")

n100_stedsnavn <- st_transform(n100_stedsnavn, 4326)
n100_stedsnavn = st_intersection(n100_stedsnavn, df_sf)
n100_stedsnavn_x <- n100_stedsnavn[!duplicated(n100_stedsnavn$fulltekst),]

n100_navn_kommune <- n100_stedsnavn[n100_stedsnavn$navneobjekttype %in% c("tettsted","by", "kommune", "tettsted"),]

# plot -------------------------------------------------------------------------

q <- p + 
  geom_sf(data = n100_vei[n100_vei$typeVeg == "enkelBilveg" & n100_vei$vegkategori == c("F"),], aes(linetype = vegkategori), color = "#FB2445") +
  geom_sf(data = n100_vei[n100_vei$typeVeg == "enkelBilveg" & n100_vei$vegkategori == c("E"),], aes(linetype = vegkategori), color = "#FB2445") +
  geom_sf(data = n100_vei[n100_vei$typeVeg == "enkelBilveg" & n100_vei$vegkategori == c("R"),], aes(linetype = vegkategori), color = "#FB2445") +
  scale_linetype_manual(values = c("F" = "dotted", "R" = "dashed","E" = "solid")) +
  geom_sf(data = n100_elv, linewidth = 1, color = "black") +
  geom_image(data = mnts[mnts$høyde > 800,], aes(x = lon, y = lat, image = IMAGE), size = 0.08) +
  geom_image(data = mnts[mnts$høyde > 100,], aes(x = lon, y = lat, image = IMAGE), size = 0.06) +
  geom_image(data = mnts[mnts$høyde < 100,], aes(x = lon, y = lat, image = IMAGE), size = 0.03) +
  # geom_sf_text(data = mnts, label = "^", size = (mnts$høyde)/100, family="xkcd Script", color = "black") + # scale factor 100
  geom_sf(data = st_centroid(n100_navn_kommune), shape = 0,) + 
  geom_sf_text(data = n100_stedsnavn_x, aes(label = fulltekst), color = "#FB2445", family = "Aniron", fontface = "bold") +
  geom_sf_text(data = n100_navn_kommune, aes(label = fulltekst), color = "#FB2445", family = "Aniron", fontface = "bold") +
  annotation_scale(aes(style = "ticks", location = "br"), text_family = "AnironC") # scale bar 
q

