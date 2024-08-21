library(sf)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(grid)
library(extrafont)
library(mapview)
library(xkcd)

# Note: xkcd font is downloaded from https://github.com/ipython/xkcd-font.
#       import font with the extrafont package, and set working directory to where you 
#       downloaded the font, then import font with this code: 
#       font_import(paths = getwd(), pattern = "[X/x]kcd-script", prompt=FALSE)
#       I could not import it any other way, and this seems to be a known issue!

# Global variables: ------------------------------------------------------------

# Sampling area
box = c(xmin=14, xmax=16, ymin=67.1, ymax=67.31)


# Sampling data: ----------------------------------------------------------------

lat <- c(67.239483, 67.245283, 67.25145, 67.26085, 67.247767, 67.248683, 67.249683, 67.24805, 67.24795, 67.262517, 67.24495, 67.252383, 67.251917, 67.234674, 67.253265)
long <- c(14.6665, 14.552833, 14.470883, 14.434667, 14.899767, 14.894567, 14.882133, 14.8616, 14.84895, 14.601917, 14.654117, 14.603917, 14.63705, 14.734726, 14.470267)
loc <- c("Saltenfjorden", "Saltenfjorden", "Saltenfjorden", "Saltenfjorden", "Skjerstadfjorden", "Skjerstadfjorden", "Skjerstadfjorden", "Skjerstadfjorden", "Skjerstadfjorden", "Saltenfjorden", "Saltenfjorden", "Saltenfjorden", "Saltenfjorden", "Skjerstadfjorden", "Saltenfjorden")
station_id <- c("SAL30", "SAL100", "SAL200", "SAL350", "SKJ30", "SKJ100", "SKJ200", "SKJ350", "SKJ500", "SAL350", "SAL100", "SAL200", "SAL200", "SKJ500", "SAL350")
label_id <- c("SAL 30m", "SAL 100m", "SAL 200m", "SAL 350m", "SKJ 30m", "SKJ 100m", "SKJ 200m", "SKJ 350m", "SKJ 500m", "SAL 350m", "SAL 100m", "SAL 200m", "SAL 200m", "SKJ 500m", "SAL 350m")
season <- c("Summer", "Summer", "Summer", "Summer", "Summer", "Summer", "Summer", "Summer", "Summer", "Winter", "Winter", "Winter", "Winter", "Summer", "Summer")
type <- c("ecology", "ecology","ecology","ecology","ecology","ecology","ecology","ecology","genetics","genetics","ecology","ecology","ecology","genetics","genetics")
samples <- data.frame(lat, long, station_id, label_id, season, type)

samples <- samples[-c(2:4, 12),]
rownames(samples) <- 1:nrow(samples)


# Nordland depth data: --------------------------------------------------------------

# Data link:
# https://kartkatalog.geonorge.no/metadata/sjoekart-dybdedata/2751aacf-5472-4850-a208-3532a51c529a

# Product specification
# https://register.geonorge.no/data/documents/Produktspesifikasjoner_sjokart-dybdedata_v2_produktspesifikasjon-kartverket-dybdedata-20201001_.pdf

# Load depth area
depth_nrdl <- read_sf("Basisdata_18_Nordland_25833_Dybdedata_GML.gml", "Dybdeareal")

# Format data (naming convention and WGS82 projection)
st_geometry(depth_nrdl) <- "geometry"
depth_nrdl <- st_transform(depth_nrdl, 4326)


# Subsample sf object
sf_use_s2(FALSE) # Can't subset without this function turned off

box = c(xmin=14, xmax=16, ymin=67.1, ymax=67.31)
depth_nrdl_crop <- st_crop(depth_nrdl, box)

# Plot depth map
ggplot() +
  geom_sf(data = depth_nrdl_crop, mapping = aes(fill = as.factor(-minimumsdybde)), colour = NA) +
  # geom_sf(data = nor_sf[nor_sf$kommunenummer == "1804",],) + # map layer +
  scale_fill_grey(aesthetics = "fill") + 
  theme(legend.position = "none", 
        # panel.background = element_rect(fill = "white"),
        # panel.border = element_blank(),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank()
  )

# Load shoreline maps: ---------------------------------------------------

# Link to geojson file
# https://github.com/robhop/fylker-og-kommuner/blob/main/Kommuner-L.geojson

library(sf)
nor_sf <- read_sf("Kommuner-L.geojson") # Sampling map shoreline

nor_sf_crop <- st_crop(nor_sf, box) # Sampling area
nor_sf_fylker <- read_sf("Fylker-L.geojson") # Inset map shoreline

# Sampling map -----------------------------------------------------------------
p <- ggplot() +
  
  # Depth data
  # geom_sf(data = depth_bodo, mapping = aes(fill = -minimumsdybde), colour = NA) + # depth bod? kommune
  geom_sf(data = depth_nrdl_crop, mapping = aes(fill = as.factor(-minimumsdybde)), colour = NA) + # Depth Nordland fylke
  scale_fill_grey(aesthetics = "fill", guide = "none") + 
  
  # Shoreline data
  # geom_sf(data = nor_sf[nor_sf$kommunenummer == "1804",],) + # map layer Bod? kommune subset
  geom_sf(data = nor_sf_crop, color = "azure4") + # Map area subset
  scale_size_manual(values=c(2,3.5)) +
  xlim(14.3, 15) +
  ylim(67.2, 67.3) +
  
  #grid and legend manipulation
  theme(legend.position = "bottom", 
        legend.key = element_rect(fill = NA),
        panel.background = element_rect(fill = "lightgray"),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major = element_blank(),
        axis.title=element_blank(),
        text=element_text(size=14, family="xkcd"),
        panel.grid.minor = element_blank()) +
  labs(shape = "Sample type:") + 
  
  # Datapoints
  geom_point(data = samples, aes(x = long, y = lat, shape = type),color = "black") + # All sample layers
  geom_point(data = samples[samples$type == "genetics",], aes(x = long, y = lat), pch = 1, size = 5, color = "black") + # Gen samples
  geom_point(data = samples, aes(x = long, y = lat)) +
  scale_shape_manual(values=c(19,1)) +
  # geom_text_repel(data = samples, label = station_id, x = long, y = lat, label.size = 0, fill = "transparent")
  # annotate("text", family="xkcd", label = samples$station_id, x = samples$long + 0.015, y = samples$lat, size = 3, color = "black") # static lables
  geom_text_repel(data = samples, x = samples$long, y = samples$lat, label = samples$label_id, family="xkcd", size = 3, color = "black")+ # sample labels
  
  # Location labels
  geom_text(family = "xkcd")+
  annotate("text", label = "Bodo", size = 3, x = 14.34, y = 67.262, family="xkcd") + # text layer for Bod? city
  annotate("text", label = "Saltstrumen", size = 3, x = 14.6233, y = 67.2350, family="xkcd") + # text layer for Saltstraumen
  annotate("text", label = "Saltenfjorden", size = 4, x = 14.54, y = 67.265, family="xkcd") + # text layer for 
  annotate("text", label = "Skjerstadfjorden", size = 4, x = 14.85, y = 67.26, family="xkcd") + # text layer for 
  
  annotation_scale(aes(style = "ticks", location = "br"), text_family = "xkcd") # scale bar 

p

# Add mountains: ---------------------------------------------------------------

mnts <- read_sf("Basisdata_1804_Bodo_25833_N50Hoyde_GML.gml", 
                "Terrengpunkt")

st_geometry(mnts) <- "geometry"
mnts <- st_transform(mnts, 4326)


# Subsample sf object
sf_use_s2(FALSE) # Can't subset without this function turned off

box = c(xmin=14, xmax=16, ymin=67.1, ymax=67.31)
mnts <- st_crop(mnts, box)

p <- p + 
  # geom_sf(data = mnts)
  # geom_sf_text(data = mnts[mnts$høyde > 200,], label = "^", size = log(mnts[mnts$høyde > 200,]$høyde), family="xkcd Script", color = "gray40") # scale factor 100
  geom_sf_text(data = mnts, label = "^", size = (mnts$høyde)/50, family="xkcd Script", color = "gray40") # scale factor 100


# Add funny bits: ------------------------------------------------------------
b <- p +  
  annotate("text", label = "Salty twirly\nnarrow\nwater passage", size = 3, x = 14.63, y = 67.21, family="xkcd", color = "gray30") + # text layer for Bod? city
  annotate("text", label = "the abyss!", size = 3, x = 14.975, y = 67.276, family="xkcd", color = "gray30", angle = -30) +
  annotate("text", label = "Cool\namphipods", size = 3, x = 14.31, y = 67.23, family="xkcd", color = "gray30", angle = 15) +
  annotate("text", label = "My\nhouse", size = 3, x = 14.75, y = 67.205, family="xkcd", color = "gray30") 

# add pointer lines
pointer_1 <- data.frame(x1 = 14.62, x2 = 14.611, y1 = 67.219, y2 =67.23) # saltstraumen
pointer_2 <- data.frame(x1 = 14.73, x2 = 14.68, y1 = 67.205, y2 =67.20392663906014)

xkcd_map <-  b + geom_curve(mapping = aes(x = x1, y = y1, xend = x2, yend = y2), data = pointer_1, 
                            arrow = arrow(length = unit(0.03, "npc")), size = 0.9, curvature = 0.3, colour = "gray30") +
  geom_point(aes(x=14.672876683009468, y=67.20392663906014), size = 0.5) +
  geom_curve(mapping = aes(x = x1, y = y1, xend = x2, yend = y2), data = pointer_2, 
             arrow = arrow(length = unit(0.03, "npc")), size = 0.9, curvature = -0.2, colour = "gray40")



# Inset for sampling map: ------------------------------------------------------
inset <- ggplot() + 
  geom_sf(data = nor_sf_fylker[nor_sf_fylker$name == "Nordland",], ) +
  xlim(11, 18) +
  ylim(65.5, 68.3) +
  # Add guide lines to plot
  # guides(size = FALSE)
  # geom_hline(yintercept = 66, lty = 2, colour = "red") +
  # geom_hline(yintercept = 68, lty = 2, colour = "red") +
  # geom_vline(xintercept = 12, lty = 2, colour = "red") +
  # geom_vline(xintercept = 14, lty = 2, colour = "red") +
  
  # Rectangle to show area if interest
  geom_rect(aes(xmin = 14.3, xmax = 15, ymin = 67.15, ymax = 67.4), size = 1, color = "red", fill = NA) + 
  annotate("text", label = "Nordland", size = 4, x = 14, y = 66.3, family="xkcd", angle = 60) +
  
  # grid and legend manipulation
  theme(rect=element_rect(fill="transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        axis.ticks = element_blank(),
        plot.margin = margin(-1, -1, -1, -1, "cm"),
        panel.background= element_rect(color = "black", fill = "white")) 
inset

# Add inset to sampling map
library(grid)

xkcd_map
print(inset, vp = viewport(0.2, 0.8, width = 0.1, height = 0.1))

