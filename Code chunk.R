knitr::opts_chunk$set

####load libraries, etc
library(tidyverse)
library(sf)
library(spdep)
library(caret)
library(ckanr)
library(FNN)
library(grid)
library(knitr)
library(gridExtra)
library(ggcorrplot)
library(stargazer)
library(mapview)
library(osmdata)
library(tidycensus)
library(tidygeocoder)
library(raster)
library(rnaturalearth)
library(RColorBrewer)
library(rnaturalearthdata)
library(geosphere)

mapTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 14,colour = "black"),
    plot.subtitle=element_text(face="italic"),
    plot.caption=element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),axis.title = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2)
  )
}

plotTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 14,colour = "black"),
    plot.subtitle = element_text(face="italic"),
    plot.caption = element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_line("grey80", size = 0.1),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    strip.background = element_rect(fill = "grey80", color = "white"),
    strip.text = element_text(size=12),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(colour = "black", face = "italic"),
    legend.text = element_text(colour = "black", face = "italic"),
    strip.text.x = element_text(size = 14)
  )
}

palette5 <- c("#c8ddfa", "#8cb5ed", "#5890db",   "#2868bd", "#023578")
paletteMap <- c("gray90","gray70","gray50","gray30","gray10")

qBr <- function(df, variable, rnd) {
  if (missing(rnd)) {
    as.character(quantile(round(df[[variable]],0),
                          c(.01,.2,.4,.6,.8), na.rm=T))
  } else if (rnd == FALSE | rnd == F) {
    as.character(formatC(quantile(df[[variable]]), digits = 3),
                 c(.01,.2,.4,.6,.8), na.rm=T)
  }
}

q5 <- function(variable) {as.factor(ntile(variable, 5))}

#nearest neighbor function
nn_function <- function(measureFrom,measureTo,k) {
  measureFrom_Matrix <- as.matrix(measureFrom)
  measureTo_Matrix <- as.matrix(measureTo)
  nn <-   
    get.knnx(measureTo, measureFrom, k)$nn.dist
  output <-
    as.data.frame(nn) %>%
    rownames_to_column(var = "thisPoint") %>%
    gather(points, point_distance, V1:ncol(.)) %>%
    arrange(as.numeric(thisPoint)) %>%
    group_by(thisPoint) %>%
    summarize(pointDistance = mean(point_distance)) %>%
    arrange(as.numeric(thisPoint)) %>%
    dplyr::select(-thisPoint) %>%
    pull()
  
  return(output)  
}
#projected to NAD 1983 StatePlane Florida East FIPS 0901 Feet


#STUDY AREA
#miamiBound <- st_read("E:/Upenn/CPLN508/miami/2_Miami-Prediction/Raw Data/Municipal_Boundary.geojson") %>%
miamiBound <- st_read("/Users/annaduan/Documents/GitHub/2_Miami\ Prediction/Raw\ Data/Municipal_Boundary.geojson") %>%
  filter(NAME == "MIAMI BEACH" | NAME == "MIAMI") %>%
  st_union() %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')



#STUDY AREA OSM (not projected so that it works)
#miamiBoundOSM <- st_read("E:/Upenn/CPLN508/miami/2_Miami-Prediction/Raw Data/Municipal_Boundary.geojson") %>%
miamiBoundOSM <- st_read("/Users/annaduan/Documents/GitHub/2_Miami\ Prediction/Raw\ Data/Municipal_Boundary.geojson") %>%
  filter(NAME == "MIAMI BEACH" | NAME == "MIAMI") %>%
  st_union()



#HOUSE DATA
#houses <- st_read("E:/Upenn/CPLN508/miami/2_Miami-Prediction/Raw Data/studentsData.geojson") %>%
houses <- st_read("/Users/annaduan/Documents/GitHub/2_Miami\ Prediction/Raw\ Data/studentsData.geojson") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658') %>%
  st_centroid()




#HOUSE DATA OSM (Not projected)
#housesOSM <- st_read("E:/Upenn/CPLN508/miami/2_Miami-Prediction/Raw Data/studentsData.geojson")
housesOSM <- st_read("/Users/annaduan/Documents/GitHub/2_Miami\ Prediction/Raw\ Data/studentsData.geojson")



#CENSUS
census_api_key("d9ebfd04caa0138647fbacd94c657cdecbf705e9", install = TRUE, overwrite = TRUE)
#read in: vacant property, total housing units, mhhinc, white, population, owner occ, renter occ, travel time to work
acs <-
  get_acs(geography = "tract", variables = c("B25002_003E", "B25001_001E", "B19013_001E", "B01001A_001E", "B01003_001E", "B07013_002E", "B07013_003E", "B08012_001E", "B25104_001E"), year=2018, state=12, county=086, geometry=T) %>%
  st_transform('ESRI:102658')
#filter for Miami/Miami beach tracts
acs <-
  rbind(
    st_centroid(acs)[miamiBound,] %>%
      st_drop_geometry() %>%
      left_join(acs) %>%
      st_sf() %>%
      mutate(inMiami = "YES"),
    st_centroid(acs)[miamiBound, op = st_disjoint] %>%
      st_drop_geometry() %>%
      left_join(acs) %>%
      st_sf() %>%
      mutate(inMiami = "NO")) %>%
  filter(inMiami == "YES") %>%
  dplyr::select(-inMiami)
#long to wide form
acs <-
  acs %>%
  dplyr::select(-moe, -GEOID) %>%
  spread(variable, estimate) %>%
  dplyr::select(-geometry) %>%
  rename(vacantUnits = B25002_003,
         totalUnits = B25001_001,
         medHHInc = B19013_001,
         white = B01001A_001,
         population = B01003_001,
         ownerOcc = B07013_002,
         renterOcc = B07013_003,
         timeToWork = B08012_001,
         monthhousingcost = B25104_001)
acs["119", "medHHInc"] = 33194   #input value from nearby tract in same neighborhod because NA was messing up MAE
acs %>% na.omit()


#mutate
acs <-
  acs %>%
  mutate(pctVacant = ifelse(totalUnits > 0, vacantUnits / totalUnits, 0),
         pctWhite = ifelse(population > 0, white / population, 0),
         totalOcc = ownerOcc + renterOcc,
         pctRenterOcc = ifelse(totalOcc > 0, renterOcc / totalOcc, 0)) %>%
  dplyr::select(-totalUnits,-vacantUnits,-totalUnits,-population,-white, -ownerOcc, -renterOcc, -totalOcc)


#OSM BBOX (uses the non-projected base)
xmin = st_bbox(miamiBoundOSM)[[1]]
ymin = st_bbox(miamiBoundOSM)[[2]]
xmax = st_bbox(miamiBoundOSM)[[3]]  
ymax = st_bbox(miamiBoundOSM)[[4]]



#FOOD AND BEVERAGE SPOTS
foodBev <- opq(bbox = c(xmin, ymin, xmax, ymax)) %>%
  add_osm_feature(key = 'amenity', value = c("bar","pub","restaurant","cafe")) %>%
  osmdata_xml(filename = 'foodBev.osm')
#project
foodBev <- sf::st_read('foodBev.osm', layer = 'points') %>%
  st_as_sf(coords = c("LON", "LAT"), crs = EPSG:3857, agr = "constant") %>%
  st_transform('ESRI:102658')
#filter for facilities in study area
foodBev <- rbind(
  st_centroid(foodBev)[miamiBound,] %>%
    st_drop_geometry() %>%
    left_join(foodBev) %>%
    st_sf() %>%
    mutate(inMiami = "YES"),
  foodBev[miamiBound, op = st_disjoint] %>%
    st_drop_geometry() %>%
    left_join(foodBev) %>%
    st_sf() %>%
    mutate(inMiami = "NO")) %>%
  filter(inMiami == "YES") %>%
  dplyr::select(name)



#COASTLINE
Coastline<-opq(bbox = c(xmin, ymin, xmax, ymax)) %>%
  add_osm_feature("natural", "coastline") %>%
  osmdata_sf()
#add to housesOSM and convert to miles, then add to houses
housesOSM <-
  housesOSM %>%  
  mutate(CoastDist=(geosphere::dist2Line(p=st_coordinates(st_centroid(housesOSM)),
                                         line=st_coordinates(Coastline$osm_lines)[,1:2])*0.00062137)[,1])
houses <-
  houses %>%
  mutate(distWater = housesOSM$CoastDist,
         SPSqFt = ifelse(!is.na(ActualSqFt)&!is.na(SalePrice), SalePrice / ActualSqFt, 0))



#PARKS
muniParks <- st_read("https://opendata.arcgis.com/datasets/16fe02a1defa45b28bf14a29fb5f0428_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658') %>%
  dplyr::select(NAME, ADDRESS, CITY, CLASS, Shape__Area)

countyParks <- st_read("https://opendata.arcgis.com/datasets/aca1e6ff0f634be282d50cc2d534a832_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658') %>%
  dplyr::select(NAME, ADDRESS, CITY, CLASS, Shape__Area)
parks <- bind_rows(muniParks, countyParks) %>%
  filter(CITY == "Miami" | CITY == "Miami Beach") %>%
  mutate(counter = 1)



#SCHOOL DISTRICT
schoolDist <- st_read("https://opendata.arcgis.com/datasets/bc16a5ebcdcd4f3e83b55c5d697a0317_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658') %>%
  dplyr::select(ID)



#PUBLIC SCHOOL CATCHMENT/ATTENDANCE ZONES
#elementary
elementary <- st_read("https://opendata.arcgis.com/datasets/19f5d8dcd9714e6fbd9043ac7a50c6f6_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')
elementary <- rbind(
  st_centroid(elementary)[miamiBound,] %>%
    st_drop_geometry() %>%
    left_join(elementary) %>%
    st_sf() %>%
    mutate(inMiami = "YES"),
  st_centroid(elementary)[miamiBound, op = st_disjoint] %>%
    st_drop_geometry() %>%
    left_join(elementary) %>%
    st_sf() %>%
    mutate(inMiami = "NO")) %>%
  filter(inMiami == "YES") %>%
  dplyr::select(NAME)
#middle
middle <- st_read("https://opendata.arcgis.com/datasets/dd2719ff6105463187197165a9c8dd5c_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')
middle <- rbind(
  st_centroid(middle)[miamiBound,] %>%
    st_drop_geometry() %>%
    left_join(middle) %>%
    st_sf() %>%
    mutate(inMiami = "YES"),
  st_centroid(middle)[miamiBound, op = st_disjoint] %>%
    st_drop_geometry() %>%
    left_join(middle) %>%
    st_sf() %>%
    mutate(inMiami = "NO")) %>%
  filter(inMiami == "YES") %>%
  dplyr::select(NAME)
#high
high <- st_read("https://opendata.arcgis.com/datasets/9004dbf5f7f645d493bfb6b875a43dc1_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')
high <- rbind(
  st_centroid(high)[miamiBound,] %>%
    st_drop_geometry() %>%
    left_join(high) %>%
    st_sf() %>%
    mutate(inMiami = "YES"),
  st_centroid(high)[miamiBound, op = st_disjoint] %>%
    st_drop_geometry() %>%
    left_join(high) %>%
    st_sf() %>%
    mutate(inMiami = "NO")) %>%
  filter(inMiami == "YES") %>%
  dplyr::select(NAME)



#PUBLIC TRANSPORTATION
#bus
bus <- st_read("https://opendata.arcgis.com/datasets/021adadcf6854f59852ff4652ad90c11_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant")  %>%
  st_transform('ESRI:102658')
bus <- rbind(
  bus[miamiBound,] %>%
    st_drop_geometry() %>%
    left_join(bus) %>%
    st_sf() %>%
    mutate(inMiami = "YES"),
  bus[miamiBound, op = st_disjoint] %>%
    st_drop_geometry() %>%
    left_join(bus) %>%
    st_sf() %>%
    mutate(inMiami = "NO")) %>%
  filter(inMiami == "YES")
#metro mover
metromover <- st_read("https://opendata.arcgis.com/datasets/aec76104165c4e879b9b0203fa436dab_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')
metromover <- rbind(
  metromover[miamiBound,] %>%
    st_drop_geometry() %>%
    left_join(metromover) %>%
    st_sf() %>%
    mutate(inMiami = "YES"),
  metromover[miamiBound, op = st_disjoint] %>%
    st_drop_geometry() %>%
    left_join(metromover) %>%
    st_sf() %>%
    mutate(inMiami = "NO")) %>%
  filter(inMiami == "YES")
#metro rail
metrorail <- st_read("https://opendata.arcgis.com/datasets/ee3e2c45427e4c85b751d8ad57dd7b16_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')
metrorail <- rbind(
  metrorail[miamiBound,] %>%
    st_drop_geometry() %>%
    left_join(metrorail) %>%
    st_sf() %>%
    mutate(inMiami = "YES"),
  metrorail[miamiBound, op = st_disjoint] %>%
    st_drop_geometry() %>%
    left_join(metrorail) %>%
    st_sf() %>%
    mutate(inMiami = "NO")) %>%
  filter(inMiami == "YES")



#CULTURE SPOTS
culture <- st_read("https://opendata.arcgis.com/datasets/70c48f0eb067448c8a787cfa1c1c3bb9_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')
culture <- rbind(
  culture[miamiBound,] %>%
    st_drop_geometry() %>%
    left_join(culture) %>%
    st_sf() %>%
    mutate(inMiami = "YES"),
  culture[miamiBound, op = st_disjoint] %>%
    st_drop_geometry() %>%
    left_join(culture) %>%
    st_sf() %>%
    mutate(inMiami = "NO")) %>%
  filter(inMiami == "YES")



#COMMERCIAL PROPERTIES
#read, project
commercial <- st_read("https://opendata.arcgis.com/datasets/fb8303c577c24ea386a91be7329842be_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')
#filter
commercial <- rbind(
  commercial[miamiBound,] %>%
    st_drop_geometry() %>%
    left_join(commercial) %>%
    st_sf() %>%
    mutate(inMiami = "YES"),
  commercial[miamiBound, op = st_disjoint] %>%
    st_drop_geometry() %>%
    left_join(commercial) %>%
    st_sf() %>%
    mutate(inMiami = "NO")) %>%
  filter(inMiami == "YES")



#FLOOD RISK ZONES
floodRisk <- st_read("https://opendata.arcgis.com/datasets/ef3bdd041b2e424695eb4dfe965966c4_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')
#filter
floodRisk <-
  rbind(
    st_centroid(floodRisk)[miamiBound,] %>%
      st_drop_geometry() %>%
      left_join(floodRisk) %>%
      st_sf() %>%
      mutate(inMiami = "YES"),
    st_centroid(floodRisk)[miamiBound, op = st_disjoint] %>%
      st_drop_geometry() %>%
      left_join(floodRisk) %>%
      st_sf() %>%
      mutate(inMiami = "NO")) %>%
  filter(inMiami == "YES") %>%
  dplyr::select(-inMiami, -SHAPE_Length, -ELEV, -FID) %>%
  dplyr::rename(FloodZone = FZONE, FloodHazard = ZONESUBTY)


floodInsure <- st_read("https://opendata.arcgis.com/datasets/f589473ddada46e78d437aaf09205b04_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')
filter
floodInsure <-
  rbind(
    st_centroid(floodInsure)[miamiBound,] %>%
      st_drop_geometry() %>%
      left_join(floodInsure) %>%
      st_sf() %>%
      mutate(inMiami = "YES"),
    st_centroid(floodInsure)[miamiBound, op = st_disjoint] %>%
      st_drop_geometry() %>%
      left_join(floodInsure) %>%
      st_sf() %>%
      mutate(inMiami = "NO")) %>%
  filter(inMiami == "YES") %>%
  mutate(floodInsureType = PANELID)



#CONTAMINATED SITES
contaminated <- st_read("https://opendata.arcgis.com/datasets/43750f842b1e451aa0347a2ca34a61d7_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')
contaminated <-
  rbind(
    st_centroid(contaminated)[miamiBound,] %>%
      st_drop_geometry() %>%
      left_join(contaminated) %>%
      st_sf() %>%
      mutate(inMiami = "YES"),
    st_centroid(contaminated)[miamiBound, op = st_disjoint] %>%
      st_drop_geometry() %>%
      left_join(contaminated) %>%
      st_sf() %>%
      mutate(inMiami = "NO")) %>%
  filter(inMiami == "YES")
#CONTAMINATION BUFFER

contamBuffer <- contaminated %>%
  st_buffer(800) %>%
  st_union() %>%
  st_as_sf() %>%
  mutate(contam = 1)
houses$contaminated <- houses %>%
  st_join(contamBuffer) %>%
  mutate(contam = ifelse(is.na(contam), 0, 1)) %>%
  pull(contam)



#NEAREST NEIGHBOR (some are used for testing, to determine feature buffer distances)
st_c <- st_coordinates
houses <-
  houses %>%
  mutate(
    #commercial properties NN
    commNN1 = nn_function(st_c(st_centroid(houses)), st_c(st_centroid(commercial)), 1),
    commNN5 = nn_function(st_c(st_centroid(houses)), st_c(st_centroid(commercial)), 5),
    #metro mover stations
    metroMNN1 = nn_function(st_c(st_centroid(houses)), st_c(metromover), 1),
    metroMNN5 = nn_function(st_c(st_centroid(houses)), st_c(metromover), 5),
    #metro rail stations
    metroRNN1 = nn_function(st_c(st_centroid(houses)), st_c(metrorail), 1),
    metroRNN5 = nn_function(st_c(st_centroid(houses)), st_c(metrorail), 5),
    #food/drinks
    foodBevNN1 = nn_function(st_c(st_centroid(houses)), st_c(foodBev), 1),
    foodBevNN5 = nn_function(st_c(st_centroid(houses)), st_c(foodBev), 5)
  ) 


#COMMERCIAL BUFFER
commercial <- commercial %>%
  mutate(counter = 1) %>%
  dplyr::select(counter)
#count properties within each buffer
houses$commercialProperties <-
  st_buffer(houses, 846) %>%
  aggregate(commercial, ., sum) %>%
  st_drop_geometry() %>%
  mutate(counter = ifelse(is.na(counter), 0, counter)) %>%
  pull(counter)



#FOOD AND BEV BUFFER
foodBev <- foodBev %>%
  mutate(counter = 1) %>%
  dplyr::select(counter)
#count parks within each buffer
houses$foodEstablishments <-
  st_buffer(houses, 2774) %>%
  aggregate(foodBev, ., sum) %>%
  st_drop_geometry() %>%
  mutate(counter = ifelse(is.na(counter), 0, counter)) %>%
  pull(counter)



#CULTURE BUFFER
culture <- culture %>%
  mutate(counter = 1) %>%
  dplyr::select(counter)
#count culture within each buffer
houses$cultureSpots <-
  st_buffer(houses, 774) %>%
  aggregate(culture, ., sum) %>%
  st_drop_geometry() %>%
  mutate(counter = ifelse(is.na(counter), 0, counter)) %>%
  pull(counter)




#METRORAIL BUFFER
metrorail <- metrorail %>%
  mutate(counter = 1) %>%
  dplyr::select(counter)
#count stops within each buffer
houses$metrorailStops <-
  st_buffer(houses, 12076.7) %>%
  aggregate(metrorail, ., sum) %>%
  st_drop_geometry() %>%
  mutate(counter = ifelse(is.na(counter), 0, counter)) %>%
  pull(counter)



#METROMOVER BUFFER
metromover <- metromover %>%
  mutate(counter = 1) %>%
  dplyr::select(counter)
#count metroM stops within each buffer
houses$metromoverStops <-
  st_buffer(houses, 18845) %>%
  aggregate(metromover, ., sum) %>%
  st_drop_geometry() %>%
  mutate(counter = ifelse(is.na(counter), 0, counter)) %>%
  pull(counter)



#BUS BUFFER
bus <- bus %>%
  mutate(counter = 1) %>%
  dplyr::select(counter)
#count bus within each buffer
houses$busStops <-
  st_buffer(houses, 775) %>%
  aggregate(bus, ., sum) %>%
  st_drop_geometry() %>%
  mutate(counter = ifelse(is.na(counter), 0, counter)) %>%
  pull(counter)



#PARKS BUFFER + AREA CALCULATION (using 1600ft buffer distance because the mean NN1 = 1600)
#get centroids
parkCentroids <- parks %>%
  st_centroid(parks) %>%    #get centroids of park layer
  dplyr::select(counter)
#count parks within each buffer
houses$parkCount <-
  st_buffer(houses, 1600) %>%
  aggregate(parkCentroids, ., sum) %>%
  st_drop_geometry() %>%
  mutate(counter = ifelse(is.na(counter), 0, counter)) %>%
  pull(counter)
#make buffer for each house
parkBuffer <- st_buffer(houses, 1600) %>%
  dplyr::select(Property.Address) %>%
  st_as_sf()
#calculate area of park space in each buffer
bufferedParks <- st_intersection(parkBuffer, parks) %>%
  group_by(Property.Address) %>%
  summarise() %>%
  mutate(parkArea = units::drop_units(st_area(.))) %>%
  st_drop_geometry()
#add park area back to houses file
houses <-
  left_join(houses, bufferedParks)



#SCHOOL CATCHMENT CATEGORIES
houses <-
  st_join(houses, elementary) %>%
  rename(elemCatch = 'NAME')

houses <-
  st_join(houses, middle) %>%
  rename(middleCatch = 'NAME')

houses <-
  st_join(houses, high) %>%
  rename(highCatch = 'NAME')



#SCHOOL DISTRICT CATEGORIES
houses <-
  st_join(houses, schoolDist) %>%
  rename(schoolDist = ID)


#FLOOD INSURANCE CATEGORIES
floodInsure <- floodInsure %>%
  dplyr::select(floodInsureType)
houses <- houses %>%
  st_join(., floodInsure) %>%
  mutate(floodInsureType = ifelse(is.na(floodInsureType), "other", floodInsureType))



#ADD ACS DATA
houses <-
  st_join(houses, acs)




#HOUSE AGE
houses <-
  houses %>%
  mutate(age = ifelse(is.na(YearBuilt), 0, (2020 - YearBuilt)))



#MAKE CATEGORICAL VARIABLES
houses <-
  houses %>%
  mutate(Bed.cat = case_when(
    Bed >= 0 & Bed < 3  ~ "Up to 2 Beds",
    Bed >= 3 & Bed < 4  ~ "3 Beds",
    Bed >= 4                    ~ "4+ Beds"))


houses <-
  houses %>%
  mutate(Bath.cat = case_when(
    Bath < 2  ~ "Up to 1 Bathroom",
    Bath == 2  ~ "2 Bathrooms",
    Bath >= 3                    ~ "3+ Bathrooms"))


houses <-
  houses %>%
  mutate(Stories.cat = case_when(
    Stories < 2  ~ "Up to 1 Stories",
    Stories == 2  ~ "2 Stories",
    Stories >= 3                    ~ "3+ Stories"))



houses <-
  houses %>%
  mutate(distWater.cat = case_when(
    distWater < 0.25  ~ "Less than 1/4 Mile",
    distWater >= 0.25   ~ "More than 1/4 Mile"))


houses <-
  houses %>%
  mutate(parkArea.cat = case_when(
    parkArea < 57000  ~ "Less than 57,000 SqFt",
    parkArea >= 57000 & parkArea < 320000  ~ "57,000 - 320,000 SqFt",
    parkArea >= 320000              ~ "More than 320,000 SqFt"))



#OTHER HOUSE FEATURES
houses <- houses %>%
  mutate(Pool = ifelse(str_detect(XF1, "Pool") | str_detect(XF2, "Pool") | str_detect(XF3, "Pool") | str_detect(XF1, "Whirlpool") | str_detect(XF2, "Whirlpool") | str_detect(XF3, "Whirlpool") | str_detect(XF1, "Jacuzzi") | str_detect(XF2, "Jacuzzi") | str_detect(XF3, "Jacuzzi"), "Pool", "No Pool"),
         Patio = ifelse(str_detect(XF1, "Patio") | str_detect(XF2, "Patio") | str_detect(XF3, "Patio"), "Patio", "No Patio"),
         Fence = ifelse(str_detect(XF1, "Fence") | str_detect(XF2, "Fence") | str_detect(XF3, "Fence"), "Fence", "No Fence"),
         Gazebo = ifelse(str_detect(XF1, "Gazebo") | str_detect(XF2, "Gazebo") | str_detect(XF3, "Gazebo"), "Gazebo", "No Gazebo"),
         Carport = ifelse(str_detect(XF1, "Carport") | str_detect(XF2, "Carport") | str_detect(XF3, "Carport"), "Carport", "No Carport"),
         Wall = ifelse(str_detect(XF1, "Wall") | str_detect(XF2, "Wall") | str_detect(XF3, "Wall"), "Wall", "No Wall"),
         Dock = ifelse(str_detect(XF1, "Dock") | str_detect(XF2, "Dock") | str_detect(XF3, "Dock"), "Dock", "No Dock"),
  )

#FIX NA VALUES
#zip
houses <-
  houses %>%
  mutate(Mailing.Zip = as.numeric(Mailing.Zip),
         Mailing.Zip = ifelse(is.na(Mailing.Zip), 0, Mailing.Zip))

#School dist
houses <-
  houses %>%
  mutate(elemCatch = ifelse(is.na(elemCatch), "other", elemCatch),
         middleCatch = ifelse(is.na(middleCatch), "other", middleCatch),
         highCatch = ifelse(is.na(highCatch), "other", highCatch),
  )



#park
houses <- houses %>%
  mutate(parkArea = ifelse(is.na(parkArea), 0, parkArea))

#acs - here I found the location of NA values and gave these houses the ACS values of houses nearby/in same neighborhood
houses["3487", "NAME"] = houses["1099", "NAME"]
houses["3487","timeToWork"] =  houses["1099","timeToWork"]
houses["3487", "medHHInc"] = houses["1099", "medHHInc"]
houses["3487", "monthhousingcost"] = houses["1099", "monthhousingcost"]
houses["3487", "pctVacant"] = houses["1099", "pctVacant"]
houses["3487", "pctWhite"] = houses["1099", "pctWhite"]
houses["3487", "pctRenterOcc"] = houses["1099", "pctRenterOcc"]


houses[c("579","580","581","1016","1372","1557","1853","2140","2557","2563","2571","2603","2786","2981","3050","3361"), "NAME"] = houses["3458", "NAME"]
houses[c("579","580","581","1016","1372","1557","1853","2140","2557","2563","2571","2603","2786","2981","3050","3361"), "timeToWork"] = houses["3458", "timeToWork"]
houses[c("579","580","581","1016","1372","1557","1853","2140","2557","2563","2571","2603","2786","2981","3050","3361"), "medHHInc"] = houses["3458", "medHHInc"]
houses[c("579","580","581","1016","1372","1557","1853","2140","2557","2563","2571","2603","2786","2981","3050","3361"), "monthhousingcost"] = houses["3458", "monthhousingcost"]
houses[c("579","580","581","1016","1372","1557","1853","2140","2557","2563","2571","2603","2786","2981","3050","3361"), "pctVacant"] = houses["3458", "pctVacant"]
houses[c("579","580","581","1016","1372","1557","1853","2140","2557","2563","2571","2603","2786","2981","3050","3361"), "pctWhite"] = houses["3458", "pctWhite"]
houses[c("579","580","581","1016","1372","1557","1853","2140","2557","2563","2571","2603","2786","2981","3050","3361"), "pctRenterOcc"] = houses["3458", "pctRenterOcc"]

houses["441", c("NAME","timeToWork","medHHInc","monthhousingcost","pctVacant","pctWhite","pctRenterOcc")] = houses["2090", c("NAME","timeToWork","medHHInc","monthhousingcost","pctVacant","pctWhite","pctRenterOcc")]

houses %>%
  dplyr::select(-Property.Zip)


#table of summary statistics with variable descriptions, sorted by category

housesSub <- houses %>%
  dplyr::select("AdjustedSqFt", "LotSize", "Bed", "Bath", "Stories", "commercialProperties", "age", "distWater", "foodEstablishments", "cultureSpots", "busStops", "parkArea", "timeToWork","monthhousingcost","pctVacant")

housesSub <- st_drop_geometry(housesSub)
stargazer(as.data.frame(housesSub), type="text", digits=1, title="Table 1: Descriptive Statistics for Miami Houses", out = "Miami Data.txt")
corrPlotVars <- houses %>%
  dplyr::select(-saleDate, -saleType, -saleYear, -Bldg, -Land, -Assessed, -WVDB, -HEX, -GPAR, -County.2nd.HEX, -County.Senior, -County.LongTermSenior, -County.Other.Exempt, -County.Taxable, -City.2nd.HEX, -City.Senior, -City.LongTermSenior, -City.Other.Exempt, -City.Taxable, -MillCode, -Owner1, -Owner2, -Mailing.State, -Mailing.Country, -Legal1, -Legal2, -Legal3, -Legal4, -Legal5, -Legal6, -YearBuilt, -EffectiveYearBuilt, -toPredict)

numericVars <-
  select_if(st_drop_geometry(corrPlotVars), is.numeric) %>% na.omit()

ggcorrplot(
  round(cor(numericVars), 1),
  p.mat = cor_pmat(numericVars),
  colors = c("#25CB10", "white", "#FA7800"),
  type="lower",
  insig = "blank") +  
  labs(title = "Figure 1: Sale Price Correlation with Numeric Variables") +
  plotTheme()
housesKnown <- houses %>%    
  filter(.,toPredict == 0)
housesUnknown <- houses %>%
  filter(.,toPredict ==1)

#1: distWater
ggplot(data = housesKnown, aes(x = distWater, y = SalePrice)) +
  geom_point(size=2, shape=20)  +
  labs(title = "Figure 2.1: Distance to Shoreline and Sale Price", subtitle = "Miami and Miami Beach, FL") +
  geom_smooth(method = "lm", se=F, colour = "green") +
  plotTheme()
#2: parkArea
ggplot(data = housesKnown, aes(x = parkArea, y = SalePrice)) +
  geom_point(size=2, shape=20) +
  labs(title = "Figure 2.2: Nearby Park Space and Sale Price", subtitle = "Miami and Miami Beach, FL") +
  geom_smooth(method = "lm", se=F, colour = "green") +
  plotTheme()
#3: vacant
ggplot(data = housesKnown, aes(x = pctVacant, y = SalePrice)) +
  geom_point(size=2, shape=20) +
  labs(title = "Figure 2.3: Share of Vacant Properties and Sale Price", subtitle = "Miami and Miami Beach, FL") +
  geom_smooth(method = "lm", se=F, colour = "green") +
  plotTheme()
#4: bus stops
ggplot(data = housesKnown, aes(x = busStops, y = SalePrice)) +
  geom_point(size=2, shape=20) +
  labs(title = "Figure 2.4: Bus Stops and Sale Price", subtitle = "Miami and Miami Beach, FL") +
  geom_smooth(method = "lm", se=F, colour = "green") +
  plotTheme()
#map of your dependent variable (sale price)
#water is just for mapping visuals
water <- st_read("https://opendata.arcgis.com/datasets/bf9de3192c9c4e458d1453f6d4c88d6c_0.geojson") %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658') %>%
  st_union() %>%
  st_intersection(.,miamiBound)

ggplot() +
  geom_sf(data = acs, fill = "gray80", colour = "white") +
  geom_sf(data = water, fill = "light blue", colour = "light blue") +
  geom_sf(data = housesKnown, aes(colour = q5(SalePrice))) +
  scale_colour_manual(values = paletteMap) +
  labs(title = "Figure 3: Home Sale Price", subtitle = "Miami and Miami Beach, FL") +
  mapTheme()
#1: distWater
ggplot() +
  geom_sf(data = acs, fill = "gray90", colour = "white") +
  geom_sf(data = water, fill = "light blue", colour = "light blue") +
  geom_sf(data = houses, aes(colour = q5(distWater))) +
  labs(title = "Figure 4.1: Distance to Shore", subtitle = "Miami and Miami Beach, FL") +
  scale_colour_manual(values = palette5) +
  mapTheme()
#2: parkArea
ggplot() +
  geom_sf(data = acs, fill = "gray90", colour = "white") +
  geom_sf(data = water, fill = "light blue", colour = "light blue") +
  geom_sf(data = housesKnown, aes(colour = q5(parkArea))) +
  scale_colour_manual(values = palette5) +
  labs(title = "Figure 4.2: Nearby Park Area", subtitle = "Miami, FL") +
  mapTheme()
#3: Food and bev
ggplot() +
  geom_sf(data = acs, fill = "gray90", colour = "white") +
  geom_sf(data = water, fill = "light blue", colour = "light blue") +
  geom_sf(data = housesKnown, aes(colour = q5(foodEstablishments))) +
  scale_colour_manual(values = palette5) +
  labs(title = "Figure 4.3: Number of Nearby Food/Bev Places", subtitle = "Miami, FL") +
  mapTheme()
#4: AdjustedSqFt
ggplot() +
  geom_sf(data = acs, fill = "gray90", colour = "white") +
  geom_sf(data = water, fill = "light blue", colour = "light blue") +
  geom_sf(data = houses, aes(colour = q5(AdjustedSqFt))) +
  labs(title = "Figure 4.4: Adjusted Square Feet", subtitle = "Miami and Miami Beach, FL") +
  scale_colour_manual(values = palette5) +
  mapTheme()
#5: MedHHInc
ggplot() +
  geom_sf(data = acs, fill = "gray90", colour = "white") +
  geom_sf(data = water, fill = "light blue", colour = "light blue") +
  geom_sf(data = houses, aes(colour = q5(medHHInc))) +
  labs(title = "Figure 4.5: Tract Median Household Income", subtitle = "Miami and Miami Beach, FL") +
  scale_colour_manual(values = palette5) +
  mapTheme()
#0.72 R2
reg2 <- lm(SalePrice ~ ., data = st_drop_geometry(housesKnown) %>%
             dplyr::select(SalePrice, ActualSqFt, LotSize, Zoning, Stories.cat, Bath.cat, Pool, medHHInc, Dock, Bed.cat, middleCatch, age, pctVacant, pctRenterOcc, monthhousingcost, Patio, foodEstablishments, timeToWork, metromoverStops, metrorailStops, parkArea, floodInsureType))
#read FL neighborhoods
#nhoods_fl <- aoi_boundary_HARV <- st_read("E:/Upenn/CPLN508/miami/zillow_nghbrhd_feb17/zillow_nghbrhd_feb17.shp")
nhoods_fl <- aoi_boundary_HARV <- st_read("/Users/annaduan/Documents/GitHub/Miami-Oct12/Raw\ Data/zillow_nghbrhd_feb17/zillow_nghbrhd_feb17.shp")
nhoods_mb <- subset(nhoods_fl, CITY == "MIAMI BEACH")%>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')
nhoods_m <- subset(nhoods_fl, CITY == "MIAMI")%>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')

#Join neighborhoods
nhoods <- rbind(nhoods_mb, nhoods_m)
nhoods <- nhoods %>%
  dplyr::select(NAME) %>%
  rename(neighborhood = NAME)
housesKnown <- housesKnown %>% st_join(., nhoods, join = st_within)
houses <- houses %>% st_join(., nhoods, join = st_within)
housesUnknown <- housesUnknown %>% st_join(., nhoods, join = st_within)

#Separate test/train sets

inTrain <- createDataPartition(
  y = paste(housesKnown$Zoning, housesKnown$floodInsureType, housesKnown$neighborhood),
  p = .60, list = FALSE)
miami.training <- housesKnown[inTrain,]
miami.test <- housesKnown[-inTrain,]  


#Training regression 
reg.training <- 
  lm(SalePrice ~ ., data = st_drop_geometry(miami.training) %>% 
       dplyr::select(SalePrice, ActualSqFt, LotSize, Zoning, Stories.cat, Bath.cat, Pool, medHHInc, Dock, Bed.cat, middleCatch, age, pctVacant, pctRenterOcc, monthhousingcost, Patio, foodEstablishments, timeToWork, metromoverStops, metrorailStops, parkArea, floodInsureType))
#polished table of  (training set) lm summary results (coefficients, R2 etc)
stargazer(reg.training, type="text", digits=1, title="Table 2: LM of Training Data", out = "Training LM.txt")
#Test regression on miami.test
miami.test <-
  miami.test %>%
  mutate(Regression = "Baseline Regression",
         SalePrice.Predict = predict(reg.training, miami.test), #751571.5
         SalePrice.Error = SalePrice.Predict - SalePrice, #60608.35
         SalePrice.AbsError = abs(SalePrice.Predict - SalePrice), #363363
         SalePrice.APE = SalePrice.AbsError / SalePrice) %>% #0.05589877  #corrected 
  filter(SalePrice < 5000000) 

#Mean error and APE 432255.1   0.3407089   443611.7 0.1643853 467809 -0.0241117
mean(miami.test$SalePrice.AbsError, na.rm = T)#[1] 351280.6
mean(miami.test$SalePrice.APE, na.rm = T)#[1] 0.7101703
mean(miami.test$SalePrice.Predict, na.rm = T)#[1] 770165.4

ggplot(data = miami.test) +
  geom_point(aes(x = SalePrice, y = SalePrice.AbsError)) +
  labs(title = "Figure 5.1: Observed Sale Price and Absolute Error") +
  plotTheme()
ggplot(data = miami.test) +
  geom_point(aes(x = SalePrice, y = SalePrice.APE)) +
  labs(title = "Figure 5.2: Observed Sale Price and Absolute Percent Error") +
  plotTheme()
#K Folds test
fitControl <- trainControl(method = "cv", number = 100)
set.seed(825)

reg.cv <- 
  train(SalePrice ~ ., data = st_drop_geometry(housesKnown) %>% 
          dplyr::select(SalePrice, ActualSqFt, LotSize, Zoning, Stories.cat, Bath.cat, Pool, medHHInc, Dock, Bed.cat, middleCatch, age, pctVacant, pctRenterOcc, monthhousingcost, Patio, foodEstablishments, timeToWork, metromoverStops, metrorailStops, parkArea, floodInsureType), 
        method = "lm", trControl = fitControl, na.action = na.pass)


#k-fold function online
kfold.MLR = function(fit,k=10,data=fit$model) {    
  sum.sqerr = rep(0,k)
  sum.abserr = rep(0,k)
  sum.pererr = rep(0,k)
  y = fit$model[,1]
  x = fit$model[,-1]
  n = nrow(data)
  folds = sample(1:k,nrow(data),replace=T)
  for (i in 1:k) {
    fit2 <- lm(formula(fit),data=data[folds!=i,])
    ypred = predict(fit2,newdata=data[folds==i,])
    sum.sqerr[i] = sum((y[folds==i]-ypred)^2)
    sum.abserr[i] = sum(abs(y[folds==i]-ypred))
    sum.pererr[i] = sum(abs(y[folds==i]-ypred)/y[folds==i])
  }
  cv = return(data.frame(RMSEP=sqrt(sum(sum.sqerr)/n),
                         MAE=sum(sum.abserr)/n,
                         MAPE=(sum(sum.pererr)/n)*100))
}

#ADDED TODAY
reg.cv.known <-
  housesKnown %>%
  mutate(Regression = "Baseline Regression",
         SalePrice.Predict = predict(reg.cv, housesKnown), #751571.5
         SalePrice.Error = SalePrice.Predict - SalePrice, #60608.35
         SalePrice.AbsError = abs(SalePrice.Predict - SalePrice), #363363
         SalePrice.APE = SalePrice.AbsError / SalePrice) %>% #0.05589877  #corrected 
  filter(SalePrice < 5000000) 
fold75 <- reg.cv$control$indexOut$Resample075

reg75 <- reg.cv.known[fold75,c("SalePrice", "SalePrice.Predict")]
reg75.test <-
  reg75 %>%
  mutate(SalePrice.Error = SalePrice.Predict - SalePrice, 
         SalePrice.AbsError = abs(SalePrice.Predict - SalePrice), 
         SalePrice.APE = SalePrice.AbsError / SalePrice) %>% 
  filter(SalePrice < 5000000) 
reg.cv.rs.min <- reg.cv$resample[75,]
reg.cv.rs.min$MAPE <- mean(reg75.test$SalePrice.APE)

round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}
reg.cv.rs.min <- round_df(reg.cv.rs.min, 2)

library(kableExtra)


reg.cv.rs.min <- reg.cv$resample[75,]


reg.cv.rs.min %>%                     #AD: knitting stops here
  gather(Variable, Value) %>%
  #filter(Variable == "MAE" | Variable == "RMSE") %>%
  group_by(Variable) %>%
  spread(Variable, Value) %>%
  kable(caption = "Table 3: Regression Results of One Test Set") %>%
  kable_classic(full_width = F, html_font = "Cambria")
ggplot(reg.cv$resample, aes(x=MAE)) +
  geom_histogram() +
  labs(title = "Figure 6.1: Mean Average Error in Cross Validation Tests") +
  plotTheme()
library(knitr)
library(kableExtra)
library(scales)

coords <- st_coordinates(housesKnown)
neighborList <- knn2nb(knearneigh(coords, 5)) #5 nearest neighborhoods

spatialWeights <- nb2listw(neighborList, style="W") #not sure what is W here
housesKnown$lagPrice <- lag.listw(spatialWeights, housesKnown$SalePrice)
#plot(housesKnown$SalePrice, housesKnown$lagPrice)

coords.test <-  st_coordinates(miami.test)
neighborList.test <- knn2nb(knearneigh(coords.test, 5))
spatialWeights.test <- nb2listw(neighborList.test, style="W")
#10.1 Map of test set residuals
library(modelr)

miami.test$resid <- 
  miami.test %>%
  as_data_frame() %>%
  add_residuals(., reg.training, var = "resid") %>%
  dplyr::select(resid, Folio) %>%
  pull(resid)

ggplot() +
  geom_sf(data = acs, fill = "gray90", colour = "white") +
  geom_sf(data = water, fill = "light blue", colour = "light blue") +
  geom_sf(data = miami.test, aes(colour = q5(resid))) +
  scale_colour_manual(values = palette5) +
  labs(title = "Figure 7.1: Test Set Residual Errors", subtitle = "Miami and Miami Beach, FL") +
  mapTheme()
miami.test %>%                
  mutate(lagPriceError = lag.listw(spatialWeights.test, SalePrice.Error)) %>%  
  ggplot(aes(lagPriceError, SalePrice)) +
  geom_point() +
  stat_smooth(aes(lagPriceError, SalePrice), 
              method = "lm", se = FALSE, size = 1, colour="#FA7800")+
  labs(title = "Figure 8.1: Spatial Lag of Price Errors") +
  plotTheme() + theme(plot.title = element_text(size = 18, colour = "black")) 
moranTest <- moran.mc(miami.test$SalePrice.Error,
                      spatialWeights.test, nsim = 999)

ggplot(as.data.frame(moranTest$res[c(1:999)]), aes(moranTest$res[c(1:999)])) +
  geom_histogram(binwidth = 0.01) +
  geom_vline(aes(xintercept = moranTest$statistic), colour = "#FA7800",size=1) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(title="Figure 9.1: Observed and Permuted Moran's I",
       subtitle= "Observed Moran's I in orange",
       x="Moran's I",
       y="Count") +
  plotTheme()
#Make new neighborhood regression
reg.nhood <- lm(SalePrice ~ ., data = as.data.frame(miami.training) %>% 
                  dplyr::select(neighborhood, SalePrice, ActualSqFt, LotSize, Zoning, Stories.cat, Bath.cat, Pool, medHHInc, Dock, Bed.cat, middleCatch, age, pctVacant, pctRenterOcc, monthhousingcost, Patio, foodEstablishments, timeToWork, metromoverStops, metrorailStops, parkArea, floodInsureType))

#Outcomes
miami.test.nhood <-
  miami.test %>%
  mutate(Regression = "Neighborhood Effects",
         SalePrice.Predict = predict(reg.nhood, miami.test), #613237.3
         SalePrice.Error = SalePrice - SalePrice.Predict, # -108442.4      
         SalePrice.AbsError = abs(SalePrice - SalePrice.Predict), # 491238.7
         SalePrice.APE = (abs(SalePrice - SalePrice.Predict)) / SalePrice)%>% #0.7109973
  filter(SalePrice < 5000000)

#Check accuracy
bothRegressions <-
  rbind(
    dplyr::select(miami.test, starts_with("SalePrice"), Regression, neighborhood) %>%
      mutate(lagPriceError = lag.listw(spatialWeights.test, SalePrice.Error)),
    dplyr::select(miami.test.nhood, starts_with("SalePrice"), Regression, neighborhood) %>%
      mutate(lagPriceError = lag.listw(spatialWeights.test, SalePrice.Error)))   

st_drop_geometry(bothRegressions) %>%
  gather(Variable, Value, -Regression, -neighborhood) %>%
  filter(Variable == "SalePrice.AbsError" | Variable == "SalePrice.APE") %>%
  group_by(Regression, Variable) %>%
  summarize(meanValue = mean(Value, na.rm = T)) %>%
  spread(Variable, meanValue) %>%
  kable(caption = "Table 4: Neighborhood Effect on Error")
bothRegressions %>%
  dplyr::select(SalePrice.Predict, SalePrice, Regression) %>%
  ggplot(aes(SalePrice, SalePrice.Predict)) +
  geom_point() +
  stat_smooth(aes(SalePrice, SalePrice),
              method = "lm", se = FALSE, size = 1, colour="#FA7800") +
  stat_smooth(aes(SalePrice.Predict, SalePrice),
              method = "lm", se = FALSE, size = 1, colour="#25CB10") +
  facet_wrap(~Regression) +
  labs(title="Figure 10.1: Predicted Sale Price and Observed Price",
       subtitle="Orange line represents a perfect prediction; Green line represents prediction") +
  plotTheme() + theme(plot.title = element_text(size = 18, colour = "black"))
#876 rows
#filter by toPredict = 1
housesPredictions <-                   
  houses %>%
  mutate(prediction = predict(reg.nhood, houses),
         team_name = 'Panda')

predictions <- housesPredictions[,c("Folio", "prediction", "team_name", "toPredict")] %>%
  st_drop_geometry() %>%
  filter(toPredict == 1) %>%
  dplyr::select(-toPredict)

# write.csv(predictions, "PANDA.csv") #"The column names MUST to be "prediction", "Folio", and "team_name"" - piazza

#Map values
ggplot() +
  geom_sf(data = acs, fill = "gray90", colour = "white") +
  geom_sf(data = water, fill = "light blue", colour = "light blue") +
  geom_sf(data = housesPredictions, aes(colour = q5(prediction))) +
  scale_colour_manual(values = palette5) +
  labs(title = "Figure 11.1: Predicted Sale Price Values", subtitle = "Miami and Miami Beach, FL") +
  # facet_wrap(~toPredict) +
  mapTheme()
#Using the test set predictions, provide a map of mean absolute percentage error (MAPE) by neighborhood
names(bothRegressions)[names(bothRegressions) == "neighborhood"] <- "neighborhood"
st_drop_geometry(bothRegressions) %>%
  group_by(Regression, neighborhood) %>%
  summarise(mean.MAPE = mean(SalePrice.APE, na.rm = T)) %>%
  ungroup() %>%
  left_join(nhoods) %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(data = water, fill = "light blue", colour = "light blue") +
  geom_sf(colour = "gray", aes(fill = q5(mean.MAPE))) +
  scale_fill_manual(values = paletteMap) +
  labs(title = "Figure 12.1: Mean Average Percentage Error by Neighborhood") +
  mapTheme()
scatter_hood <-
  miami.test.nhood %>%
  group_by(neighborhood) %>%
  dplyr::select(neighborhood, SalePrice.APE, SalePrice.Predict)

mean_sca_hd <-
  scatter_hood %>%
  group_by(neighborhood) %>%
  summarise_at(vars("SalePrice.APE", "SalePrice.Predict"), mean)

plot(mean_sca_hd$SalePrice.Predict, mean_sca_hd$SalePrice.APE, main="Figure 13.1: MAPE by Neighborhood and Mean Price by Neighborhood", xlab="Mean Price by Neighborhood", ylab="MAPE by neighborhood") +
  plotTheme()
#https://r-graphics.org/recipe-scatter-labels or we can use the method below:

#scatter_mae_mean <- ggplot(mean_sca_hd, aes(x = SalePrice.Predict, y = SalePrice.APE)) +
#    geom_point() +
#  plotTheme()

#scatter_mae_mean +    #AD: this is cool but we can't really see which dot is which neighborhood
#  annotate("text", x = 820000, y = -11.7, label = "UPPER EASTSIDE") +
#  annotate("text", x = 330000, y = -6.5, label = "LITTLE HAITI")+
#  annotate("text", x = 1300000.82, y = 6.85325579, label = "NAUTILUS")+
#  annotate("text", x = 631682.30, y = 3.57218432, label = "BISCAYNE POINT")+
#  annotate("text", x = 1557855.78, y = 2.81996319, label = "LA GORCE")+
#  annotate("text", x = 2031665.51, y = 1.95167248, label = "WYNWOOD - EDGEWATER")+
#  annotate("text", x = 890031.19, y = -2.16394539, label = "OVERTOWN")
#RACE & INCOME
#make new layer
acsRaceIncome <-
  acs %>%
  mutate(raceContext= ifelse(pctWhite > .5, "Majority White", "Majority Non-White"),
         incomeContext = ifelse(medHHInc > 49256.56, "High Income", "Low Income"))
#context
grid.arrange(ncol = 2,
             ggplot() + geom_sf(data = na.omit(acsRaceIncome), aes(fill = raceContext)) +
               scale_fill_manual(values = c("#25CB10", "#FA7800"), name="Race Context") +
               labs(title = "Figure 14.1: Race Context of Miami and Miami Beach") +
               mapTheme() + theme(legend.position="bottom"),
             ggplot() + geom_sf(data = na.omit(acsRaceIncome), aes(fill = incomeContext)) +
               scale_fill_manual(values = c("#25CB10", "#FA7800"), name="Income Context") +
               labs(title = "Income Context of Miami and Miami Beach") +
               mapTheme() + theme(legend.position="bottom"))
#tables
st_join(bothRegressions, acsRaceIncome) %>%
  filter(!is.na(raceContext)) %>%
  group_by(Regression, raceContext) %>%
  summarize(mean.MAPE = scales::percent(mean(SalePrice.APE, na.rm = T))) %>%
  st_drop_geometry() %>%
  spread(raceContext, mean.MAPE) %>%
  kable(caption = "Table 5: MAPE by Neighborhood Racial Context")
st_join(bothRegressions, acsRaceIncome) %>%
  filter(!is.na(incomeContext)) %>%
  group_by(Regression, incomeContext) %>%
  summarize(mean.MAPE = scales::percent(mean(SalePrice.APE, na.rm = T))) %>%
  st_drop_geometry() %>%
  spread(incomeContext, mean.MAPE) %>%
  kable(caption = "Table 6: MAPE by Neighborhood Income Context")