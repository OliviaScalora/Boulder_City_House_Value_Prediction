
#Setup Install Libraries
install.packages('confidence')
library(tidycensus)
library(tigris)
library(kableExtra)
library(ggplot2)
library(tidyverse)
library(sf)
library(spdep)
library(caret)
library(ckanr)
library(FNN)
library(grid)
library(gridExtra)
library(ggcorrplot)
library(ggmap)
library(stargazer)
library(geosphere)
library(scales)
library(broom)
library(hrbrthemes)
library(confidence)


root.dir = "https://raw.githubusercontent.com/urbanSpatial/Public-Policy-Analytics-Landing/master/DATA/"

source("https://raw.githubusercontent.com/urbanSpatial/Public-Policy-Analytics-Landing/master/functions.r")

options(scipen=999)
options(tigris_class = "sf")
options(tigris_use_cache = TRUE)

qBr <- function(df, variable, rnd) {
  if (missing(rnd)) 
    as.character(quantile(round(df[[variable]],0),
                          c(.01,.2,.4,.6,.8), na.rm=T))
  else if (rnd == FALSE | rnd == F) 
    as.character(formatC(quantile(df[[variable]],
                                  c(.01,.2,.4,.6,.8), na.rm=T), digits = 3))}


q5 <- function(variable) {
  as.factor(ntile(variable, 5))}

mapTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 16,colour = "black"),
    plot.subtitle=element_text(face="italic"),
    plot.caption=element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),axis.title = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    strip.text.x = element_text(size = 14)
  )
}


plotTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 16,colour = "black"),
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

paletteMap <- c("#003A6B","#5880A2","#83A3BE","#AFC6D9","#DBE9F5")
palette5 <- c("#25CB10", "#5AB60C", "#8FA108",   "#C48C04", "#FA7800")

#Load Student Data
StudentData <- st_read("studentData.geojson",crs = 'ESRI:102254')%>%
  st_transform('EPSG:26913')

#CENSUS
census_api_key("3909c4b14d693e1a81c3e86c9cfc3b22fc112324", overwrite = TRUE, install = TRUE)

#Call Census Variables
#Median Income - B19013_001 
#Median Rent - B25058_001 
#Total Below 100 percent poverty level - B06012_002 
#Total Married Couple Households - B11012_002
#Total Unmarried Female Householder - B11012_008
#Total Unmarried Male Householder - B11012_013
#Total With a Bachelors Degree - B15003_022
#Aggregate Travel Time to Work per tracts -B08013_001
#Total Population - B25026_001 
#Total Race - White - B02001_002, Black - B02001_003

tracts19 <- { 
  get_acs(geography = "tract", variables = c("B19013_001E", "B25058_001E", "B06012_002E",
                                             "B15003_022E", "B08013_001E", "B25026_001E",
                                             "B02001_002E", "B02001_003E", "B11012_002E",
                                             "B11012_008E", "B11012_013E"), 
          state='CO', county='Boulder County', geometry=T, output = "wide") %>% 
    st_transform('EPSG:26913')%>%
    dplyr::select( -NAME,-B19013_001M,-B25058_001M,-B06012_002M,
                   -B15003_022M,-B08013_001M,-B25026_001M,
                   -B02001_002M, -B02001_003M, -B11012_002M,
                   -B11012_008M, -B11012_013M)%>%
    rename(MedHHInc = B19013_001E, 
           MedRent = B25058_001E,
           BelPov100 = B06012_002E,
           TotalBach = B15003_022E,
           TravelTime = B08013_001E,
           TotalPop = B25026_001E, 
           TotalWhite = B02001_002E)%>%
    mutate(PctWhite = round((TotalWhite /TotalPop),2),
           PctMarried = round(B11012_002E/(B11012_008E + B11012_013E),2))%>%
    dplyr::select(-B11012_008E,-B11012_013E, -B02001_003E, -B11012_002E)}

#TRAILHEAD DATA
Trailhead <- st_read("Boulder_Area_Trailheads.geojson",crs = 'EPSG:4326')%>%
  st_transform('EPSG:26913')%>%
  dplyr::select(OBJECTID, THNAME, geometry)
#Create column which counts how many trailheads are within 1/2 mi from each home sale
StudentData <- StudentData%>%
  mutate(NumTrailHeads =(st_buffer(StudentData, 2640) %>% 
                        aggregate(mutate((Trailhead%>%select(geometry)), counter = 1),., sum)%>%
                          pull(counter)))
StudentData$NumTrailHeads[is.na(StudentData$NumTrailHeads)] <- 0

#WILDFIRE HISTORY 
Wildfire_History <-  st_read("Wildfire_History.geojson",crs = 'EPSG:4326')%>%
  st_transform('EPSG:26913')%>%
  select(OBJECTID, NAME, geometry)

#create 1/4 mile buffer intervals
 wf_buffer_25 <- 
  st_union(st_buffer(Wildfire_History, 1320)) %>% st_sf() %>%
  mutate(Legend = "Wildfire Buffer")
 wf_buffer_5 <- 
   st_union(st_buffer(Wildfire_History, 2640)) %>% st_sf() %>%
   mutate(Legend = "Wildfire Buffer")
 wf_buffer_75 <- 
   st_union(st_buffer(Wildfire_History, 3690)) %>% st_sf() %>%
   mutate(Legend = "Wildfire Buffer")
 wf_buffer_1 <- 
   st_union(st_buffer(Wildfire_History, 5280)) %>% st_sf() %>%
   mutate(Legend = "Wildfire Buffer")
#create wildfire_risk column 1 if within 1/2 mi. 0 if not
#StudentData <- {
#  StudentData[wildfire_buffer,] %>%
#    st_drop_geometry() %>%
#    mutate(wildfire_risk = 1,
#           wildfire_risk = replace_na(wildfire_risk,0))%>%
#    full_join(StudentData)%>%
#    st_sf()}
#StudentData$wildfire_risk[is.na(StudentData$wildfire_risk)] <- 0

#GEOHAZARD RISK
GeoHazard <-  st_read("Geological_Hazard_Corridor.geojson", crs = 'EPSG:4326')%>%
  st_transform('EPSG:26913')%>%
  select(OBJECTID, geometry)

#Identify home sales that are within a quarter mile of 
#geological hazard reports 
#geohazard_buffer <-   
#  st_union(st_buffer(GeoHazard, 1320)) %>%
#  st_sf() %>%
#  mutate(Legend = "GeoHazard Buffer")

geo_buffer_25 <- 
  st_union(st_buffer(GeoHazard, 1320)) %>% st_sf() %>%
  mutate(Legend = "GeoHazardBuffer")
geo_buffer_5 <- 
  st_union(st_buffer(GeoHazard, 2640)) %>% st_sf() %>%
  mutate(Legend = "GeoHazardBuffer")
geo_buffer_75 <- 
  st_union(st_buffer(GeoHazard, 3690)) %>% st_sf() %>%
  mutate(Legend = "GeoHazard Buffer")
geo_buffer_1 <- 
  st_union(st_buffer(GeoHazard, 5280)) %>% st_sf() %>%
  mutate(Legend = "GeoHazardBuffer")
#create geohazard_risk column 1 if within 1/4 mi. 0 if not
#StudentData <- {
#  StudentData[geohazard_buffer,] %>%
#    st_drop_geometry() %>%
#    mutate(geohazard_risk = 1,
#           geohazard_risk = replace_na(geohazard_risk,0))%>%
#    full_join(StudentData)%>%
#    st_sf()}

#STREAMS 
#Streams <- st_read("onlynaturalstreams.geojson",crs = 'EPSG:4326')%>%
#  st_transform('EPSG:26913')%>% drop_na()

#Identify home sales that are within a quarter of streams
#Streams_buffer <-   
#  st_union(st_buffer(Streams, 1320)) %>%
#  st_sf() %>%
#  mutate(Legend = "Streams Buffer")

#create stream_proximity column 1 if within 1/4 mi. 0 if not
#StudentData <- {
#  StudentData[Streams_buffer,] %>%
#    st_drop_geometry() %>%
#    mutate(stream_proximity = 1,
#    stream_proximity = replace_na(stream_proximity,0))%>%
#    full_join(StudentData)%>%
#    st_sf()}

#LAKES AND RESERVOIRS 
Lakes_Reservoir <- st_read("Lakes_and_Reservoirs.geojson",crs = 'EPSG:4326')%>%
  st_transform('EPSG:26913')%>%
  dplyr::select(FeatureType, Acreage, geometry)%>%
  filter(Acreage>=1)

#identify home sales within quarter mile of a lake or reservoir that is 
#larger than 1 acre
lake_buffer_25 <- 
  st_union(st_buffer(Lakes_Reservoir, 1320)) %>% st_sf() %>%
  mutate(Legend = "LakeBuffer")
lake_buffer_5 <- 
  st_union(st_buffer(Lakes_Reservoir, 2640)) %>% st_sf() %>%
  mutate(Legend = "LakeBuffer")
lake_buffer_75 <- 
  st_union(st_buffer(Lakes_Reservoir, 3690)) %>% st_sf() %>%
  mutate(Legend = "LakeBuffer")
lake_buffer_1 <- 
  st_union(st_buffer(Lakes_Reservoir, 5280)) %>% st_sf() %>%
  mutate(Legend = "LakeBuffer")
#create lake proximity column
#StudentData <- {
#  StudentData[lake_buffer,] %>%
#    st_drop_geometry() %>%
#    mutate(lake_proximity = 1,
#           lake_proximity = replace_na(lake_proximity,0))%>%
#    full_join(StudentData)%>%
#    st_sf()}

#FLOOD ZONE
FloodZone <- st_read("Flood_Areas.geojson",crs = 'EPSG:4326')%>%
  st_transform('EPSG:26913')%>%
  dplyr::select(OBJECTID, geometry)


#identify home sales within quarter mile of a flood event
#FloodZone_buffer <-
#  st_union(st_buffer(FloodZone, 1320)) %>%
#  st_sf() %>%
#  mutate(Legend = "Flood Buffer")
fl_buffer_25 <- 
  st_union(st_buffer(FloodZone, 1320)) %>% st_sf() %>%
  mutate(Legend = "FloodZone")
fl_buffer_5 <- 
  st_union(st_buffer(FloodZone, 2640)) %>% st_sf() %>%
  mutate(Legend = "FloodZone")
fl_buffer_75 <- 
  st_union(st_buffer(FloodZone, 3690)) %>% st_sf() %>%
  mutate(Legend = "FloodZone")
fl_buffer_1 <- 
  st_union(st_buffer(FloodZone, 5280)) %>% st_sf() %>%
  mutate(Legend = "FloodZone")

#create flood proximity column
#StudentData <- {
#  StudentData[FloodZone_buffer,] %>%
#    st_drop_geometry() %>%
#    mutate(flood_risk = 1,
#           flood_risk = replace_na(flood_risk,0))%>%
#    full_join(StudentData)%>%
#    st_sf()}

# Business + Commercial 
Business <- st_read("Zoning_Districts.geojson",crs = 'EPSG:4326')%>%
  st_transform('EPSG:26913')%>%
  filter(ZONEDESC =='Business')%>%
  select(geometry, ZONEDESC)
Commercial <- st_read("Zoning_Districts.geojson",crs = 'EPSG:4326')%>%
  st_transform('EPSG:26913')%>%
  filter(ZONEDESC =='Commercial')%>%
  select(geometry, ZONEDESC)
BusinessCommercial <- rbind(Business, Commercial)

BC_buffer_25 <- 
  st_union(st_buffer(BusinessCommercial, 1320)) %>% st_sf() %>%+
  mutate(Legend = "Business Buffer")
BC_buffer_5 <- 
  st_union(st_buffer(BusinessCommercial, 2640)) %>% st_sf() %>%
  mutate(Legend = "Business Buffer")
BC_buffer_75 <- 
  st_union(st_buffer(BusinessCommercial, 3690)) %>% st_sf() %>%
  mutate(Legend = "Business Buffer")
BC_buffer_1 <- 
  st_union(st_buffer(BusinessCommercial, 5280)) %>% st_sf() %>%
  mutate(Legend = "Business Buffer")

#general Industrial
Industrial1 <- st_read("Zoning_Districts.geojson",crs = 'EPSG:4326')%>%
  st_transform('EPSG:26913')%>%
  filter(ZONEDESC =='General Industrial')%>%
  select(geometry, ZONEDESC)
Industrial2 <- st_read("Zoning_Districts.geojson",crs = 'EPSG:4326')%>%
  st_transform('EPSG:26913')%>%
  filter(ZONEDESC =='Light Industrial')%>%
  select(geometry, ZONEDESC)
Industrial <- rbind(Industrial1, Industrial2)

I_buffer_25 <- 
  st_union(st_buffer(Industrial, 1320)) %>% st_sf() %>%
  mutate(Legend = "Industrial Buffer")
I_buffer_5 <- 
  st_union(st_buffer(Industrial, 2640)) %>% st_sf() %>%
  mutate(Legend = "Industrial Buffer")
I_buffer_75 <- 
  st_union(st_buffer(Industrial, 3690)) %>% st_sf() %>%
  mutate(Legend = "Industrial Buffer")
I_buffer_1 <- 
  st_union(st_buffer(Industrial, 5280)) %>% st_sf() %>%
  mutate(Legend = "Industrial Buffer")

#Agriculture
Ag <-  st_read("Zoning_Districts.geojson",crs = 'EPSG:4326')%>%
  st_transform('EPSG:26913')%>%
  filter(ZONEDESC =='Agricultural')%>%
  select(geometry, ZONEDESC)

AG_buffer_25 <- 
  st_union(st_buffer(Ag, 1320)) %>% st_sf() %>%
  mutate(Legend = "AG Buffer")
AG_buffer_5 <- 
  st_union(st_buffer(Ag, 2640)) %>% st_sf() %>%
  mutate(Legend = "AG Buffer")
AG_buffer_75 <- 
  st_union(st_buffer(Ag, 3690)) %>% st_sf() %>%
  mutate(Legend = "AG Buffer")
AG_buffer_1 <- 
  st_union(st_buffer(Ag, 5280)) %>% st_sf() %>%
  mutate(Legend = "AG Buffer")

Zoning <- st_read("Zoning_Districts.geojson",crs = 'EPSG:4326')%>%
  st_transform('EPSG:26913')

st_c <- st_coordinates
StudentData <-
  StudentData %>%
  mutate(
    wildfire.nn = nn_function(st_c(st_centroid(StudentData)),st_c(st_centroid(Wildfire_History)), 5),
    geohazard.nn = nn_function(st_c(st_centroid(StudentData)),st_c(st_centroid(GeoHazard)), 5),
    lakes.nn = nn_function(st_c(st_centroid(StudentData)),st_c(st_centroid(Lakes_Reservoir)),5),
    floodzone.nn = nn_function(st_c(st_centroid(StudentData)),st_c(st_centroid(FloodZone)),5),
    Ag.nn = nn_function(st_c(st_centroid(StudentData)),st_c(st_centroid(Ag)),5),
    Industrial.nn = nn_function(st_c(st_centroid(StudentData)),st_c(st_centroid(Industrial)),5), 
    Business.nn = nn_function(st_c(st_centroid(StudentData)),st_c(st_centroid(BusinessCommercial)),5))

#Join ACS data to StudentData
StudentData <-
  st_join(StudentData, tracts19)

#Air Condition Categories
StudentData$AcDscr[is.na(StudentData$Ac)] <- 'No AC'

#age 
StudentData <-
  StudentData %>%
  mutate(HouseAge = ifelse(is.na(builtYear), 0, (2021 - builtYear)))

#car storage space
StudentData$carStorageSF[is.na(StudentData$carStorageSF)] <- 0
StudentData <-
  StudentData %>%
  mutate(carStorage = case_when(
    carStorageSF == 0     ~ "No car Storage",
    carStorageSF> 0 & carStorageSF< 1000 ~ "Small",
    carStorageSF >= 2000 ~ "Large"))


#home sale has $3 million with low square footage - determine to be data entry mistake
#omit outlier data observation 
StudentData <-StudentData %>%
  filter(!MUSA_ID == '8735')

#Clean Up data set to include only variables of interest
HouseData<- StudentData %>%
  select(-year_quarter,-address,-bld_num,-section_num,-bldgClass,-bldgClassDscr,-ConstCode,-ConstCodeDscr,
         -builtYear,-CompCode,-EffectiveYear,-bsmtType,-bsmtTypeDscr,-carStorageType,-carStorageTypeDscr,-Ac,
         -ExtWallDscrSec,-ExtWallSec, -IntWall, -Roof_Cover,-Roof_CoverDscr,-Stories,-UnitCount,-status_cd)

HouseDataDescr <- HouseData %>%
  dplyr::select("designCodeDscr", "qualityCodeDscr", "bsmtSF", "carStorageSF", "nbrBedRoom", 
                "TotalFinishedSF","AcDscr", "HeatingDscr", "NumTrailHeads","wildfire.nn","geohazard.nn",
                "lakes.nn","floodzone.nn","Ag.nn","Industrial.nn","Business.nn", "MedHHInc", "MedRent",
                "BelPov100", "TotalBach","TravelTime","PctWhite","PctMarried","HouseAge")%>%
  st_drop_geometry()

stargazer(as.data.frame(HouseDataDescr), 
          type="text", 
          digits=1, 
          title="Descriptive Statistics for Boulder Homes", 
          caption = 'table 1', 
          out = "DescriptiveTable.txt")

#correlation  Plot
corrPlotVars <- HouseData %>%
  dplyr::select(-saleDate, -year, -designCode, -qualityCode,-HeatingDscr,
                -ExtWallDscrPrim, -IntWallDscr, -MUSA_ID, -toPredict,-GEOID,
                -TotalPop, -TotalWhite)

numericVars <-
  select_if(st_drop_geometry(corrPlotVars), is.numeric) %>% 
  na.omit()

ggcorrplot(
  round(cor(numericVars), 1),
  p.mat = cor_pmat(numericVars),
  colors = c("#25CB10", "white", "#FA7800"),
  type="lower",
  insig = "blank") +  
  labs(title = "Sale Price Correlation with Numeric Variables", caption = "figure #") +
  plotTheme()

correlation.matrix <- cor(numericVars)
stargazer(correlation.matrix, 
          type="text", 
          title="Correlation Matrix of Predictive Variables for Boulder Homes", 
          caption = 'table 2', 
          out = "correlation.matrix.txt")

#Plot Interpretation
  # many internal features ie; bsmt sf, car storage sf, total sf, and number of bed, bath etc. have
  # a positive correlation with house price as we would expect. As the sf or number of 
  # bed/bathrooms increase, so does the price of a home. the only numeric variables with a negative
  # correlation are travel time to work and percent married. We should look at these variables.
  # Longer travel times indicate homes in suburban census tracts. Should proximity to downtown indicate
  # higher or lower home price? is travel time to work a strong predictor of home price?
  # we should also investigate the below poverty line variable. Based on what we know, census tracts that 
  # have a higher percentage of the population living in poverty should impact home sales negatively. the
  # correlation plot should show a negative correlation. We will also look at total population with a bachelors
  # degree as well. We would expect higher education levels to have positive impact on the price of homes
  # in those census tracts. And we will look at House Age. Do newer homes sell for higher prices?
#look at histograms of our dependent and predictor variables 

hist(HouseData$price)
HouseData$LNprice = log(HouseData$price+1)
hist(HouseData$LNprice)

hist(HouseData$lakes.nn)
HouseData$LNlakes.nn = log(HouseData$lakes.nn+1)
hist(HouseData$LNlakes.nn)

hist(HouseData$TotalFinishedSF)
HouseData$LNTotalSF = log(HouseData$TotalFinishedSF+1)
hist(HouseData$LNTotalSF)

hist(HouseData$Ag.nn)
HouseData$LNAg = log(HouseData$Ag.nn+1)
hist(HouseData$LNAg)

hist(HouseData$Industrial.nn)
HouseData$LNIndustry = log(HouseData$Industrial.nn+1)
hist(HouseData$LNIndustry)

#separating data between test set and training set
knownSales <- HouseData %>%    
  filter(.,toPredict == 0)
toPredictSales <- HouseData %>%
  filter(.,toPredict ==1)

#SCATTER PLOTS
#remove data observation that is out liar
knownSales <- knownSales%>%
  filter(MUSA_ID != 8197)

#Pearson Correlation Tests and Linear Relationship Plots

#bsmtSF
cor.test(knownSales$bsmtSF, knownSales$price, method = "pearson")
#marginal pearsons r of .21 but statistically significant

bsmtSFplot <- ggplot(data = knownSales, aes(x = bsmtSF, y = price)) +
  geom_point(size=2, shape=20)  +
  labs(x = 'Basement Square Footage',
       y = 'House Price',
       title = "Basement Square Footage \nand Sale Price", 
       subtitle = "Boulder County, CO") +
  geom_smooth(method = "lm", se=F, colour = "orange") +
  theme_ipsum_rc()

#carStorageSF
cor.test(knownSales$carStorageSF, knownSales$price, method = "pearson")
#marginal correlation value of .25 but statistically significant

carSFplot <- ggplot(data = knownSales, aes(x = carStorageSF, y = price)) +
  geom_point(size=2, shape=20)  +
  labs(x = 'Car Storage Square Footage',
       y = 'House Price',
       title = "Car Storage Square Footage \nand Sale Price", 
       subtitle = "Boulder County, CO") +
  geom_smooth(method = "lm", se=F, colour = "orange") +
  theme_ipsum_rc()

#total finished square footage 
cor.test(knownSales$TotalFinishedSF, knownSales$price, method = "pearson")
#pearsons r of .49 and statistically significant 

TSFplot <- ggplot(data = knownSales, aes(x = TotalFinishedSF, y = price)) +
  geom_point(size=2, shape=20)  +
  labs(x = 'Total Finished Square Footage',
       y = 'House Price',
       title = "Total Finished Square Footage \nand Sale Price", 
       subtitle = "Boulder County, CO") +
  geom_smooth(method = "lm", se=F, colour = "green") +
  plotTheme()

#travel time and sale price correlation 
cor.test(knownSales$TravelTime, knownSales$price, method = "pearson")
#negative pearsons r - weak correlation but statistically significant

TTplot <- ggplot(data = knownSales, aes(x = TravelTime, y = price)) +
  geom_point(size=2, shape=20)  +
  labs(x = 'Aggregate Travel Time to Work',
       y = 'House Price',
       title = "Travel Time to Work \nand Sale Price", 
       subtitle = "Boulder County, CO") +
  geom_smooth(method = "lm", se=F, colour = "green") +
  plotTheme()
#correlation has a slightly negative relationship. 

#Median rent and home price
cor.test(knownSales$MedRent, knownSales$price, method = "pearson")
#weak correlation but statistically significant

RentPlot <- ggplot(data = knownSales, aes(x = MedRent, y = price)) +
  geom_point(size=2, shape=20)  +
  labs(x = 'Median Rent',
       y = 'House Price',
       title = "Median Rent and Sale Price", 
       subtitle = "Boulder County, CO") +
  geom_smooth(method = "lm", se=F, colour = "green") +
  plotTheme()

#household income and home price
cor.test(knownSales$MedHHInc, knownSales$price, method = "pearson")
#marginal pearsons r of .18 and statistically significant 

HHIncplot<- ggplot(data = knownSales, aes(x = MedHHInc, y = price)) +
  geom_point(size=2, shape=20)  +
  labs(x = 'Median Household Income',
       y = 'House Price',
       title = "Median Household Income \nand Sale Price", 
       subtitle = "Boulder County, CO") +
  geom_smooth(method = "lm", se=F, colour = "green") +
  plotTheme()

#percent married folks and home price
cor.test(knownSales$PctMarried, knownSales$price, method = "pearson")
#very weak correlation, almost none

MarriedPlot <- ggplot(data = knownSales, aes(x = PctMarried, y = price)) +
  geom_point(size=2, shape=20)  +
  labs(x = 'Percent Married',
       y = 'House Price',
       title = "Percent of Homeowners with \nMarried Status and Sale Price", 
       subtitle = "Boulder County, CO") +
  geom_smooth(method = "lm", se=F, colour = "green") +
  plotTheme()

#plot percent white and home price
cor.test(knownSales$PctWhite, knownSales$price, method = "pearson")
#very weak correlation

Whiteplot<-  ggplot(data = knownSales, aes(x = PctWhite, y = price)) +
  geom_point(size=2, shape=20)  +
  labs(x = 'Percent White',
       y = 'House Price',
       title = "Percent White Population \nand Sale Price", 
       subtitle = "Boulder County, CO") +
  geom_smooth(method = "lm", se=F, colour = "green") +
  plotTheme()

#show grid plot in markdown
grid.arrange(carSFplot, bsmtSFplot, TSFplot, TTplot, RentPlot, HHIncplot, MarriedPlot, Whiteplot, ncol=2, nrow =4)


#MAPPING VARIABLE RELATIONSHIPS

#Wildfire Proximity
wfplot<- ggplot() +
  geom_sf(data= tracts19, fill = 'grey99', color = 'black')+
  geom_sf(data = Wildfire_History,fill = 'red',alpha = .75) +  
  geom_sf(data = wf_buffer_25,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = wf_buffer_5,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = wf_buffer_75,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = wf_buffer_1,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = knownSales, aes(colour = q5(price)),
          size = .4) +
  scale_colour_manual(values = paletteMap)+
  labs(title = "House Proximity to Wildfire Risk and Sale Price", 
       subtitle = "Boulder County, CO",
       caption = 'figure 1.1') +
  mapTheme()+theme_ipsum_rc()
cor.test(knownSales$wildfire.nn, knownSales$price, method = "pearson")

#Lakes Map
lakeplot<- ggplot() +
  geom_sf(data= tracts19, fill = 'grey99', color = 'black')+
  geom_sf(data = Lakes_Reservoir, fill = 'red', alpha = .85) +  
  geom_sf(data = lake_buffer_25,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = lake_buffer_5,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = lake_buffer_75,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = lake_buffer_1,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = knownSales, aes(colour = q5(price)),
          size = .4) +
  scale_colour_manual(values = paletteMap)+
  labs(title = "House Proximity to Lakes and Reservoirs \nand Sale Price" ,
       subtitle = "Boulder County, CO",
       caption = 'figure 1.2') +
  mapTheme()+theme_ipsum_rc()
cor.test(knownSales$lakes.nn, knownSales$price, method = "pearson")

#Geohazard Plot
geoplot<- ggplot() +
  geom_sf(data= tracts19, fill = 'grey99', color = 'black')+
  geom_sf(data = GeoHazard, color= 'red', fill = 'red', alpha = .85) +  
  geom_sf(data = geo_buffer_25,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = geo_buffer_5,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = geo_buffer_75,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = geo_buffer_1,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = knownSales, aes(colour = q5(price)),
          size = .4) +
  scale_colour_manual(values = paletteMap)+
  labs(title = "House Proximity to Geological Hazard Sites \nand Sale Price", 
       subtitle = "Boulder County, CO",
       caption = 'figure 1.3') +
  mapTheme()+theme_ipsum_rc()
cor.test(knownSales$geohazard.nn, knownSales$price, method = "pearson")


#FloodRisk Plot
floodplot<- ggplot() +
  geom_sf(data= tracts19, fill = 'grey99', color = 'black')+
  geom_sf(data = FloodZone,color = 'red',fill = 'red', alpha = .85) +  
  geom_sf(data = fl_buffer_25,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = fl_buffer_5,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = fl_buffer_75,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = fl_buffer_1,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = knownSales, aes(colour = q5(price)),
          size = .4) +
  scale_colour_manual(values = paletteMap)+
  labs(title = "House Proximity to Flood Zone and Sale Price", 
       subtitle = "Boulder County, CO",
       caption = 'figure 1.4') +
  mapTheme()+theme_ipsum_rc()
cor.test(knownSales$floodzone.nn, knownSales$price, method = "pearson")

#Business and Commercial Zoning
Businessplot<- ggplot() +
  geom_sf(data= tracts19, fill = 'grey99', color = 'black')+
  geom_sf(data = BusinessCommercial,color = 'red',fill = 'red', alpha = .85) +  
  geom_sf(data = BC_buffer_25,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = BC_buffer_5,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = BC_buffer_75,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = BC_buffer_1,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = knownSales, aes(colour = q5(price)),
          size = .4) +
  scale_colour_manual(values = paletteMap)+
  labs(title = "House Proximity to Business Zoning \nand Sale Price", 
       subtitle = "Boulder County, CO",
       caption = 'figure 1.5') +
  mapTheme()+theme_ipsum_rc()
cor.test(knownSales$Business.nn, knownSales$price, method = "pearson")

#Industrial Zoning
Industryplot<- ggplot() +
  geom_sf(data= tracts19, fill = 'grey99', color = 'black')+
  geom_sf(data = Industrial,color = 'red',fill = 'red', alpha = .85) +  
  geom_sf(data = I_buffer_25,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = I_buffer_5,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = I_buffer_75,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = I_buffer_1,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = knownSales, aes(colour = q5(price)),
          size = .4) +
  scale_colour_manual(values = paletteMap)+
  labs(title = "House Proximity to Industrial Zoning \nand Sale Price", 
       subtitle = "Boulder County, CO",
       caption = 'figure 1.6') +
  mapTheme()+theme_ipsum_rc()
cor.test(knownSales$Industrial.nn, knownSales$price, method = "pearson")

#Agricultural zoning
Agplot<- ggplot() +
  geom_sf(data= tracts19, fill = 'grey99', color = 'black')+
  geom_sf(data = Ag,color = 'red',fill = 'red', alpha = .85) +  
  geom_sf(data = AG_buffer_25,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = AG_buffer_5,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = AG_buffer_75,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = AG_buffer_1,fill = 'red',color = NA,alpha = .1) +
  geom_sf(data = knownSales, aes(colour = q5(price)),
          size = .4) +
  scale_colour_manual(values = paletteMap)+
  labs(title = "House Proximity to Agriculture\nand Sale Price", 
       subtitle = "Boulder County, CO",
       caption = 'figure 1.7') +
  mapTheme()+theme_ipsum_rc()
cor.test(knownSales$Ag.nn, knownSales$price, method = "pearson")

#Streams Plot
#streamplot<- ggplot() +
#  geom_sf(data= tracts19, fill = 'grey99', color = 'black')+
#  geom_sf(data = Streams,
#          fill = 'red',
#          alpha = .85) +  
#  geom_sf(data = Streams_buffer,
#          fill = 'red',
#          color = 'red',
#          alpha = .15) +
#  geom_sf(data = knownSales, aes(colour = q5(price)),
#          size = .4) +
#  scale_colour_manual(values = paletteMap)+
#  labs(title = "1/4 mi Stream Buffer and Home Sale Price", 
#       subtitle = "Boulder County, CO",
#       caption = 'figure #') +
#  mapTheme()+plotTheme()

#show grid maps in markdown
grid.arrange(wfplot, geoplot, 
             floodplot,lakeplot, ncol=2, nrow =2)
grid.arrange(Businessplot, Industryplot, 
             Agplot, ncol = 2, nrow = 2)


#Building Regression Model
knownSales_nog <- st_drop_geometry(knownSales)

reg <- lm(LNprice ~ MedRent+ HouseAge+ PctMarried+BelPov100+wildfire.nn+geohazard.nn
          +LNlakes.nn+floodzone.nn+Ag.nn+Business.nn+Industrial.nn+NumTrailHeads+bsmtSF+HeatingDscr+TotalBach+
            MedHHInc+TravelTime+nbrBedRoom+qualityCodeDscr+
            TotalFinishedSF+designCodeDscr, data = knownSales_nog)

summary(reg)
# Our adjusted R sq is 65.25% - not super great
broom::tidy(reg)%>%
  kable()%>%
  kable_styling()%>%
  footnote(general_title = "Regression Output",
           general = "Table 1")

#separate the test and training sets
inTrain <- createDataPartition(
  y = paste(knownSales$qualityCodeDscr, knownSales$designCodeDscr, knownSales$HeatingDscr),
  p = .60, list = FALSE)
boulder.training <- knownSales[inTrain,]
boulder.test <- knownSales[-inTrain,]  

#training regression
reg.training <-lm(LNprice ~ MedRent+ HouseAge+ PctMarried+BelPov100+wildfire.nn+geohazard.nn
                  +LNlakes.nn+floodzone.nn+Ag.nn+Business.nn+Industrial.nn+NumTrailHeads+bsmtSF+HeatingDscr+TotalBach+
                    MedHHInc+TravelTime+nbrBedRoom+qualityCodeDscr+
                    TotalFinishedSF+designCodeDscr, data = boulder.training)

summary(reg.training)
#stargazer summary
training.summary<- stargazer(reg.training, type="text", digits=1, title="Table 2: LM of Training Data", out = "Training LM.txt")

broom::tidy(reg.training)%>%
      kable()%>%
      kable_styling()%>%
      footnote(general_title = "Training Regression Output",
                general = "Table 2")

#test regression on boulder.test
boulder.test <-
  boulder.test %>%
  mutate(Regression = "Baseline Regression",
         SalePrice.Predict = exp(predict(reg.training, boulder.test)),
         SalePrice.Error = SalePrice.Predict - price,
         SalePrice.AbsError = abs(SalePrice.Predict - price),
         SalePrice.APE = (abs(SalePrice.Predict - price)) / SalePrice.Predict)%>%
  filter(price < 4000000)

#Testing Accuracy 

mean(boulder.test$SalePrice.Predict, na.rm = T)
#[1] 725547.6 average home sale price in the test set is $725,547.6 

mean(boulder.test$SalePrice.AbsError, na.rm = T)
#[1] 160063.7 - average absolute error is $160,063.7

mean(boulder.test$SalePrice.APE, na.rm = T)
#[1] 0.242315 - average absolute percent error is 24.2 % - model errs by 24.2%

ggplot(data = boulder.test, aes(x = SalePrice.Predict, y = price)) +
  geom_point(size=2, shape=20)  +
  labs(title = "Predicted Sale Price As A \nFunction of Observed Price", 
       subtitle = "Orange line represents a perfect prediction\nGreen line represents average prediction \nBoulder County, CO") +
  geom_smooth(method = "lm", se=F, colour = "green")+
  geom_smooth(data = boulder.test, aes(price, price), color = 'orange', se=F)+
  theme_ipsum_rc()

#Test Generalizability
fitControl <- trainControl(method = "cv", number = 100)
set.seed(825)

reg.cv <- 
  train(LNprice ~ MedRent+ HouseAge+ PctMarried+BelPov100+wildfire.nn+geohazard.nn
        +LNlakes.nn+floodzone.nn+Ag.nn+Business.nn+Industrial.nn+NumTrailHeads+bsmtSF+HeatingDscr+TotalBach+
          MedHHInc+TravelTime+nbrBedRoom+qualityCodeDscr+
          TotalFinishedSF+designCodeDscr, data = st_drop_geometry(knownSales),
        method = "lm", trControl = fitControl, na.action = na.pass)

reg.cv

reg.cv$resample[1:5,]

reg.cv.rs.min <- reg.cv$resample[80,]

reg.cv.rs.min %>%                     
  gather(Variable, Value) %>%
  group_by(Variable) %>%
  spread(Variable, Value) %>%
  kable(caption = "Regression Results of One Test Set") %>%
  kable_classic(full_width = F, html_font = "Helvetica")

ggplot(reg.cv$resample, aes(x=MAE)) +
  geom_histogram() +
  labs(title = "Mean Average Error in Cross Validation Tests",
       caption = "figue #") +
  theme_ipsum_rc()

#Test Spatial Clustering

coords <- st_coordinates(knownSales)
neighborList <- knn2nb(knearneigh(coords, 5)) #5 nearest neighborhoods
spatialWeights <- nb2listw(neighborList, style="W") #not sure what is W here
knownSales$lagPrice <- lag.listw(spatialWeights, knownSales$price)


coords.test <-  st_c(boulder.test) 
neighborList.test <- knn2nb(knearneigh(coords.test, 5))
spatialWeights.test <- nb2listw(neighborList.test, style="W")

boulder.test$lagPriceError <- boulder.test%>% 
  mutate(lagPriceError = lag.listw(spatialWeights.test, SalePrice.Error))%>%
  pull(price)


#interpretation of price as function of spatial lag of price -
#as home sale price increases, so does the average price of a home's 
# 5 nearest neighbors - the correlation of this spatial lag relationship
# is .79 and is highly statistically significant - substantial evidence
# fpr clustering of home prices.

ggplot(data = boulder.test, aes(x = lagPriceError, y = SalePrice.AbsError)) +
  geom_point(size=2, shape=20, color = 'orange')  +
  labs(x = 'Spatial lag of errors (Mean error of 5 nearest neighors)',
       y = 'Sale Price Error',
    title = "Error as a function of the spatial lag of price") +
  geom_smooth(aes(lagPriceError,SalePrice.AbsError),method = "lm", se=F, colour = "green")+
  theme_ipsum_rc()

cor.test(boulder.test$lagPriceError, boulder.test$SalePrice.AbsError, method = "pearson")

#interpret

moranTest <- moran.mc(boulder.test$SalePrice.Error, 
                      spatialWeights.test, nsim = 999)

ggplot(as.data.frame(moranTest$res[c(1:999)]), aes(moranTest$res[c(1:999)])) +
  geom_histogram(binwidth = 0.01) +
  geom_vline(aes(xintercept = moranTest$statistic), colour = "#FA7800",size=1) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(title="Observed and permuted Moran's I",
       subtitle= "Observed Moran's I in orange",
       x="Moran's I",
       y="Count") +
  theme_ipsum_rc()

#histogram to show frequency of all 999 randomly permutated I 
#the observed I is higher than all 999 randomly generated I's suggests 
#strong evidence for spatial autocorrelation. 
#both spatial lag test and moran's I confirm spatial autocorrelation - this
#suggests that variation in price likely related to the spatial process
#has been omitted from our current model

#Neighorhood Fixed Effect - accounting for the neighborhood
#Join bytract since we dont have neighborhood data

left_join(
  st_drop_geometry(boulder.test) %>%
    group_by(GEOID) %>%
    summarize(meanPrice = mean(price, na.rm = T)),
  mutate(boulder.test, predict.fe = 
           predict(lm(price ~ GEOID, data = boulder.test), 
                   boulder.test)) %>%
    st_drop_geometry %>%
    group_by(GEOID) %>%
    summarize(meanPrediction = mean(predict.fe))) %>%
  kable() %>% kable_styling()

reg.tract <- lm(LNprice ~ ., data = as.data.frame(boulder.training) %>% 
                  dplyr::select(LNprice, GEOID, MedRent, HouseAge, PctMarried, BelPov100, wildfire.nn, geohazard.nn,
                                LNlakes.nn, floodzone.nn, Ag.nn, Business.nn, Industrial.nn, NumTrailHeads, 
                                bsmtSF, HeatingDscr, TotalBach,
                                MedHHInc, TravelTime, nbrBedRoom, qualityCodeDscr,
                                TotalFinishedSF, designCodeDscr))
  
summary(reg.tract)
boulder.test.tract <-
  boulder.test %>%
  mutate(Regression = "Neighborhood Effects",
         SalePrice.Predict = exp(predict(reg.tract, boulder.test)),
         SalePrice.Error = SalePrice.Predict - price,
         SalePrice.AbsError = abs(SalePrice.Predict - price),
         SalePrice.APE = (abs(SalePrice.Predict - price)) / SalePrice.Predict)%>%
  filter(price < 4000000)

mean(boulder.test.tract$SalePrice.Predict, na.rm = T)
#[1] 725547.6 average home sale price in the test set is $725,547.6 

mean(boulder.test.tract$SalePrice.AbsError, na.rm = T)
#[1] 160063.7 - average absolute error is $160,063.7

mean(boulder.test.tract$SalePrice.APE, na.rm = T)
#[1] 0.242315 - average absolute percent error is 24.2 % - model errs by 24.2%

bothRegressions <- 
  rbind(
    dplyr::select(boulder.test, starts_with("SalePrice"), price, Regression, GEOID) %>%
      mutate(lagPriceError = lag.listw(spatialWeights.test, SalePrice.Error)),
    dplyr::select(boulder.test.tract, starts_with("SalePrice"), price, Regression, GEOID) %>%
      mutate(lagPriceError = lag.listw(spatialWeights.test, SalePrice.Error)))

both.reg <- st_drop_geometry(bothRegressions) %>%
  gather(Variable, Value, -Regression, -GEOID) %>%
  filter(Variable == "SalePrice.AbsError" | Variable == "SalePrice.APE") %>%
  group_by(Regression, Variable) %>%
  summarize(meanValue = mean(Value, na.rm = T)) %>%
  spread(Variable, meanValue) %>%
  kable()%>%kable_styling()

both.reg

bothRegressions %>%
  dplyr::select(SalePrice.Predict,price, Regression) %>%
  ggplot(aes(price, SalePrice.Predict)) +
  geom_point() +
  stat_smooth(aes(price, price), 
              method = "lm", se = FALSE, size = 1, colour="#FA7800") + 
  stat_smooth(aes(SalePrice.Predict, price), 
              method = "lm", se = FALSE, size = 1, colour="#25CB10") +
  facet_wrap(~Regression) +
  labs(title="Predicted sale price as a function of observed price",
       subtitle="Orange line represents a perfect prediction; Green line represents prediction") +
  theme_ipsum_rc()

housesPredictions <-                   
  HouseData %>%
  mutate(prediction = exp(predict(reg.tract, HouseData)),
         team_name = 'HROS')

predictions <- housesPredictions[,c("prediction", "team_name", "toPredict")] %>%
  st_drop_geometry() %>%
  filter(toPredict == 1) %>%
  dplyr::select(-toPredict)

predictplot<- ggplot() +
  geom_sf(data= tracts19, fill = 'grey99', color = 'black')+
  geom_sf(data = housesPredictions, aes(colour = q5(prediction)))+ 
  scale_colour_manual(values = paletteMap)+
  labs(title = "Predicted Sale Price Values", 
       subtitle = "Boulder County, CO",
       caption = 'figure #') +
  mapTheme()+theme_ipsum_rc()

predictplot

st_drop_geometry(bothRegressions) %>%
  group_by(Regression, GEOID)%>%
  summarise(mean.MAPE = mean(SalePrice.APE, na.rm = T)) %>%
  ungroup() %>%
  left_join(tracts19) %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(colour = "gray", aes(fill = q5(mean.MAPE))) +
  scale_fill_manual(values = paletteMap) +
  labs(title = "Figure 12.1: Mean Average Percentage Error by Neighborhood") +
  mapTheme()

scatter_tracts <-
  boulder.test.tract %>%
  group_by(GEOID) %>%
  dplyr::select(GEOID, SalePrice.APE, SalePrice.Predict)

mean_sca_t <-
  scatter_tracts %>%
  group_by(GEOID) %>%
  summarise_at(vars("SalePrice.APE", "SalePrice.Predict"), mean)

ggplot()+
  geom_point(data = mean_sca_t, 
             aes(SalePrice.Predict,SalePrice.APE),
             size = 2,
             color = "#003A6B" )+
  labs(x = 'Mean Price by Census Tract',
       y = 'MAPE by Census Tract',
       title = "MAPE by Census Tract and Mean Price by Census Tract", 
       subtitle = "Boulder County, CO",
       caption = "figure #") + 
  theme_ipsum_rc()

#Race Context

tracts19<- tracts19%>%
  mutate(raceContext = ifelse(PctWhite > .5, "Majority White", "Majority Non-White"),
       incomeContext = ifelse(MedHHInc > 65000, "High Income", "Low Income"))

grid.arrange(ncol = 2,
             ggplot() + geom_sf(data = na.omit(tracts19), aes(fill = raceContext)) +
               scale_fill_manual(values = c("#25CB10", "#FA7800"), name="Race Context") +
               labs(title = "Race Context") +
               mapTheme() + theme(legend.position="bottom"), 
             ggplot() + geom_sf(data = na.omit(tracts19), aes(fill = incomeContext)) +
               scale_fill_manual(values = c("#25CB10", "#FA7800"), name="Income Context") +
               labs(title = "Income Context") +
               mapTheme() + theme(legend.position="bottom"))

st_join(bothRegressions, tracts19) %>% 
  group_by(Regression, raceContext) %>%
  summarize(mean.MAPE = scales::percent(mean(SalePrice.APE, na.rm = T))) %>%
  st_drop_geometry() %>%
  spread(raceContext, mean.MAPE) %>%
  kable(caption = "Test set MAPE by neighborhood racial context")%>%kable_styling()

st_join(bothRegressions, tracts19) %>% 
  filter(!is.na(incomeContext)) %>%
  group_by(Regression, incomeContext) %>%
  summarize(mean.MAPE = scales::percent(mean(SalePrice.APE, na.rm = T))) %>%
  st_drop_geometry() %>%
  spread(incomeContext, mean.MAPE) %>%
  kable(caption = "Test set MAPE by neighborhood income context")%>%kable_styling()


