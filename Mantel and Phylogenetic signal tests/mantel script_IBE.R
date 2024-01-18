#Script for testing the correlations between environmental variables and genetic distance of Microcoleus strains (Mantel)
library(sp)
library(maps)
library(raster)
library(dismo)
library(maptools)
library(rgeos)
library(rgdal)
library(vegan)
library(geosphere)

coords <- read.csv("coords.csv")
#download tif files for all bioclimatic variables from https://www.worldclim.org/data/bioclim.html
#rename tif files to go from a, not from 1 (R takes order 1,10,11,2,3,etc.)
#Extract the bioclimatic variables
raster_files_bioclimatic <- list.files(path = "./wc2.1_2.5m_bio/", pattern = 'tif$', full.names = TRUE)
bioclimatic <- stack(raster_files_bioclimatic)
bioclim_variables <- extract(bioclimatic, coords[,2:3])

bioclim_names <- c("1_Ann_T", "2_Diurnal_Range", "3_Isothermality", "4_T_Seasonality", "5_MaxT_Wrmst_Month", "6_MinT_Cldst_Month", "7_T_Ann_Range", "8_T_Wettest_Qtr", "9_T_Driest_Qtr", "10_T_Wrmst_Qtr", "11_T_Cldst_Qtr", "12_Ann_Precip", "13_P_Wettest_Month", "14_P_Driest_Month", "15_P_Seasonality", "16_P_Wettest_Qtr", "17_P_Driest_Qtr", "18_P_Wrmst_Qtr", "19_P_Cldst_Qtr")
colnames(bioclim_variables) <- bioclim_names

bioclim_variables_df <- as.data.frame(bioclim_variables)
write.csv(bioclim_variables_df, "./bioclimatic variables.csv")

#test which variables are autocorrelated to remove them
library(caret)
bioclimatic_variables_cor <- cor(bioclim_variables_df)
hc <- findCorrelation(bioclimatic_variables_cor, cutoff = 0.9)
hc <- sort(hc)
reduced_variables <- bioclimatic_variables_cor[,-c(hc)]

#autocorrelated variables were 2 Diurnal range,3 Isothermality,4 T Seasonality
#5 MaxTWarmest Month, 7 T annual range, 8 T wettest Q, 9 T Driest Q, 15 P Seasonality
#16 P wettest Q, 17 P Driest Q, 18 P warmest Q, 19 P Coldest Quarter
#remove them and load that new set of variables with strains names
reduced.bioclimatic.variables <- read.csv("reduced bioclimatic variables.csv", row.names=1)

#Extract other variables downloaded from https://www.worldclim.org/data/bioclim.html
files_wind <- list.files("wc2.1_2.5m_wind/", pattern = 'tif', full.names = TRUE)
biowind <- stack(files_wind)
wind_speed <- extract(biowind, microcoleus_gps[,2:3])

files_solrad <- list.files("wc2.1_2.5m_srad/", pattern = 'tif', full.names = TRUE)
biosolrad <- stack(files_solrad)
sol_rad <- extract(biosolrad, microcoleus_gps[,2:3])

files_water_vapor_pressure <- list.files("wc2.1_2.5m_vapr/", pattern = 'tif', full.names = TRUE)
biovapor <- stack(files_water_vapor_pressure)
vapor <- extract(biovapor, microcoleus_gps[,2:3])

#Variables for soil propery, UV-B light, and HANPP were extracted using QGIS
#After extracting all variables, put them in one joint file All.variables.mantel

#load genetic distance matrix (ANI from 0-1)
#genetic distance matrix is the matrix of ANI values turned into binary matrix (0-1).
library(jvamisc)
ANI <- read_excel("ANI.xlsx", col_names = FALSE)
min(ANI, na.rm = T) #minimal ANI is 82.42258
ANI <- ANI[,-1]
#Transform lower triangle matrix into a symmetric one
ANI$...2 <- as.numeric(ANI$...2)
ANI$...292 <- as.numeric(ANI$...292)
ANI_n <- upper2full(t(ANI), diagval = 100)
ANI_b <- ANI_n / 100
genetic_distance <- as.dist(1 - ANI_b)

#perform Mantel tests
AnnT <- dist(reduced.bioclimatic.variables$Ann_T, method = "euclidean")
MinTCld <- dist(reduced.bioclimatic.variables$MinT_Cldst_Month, method = "euclidean")
TWarmQ <- dist(reduced.bioclimatic.variables$T_Wrmst_Qtr, method = "euclidean")
TCldQ <- dist(reduced.bioclimatic.variables$T_Cldst_Qtr, method = "euclidean")
AnnPrec <- dist(reduced.bioclimatic.variables$Ann_Precip, method = "euclidean")
PrecWetMon <- dist(reduced.bioclimatic.variables$P_Wettest_Month, method = "euclidean")
PrecDriestMon <- dist(reduced.bioclimatic.variables$P_Driest_Month, method = "euclidean")

mantel(genetic_distance, AnnT, permutations = 9999, na.rm = TRUE)
mantel(genetic_distance, MinTCld, permutations = 9999, na.rm = TRUE)
mantel(genetic_distance, TWarmQ, permutations = 9999, na.rm = TRUE)
mantel(genetic_distance, TCldQ, permutations = 9999, na.rm = TRUE)
mantel(genetic_distance, AnnPrec, permutations = 9999, na.rm = TRUE)
mantel(genetic_distance, PrecWetMon, permutations = 9999, na.rm = TRUE)
mantel(genetic_distance, PrecDriestMon, permutations = 9999, na.rm = TRUE)

#load UV-B variables and remove highly correlated ones
UV.radiation.mantel <- read.csv("E:/Mantel_variables_coords_new/UV radiation mantel.csv", row.names=1)
uv_cor <- cor(UV.radiation.mantel)
uv_hc <- findCorrelation(uv_cor, cutoff = 0.9)
uv_hc <- sort(uv_hc)
reduced_uv <- sort(uv_hc)
#kept 1,4,and 5th, removed mean UV-Bs
reduced.UV.radiation.mantel <- read.csv("E:/Mantel_variables_coords_new/reduced UV radiation mantel.csv", row.names=1)
AnnMeanUvB <- dist(reduced.UV.radiation.mantel$Annual.Mean.UV.B, method = "euclidean")
MonMeanUvHighQ <- dist(reduced.UV.radiation.mantel$Sum.of.Monthly.Mean.UV.B.during.Highest.Quarter, method = "euclidean")
MonMeanUvLowQ <- dist(reduced.UV.radiation.mantel$Sum.of.Monthly.Mean.UV.B.during.Lowest.Quarter, method = "euclidean")

mantel(genetic_distance, AnnMeanUvB, permutations = 9999, na.rm = TRUE)
mantel(genetic_distance, MonMeanUvHighQ, permutations = 9999, na.rm = TRUE)
mantel(genetic_distance, MonMeanUvLowQ, permutations = 9999, na.rm = TRUE)

#The rest variables
All.variables.mantel  <- read.csv("Environmental_variables_all.csv", row.names=1)
elevation <- dist(All.variables.mantel $Elevation, method = "euclidean")
mantel(genetic_distance, elevation, permutations = 9999, na.rm = TRUE)
Nfertilizer <- dist(All.variables.mantel $Fertilizer.global, method = "euclidean")
mantel(genetic_distance, Nfertilizer, permutations = 9999, na.rm = TRUE)
tdnpplcpal <- dist(All.variables.mantel $tdnpplcpal, method = "euclidean")
mantel(genetic_distance, tdnpplcpal, permutations = 9999, na.rm = TRUE)
thanpppall <- dist(All.variables.mantel $thanpppall, method = "euclidean")
mantel(genetic_distance, thanpppall, permutations = 9999, na.rm = TRUE)
thanppallg <- dist(All.variables.mantel $thanppallg, method = "euclidean")
mantel(genetic_distance, thanppallg, permutations = 9999, na.rm = TRUE)
tn0 <- dist(All.variables.mantel $tn0_all_gc, method = "euclidean")
mantel(genetic_distance, tn0, permutations = 9999, na.rm = TRUE)
tnapallg <- dist(All.variables.mantel $tnap_all_g, method = "euclidean")
mantel(genetic_distance, tnapallg, permutations = 9999, na.rm = TRUE)
tntpallg <- dist(All.variables.mantel $tntp_all_g, method = "euclidean")
mantel(genetic_distance, tntpallg, permutations = 9999, na.rm = TRUE)
maxT <- dist(All.variables.mantel $Maximum.Temperature..average., method = "euclidean")
mantel(genetic_distance, maxT, permutations = 9999, na.rm = TRUE)
minT <- dist(All.variables.mantel $Minimum.temperature..average., method = "euclidean")
mantel(genetic_distance, minT, permutations = 9999, na.rm = TRUE)
popcount <- dist(All.variables.mantel $Population.count..average., method = "euclidean")
mantel(genetic_distance, popcount, permutations = 9999, na.rm = TRUE)
popdens <- dist(All.variables.mantel $Population.density..average., method = "euclidean")
mantel(genetic_distance, popdens, permutations = 9999, na.rm = TRUE)
rural <- dist(All.variables.mantel $Rural.population.count, method = "euclidean")
mantel(genetic_distance, rural, permutations = 9999, na.rm = TRUE)
urban <- dist(All.variables.mantel $Urban.area, method = "euclidean")
mantel(genetic_distance, urban, permutations = 9999, na.rm = TRUE)
rainfed <- dist(All.variables.mantel $Total.rainfed, method = "euclidean")
mantel(genetic_distance, rainfed, permutations = 9999, na.rm = TRUE)
avgPrec <- dist(All.variables.mantel $Precipitation..average., method = "euclidean")
mantel(genetic_distance, avgPrec, permutations = 9999, na.rm = TRUE)
ocd <- dist(All.variables.mantel $ocd_0.5cm_, method = "euclidean")
mantel(genetic_distance, ocd, permutations = 9999, na.rm = TRUE)
ocarbons <- dist(All.variables.mantel $ocs_0.30cm, method = "euclidean")
mantel(genetic_distance, ocarbons, permutations = 9999, na.rm = TRUE)
soc <- dist(All.variables.mantel $soc_0.5cm_, method = "euclidean")
mantel(genetic_distance, soc, permutations = 9999, na.rm = TRUE)
silt <- dist(All.variables.mantel $silt_0.5cm, method = "euclidean")
mantel(genetic_distance, silt, permutations = 9999, na.rm = TRUE)
sand <- dist(All.variables.mantel $sand_0.5cm, method = "euclidean")
mantel(genetic_distance, sand, permutations = 9999, na.rm = TRUE)
phh2o <- dist(All.variables.mantel $phh2o_0.5c, method = "euclidean")
mantel(genetic_distance, phh2o, permutations = 9999, na.rm = TRUE)
nitrogen <- dist(All.variables.mantel $nitrogen_0, method = "euclidean")
mantel(genetic_distance, nitrogen, permutations = 9999, na.rm = TRUE)
clay <- dist(All.variables.mantel $clay_0.5cm, method = "euclidean")
mantel(genetic_distance, clay, permutations = 9999, na.rm = TRUE)
coarse <- dist(All.variables.mantel $cfvo_0.5cm, method = "euclidean")
mantel(genetic_distance, coarse, permutations = 9999, na.rm = TRUE)
cec <- dist(All.variables.mantel $cec_0.5cm_, method = "euclidean")
mantel(genetic_distance, cec, permutations = 9999, na.rm = TRUE)
bulk <- dist(All.variables.mantel $bdod_0.5cm, method = "euclidean")
mantel(genetic_distance, bulk, permutations = 9999, na.rm = TRUE)
solarrad <- dist(All.variables.mantel $Solar.radiation..average., method = "euclidean")
mantel(genetic_distance, solarrad, permutations = 9999, na.rm = TRUE)
vapor <- dist(All.variables.mantel $Vapor..average., method = "euclidean")
mantel(genetic_distance, vapor, permutations = 9999, na.rm = TRUE)
wind <- dist(All.variables.mantel $Wind.speed..average., method = "euclidean")
mantel(genetic_distance, wind, permutations = 9999, na.rm = TRUE)




