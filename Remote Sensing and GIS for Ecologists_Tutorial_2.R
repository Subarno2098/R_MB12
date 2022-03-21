library(raster)
# p224r63_2011_B1 <- raster("D:/EAGLE/MB2/Data/R_MB12/LT05_L1TP_224063_20110729_20200820_02_T1/LT05_L1TP_224063_20110729_20200820_02_T1_B1.TIF")
# p224r63_2011_B1

r1 <- raster("D:/EAGLE/MB2/Data/R_MB12/LT05_L1TP_224063_20110729_20200820_02_T1/LT05_L1TP_224063_20110729_20200820_02_T1_B1.TIF")
r2 <- raster("D:/EAGLE/MB2/Data/R_MB12/LT05_L1TP_224063_20110729_20200820_02_T1/LT05_L1TP_224063_20110729_20200820_02_T1_B2.TIF")
r3 <- raster("D:/EAGLE/MB2/Data/R_MB12/LT05_L1TP_224063_20110729_20200820_02_T1/LT05_L1TP_224063_20110729_20200820_02_T1_B3.TIF")
r4 <- raster("D:/EAGLE/MB2/Data/R_MB12/LT05_L1TP_224063_20110729_20200820_02_T1/LT05_L1TP_224063_20110729_20200820_02_T1_B4.TIF")
r5 <- raster("D:/EAGLE/MB2/Data/R_MB12/LT05_L1TP_224063_20110729_20200820_02_T1/LT05_L1TP_224063_20110729_20200820_02_T1_B5.TIF")
r6 <- raster("D:/EAGLE/MB2/Data/R_MB12/LT05_L1TP_224063_20110729_20200820_02_T1/LT05_L1TP_224063_20110729_20200820_02_T1_B6.TIF")
r7 <- raster("D:/EAGLE/MB2/Data/R_MB12/LT05_L1TP_224063_20110729_20200820_02_T1/LT05_L1TP_224063_20110729_20200820_02_T1_B7.TIF")

p224r63_2011 <- stack(c(r1, r2, r3, r4, r5, r6, r7))
p224r63_2011

# changing the names of the band (since the names of the bands are too long)
names(p224r63_2011) <- gsub(pattern = "LT05_L1TP_224063_20110729_20200820_02_T1_", replacement = "", x = names(p224r63_2011))
names(p224r63_2011)
p224r63_2011

# writing raster in the working directory (saving it)
writeRaster(p224r63_2011, filename = "D:/EAGLE/MB2/Data/R_MB12/raster_1.tif" )

# Saving as raster brick
p224r63_2011 <- brick(p224r63_2011)
p224r63_2011


# converting a raster brick into a SpatialPixelsDataFrame

p224r63_2011_sp <- as(p224r63_2011, "SpatialPixelsDataFrame")
class(p224r63_2011_sp)


# WORKING WITH VECTOR DATA in R 
library(maptools)

# Reading vector data 

# spatialpolygondataframe
studyArea <-  readOGR("D:/EAGLE/MB2/Data/R_MB12/Remote Sensing and GIS Data/vector_data", "study_area_UTMnorth")
studyArea

# spatialpolygondataframe
pa <- readOGR("D:/EAGLE/MB2/Data/R_MB12/Remote Sensing and GIS Data/vector_data", "PAs_UNEP_WCMC_p224r63")
pa

#  spatiallinesdataframe
roads <- readOGR("D:/EAGLE/MB2/Data/R_MB12/Remote Sensing and GIS Data/vector_data", "roads_p224r63_UTMsouth")
roads 


# importing GPS data
gpsPoints <- readOGR(dsn = "D:/EAGLE/MB2/Data/R_MB12/Remote Sensing and GIS Data/vector_data/field_measuerements.gpx", layer = "waypoints")


# loading sp library 
library(sp)

# reading a csv file into R 
csv_file <- read.csv("D:/EAGLE/MB2/Data/R_MB12/csv_file.csv")

class(csv_file)

# converting the csv into point data using SpatialPoints
csv.sp <- SpatialPoints(coords = csv_file[c("X", "Y")])
csv.sp

# converting spatial points into spatial points data frame. (since object does not have any values associated with it)
csv.spdf <- SpatialPointsDataFrame(coords = csv_file[c("X", "Y")], data = csv_file[3:5])
csv.spdf

# Alternative workflow for converting data frame into spatial points data frame 
coordinates(csv_file) <- c("X", "Y")
class(csv_file)

# Defining the projection of the spatial data frame same as p224r63_2011 (without having the same projection, points cannot be plotted in the raster)
projection(csv.spdf) <- projection(p224r63_2011)
csv.spdf

# writing the shape files (csv.spdf) to disk using write OGR. 
writeOGR(obj = csv.spdf, "D:/EAGLE/MB2/Data/R_MB12", "csv_spdf_as_shp", driver = "ESRI Shapefile")


# PLOTTING OF RASTER DATA
plot(p224r63_2011)

# plotting only a single layer of raster data (b5 in this case)
plot(p224r63_2011, 5)
# or
plot(p224r63_2011, "B5")

# changing the colour pallete of the plot
plot(p224r63_2011, "B5", col = grey.colors(100))


# set layer transparency to 50% 
plot(p224r63_2011, "B5", alpha = 0.5)

# if plots are taking too much time to process maxpixels = 2e+05 can be used as it reduces the number of pixels to process.
plot(p224r63_2011, "B5", maxpixels = 2e+05)

# However, when generating the final plots, it should be avoided as it will make the scene pixelated. 

#plotting using spplot (lattice plotting system)
spplot(p224r63_2011, "B5")

# plotting multiple bands together using spplot 
spplot(p224r63_2011, 1:4)

# plotting using ggR (ggplot2)
ggR(p224r63_2011, "B5")


# ggR with legend 
ggR(p224r63_2011, "B5", geom_raster = TRUE)

# with a custom legend
# ggR(p224r63_2011, 5, geom_raster = TRUE) + scale_fill_gradient(colours = rainbow(100))



# plotting the raster using plotRGB (helps to calculate NDVI etc)
plotRGB(p224r63_2011, r = "B3", g = "B2", b = "B1", scale = 1) #does not work
plotRGB(p224r63_2011, r = "B3", g = "B2", b = "B1", stretch = "lin") 

# plotting false composite color (generally used to highlight the differences in vegetation)
plotRGB(p224r63_2011, r = "B4", g = "B3", b = "B2", stretch = "lin")

# same plot as above but using ggplot 
ggRGB(p224r63_2011, r = "B4", g = "B3", b = "B2", stretch = "lin")


# adding spatial points on top of raster data
plot(p224r63_2011, "B5")
points(csv.spdf, col = "blue")


# plotting points on the raster using ggplot
# first the spatial point data frame has to be converted into a data frame 
df_pts <- as.data.frame(csv.spdf)

# plotting df_pts using ggplot on rgb
ggRGB(p224r63_2011, 4, 3, 2, stretch = "lin") + geom_point(data = df_pts, aes(x = X, y = Y), size = 5, col = "yellow")


# BASIC SPATIAL DATA HANDLING in R
# Indexing Raster Data 

# ex. indexing single raster layer of p224r63_2011
values <- p224r63_2011_B1[1, ] #first row
head(values)

values <- p224r63_2011_B1[, 2] #second coloumn
values <- p224r63_2011[, 7:300] #coloumns from 7 to 300
values <- p224r63_2011[1, 2] #single pixel: 1st row and 2nd Coloumn

# extracting the first number of rows and coloumns 

nr <- nrow(p224r63_2011_B1)
nr

nc <- ncol(p224r63_2011_B1)
nc

p224r63_2011_B1[c(1, nr), c(1, nc)]


# accessing different layers of a raster 
# Band 2 
p224r63_2011[[2]]
p224r63_2011[["B2"]]

# band 5
p224r63_2011[[5]]
p224r63_2011[["B5"]]


# Querying multiple layers of rasters 
# for band 1 to 5
p224r63_2011[[1:5]]
# or
p224r63_2011[[c("B1", "B5")]]

# for band 4 to 7
p224r63_2011[[4:7]]
# or
p224r63_2011[[c("B4", "B7")]]

# eg
# p224r63_2011[c(1, 3, 7)]

# removing layers from raster 
p224r63_2011.drop1and5 <- p224r63_2011[[-c(1, 5)]]
p224r63_2011.drop1and5
# or
p224r63_2011.drop1and5 <- dropLayer(p224r63_2011, c(1, 5))
p224r63_2011.drop1and5
# or
p224r63_2011.drop1and5 <- dropLayer(p224r63_2011, c("B1", "B5"))
p224r63_2011.drop1and5


# adding a new layer to a raster 
newLayer <- p224r63_2011[[1]]/2
names(newLayer) <- "new"
p224r63_2011.add <- addLayer(p224r63_2011, newLayer)
p224r63_2011.add
# or
p224r63_2011.add[["myNewLayer"]] <- newLayer
p224r63_2011.add


# dropping the added layers (Test for myself)
p224r63_2011.rmls <- dropLayer(p224r63_2011.add, c("new", "myNewLayer"))
p224r63_2011.rmls


# get layer 5, rows 3– 7 and column 90
p224r63_2011[[5]][3:7, 90]

# modifying the existing values 
layerCopy <- p224r63_2011[[5]]
layerCopy[70:71, 90:91]

# replacing the values for row 70 and coloumn 90 to 91 as 3
layerCopy[70, 90:91] <- 3

layerCopy[70:71, 90:91]

# changing all the values of layerCopy
layerCopy[] <- 0
layerCopy

# logical raster 
# pixels with < 0.2 
queryRaster <- p224r63_2011[[5]] < 0.2
dataType(queryRaster)

# returing the values which are true 
values <- p224r63_2011[queryRaster]
values

# setting all the values less than 0.2 to NA
layerCopy[queryRaster] <- NA
layerCopy



# INDEXING VECTOR DATA 

# retrieving the attribute table of csv.spdf by using @data
csv.spdf@data

# accessing the coordinates of spatialpointsdataframe 
coordinates(csv.spdf) 

# indexing the first row of csv.spdf
subs <- csv.spdf[1, ]
subs

# keeping coloumns 1 and 3 except the first row
subs <- csv.spdf[-1, c(1, 3)]
subs

# returning just the values of coloumn without the whole spatial object
subs <-  csv.spdf[["plot_name"]]
subs
# or
subs <- csv.spdf$plot_name
subs

# subsetting csv.spdf to return values only less than 20

subs <- csv.spdf[csv.spdf[["value"]] < 20, ]
subs
# or 
subs <- csv.spdf[csv.spdf$value < 20, ]
subs


# modifying the values inside csv.spdf (changing the value of 1st coloumn to 77)
csv.spdf[1, "value"] <- 77
csv.spdf@data

# ADDING INFORMATION TO VECTOR DATA
# adding a new coloumn (however it should always match the order of the existing spatial object)
csv.spdf$myNewCol <- c(7, 19, 1, NA)
csv.spdf
# or
csv.spdf[["myNewCol"]] <- c(7, 19, 1, NA)
csv.spdf

# adding more than one coloumn 
newData <- data.frame(d1 = c("a", "b", "c", "d"), d2 = 1:4, name = csv.spdf$plot_name)
newData

# merging csv.spdf and newData using cbind
csv.spdf@data <- cbind(csv.spdf@data, newData)
csv.spdf@data


# merging using common ID (IMP) 
csv.spdf <- csv.spdf[, 1:4] #reset csv.spdf to initial coloumns
csv.spdf@data

newData <- newData[c(4, 2, 1, 3), ] #shuffling the data rowwise 
newData 

# merging  using sp::merge()
library(sp)

csv.spdf <- merge(x = csv.spdf, y = newData, by.x = "plot_name", by.y = "name")
csv.spdf@data


# CRS AND DATA PROJECTIONS 
