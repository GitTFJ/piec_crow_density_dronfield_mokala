library(sp)
library(sf)
library(dplyr)
library(raster)
library(geosphere)
library(Distance)
library(dsm)
library(ggplot2)
library(ggeffects)
library(ggpubr)

segments_df = NULL
tran_list = list()
for(a in 1:8){
  gpx = maptools::getKMLcoordinates(paste0("tran",a,".kml"))
  sfc <- st_as_sfc(
    lapply(1, function(x) {
      st_linestring(gpx[[1]][,c(1,2)], dim = "XY")
    }),
    crs = 4326
  )
  tran_list[[a]] = sfc
  lin <- st_transform(sfc, CRS("+proj=utm +zone=35 +datum=WGS84 +units=m +no_defs"))
  segs <- st_segmentize(lin, dfMaxLength=units::set_units(100, "metres"))
  
  stdh_cast_substring <- function(x, to = "MULTILINESTRING") {
    ggg <- st_geometry(x)
    
    if (!unique(st_geometry_type(ggg)) %in% c("POLYGON", "LINESTRING")) {
      stop("Input should be  LINESTRING or POLYGON")
    }
    for (k in 1:length(st_geometry(ggg))) {
      sub <- ggg[k]
      geom <- lapply(
        1:(length(st_coordinates(sub)[, 1]) - 1),
        function(i)
          rbind(
            as.numeric(st_coordinates(sub)[i, 1:2]),
            as.numeric(st_coordinates(sub)[i + 1, 1:2])
          )
      ) %>%
        st_multilinestring() %>%
        st_sfc()
      
      if (k == 1) {
        endgeom <- geom
      }
      else {
        endgeom <- rbind(endgeom, geom)
      }
    }
    endgeom <- endgeom %>% st_sfc(crs = st_crs(x))
    if (class(x)[1] == "sf") {
      endgeom <- st_set_geometry(x, endgeom)
    }
    if (to == "LINESTRING") {
      endgeom <- endgeom %>% st_cast("LINESTRING")
    }
    return(endgeom)
  }
  segs2 <- stdh_cast_substring(segs, to="LINESTRING")
  segments_tmp = data.frame(
    tran = a,
    lat = st_transform(segs, 4326)[[1]][-1,2],
    lon = st_transform(segs, 4326)[[1]][-1,1],
    sou = segs[[1]][-1,2],
    eas = segs[[1]][-1,1],
    eff = st_length(segs2),
    lab = paste0("tran", a, " - ", 1:length(segs2))
  )
  segments_df = rbind(segments_df, segments_tmp)
}
hist(segments_df$eff)


tran_report = read.csv("C:/Users/mn826766/OneDrive - University of Reading/Personal/White-backedVulture/White-backedVulture/Publications/pied_crow_abundance/transect_report.csv")
tran_report = tran_report %>%
  group_by(Transect) %>%
  dplyr::summarise(reps = max(Repeat))

segments_df = left_join(segments_df, tran_report, by = c("tran" = "Transect"))
segments_df$eff_adj = segments_df$eff*segments_df$reps

NIR = raster("landsat8/LC08_L1TP_172080_20150610_20180526_01_T1_B5.tif")
RED = raster("landsat8/LC08_L1TP_172080_20150610_20180526_01_T1_B4.tif")
NDVI = (NIR - RED)/(NIR + RED)
NDVI_points <- SpatialPoints(segments_df[,c("eas", "sou")], proj4string = CRS("+proj=utm +zone=35 +datum=WGS84 +units=m +no_defs"))
segments_df$ndvi = extract(NDVI, NDVI_points, method='bilinear')

dron = Polygon(maptools::getKMLcoordinates("poly_dronfield.kml")[[1]][,c(1,2)])
dron = Polygons(list(dron),1)
dron = SpatialPolygons(list(dron), proj4string = crs("+proj=longlat +datum=WGS84 +no_defs"))
dron_utm = spTransform(dron, CRS("+proj=utm +zone=35 +datum=WGS84 +units=m +no_defs"))
dron_ndvi_utm = crop(NDVI, dron_utm)
dron_ndvi_utm = mask(dron_ndvi_utm, dron_utm)
dron_ndvi_ll = projectRaster(dron_ndvi_utm, crs = crs("+proj=longlat +datum=WGS84 +no_defs"))
dron_ndvi_ll_df = as.data.frame(rasterToPoints(dron_ndvi_ll))
colnames(dron_ndvi_ll_df) = c("lon", "lat", "ndvi")

mok = Polygon(maptools::getKMLcoordinates("poly_mokala.kml")[[1]][,c(1,2)])
mok = Polygons(list(mok),1)
mok = SpatialPolygons(list(mok), proj4string = crs("+proj=longlat +datum=WGS84 +no_defs"))
mok_utm = spTransform(mok, CRS("+proj=utm +zone=35 +datum=WGS84 +units=m +no_defs"))
mok_ndvi_utm = crop(NDVI, mok_utm)
mok_ndvi_utm = mask(mok_ndvi_utm, mok_utm)
mok_ndvi_ll = projectRaster(mok_ndvi_utm, crs = crs("+proj=longlat +datum=WGS84 +no_defs"))
mok_ndvi_ll_df = as.data.frame(rasterToPoints(mok_ndvi_ll))
colnames(mok_ndvi_ll_df) = c("lon", "lat", "ndvi")

obs_df = read.csv("obs.csv")
obs_df$lab = NA
obs_df$ndvi = NA
for(a in 1:nrow(obs_df)){
  tran = obs_df[a,"Transect.No"]
  sg_df = segments_df[which(segments_df$tran == tran ),]
  pnt = which.min(distm(sg_df[,c("lon","lat")], obs_df[a,c("Longitude","Latitude")]))
  obs_df$lab[a] = paste0("tran", tran, " - ",pnt)
  obs_df$ndvi[a] = sg_df[pnt,"ndvi"]
}
obs_df$object = 1:nrow(obs_df)
segments_df$site = ifelse(segments_df$tran < 5, "Dronfield", "Mokala")

dron_obs = obs_df[which(obs_df$Location == "Dronfield" & obs_df$Poss_psr == "No"),]
dron_obs$Transect.No = factor(dron_obs$Transect.No, levels = c("1", "2", "3", "4"))
plt_dron_tran = ggplot() +
  geom_tile(data = dron_ndvi_ll_df, aes(x = lon, y = lat, fill = ndvi)) +
  geom_sf(data = tran_list[[1]], colour = "red", size = 1.2, alpha = 0.5) +
  geom_sf(data = tran_list[[2]], colour = "green", size = 1.2, alpha = 0.5) +
  geom_sf(data = tran_list[[3]], colour = "blue", size = 1.2, alpha = 0.5) +
  geom_sf(data = tran_list[[4]], colour = "purple", size = 1.2, alpha = 0.5) +
  geom_jitter(data = dron_obs,aes(x = Longitude, y = Latitude), width = 0.004, height = 0.004, size = 2, alpha = 0.55) +
  scale_fill_gradient(low = "white", high = "black", guide = F) +
  scale_x_continuous(breaks = c(24.78, 24.83, 24.88)) +
  scale_y_continuous(breaks = c(-28.7, -28.65, -28.6)) +
  #scale_colour_discrete(guide = F) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic()  

ggsave(plot = plt_dron_tran,"transect_routes_dron.png", width = 4, height = 4, dpi = 300)

mok_obs = obs_df[which(obs_df$Location == "Mokala" & obs_df$Poss_psr == "No"),]
mok_obs$Transect.No = factor(mok_obs$Transect.No, levels = c("5", "6", "7", "8"))
plt_mok_tran = ggplot() +
  geom_tile(data = mok_ndvi_ll_df, aes(x = lon, y = lat, fill = ndvi)) +
  geom_sf(data = tran_list[[5]], colour = "red", size = 1.2, alpha = 0.5) +
  geom_sf(data = tran_list[[6]], colour = "green", size = 1.2, alpha = 0.5) +
  geom_sf(data = tran_list[[7]], colour = "blue", size = 1.2, alpha = 0.5) +
  geom_sf(data = tran_list[[8]], colour = "purple", size = 1.2, alpha = 0.5) +
  geom_jitter(data = mok_obs,aes(x = Longitude, y = Latitude), width = 0.004, height = 0.004, size = 2, alpha = 0.55) +
  scale_fill_gradient(low = "white", high = "black", guide = F) +
  scale_x_continuous(breaks = c(24.3, 24.4, 24.5)) +
  scale_y_continuous(breaks = c(-29.2, -29.1, -29)) +
  #scale_colour_discrete(guide = F) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() 

ggsave(plot = plt_mok_tran,"transect_routes_mok.png", width = 4.6, height = 4.6, dpi = 300)

segdata = segments_df[,c(7,1,3,2,5,4,9,10,11)]
obs_df = subset(obs_df, Poss_psr == "No")
distdata = obs_df[,c(17,15,6,8,9,16,5)]


colnames(segdata) = c("Sample.Label", "Transect", "lon", "lat", "eas", "sou", "Effort","ndvi","site")
colnames(distdata) = c("object", "Sample.Label", "Transect", "size", "distance", "ndvi","site")

distdata$distance2 = ifelse(distdata$distance == "0-25m",
                            12.5,
                            NA)
distdata$distance2 = ifelse(distdata$distance == "25-50m",
                            37.5,
                            distdata$distance2)
distdata$distance2 = ifelse(distdata$distance == "50-100m",
                            75,
                            distdata$distance2)
distdata$distance2 = ifelse(distdata$distance == "100-200m",
                            150,
                            distdata$distance2)
distdata$distance2 = ifelse(distdata$distance == "+200m",
                            200,
                            distdata$distance2)
distdata$distance = distdata$distance2


distdata$site = as.factor(distdata$site)
segdata$Effort = as.numeric(segdata$Effort)
segdata$Effort = ifelse(segdata$Effort < 100, 100 , segdata$Effort)
segdata$site = as.factor(segdata$site)

detfc_ndvi_site = ds(data = distdata, 
           truncation = list(left = 1, right = 250), 
           transect = "line",
           formula=~ ndvi + site,
           cutpoints = c(1, 25, 50, 100, 200, 250))
#AIC 128.537

detfc_ndvi = ds(data = distdata, 
                     truncation = list(left = 1, right = 250), 
                     transect = "line",
                     formula=~ ndvi,
                     cutpoints = c(1, 25, 50, 100, 200, 250))
#AIC 126.722

detfc_site = ds(data = distdata, 
                     truncation = list(left = 1, right = 250), 
                     transect = "line",
                     formula=~ site,
                     cutpoints = c(1, 25, 50, 100, 200, 250))
#AIC 127.578

detfc = ds(data = distdata, 
                     truncation = list(left = 1, right = 250), 
                     transect = "line",
                     cutpoints = c(1, 25, 50, 100, 200, 250))
#AIC 121.497
summary(detfc)
png("detection.png", res = 300, width = 1200, height = 1200)
plot(detfc)
dev.off()


dron_segdata = subset(segdata,site == "Dronfield") 
dron_distdata = subset(distdata,site == "Dronfield")
dron_obsdata = dron_distdata[,c(1,2,4,5)]
dron_dsm.xy <- dsm(count ~ s(eas,sou, k = 10) + ndvi, detfc, dron_segdata, dron_obsdata, method="REML", family = tw(), select = T)
summary(dron_dsm.xy)
plot(dron_dsm.xy)
gam.check(dron_dsm.xy)
rqgam.check(dron_dsm.xy)

dron_ndvi_utm_df = as.data.frame(rasterToPoints(dron_ndvi_utm))
colnames(dron_ndvi_utm_df) = c("eas", "sou", "ndvi")
dron_ndvi_utm_df$site = "Dronfield"
dron_ndvi_utm_df = dron_ndvi_utm_df[which(
  dron_ndvi_utm_df$eas < 
    max(segdata[which(
      segdata$site == "Dronfield"
    ),"eas"]) 
  &
    dron_ndvi_utm_df$eas > 
    min(segdata[which(
      segdata$site == "Dronfield"
    ),"eas"]) 
  &
    dron_ndvi_utm_df$sou < 
    max(segdata[which(
      segdata$site == "Dronfield"
    ),"sou"]) 
  &
    dron_ndvi_utm_df$sou > 
    min(segdata[which(
      segdata$site == "Dronfield"
    ),"sou"]) 
),]
dron_dsm.xy.pred <- predict(dron_dsm.xy, dron_ndvi_utm_df, off.set = 900)
dsm.var.gam(dron_dsm.xy, dron_ndvi_utm_df, off.set = 900)
dron_ndvi_utm_df$pred = dron_dsm.xy.pred
dron_spdf <- SpatialPointsDataFrame(coords = dron_ndvi_utm_df[,c("eas", "sou")], 
                                    data = dron_ndvi_utm_df[,c("ndvi","pred")],
                                    proj4string = CRS("+proj=utm +zone=35 +datum=WGS84 +units=m +no_defs"))

dron_spdf = spTransform(dron_spdf,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
dron_spdf = data.frame(
  lon = dron_spdf@coords[,1],
  lat = dron_spdf@coords[,2],
  pred = dron_spdf@data[,2]
)

dron_boundary = data.frame(
  lon = dron@polygons[[1]]@Polygons[[1]]@coords [,1],
  lat = dron@polygons[[1]]@Polygons[[1]]@coords [,2]
)

8.1/(0.0009*nrow(dron_spdf))
16.6/(0.0009*nrow(dron_spdf))
34.3/(0.0009*nrow(dron_spdf))

plot_den_dron = ggplot() + 
  geom_tile(data = dron_spdf, aes(x = lon, y = lat, fill = pred/0.0009),width=0.00031,height=0.00031) +
  scale_fill_gradient(low = "grey90", high = "grey10", name = expression("Individuals/"~km^2)) +
  geom_polygon(data = dron_boundary, aes(x = lon, y = lat), alpha = 0.1, fill = NA, colour = "black") +
  labs(x = "Longitude", y = "Latitude", title = paste0(
    "Dronfield\n\nAbundance: 16.6 [8.1, 34.3]\nDensity: 0.24 [0.12, 0.50]\nCV: 0.38")) +
  scale_x_continuous(breaks = c(24.78, 24.83, 24.88)) +
  scale_y_continuous(breaks = c(-28.7, -28.65, -28.6)) +
  theme_classic() +
  theme(
    legend.position = "bottom"
  )



mok_segdata = subset(segdata,site == "Mokala") 
mok_distdata = subset(distdata,site == "Mokala")
mok_obsdata = mok_distdata[,c(1,2,4,5)]
mok_dsm.xy <- dsm(count ~ s(eas,sou, k = 10) + ndvi, detfc, mok_segdata, mok_obsdata, method="REML", family = tw(), select = T)
summary(mok_dsm.xy)
plot(mok_dsm.xy)
gam.check(mok_dsm.xy)
rqgam.check(mok_dsm.xy)


mok_ndvi_utm_df = as.data.frame(rasterToPoints(mok_ndvi_utm))
colnames(mok_ndvi_utm_df) = c("eas", "sou", "ndvi")
mok_ndvi_utm_df$site = "Mokala"
mok_ndvi_utm_df = mok_ndvi_utm_df[which(
  mok_ndvi_utm_df$eas < 
    max(segdata[which(
      segdata$site == "Mokala"
    ),"eas"]) 
  &
    mok_ndvi_utm_df$eas > 
    min(segdata[which(
      segdata$site == "Mokala"
    ),"eas"]) 
  &
    mok_ndvi_utm_df$sou < 
    max(segdata[which(
      segdata$site == "Mokala"
    ),"sou"]) 
  &
    mok_ndvi_utm_df$sou > 
    min(segdata[which(
      segdata$site == "Mokala"
    ),"sou"]) 
),]
mok_dsm.xy.pred <- predict(mok_dsm.xy, mok_ndvi_utm_df, off.set = 900)
dsm.var.gam(mok_dsm.xy, mok_ndvi_utm_df, off.set = 900)
mok_ndvi_utm_df$pred = mok_dsm.xy.pred
mok_ndvi_utm_df$pred = ifelse(mok_ndvi_utm_df$pred > 0.0027, 0.0027, mok_ndvi_utm_df$pred)
mok_spdf <- SpatialPointsDataFrame(coords = mok_ndvi_utm_df[,c("eas", "sou")], 
                                    data = mok_ndvi_utm_df[,c("ndvi","pred")],
                                    proj4string = CRS("+proj=utm +zone=35 +datum=WGS84 +units=m +no_defs"))

mok_spdf = spTransform(mok_spdf,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
mok_spdf = data.frame(
  lon = mok_spdf@coords[,1],
  lat = mok_spdf@coords[,2],
  pred = mok_spdf@data[,2]
)

mok_boundary = data.frame(
  lon = mok@polygons[[1]]@Polygons[[1]]@coords [,1],
  lat = mok@polygons[[1]]@Polygons[[1]]@coords [,2]
)

11.1/(0.0009*nrow(mok_spdf))
28.1/(0.0009*nrow(mok_spdf))
71.1/(0.0009*nrow(mok_spdf))

plot_den_mok = ggplot() + 
  geom_tile(data = mok_spdf, aes(x = lon, y = lat, fill = pred/0.0009),width=0.00031,height=0.00031) +
  scale_fill_gradient(low = "grey90", high = "grey10", name = expression("Individuals/"~km^2)) +
  geom_polygon(data = mok_boundary, aes(x = lon, y = lat), alpha = 0.1, fill = NA, colour = "black") +
  labs(x = "Longitude", y = " ", title = paste0(
    "Mokala\n\nAbundance: 28.1 [11.1, 71.1]\nDensity: 0.18 [0.07, 0.45]\nCV: 0.50")) +
  scale_x_continuous(breaks = c(24.3, 24.4, 24.5)) +
  scale_y_continuous(breaks = c(-29.2, -29.1, -29)) +
  theme_classic() +
  theme(
    legend.position = "bottom"
  )

ggarrange(plot_den_dron, plot_den_mok, ncol = 2, common.legend = T, legend = "bottom")
ggsave("density_map.png", width = 8, height = 5, dpi = 300)

ndvi_dron = ggpredict(dron_dsm.xy, terms = "ndvi")
ndvi_mok = ggpredict(mok_dsm.xy, terms = "ndvi")

ggplot() +
  geom_line(data = ndvi_dron, aes(x = x, y = predicted/0.0009), colour = "red") +
  geom_ribbon(data = ndvi_dron, aes(x = x, ymin = conf.low/0.0009, ymax = conf.high/0.0009), fill = "red", alpha = 0.1) +
  geom_line(data = ndvi_mok, aes(x = x, y = predicted/0.0009), colour = "blue") +
  geom_ribbon(data = ndvi_mok, aes(x = x, ymin = conf.low/0.0009, ymax = conf.high/0.0009), fill = "blue", alpha = 0.1) +
  labs(x = "NDVI", y = expression(Individuals/km^2)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()
ggsave("ndvi.png", width = 4, height = 4, dpi = 300)
