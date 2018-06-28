
library(rgdal)
library(pgrid)
library(cshapes)
library(raster)
library(velox)

TARGET <- 'Uganda'
YEARS <- 1997:2017

####################################
# Load GAUL
####################################

admin.spdf <- readOGR('/home/hunzikp/Data/gaul/g2015_1997_1/', 'g2015_1997_1')
units.spdf <- admin.spdf[admin.spdf$ADM0_NAME == TARGET,]


####################################
# Add 1990 population data
####################################

pop.ras <- raster('/home/hunzikp/Data/grump/africa/afup90ag.bil')
pop.vx <- velox(pop.ras)
pop.vec <- pop.vx$extract(sp = units.spdf, fun = function(x) sum(x, na.rm = TRUE), legacy = TRUE)
units.spdf$pop <- pop.vec


####################################
# Add time-invariant prio grid data
####################################

variables <- c('gwno', 'ttime_mean', 'mountains_mean')
prio.ls <- getPrioRaster(names = c('gid', variables), years = 1997)
prio.meta.df <- prio.ls$meta
prio.vx <- prio.ls$raster
prio.vx$crop(units.spdf)

pdata.mat <- prio.vx$extract(sp = units.spdf, fun = function(x) mean(x, na.rm = TRUE), small = TRUE)
units.spdf$ttime <- pdata.mat[,3]
units.spdf$mountains <- pdata.mat[,4]


####################################
# Make admin units cross sections
####################################

units.spdf <- units.spdf[,c('ADM1_CODE', 'Shape_Area', 'pop', 'ttime', 'mountains')]
units.df <- units.spdf@data

crosssection.ls <- rep(list(units.df), length(YEARS))
crosssection.ls <- lapply(1:length(YEARS), function(i) {
  crosssection.ls[[i]]$year <- YEARS[i]
  crosssection.ls[[i]]
})

###################################
# Add acled data
###################################

acled17.df <- read.csv("/home/hunzikp/Data/acled/2017/ACLED-Africa_1997-2018_upd-Feb20.csv", stringsAsFactors = FALSE, fileEncoding="latin1")
keep <- c("ISO", "EVENT_ID_CNTY", "EVENT_DATE", "YEAR", "TIME_PRECISION",
          "EVENT_TYPE", "ACTOR1", "ACTOR2",
          "INTERACTION",
          "COUNTRY", "LOCATION", "LATITUDE",
          "LONGITUDE","FATALITIES")
acled17.df <- acled17.df[,keep]
acled.df <- acled17.df

## Correct dates
acled.df$EVENT_DATE <- as.Date(acled.df$EVENT_DATE, "%d-%B-%Y")
names(acled.df) <- tolower(names(acled.df))

## Add event type markers
acled.df$battle_event <- grepl(acled.df$event_type, pattern = 'attle')
acled.df$protest_event <- grepl(acled.df$event_type, pattern = 'rotest')
acled.df$civilian_event <- grepl(acled.df$event_type, pattern = 'ivilian')

## Make spatial points
acled.spdf <- SpatialPointsDataFrame(coords = acled.df[,c('longitude', 'latitude')], data = acled.df)

## Add acled data to priogrid cross-sections
variables <- c('battle_event', 'protest_event', 'civilian_event')
for (t in 1:length(YEARS)) {
  year <- YEARS[t]
  thisgid.df <- crosssection.ls[[t]]
  this.acled.spdf <- acled.spdf[acled.spdf$year == year,]
  
  for (j in 1:length(variables)) {
    variable <- variables[j]
    this.event.spdf <- this.acled.spdf[,variable]
    names(this.event.spdf) <- 'event'
    this.event.spdf <- this.event.spdf[this.event.spdf$event == 1,]
    
    intrs <- gIntersects(units.spdf, this.event.spdf, byid = TRUE)
    event.df <- data.frame(count = colSums(intrs))
    names(event.df) <- variable
    
    thisgid.df <- cbind(thisgid.df, event.df)
  }
  crosssection.ls[[t]] <- thisgid.df
}

###################################
# Make panel
###################################

panel.df <- do.call(rbind, crosssection.ls)


###################################
# Save the data set
###################################

data.ls <- list(spdf = units.spdf, panel = panel.df)
saveRDS(data.ls, file = 'data/admin_acled.rds')
