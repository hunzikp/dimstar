
###################################
# Dependencies
###################################

library(pgrid)
library(cshapes)
library(raster)
library(velox)

###################################
# Get prio-grid data
###################################

variables <- c('gwno', 'ttime_mean', 'mountains_mean', 'pop_gpw_sum', 'nlights_mean', 'excluded')
years <- 1997:2014
prio.ls <- getPrioRaster(names = c('gid', variables), years = 1997:2014)
prio.meta.df <- prio.ls$meta
prio.vx <- prio.ls$raster

###################################
# Get African country borders
###################################

cshp.spdf <- cshp(as.Date('2013-01-01'))
africa.spdf <- cshp.spdf[(cshp.spdf$GWCODE >= 400 & cshp.spdf$GWCODE <= 626) | cshp.spdf$GWCODE == 651,]

###################################
# Make prio-grid base sample
###################################

# Crop all rasters to Africa
prio.vx$crop(x = africa.spdf)

# Get the base raster of all relevant raster cells
base.ras <- prio.vx$as.RasterLayer(band = prio.meta.df$band[prio.meta.df$name == 'gwno' & prio.meta.df$year == 2014])
base.ras[base.ras < 400 | base.ras > 651] <- NA
base.ras[base.ras > 626 & !(base.ras == 651)] <- NA
base.ras[!is.na(base.ras)] <- 1

# Get base data frame of GIDs
gid.ras <- prio.vx$as.RasterLayer(band = prio.meta.df$band[prio.meta.df$name == 'gid' & prio.meta.df$year == 2014])
gid.ras <- gid.ras*base.ras
gid.df <- na.omit(as.data.frame(gid.ras))
names(gid.df) <- 'gid'

# Make spdf corresponding to base data frame GIDs
gid.spdf <- rasterToPolygons(gid.ras)
names(gid.spdf) <- 'gid'

# Get list of prio-grid cross-sections
crosssection.ls <- vector('list', length(years))

for (t in 1:length(years)) {
  thisgid.df <- gid.df
  year <- years[t]
  for (j in 1:length(variables)) {
    variable <- variables[j]

    years.avail <- prio.meta.df$year[prio.meta.df$name == variable]
    if (any(is.na(years.avail))) {
      # if variable is static, just select the (only) band
      variable_band <- prio.meta.df$band[prio.meta.df$name == variable]
    } else {
      # otherwise, select the closest year
      diff <- abs(years.avail - year)
      closest_year <- (years.avail[which.min(diff)])[1]
      variable_band <- prio.meta.df$band[prio.meta.df$name == variable & prio.meta.df$year == closest_year]
    }

    variable.ras <- prio.vx$as.RasterLayer(band = variable_band)
    variable.df <- as.data.frame(variable.ras)
    base.vec <- as.data.frame(base.ras)[,1]
    variable.df <- variable.df[!is.na(base.vec),,drop=FALSE]
    names(variable.df) <- variable
    thisgid.df <- cbind(thisgid.df, variable.df)
  }
  thisgid.df$year <- year
  crosssection.ls[[t]] <- thisgid.df
}


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
for (t in 1:length(years)) {
  year <- years[t]
  thisgid.df <- crosssection.ls[[t]]
  this.acled.spdf <- acled.spdf[acled.spdf$year == year,]

  for (j in 1:length(variables)) {
    variable <- variables[j]
    this.event.spdf <- this.acled.spdf[,variable]
    names(this.event.spdf) <- 'event'
    this.event.spdf <- this.event.spdf[this.event.spdf$event == 1,]

    event.ras <- rasterize(x = this.event.spdf, base.ras, field = 'event', fun = sum)
    event.ras[is.na(event.ras)] <- 0
    event.df <- as.data.frame(event.ras)
    base.vec <- as.data.frame(base.ras)[,1]
    event.df <- event.df[!is.na(base.vec),,drop=FALSE]
    names(event.df) <- variable

    thisgid.df <- cbind(thisgid.df, event.df)
  }
  crosssection.ls[[t]] <- thisgid.df
}

###################################
# Only keep those GIDs for which we have no missings
###################################

## Replace excluded NAs with zeros
for (t in 1:length(years)) {
  this.df <- crosssection.ls[[t]]
  this.df$excluded[is.na(this.df$excluded)] <- 0
  crosssection.ls[[t]] <- this.df
}

## Remove NAs
nonmissing.mat <- matrix(NA, nrow(gid.spdf), length(years))
for (t in 1:length(years)) {
  complete <- complete.cases(crosssection.ls[[t]])
  nonmissing.mat[,t] <- complete
}
nonmissing <- apply(nonmissing.mat, 1, all)

## Subset polygons
gid.spdf <- gid.spdf[nonmissing,]

## Subset crosssections and make panel
for (t in 1:length(years)) {
  crosssection.ls[[t]] <- crosssection.ls[[t]][nonmissing,]
}
gidpanel.df <- do.call(rbind, crosssection.ls)


###################################
# Save the data set
###################################

data.ls <- list(spdf = gid.spdf, panel = gidpanel.df)
saveRDS(data.ls, file = 'data/prio_acled.rds')






