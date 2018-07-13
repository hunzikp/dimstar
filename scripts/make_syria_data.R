
library(rgdal)
library(pgrid)
library(cshapes)
library(raster)
library(velox)

TARGET <- 'Syria'
CODE <- 'ADM2_CODE'

####################################
# Load GAUL
####################################

admin.spdf <- readOGR('/home/hunzikp/Data/gaul/g2015_2014_2/', 'g2015_2014_2')
units.spdf <- admin.spdf[grepl(TARGET, admin.spdf$ADM0_NAME),]
units.spdf$gid <- units.spdf@data[,CODE]
units.spdf <- units.spdf[,c("gid", "Shape_Area")]

####################################
# Add 1990 population data
####################################

pop.ras <- raster('/home/hunzikp/Data/grump/asia/asup90ag.bil')
#pop.ras <- raster("/home/hunzikp/Data/grump/africa/afup90ag.bil")
pop.vx <- velox(pop.ras)
pop.vx$crop(units.spdf)
pop.vec <- pop.vx$extract(sp = units.spdf, fun = function(x) sum(x, na.rm = TRUE), legacy = TRUE)
units.spdf$pop <- pop.vec


####################################
# Add time-invariant prio grid data
####################################

variables <- c('gwno', 'ttime_mean', 'mountains_mean')
prio.ls <- getPrioRaster(names = c('gid', variables), years = 2014)
prio.meta.df <- prio.ls$meta
prio.vx <- prio.ls$raster
prio.vx$crop(units.spdf)

pdata.mat <- prio.vx$extract(sp = units.spdf, fun = function(x) mean(x, na.rm = TRUE), small = TRUE)
units.spdf$ttime <- pdata.mat[,3]
units.spdf$mountains <- pdata.mat[,4]


####################################
# Make admin units cross sections
####################################

units.spdf <- units.spdf[,c('gid', 'Shape_Area', 'pop', 'ttime', 'mountains')]
units.df <- units.spdf@data

# startdates <- seq(as.Date("2008-07-01"), as.Date("2018-06-30"), by = "month")
# enddates <- seq(as.Date("2008-08-01"), as.Date("2018-07-01"), by = "month")-1
startdates <- seq(as.Date("2017-01-01"), as.Date("2018-06-23"), by = "week")
enddates <- seq(as.Date("2017-01-01"), as.Date("2018-06-23"), by = "week")+6

crosssection.ls <- rep(list(units.df), length(startdates))
crosssection.ls <- lapply(1:length(startdates), function(i) {
  crosssection.ls[[i]]$startdate <- startdates[i]
  crosssection.ls[[i]]$enddate <- enddates[i]
  crosssection.ls[[i]]
})

###################################
# Add acled data
###################################

acled17.df <- read.csv("/home/hunzikp/Data/acled/2018/MiddleEast_2016-2018_upd-Jul3.csv", stringsAsFactors = FALSE, fileEncoding="latin1")
#acled17.df <- read.csv("/home/hunzikp/Data/acled/2018/Africa_1997-2018_upd-Jul2.csv", stringsAsFactors = FALSE, fileEncoding="latin1")
keep <- c("ISO", "EVENT_ID_CNTY", "EVENT_DATE", "YEAR", "TIME_PRECISION",
          "EVENT_TYPE", "ACTOR1", "ACTOR2",
          "INTERACTION",
          "COUNTRY", "LOCATION", "LATITUDE",
          "LONGITUDE","FATALITIES", "NOTES")
acled17.df <- acled17.df[,keep]
acled.df <- acled17.df
acled.df <- acled.df[acled.df$COUNTRY == TARGET,]

## Correct dates
acled.df$EVENT_DATE <- as.Date(acled.df$EVENT_DATE, "%d-%B-%Y")
names(acled.df) <- tolower(names(acled.df))

## Add event type markers
acled.df$battle_event <- grepl(acled.df$event_type, pattern = 'attle')
acled.df$protest_event <- grepl(acled.df$event_type, pattern = 'rotest')
acled.df$civilian_event <- grepl(acled.df$event_type, pattern = 'ivilian') | grepl(acled.df$actor2, pattern = 'ivilian')

## Add special event markers
acled.df$strike_event <- grepl(acled.df$notes, pattern = "strike")
acled.df$isis_event <- grepl(acled.df$actor1, pattern = "slamic State") | grepl(acled.df$actor2, pattern = "slamic State")
acled.df$isis_battle <- acled.df$isis_event*acled.df$battle_event
acled.df$isis_civilian <- acled.df$isis_event*acled.df$civilian_event
acled.df$isis_strike <- acled.df$isis_event*acled.df$strike_event

## Get rid of notes
acled.df <- acled.df[,-which(names(acled.df) == "notes")]

## Make spatial points
acled.spdf <- SpatialPointsDataFrame(coords = acled.df[,c('longitude', 'latitude')], data = acled.df)

## Add acled data to priogrid cross-sections
variables <- c('battle_event', 'protest_event', 'civilian_event', "strike_event", "isis_event", "isis_battle", "isis_civilian", "isis_strike")
periods <- 1:length(crosssection.ls)
for (t in 1:length(periods)) {
  startdate <- crosssection.ls[[t]]$startdate[1]
  enddate <- crosssection.ls[[t]]$enddate[1]
  thisgid.df <- crosssection.ls[[t]]
  this.acled.spdf <- acled.spdf[acled.spdf$event_date >= startdate & acled.spdf$event_date <= enddate,]
  
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
saveRDS(data.ls, file = 'data/syria_acled.rds')



