#' gseason_prec: Calculates the growing season period as determined by precipitation and evapotranspiration balance
#'
#' @description `gseason_prec` calculates growing season period where precipitation exceeds half the potential evapotranspiration.
#'
#' @param prec a vector of precipitation values.
#' @param evap a vector of evapotranspiration values.
#' @param tme a `POSIXlt` object representing the date and time of each `evap` value.
#' @param period an optional character string defining time period of precipitation and evapotranspiration measurements. Options include "monthly" and "bimonthly". If left unspecified, assumed 'submonthly'.
#' @param surplus default = TRUE
#'
#' @return A vector of binary values, where 1 indicates growing season and 0 indicates not growing season.
#' @export
#'
#' @details
#' The growing season is defined as the period when precipitation exceeds half the potential evapotranspiration.
#' If period is "monthly", data are aggregated by each month.
#' If period is "bimonthly", data are aggregated by each half-month.
#' If surplus = 'TRUE', generates a vector of growing season values for the specified time period.
#' Requires wrapper functions `submonthly` and `cumsumseq` to be loaded.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#' @seealso [submonthly()]
#' @seealso [cumsumseq()]
#'
#' @examples
#' prec <- (10 * sin(c(0:364) * (pi / -360)) + rnorm(365) + 15)
#' evap <- (10 * sin(c(0:364) * (pi / +360)) + rnorm(365) + 2)
#' tme <- tmecreate(2010, 24)
#' gsp <- gseason_prec(prec, evap, tme, period = "monthly", surplus = TRUE)
#' gsp <- gseason_prec(prec, evap, tme, period = "bimonthly", surplus = TRUE)
#' plot(gsp)
#'
gseason_prec <- function(prec, evap, tme, period = "monthly", surplus = TRUE) {
  if (length(unique(tme$year)) > 1) warna()
  if (is.na(sd(prec, na.rm = TRUE)) | is.na(sd(evap, na.rm = TRUE))) {
    gs <- NA
  } else {
    if (period == "monthly") {
      pmth <- aggregate(prec, by = list(tme$mon), sum, na.rm = TRUE)$x
      emth <- aggregate(evap, by = list(tme$mon), sum, na.rm = TRUE)$x
      phour <- spline(pmth, n = length(prec))$y
      ehour <- spline(emth, n = length(evap))$y
      phour <- ifelse(phour<0, 0, phour)
      ehour <- ifelse(ehour<0, 0, ehour)
      gs <- ifelse(phour > 0.5 * ehour, 1, 0)
      if (surplus & length(gs[gs < 1]) > 0) {
        gs <- cumsumseq(gs, ehour, 1)
      }
    } else if (period == "bimonthly") {
      sm <- submonthly(tme, 2)
      pmth <- aggregate(prec, by = list(sm, tme$mon), sum, na.rm = TRUE)$x
      emth <- aggregate(evap, by = list(sm, tme$mon), sum, na.rm = TRUE)$x
      phour <- spline(pmth, n = length(prec))$y
      ehour <- spline(emth, n = length(evap))$y
      phour <- ifelse(phour<0, 0, phour)
      ehour <- ifelse(ehour<0, 0, ehour)
      gs <- ifelse(phour > (0.5 * ehour), 1, 0)
      if (surplus & length(gs[gs < 1]) > 0) {
        gs <- cumsumseq(gs, ehour, 2)
      }
    } else {
      sm <- submonthly(tme, 4)
      pmth <- aggregate(prec, by = list(sm, tme$mon), sum, na.rm = TRUE)$x
      emth <- aggregate(evap, by = list(sm, tme$mon), sum, na.rm = TRUE)$x
      phour <- spline(pmth, n = length(prec))$y
      ehour <- spline(emth, n = length(evap))$y
      phour <- ifelse(phour<0, 0, phour)
      ehour <- ifelse(ehour<0, 0, ehour)
      gs <- ifelse(phour > (0.5 * ehour), 1, 0)
      if (surplus & length(gs[gs < 1]) > 0) {
        gs <- cumsumseq(gs, ehour, 4)
      }
    }
  }
  gs
}
#' gseason_temp: Calculates the growing season period as determined by temperature values
#'
#' @description `gseason_temp` calculates growing season period where temperatures are above 'lower' and below 'upper' specified limits.
#'
#' @param temp a vector of temperature values.
#' @param tme a `POSIXlt` object representing the date and time of each `temps` value.
#' @param lower integer specifying lower temperature limit of growing season.
#' @param upper integer specifying upper temperature limit of growing season.
#' @param nday specifies the number of consecutive days of temperatures above 'lower' or below 'upper' before growing season start/end is accepted.
#'
#' @return A time-series object of binary values, where 1 indicates growing season and 0 indicates not growing season.
#' @export
#'
#' @details To satisfy the requirements for `tme`, a POSIXlt object can be created using the `tmecreate` wrapper function.
#'
#' @seealso [tmecreate()] can be used to create a POSIXlt object.
#'
#' @examples
#' temps <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' tme <- tmecreate(2010, 6)
#' gseast <- gseason_temp(temps, tme)
#' plot(gseast)
#'
gseason_temp <- function(temp, tme, lower = 5, upper = 35, nday = 5) {
  dint <- (24 * 3600) / (as.numeric(tme[2]) - as.numeric(tme[1]))
  ma <- function(x, n = nday * dint)
  {
    filter(x, rep(1/n, n), sides = 2, circular = TRUE)
  }
  tma <- ma(temp)
  gs <- ifelse(tma > lower, 1, 0)
  gs <- ifelse(tma > upper, 0, gs)
  gs
}
#' gseason_day: Determine daytime hours
#'
#' @description `gseason_day` is used to calculate daytime (1) and non-daytime (0) hours.
#'
#' @param tme a `POSIXlt` object.
#' @param lat latitude of the location for which gseason_day is required.
#' @param long longitude of the location for which gseason_day is required.
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is 0.
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0). Default is 0.
#'
#' @importFrom microclima julday
#'
#' @return a vector of binary values indicating day (1) or night (0).
#' @export
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#' @seealso Requires wrapper function `sisimple` to be loaded.
#'
#' @examples
#' Daytime hours in Porthleven, Cornwall for 2010.
#' tme <- tmecreate(2010, 1)
#' gseason_day(tme, 50.08, -5.31)
#'
gseason_day <- function(tme, lat, long, merid = 0, dst =0) {
  jd <- julday(tme$year + 1900, tme$mon + 1, tme$mday)
  lt <- tme$hour + tme$min / 60 + tme$sec / 3600
  si <- sisimple(lt, lat, long, jd, merid, dst)
  gs <- ifelse(si > 0, 1, 0)
  gs
}
#'gseason: Calculates the period where plant growth is possible
#'
#' @description `gseason` calculates total annual growing hours where temperatures are within defined upper and lower limits and precipitation exceeds half potential evapotranspiration.
#'
#' @param year calendar year.
#' @param precipnc nc file containing precipitation values.
#' @param evapnc nc file containing evpotranspiration values.
#' @param tempnc nc file containing temperature values.
#' @param period an optional character string defining time period of precipitation and evapotranspiration measurements. Options include "monthly" and "bimonthly". If left unspecified, assumed 'submonthly'.
#' @param surplus either TRUE or FALSE.
#' @param lower defines lower temperature limit where plant growth ceases (degrees Celcius) and growing season is switched off.
#' @param upper defines upper temeprature limit where plant growth ceases (degrees Celcius) and growing season is switched off.
#' @param nday specifies the number of consecutive days of temperatures above 'lower' or below 'upper' before growing season start/end is accepted.
#' @param daynight if TRUE, growing season is continuous over 24 hours, if FALSE, only calculates growing season during the day.
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is 0.
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0). Default is 0.
#'
#' @import ncdf4
#' @import raster
#' @importFrom microclima latlongfromraster
#'
#' @return a three dimensional array of binary values indicating growing conditions (1) = yes, (0) = no.
#' @export
#'
#' @details
#' The growing season period is defined as the period where temperatures are >5 degree Celcius and <35 degree Celcius and precipitation > 0.5 PET.
#'
#' @seealso [gseason_prec()] calculates period where precipitation exceeds half evapotranspiration.
#' @seealso [gseason_temp()] calculates period where temperature are above lower limit for plant growth and below upper limit for plant growth.
#' @seealso [gseason_day()] calculates day/night time.
#' @seealso [tmecreate()] can be used to create a POSIXlt object.
#' @seealso [nctoarray()] can be used to create an nc file from array output.
#'
gseason <- function(year, precipnc, evapnc, tempnc, period = "monthly", surplus = TRUE, lower = 5,
                    upper = 35, nday = 1/24, daynight = TRUE, merid = 0, dst = 0) {
  prec <- nctoarray(precipnc)
  evap <- nctoarray(evapnc)
  temp <- nctoarray(tempnc)
  tme<-tmecreate(year,1)
  r <- raster(precipnc)
  if (daynight) {
    lats <- latlongfromraster(r)$lats
    lons <- latlongfromraster(r)$lons
  }
  gdn <- rep(1,length(tme))
  gshourly <- array(NA, dim = dim(prec))
  for (x in 1:73) {
    print(x)
    for (y in 1:144) {
      gprec <- gseason_prec(prec[x,y,], evapa[x,y,], tme, period, surplus)
      gtemp <- gseason_temp(temp[x,y,], tme, lower, upper, nday)
      if (daynight) gdn <- gseason_day(tme, lats[x,y], lons[x,y], merid, dst)
      gshourly[x,y,] <- gprec * gtemp * gdn
      gshourly
    }
  }
}
#' events: determines the number of events
#'
#' @description `events` is used to calculate the total number of events (e.g. number of rainfall events).
#' @param x a vector with 0, indicating non-event (e.g. no rainfall), and any other value indicating an event.
#'
#' @return an integer representing the total number of events.
#' @export
#'
#' @examples
#' prec <- sample(rep(0:10,10))
#' events(prec)
#'
events <- function(x) {
  x[x>0] <- 1
  x[x<0] <- 0
  x<-c(x,0)
  length(which(diff(x) == -1))
}
#' gssm: Soil moisture content during the growing season
#'
#' @description `gssm` calculates mean soil moisture content during the growing season period.
#'
#' @param gseason a three-dimensional array of binary values indicating growing conditions (1) = yes, (0) = no.
#' @param soilm a three-dimensional array of fractional soil moisture values.
#'
#' @return a matrix of mean soil moisture values during the growing season.
#' @export
#'
#' @seealso the [gseason()] function can be used to create an array of growing conditions (1) = yes, (0) = no accounting for temperature, precipitation and daylight hours.
#'
#' @examples
#' require(microclima)
#' tme <- as.POSIXlt(c(0:1459) * 3600 * 6, origin = "2010-01-01 00:00", tz = "GMT")
#' gs <- gseason_day(tme, 6, 21)
#' gseason <- array(gs, dim=c(1, 1, 1460))
#' soilm <- array(runif(1460, 0, 100), dim= c(1,1,1460))
#' gssm <- gssm(gseason, soilm)
#'
gssm <- function(gseason, soilm) {
  soilm[is.na(soilm)==TRUE]<-0
  ym <- gseason * soilm
  ysoil <- apply(ym, c(1,2), mean, na.rm=T)
  ysoil
}
#' gst: Mean growing season temperature
#'
#' @description `gst` calculates mean temperature during the growing season.
#'
#' @param temp a three-dimensional array of air temperature values.
#' @param gseason a three-dimensional array of binary values indicating growing conditions (1) = yes, (0) = no.
#'
#' @return a matrix of mean air temperature during the growing season.
#' @export
#'
#' @seealso the [gseason()] function can be used to create an array of growing conditions (1) = yes, (0) = no accounting for temperature, precipitation and daylight hours.
#'
#' @examples
#' tme <- as.POSIXlt(c(0:1459) * 3600 * 6, origin = "2010-01-01 00:00", tz = "GMT")
#' gs <- gseason_day(tme, 6, 21)
#' gseason <- array(gs, dim=c(1, 1, 1460))
#' temp <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' meangst <- gst(temp, gseason)
gst <- function(temp, gseason) {
  tempg <- temp * gseason
  tempg[tempg==0]<-NA
  mtgs <- apply(tempg, c(1,2), mean, na.rm=T)
  mtgs
}
#' gsp: Total precipitation during growing season
#'
#' @description `gsp` calculates total precipitation during the growing season.
#'
#' @param prec a three dimensional array of precipitation values.
#' @param gseason a three dimensional array of binary values (1 = growing season, 0 = not growing season).
#'
#' @return a matrix of values of total precipitation during the growing season.
#' @export
#'
#' @seealso the [gseason()] function can be used to create an array of growing conditions (1) = yes, (0) = no accounting for temperature, precipitation and daylight hours.
#'
#' @examples
#' prec <- array(10 * sin(c(0:1459) * (pi / -1400)) + runif(1460, 0, 10) +10, dim=c(1,1,1460))
#' tme <- tmecreate(2010, 6)
#' gs <- gseason_day(tme, 6, 21)
#' gseason <- array(gs, dim=c(1, 1, 1460))
#' gsp(prec, gseason)
#'
gsp <- function(prec, gseason) {
  precg <- prec * gseason
  spg <- apply(precg, c(1,2), sum, na.rm=T)
  spg
}
#' tsp: Total summer precipitation
#'
#' @description `tsp` calculates total precipitation during the summer season
#'
#' @param prec a three dimensional array of precipitation values.
#' @param startday assumed day of year of start of summer in northen hemisphere in non-leap year. Default is day 152 (1st June).
#' @param endday assumed day of year of start of summer in northen hemisphere in non-leap year. Default is day 243 (31st August).
#' @param r a raster object of same extent as temp coded as 1 for northern hemisphere and 0 for southern hemisphere.
#'
#' @importFrom raster getValues
#'
#' @return a matrix of total summer precipitation values for northern and southern hemisphere.
#' @export
#'
#' @details Seasons are flipped in the southern hemisphere i.e. 1st June (day 152) = day 152+365/2+0.5 = 334 = 1st Dec.
#' In leap years, 1 day is added.
#'
#' @examples
#' prec <- array(10 * sin(c(0:1459) * (pi / -1400)) + runif(1460, 0, 10) +10, dim=c(73,144,1460))
#' m <- matrix(1, 73, 144)
#' r <- raster(m, crs="+init=epsg:4326")
#' extent(r) <- c(-1.25, 358.75, -91.25, 91.25)
#' enorth<-extent(-1.25,358.75,0,91.25)
#' esouth<-extent(-1.25,358.75,-91.25,0)
#' rn<-crop(r,enorth) * 0 + 1
#' rs<-crop(r,esouth) * 0
#' r<-mosaic(rn,rs,fun=mean)
#' tsp(prec, 2010, startday = 152, endday = 243, r)
tsp <- function(prec, year, startday = 152, endday = 243, r) {
  lpfun <- function(year) {
    diy <- ifelse(year%%4 == 0, 366, 365)
    if (year%%100 == 0 & year%%400 != 0) diy<-365
    diy
  }
  diy<-lpfun(year)
  startday <- ifelse(diy>365,startday+1,startday)
  startsth <- startday +  floor(diy/2)
  endday <- ifelse(diy>365,endday+1,endday)
  endsth <- (endday + floor(diy/2))%%diy
  rid <- dim(prec)[3] / diy
  sn <- (startday - 1) * rid + 1
  ss <- (startsth - 1) * rid + 1
  en <- endday * rid
  es <- endsth * rid
  rcs <- en-sn+1
  mn <- getValues(r,format="matrix")
  ms <- mn+1
  ms[ms==2] <-0
  pnth<-prec[,,sn:en]
  psth1 <- prec[,,ss:(dim(prec)[3])]
  psth2 <- prec[,,1:es]
  pnorth <- apply(pnth,c(1,2),sum,na.rm=T)*mn
  psouth <- (apply(psth1,c(1,2),sum,na.rm=T) + apply(psth2,c(1,2),sum,na.rm=T)) * ms
  psummer <- pnorth + psouth
  psummer
}
#' gsl: Growing season length
#' @description `gsl` calculates the length of the growing season in decimal days.
#'
#' @param gseason a three dimensional array of binary values (1 = growing season, 0 = not growing season).
#' @param year calander year.
#'
#' @return a matrix of values of growing season length in decimal days.
#' @export
#'
#' @seealso the `gseason` function can be used to create an array of growing conditions (1) = yes, (0) = no accounting for temperature, precipitation and daylight hours.
#'
#' @examples
#' tme <- as.POSIXlt(c(0:1459) * 3600 * 6, origin = "2010-01-01 00:00", tz = "GMT")
#' gs <- gseason_day(tme, 6, 21)
#' gsa <- array(gs, dim=c(1, 1, 1460))
#' gsl <- gsl(gsa, 2010)
#'
gsl <- function(gseason, year) {
  lpfun <- function(year) {
    diy <- ifelse(year%%4 == 0, 366, 365)
    if (year%%100 == 0 & year%%400 != 0) diy<-365
    diy
  }
  diy<-lpfun(year)
  rid <- dim(gseason)[3] / diy
  gsl <- apply(gseason,c(1,2),sum,na.rm=T) / rid
  gsl
}
#' gsmax: Maximum temperature during the growing season
#'
#' @description Calculates maximum air temperature during the growing season.
#'
#' @param temp a three dimensional array of temperature values.
#' @param gseason a three dimensional array of binary values (1 = growing season, 0 = not growing season).
#'
#' @return a matrix of values of maximum growing season temperatures.
#' @export
#'
#' @seealso the [gseason()] function can be used to create an array of growing conditions (1) = yes, (0) = no accounting for temperature, precipitation and daylight hours.
#'
#' @examples
#' temp <- array(10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460), dim=c(73,144,1460))
#' tme <- tmecreate(2010, 6)
#' gs <- gseason_day(tme, 6, 21)
#' gseason <- array(gs, dim=c(73, 144, 1460))
#' maxgst <- gsmax(temp, gseason)
gsmax <- function(temp, gseason) {
  gseason[gseason == 0] <- NA
  gstemp<-temp*gseason
  gstemp[is.na(gstemp)==TRUE] <- 0
  gsmax<-apply(gstemp, c(1,2), max, na.rm=TRUE)
  gsmax
}
#' mst: Mean summer temperature
#'
#' @description `mst` calculates the mean temperature during summer, accounting for differences in the timing of the summer period for the northern and southern hemispheres.
#'
#' @param temp a three dimensional array of temperature values (deg C).
#' @param year calendar year.
#' @param startday assumed day of year of start of summer in northen hemisphere in non-leap year.
#' @param endday assumed end of summer. Defaults are 1st June to 31st Aug.
#' @param r a raster of same extent as temp coded as 1 for northern hemisphere and 0 for southern hemisphere.
#'
#' @return a matrix of mean summer temperature values.
#' @export
#'
#' @importFrom raster getValues
#'
#' @details
#' Seasons are flipped in the southern hemisphere. I.e. 1st June (day 152) = day 152+365/2+0.5 = 334 = 1st Dec.
#' in leap years, 1 day added. Startday and endday should be for northern hemisphere, and are calculated for southern hemisphere within the function.
#'
#' @examples
#' temp <- array(10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460), dim=c(73,144,1460))
#' m <- matrix(1, 73, 144)
#' r <- raster(m, crs="+init=epsg:4326")
#' extent(r) <- c(-1.25, 358.75, -91.25, 91.25)
#' enorth<-extent(-1.25,358.75,0,91.25)
#' esouth<-extent(-1.25,358.75,-91.25,0)
#' rn<-crop(r,enorth) * 0 + 1
#' rs<-crop(r,esouth) * 0
#' r<-mosaic(rn,rs,fun=mean)
#' mst(temp, 2010, startday = 152, endday = 243, r)
#'
mst <- function(temp, year, startday = 152, endday = 243, r) {
  lpfun <- function(year) {
    diy <- ifelse(year%%4 == 0, 366, 365)
    if (year%%100 == 0 & year%%400 != 0) diy<-365
    diy
  }
  diy<-lpfun(year)
  startday <- ifelse(diy>365,startday+1,startday)
  startsth <- startday +  floor(diy/2)
  endday <- ifelse(diy>365,endday+1,endday)
  endsth <- (endday + floor(diy/2))%%diy
  rid <- dim(temp)[3] / diy
  sn <- (startday - 1) * rid + 1
  ss <- (startsth - 1) * rid + 1
  en <- endday * rid
  es <- endsth * rid
  rcs <- en-sn+1
  mn<-getValues(r,format="matrix")
  ms <- mn+1
  ms[ms==2] <-0
  tnth<-temp[,,sn:en]
  tsth1 <- temp[,,ss:(dim(temp)[3])]
  tsth2 <- temp[,,1:es]
  tnorth <- (apply(tnth,c(1,2),sum,na.rm=T)/rcs)*mn
  tsouth <- ((apply(tsth1,c(1,2),sum,na.rm=T) + apply(tsth2,c(1,2),sum,na.rm=T))
             / rcs) * ms
  tsummer <- tnorth + tsouth
  tsummer
}
#' ssm: Summer soil moisture content
#'
#' @description `ssm` calculates average soil moisture content over the summer period.
#'
#' @param soilm a three dimensional array of soil moisture values for one year.
#' @param year calendar year.
#' @param startday assumed day of year of start of summer in northen hemisphere in non-leap year.
#' @param endday assumed end of summer. Defaults are 1st June to 31st Aug.
#' @param r a raster of same extent as temp coded as 1 for northern hemisphere and 0 for southern hemisphere.
#'
#' @importFrom raster getValues
#'
#' @return a matrix of mean summer soil moisture values.
#' @export
#'
#' @details Seasons are flipped in southern hemisphere. I.e. 1st June (day 152) = day 152+365/2+0.5 = 334 = 1st Dec.
#' in leap years, 1 day added. Startday and endday should be for northern hemisphere, and are calculated for southern hemisphere within the function.
#'
#' @examples
#' soilm <- array(runif(1460, 0, 65), dim= c(73,144,1460))
#' m <- matrix(1, 73, 144)
#' r <- raster(m, crs="+init=epsg:4326")
#' extent(r) <- c(-1.25, 358.75, -91.25, 91.25)
#' enorth<-extent(-1.25,358.75,0,91.25)
#' esouth<-extent(-1.25,358.75,-91.25,0)
#' rn<-crop(r,enorth) * 0 + 1
#' rs<-crop(r,esouth) * 0
#' r<-mosaic(rn,rs,fun=mean)
#' ssm(soilm, 2010, startday = 152, endday = 243, r)
#'
ssm <- function(soilm, year, startday = 152, endday = 243, r) {
  lpfun <- function(year) {
    diy <- ifelse(year%%4 == 0, 366, 365)
    if (year%%100 == 0 & year%%400 != 0) diy<-365
    diy
  }
  diy<-lpfun(year)
  startday <- ifelse(diy>365,startday+1,startday)
  startsth <- startday +  floor(diy/2)
  endday <- ifelse(diy>365,endday+1,endday)
  endsth <- (endday + floor(diy/2))%%diy
  rid <- dim(soilm)[3] / diy
  sn <- (startday - 1) * rid + 1
  ss <- (startsth - 1) * rid + 1
  en <- endday * rid
  es <- endsth * rid
  rcs <- en-sn+1
  mn<-getValues(r,format="matrix")
  ms <- mn+1
  ms[ms==2] <-0
  pnth<-soilm[,,sn:en]
  psth1 <- soilm[,,ss:(dim(soilm)[3])]
  psth2 <- soilm[,,1:es]
  smnorth <- apply(pnth,c(1,2),mean,na.rm=T)*mn
  smsouth <- (apply(psth1,c(1,2),mean,na.rm=T)) + (apply(psth2,c(1,2),mean,na.rm=T)) * ms
  smsummer <- smnorth + smsouth
  smsummer
}
#' gseason_soilmoist: Growing season soil moisture content
#'
#' @description `gseason_soilmoist` calculates the average soil moisture content during the growing season over specified years.
#'
#' @param startyear earliest calender year to be considered in calculations.
#' @param endyear latest calender year to be considered in calculations.
#' @param ga array of growing binary values indicating growing (1) or non-growing (0) season.
#' @param gs array of growing binary values indicating growing (1) or non-growing (0) season for one year.
#' @param sm array of soil moisture values for one year.
#'
#' @import raster
#' @import ncdf4
#'
#' @return A matrix of mean growing season soil moisture values over specified years.
#' @export
#'
#' @details
#'
#' @seealso Requires function [gssm()] to be loaded.
#' @seealso the [gseason()] function can be used to create an array of growing conditions (1) = yes, (0) = no accounting for temperature, precipitation and daylight hours.
#' @seealso [mtoraster()]
#' @seealso Requires that function [gssm()] is also loaded.
#'
#' @examples
#' require(microclima)
#' tme <- as.POSIXlt(c(0:1459) * 3600 * 6, origin = "2010-01-01 00:00", tz = "GMT")
#' gs <- gseason_day(tme, 6, 21)
#' gseason <- matrix(gs, dim=c(1, 1, 1460))
#' soilm <- array(runif(1460, 0, 100), dim= c(1,1,1460))
#' gseason_soilmoist(2010, 2010, gs, soilm)
#'
gseason_soilmoist <- function(startyear, endyear, gs, sm){
  dim3 <- endyear - startyear + 1
  ysoil <- array(NA, dim = c(dim(sm)[1:2], dim3))
  i <- 1
  for (year in startyear:endyear) {
    print(year)
    ysoil[,,i] <- gssm(gs, sm)
    i <- i +1
  }
  soilm <- apply(ysoil, c(1,2), mean, na.rm=T)
  soilm
}
#'gseastemp: Mean temperature during the growing season
#'
#' @description `gseastemp` calculates the mean air temperature during the growing season over specified years.
#'
#' @param startyear earliest calender year to be considered in calculations.
#' @param endyear latest calender year to be considered in calculations.
#' @param temp array of temperature values for one year.
#' @param gseason array of binary values indicating growing (1) or non-growing (0) season for one year.
#'
#' @import raster
#' @import ncdf4
#'
#' @return A matrix of mean growing season temperature values over the specified years.
#' @export
#'
#' @details details to be added
#'
#' @seealso the [gseason()] function can be used to create an array of growing conditions (1) = yes, (0) = no accounting for temperature, precipitation and daylight hours.
#' @seealso Requires function [gst()] to be loaded.
#' @seealso [mtoraster()] to convert output matrix to a raster.
#' @seealso [nctarray()] to create an array from nc file.
#'
#'
gseastemp <- function(startyear, endyear, temp, gseason) {
  dim3 <- endyear - startyear + 1
  styear <- array(NA, dim = c(dim(temp)[1:2], dim3))
  i <- 1
  for (year in startyear:endyear) {
    print(year)
    styear[,,i] <- gst(temp, gseason)
    i <- i+1
  }
  gstemp<-apply(styear,c(1,2),mean,na.rm=T)
  gstemp
}
#'gseasonprecip: Total precipitation during the growing season
#'
#' @description `gseasonprecip` calculates the total precipitation during the growing season.
#'
#' @param startyear earliest calender year to be considered in calculations.
#' @param endyear latest calender year to be considered in calculations.
#' @param precip 3-dimensional array of precipitation values.
#' @param gseason 3-dimensional array of values indicating growing season (1) or non-growing season (0).
#'
#' @import raster
#' @import ncdf4
#'
#' @return A matrix of mean total precipitation values during the growing season over the specified years.
#' @export
#'
#' @details details to be added
#'
#' @seealso Requires function `gsp` to be loaded.
#' @seealso the `gseason` function can be used to create an array of growing conditions (1) = yes, (0) = no accounting for temperature, precipitation and daylight hours.
#' @seealso [mtoraster()] to convert output matrix to a raster.
#' @seealso [nctarray()] to create an array from an nc file.
#'
gseasprec <- function(startyear, endyear, precip, gseason) {
  dim3 <- endyear - startyear + 1
  styear <- array(NA, dim = c(dim(precip)[1:2], dim3))
  i <- 1
  for (year in startyear:endyear) {
    print(year)
    styear[,,i] <- gsp(prec, gseason)
    i <- i+1
  }
  gsprec<-apply(styear,c(1,2),mean,na.rm=T)
  gsprec
}
#' summerprecip: Total Summer precipitation
#'
#' @description `summerprecip` calculates the total amount of precipitation falling in summer.
#'
#' @param startyear earliest calender year to be considered in calculations.
#' @param endyear latest calender year to be considered in calculations.
#' @param precipnc full path name of nc file containing precipitation values  and with data extent: -1.25, 358.75, -91.25, 91.25 when converted to raster format.
#' @param startday assumed day of year of start of summer in northen hemisphere in non-leap year* (default 1st June).
#' @param endday assumed day of year of end of summer in northen hemisphere in non-leap year* (default 31st Aug).
#'
#' @import raster
#' @import ncdf4
#'
#' @return A matrix of mean growing degree days above 10 degrees Celcius during the growing season over the specified years.
#' @export
#'
#' @details Function has been designed to run with raster data of dimensions 73 x 144.
#' Seasons are flipped in southern hemisphere. I.e. 1st June (day 152) = day 152+365/2+0.5 = 334 = 1st Dec.
#' in leap years, 1 day added. Startday and endday should be for northern hemisphere, and are calculated for southern hemisphere within the function.
#'
#' @seealso Requires function [tsp()] to be loaded.
#' @seealso [mtoraster()] to convert output matrix to a raster.
#' @seealso [nctarray()] to create an array from an nc file.
#'
summerprecip <- function(startyear, endyear, precipnc, startday = 152, endday = 243) {
  dim3 <- endyear - startyear + 1
  styear <- array(NA, dim = c(73, 144, dim3))
  r<-raster(precipnc)
  enorth<-extent(-1.25,358.75,0,91.25)
  esouth<-extent(-1.25,358.75,-91.25,0)
  rn<-crop(r,enorth) * 0 + 1
  rs<-crop(r,esouth) * 0
  r<-mosaic(rn,rs,fun=mean)
  i<-1
  for (year in startyear:endyear) {
    print(year)
    prec <- nctoarray(precipnc)
    prec[prec<0] <- 0
    styear[,,i] <- tsp(prec, year, startday = startday, endday = endday, r)
    i <- i +1
  }
  sumprec<-apply(styear,c(1,2),mean,na.rm=T)
  sumprec
}
#' gseasonlength: Growing season length
#'
#'@description `gseasonlength` calculates the average length of the growing season period (in days) across specified years.
#'
#' @param startyear earliest calender year to be considered in calculations.
#' @param endyear latest calender year to be considered in calculations.
#' @param gseason a 3-dimensional array of growing season values for each year (1 = growing season, 0 = non-growing season).
#'
#' @import raster
#' @import ncdf4
#'
#' @return A matrix of mean growing season length (number of days) over specified years.
#' @export
#'
#' @details Some details required here
#'
#' @seealso the [gseason()] function can be used to create an array of growing conditions (1) = yes, (0) = no accounting for temperature, precipitation and daylight hours.
#' @seealso [mtoraster()] to convert a matrix to a raster object.
#' @seealso [nctarray()] to create array from an nc file.
#' @seealso Requires function `gsl` to be loaded.
#'
gseasonlength <- function(startyear, endyear, gseason)  {
  dim3 <- endyear - startyear + 1
  styear <- array(NA, dim = c(dim(gseason)[1:2], dim3))
  i <- 1
  for (year in startyear:endyear) {
    print(year)
    styear[,,1] <- gsl(gseason, year)
    i <- i +1
  }
  gsl <- apply(styear,c(1,2),mean,na.rm=T)
  gsl
}
#'maxgstemp: Maximum temperature during the growing season
#'
#' @description Calculates the maximum air temperature during the growing season across specified years.
#'
#' @param startyear earliest calender year to be considered in calculations.
#' @param endyear latest calender year to be considered in calculations.
#' @param temp a three-dimensional array of temperature values for one year.
#' @param gs a three-dimensional array of binary values indicating growing (1) or non-growing (0) season for one year.
#'
#' @return A matrix of mean maximum growing season air temperature values over the specified years.
#' @export
#'
#' @import raster
#' @import ncdf4
#'
#' @details some details needed here
#'
#' @seealso [mtoraster()] to convert a matrix to a raster object.
#' @seealso [nctarray()] to create array of temperature values from an nc file.
#' @seealso Requires function [gsmax()] to be loaded.
#'
maxgstemp <- function(startyear, endyear, temp, gs) {
  dim3 <- endyear - startyear + 1
  styear <- array(NA, dim = c(dim(temp)[1:2], dim3))
  i <- 1
  for (year in startyear:endyear) {
    print(year)
    styear[,,i] <- gsmax(temp, gs)
    i <- i+1
  }
  gsmax<-apply(styear,c(1,2),mean,na.rm=T)
  gsmax
}
#' summertemp: Mean summer temperature
#'
#' @description `summertemp` calculates mean summer temperature across specified time period (years).
#'
#' @param startyear the earliest calendar year (AD) to be considered in the
#' calculation.
#' @param endyear the latest calendar year (AD) to be considered in the
#' calculation.
#' @param tempnc full path name of nc file containing temperature values for each year and with data extent: -1.25, 358.75, -91.25, 91.25 when converted to raster format.
#' @param startday Indicates assumed day of year of start of summer in northern
#' hemisphere in non-leap year. Default, 1st June.
#' @param endday Indicates assumed day of year of end of summer in northern
#' hemisphere in non-leap year. Default, 31st August.
#'
#' @import raster
#' @import ncdf4
#'
#' @return a matrix of mean summer temperature values over the specified years.
#' @export
#'
#' @details Function has been designed to run with raster data of dimensions 73 x 144.
#' Seasons are flipped in southern hemisphere. I.e. 1st June (day 152) = day 152+365/2+0.5 = 334 = 1st Dec.
#' in leap years, 1 day added. Startday and endday should be for northern hemisphere, and are calculated for southern hemisphere within the function.
#'
#' @seealso [mtoraster()]
#' @seealso [nctarray()] to create array of temperature values from an nc file.
#' @seealso Requires function [mst()] to be loaded.
summertemp <- function(startyear, endyear, tempnc, startday = 152, endday = 243) {
  dim3 <- endyear - startyear + 1
  styear <- array(NA, dim = c(73, 144, dim3))
  i <- 1
  r<-raster(tempnc)
  enorth<-extent(-1.25,358.75,0,91.25)
  esouth<-extent(-1.25,358.75,-91.25,0)
  rn<-crop(r,enorth) * 0 + 1
  rs<-crop(r,esouth) * 0
  r<-mosaic(rn,rs,fun=mean)
  for (year in startyear:endyear) {
    print(year)
    temp <- nctoarray(tempnc)
    styear[,,i] <- mst(temp, year, startday = startday, endday = endday, r)
    i <- i +1
  }
  sumtemp<-apply(styear,c(1,2),mean,na.rm=T)
  sumtemp
}
#' summersoilmoist: Soil moisture content during summer
#'
#' @description Calculates mean soil moisture content during the 6-month summer period.
#'
#' @param startyear earliest calender year to be considered in calculations.
#' @param endyear latest calender year to be considered in calculations.
#' @param soilnc full path name of nc file containing soil moisture values for each year and with data extent: -1.25, 358.75, -91.25, 91.25 when converted to raster format.
#' @param fi nc file containing soil moisture values for one year.
#' @param startday Indicates assumed day of year of start of summer in northern
#' hemisphere in non-leap year. Default, 1st June.
#' @param endday indicates assumed day of year of end of summer in northern
#' hemisphere in non-leap year. Default, 31st August.
#'
#' @import raster
#' @import ncdf4
#'
#' @return a matrix of mean summer soil moisture content values over the specified years.
#' @export
#'
#' @details Function has been designed to run with raster data of dimensions 73 x 144.
#' Seasons are flipped in southern hemisphere. I.e. 1st June (day 152) = day 152+365/2+0.5 = 334 = 1st Dec.
#' in leap years, 1 day added. Startday and endday should be for northern hemisphere, and are calculated for southern hemisphere within the function.
#'
#' @seealso [mtoraster()]
#' @seealso [nctoarray()] to create array of temperature values from an nc file.
#' @seealso Requires function [ssm()] to be loaded.
#'
summersoilmoist <- function(startyear, endyear, soilnc, fi, startday = 152, endday = 243) {
  dim3 <- endyear - startyear + 1
  styear <- array(NA, dim = c(73, 144, dim3))
  r<-raster(soilnc)
  enorth<-extent(-1.25,358.75,0,91.25)
  esouth<-extent(-1.25,358.75,-91.25,0)
  rn<-crop(r,enorth) * 0 + 1
  rs<-crop(r,esouth) * 0
  r<-mosaic(rn,rs,fun=mean)
  i<-1
  for (year in startyear:endyear) {
    print(year)
    soilm <- nctoarray(fi)
    soilm[soilm<0] <- 0
    soilm[is.na(soilm)==TRUE]<-0
    styear[,,i] <- ssm(soilm, year, startday = startday, endday = endday, r)
    i <- i +1
  }
  meanmoist<-apply(styear,c(1,2),mean,na.rm=T)
  meanmoist
}
