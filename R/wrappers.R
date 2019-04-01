#'sisimple: Calculates the solar index
#'
#'@description Calculates the proportion of direct beam radiation incident on an inclined surface at a specified time and location.
#'
#'@param locatime local time (decimal hour, 24 hour clock).
#'@param lat latitude of the location for which the solar index is required (decimal degrees, -ve south of the equator).
#'@param long longitude of the location for which the solar index is required (decimal degrees, -ve west of Greenwich meridian).
#'@param julian Julian day expressed as an integer as returned by julday().
#'@param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is 0.
#'@param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).
#'
#'@return a numeric value of the proportion of direct beam radiation on a horizontal surface at the specified latitude and longitude.
#'@export
#'
#'@seealso the microclima function `julday` can be used to derive `julian`.
#'
#'@examples
#'jd <- julday (2010, 6, 21)
#'si <- sisimple(11, 50, -5, jd, merid=0, dst=0)
#'
sisimple <- function(localtime, lat, long, julian, merid = 0, dst = 0) {
  saltitude <- solalt(localtime, lat, long, julian, merid, dst)
  alt <- saltitude * (pi / 180)
  index <- cos(pi / 2 - alt)
  index[index < 0] <- 0
  index
}
#'nctoarray: Creates an array containing data within an nc file
#'
#'@description `nctoarray` is used to create a three-dimensional array
#'containing data from a .nc file.
#'
#'@param filein A character string with the full path name of the nc file. Tilde-expansion is performed.
#'@param varid The variable to be read from the ncfile. Can be a string with the
#'name of the variable or an object of class ncvar4. If left unspecified, the
#'function will determine if there is only one variable in the file and, if
#'so, read from that. If left unspecified and there are multiple variables in
#'the file, an error is generated.
#'
#'@import ncdf4
#'
#'@return a three-dimensional array of the specified variable, permuted to dimension order 2, 1, 3.
#'@export
#'
#'@examples
#'# ========= Create nc file of random data and save to disk ============= #
#' longs <- ncdim_def(name="lon", units = "degrees_east",
#'                    vals = seq(1,10,1),longname="longitude")
#' lats <- ncdim_def(name="lat",  units = "degrees_north",vals=seq(1,10,1),
#'                   longname="latitude")
#' times <- ncdim_def(name="time",units="hours since 1983-01-01 00:00:00",
#'                    vals=c(101:110))
#' mydata <-ncvar_def(name="Random data", units="mm", dim = list(longs, lats, times))
#' ncnew <- nc_create(filename="Random.nc", mydata)
#' ncvar_put(ncnew, mydata, vals = array(rnorm(10000), dim = c(10,10,10)))
#' nc_close(ncnew)
#'
#'# ========= Read in as array ============================= #
#'arr <- nctoarray("Random.nc", varid = "Random data")
#'dim(arr)
nctoarray <- function(filein, varid = NA) {
  nc <- nc_open(filein)
  a <- aperm(ncvar_get(nc, varid = varid), c(2,1,3))
  nc_close(nc)
  a
}
#'resamplearray: resamples an array to the dimensions of another array
#'
#'@description reformats data in an array and creates a new array objects to
#'  dimensions specified by `rin`.
#'
#'@param a the array object to be resampled.
#'@param rin a raster object used to resample extent of 'a' to extent of 'rout'.
#'@param rout a raster object with extent parameters that 'a' should be resampled to.
#'
#'@import raster
#'
#'@return an array with the same latitude and longitude as 'rout'.
#'@export
#'
#'@examples
#'a <- array(rnorm(100), dim=c(94,192,1464))
#'rin <- raster(a[,,1])
#'extent(rin) <- c(-0.9375, 359.0625, -89.49406, 89.49406)
#'routm <- matrix(rnorm(10512), 73, 144)
#'rout <- raster(routm)
#'extent(rout) <- c(-1.25, 358.75, -91.25, 91.25)
#'newarray <- resamplearray(a, rin, rout)
#'dim(newarray)
#'
resamplearray <- function(a, rin, rout) {
  b <- brick(a)
  rin <- extend(rin, extent(rout))
  extent(b) <- extent(rin)
  b <- resample(b, rout)
  ao <- array(NA, dim = dim(b))
  for (i in 1:dim(b)[3])
    ao[,,i] <- getValues(raster(b, i), format = "matrix")
  ao
}
#'tmecreate: Creates a `POSIXlt` object representing calendar dates and times
#'
#'@description `tmecreate` is used to create a `POSIXlt` object representing
#'calendar dates and times for the years specified by `years` and for time
#'intervals specified by `hourint`.
#'
#'@param years a vector of years.
#'@param hourint a single numeric value representing the time interval in hours
#'   between dates and times.
#'
#'@return an object of class `POSIXlt` with calander dates and times.
#'@export
#'
#'@examples
#'head(tmecreate(2010, 24))
#'tail(tmecreate(2010, 6))
tmecreate<- function(years, hourint = 6) {
  lpfun <- function(year) {
    diy <- ifelse(year%%4 == 0, 366, 365)
    if (year%%100 == 0 & year%%400 != 0) diy<-365
    diy
  }
  diy <- sapply(years, lpfun)
  i <- sum(diy) * 24 / hourint -1
  tme <- as.POSIXlt(c(0:i) * 3600 * hourint, origin = paste0(min(years), "-01-01 00:00"), tz = "GMT")
  tme
}
#'arraytonc: Create an nc file from an array
#'
#'@description `arraytonc` creates a new netCDF file from an array object.
#'
#'@param a the array containing the data to be used to create the new netCDF file.
#'@param fileout the name and location of the nteCDF file to be created.
#'@param varname the name of the variable to be created (character string).
#'@param units The variable's units (character string). Or, pass a zero length string (") to have no units attribute.
#'@param r a raster file with extent parameters that the new netCDF file should follow.
#'@param tme an object of class `POSIXlt` representing calendar dates and times.
#'@param baseyear the calendar year that measurements began.
#'
#'@import ncdf4
#'@import raster
#'
#'@return an object of class ncdf4 that has the fields described above and is saved to the location specified by 'fileout'.
#'@export
#'
#'@seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#'temps <- array(10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460), dim=c(180,360,1460))
#'r <- raster()
#'tme <- tmecreate(2010, 6)
#'arraytonc (temps, "new.nc", varname = "6-hourly air temperature", units = "degrees Celcius", r, tme, baseyear = 1900)
#'
arraytonc <- function(a, fileout, varname, units, r, tme, baseyear = 1900) {
  yrs <- c(baseyear:1969)
  hiy <- ifelse(yrs%%4 == 0, 366,365)
  hiy <- ifelse(yrs%%100 == 0 & yrs%%400 != 0, 365, hiy)
  tout <- as.numeric(tme) / 3600 + sum(hiy)
  e <- extent(r)
  ao <- aperm(a, c(2,1,3))
  stplat <- (e@ymax - e@ymin)/ dim(r)[1]
  stplon <- (e@xmax - e@xmin)/ dim(r)[2]
  latseq <- rev(seq((e@ymin + 0.5 * stplat), (e@ymax - 0.5 * stplat), stplat))
  lonseq <- seq((e@xmin + 0.5 * stplon), (e@xmax - 0.5 * stplon), stplon)
  long <- ncdim_def(name = "lon", units = "degrees_east",
                    vals = lonseq, longname = "longitude")
  lat<-ncdim_def(name = "lat", units = "degrees_north",
                 vals = latseq, longname = "latitude")
  ntme<-ncdim_def(name = "time",
                  units = paste0("hours since ", baseyear, "-01-01 00:00:00"),
                  vals = tout)
  varn <- ncvar_def(varname, units = units, dim=list(long,lat,ntme))
  ncnew <- nc_create(filename = fileout, varn)
  ncvar_put(ncnew, varn, vals = ao)
  nc_close(ncnew)
}
#'submonthly'
#'
#'@description `submonthly` is used to create a list representing each measurement per half-month.
#'
#'@param tme a POSIXlt object representing calendar dates and times.
#'@param div numeric object specifying the value by which to divide total days of each month. Default = 2.
#'
#'@return a vector of numeric values alternating between 0 and 1.
#'@export
#'
#'@details If div=2, value changes after the last 'tme' measurement in each half-month.
#'
#'@seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#'@examples
#'tme <- tmecreate(2010, 6)
#'sm <- submonthly(tme, div = 2)
#'sm2 <- submonthly(tme, div = 4)
submonthly <- function(tme, div = 2) {
  bwm <-0
  for (i in 1:12) {
    dm <-max(tme$mday[which(tme$mon == (i - 1))]) #max days in each month (e.g. Nov = 31). will account for leap years
    bw <- floor((tme$mday[which(tme$mon == (i - 1))] - 1)/ dm * div)
    bwm <- c(bwm, bw)
  }
  bwm[-1]
}
#'cumsumseq
#'#'@description `cumsumseq` is used to create a vector of values indicating growing season conditions for the specified time period.
#'
#'@param gs a vector of binary values indicating growing season (1) or non-growing season (0).
#'@param ehour a vector of evapotranspiration values.
#'@param div numeric object specifying the value by which to divide total days of each month.
#'
#'@return a vector of numeric values alternating between 0 and 1.
#'@export
#'
#'@examples
#'gs <- c(0,0,1,1,1,1,0,0,1,1,0)
#'ehour <- c(80, 81, 82, 83, 84, 85, 84, 83, 82, 81, 80)
#'div = 2
#'cs <- cumsumseq(gs, ehour, div = 2)
#'
cumsumseq <- function(gs, ehour, div) {
  sel <- which(gs < 1)
  sq <- 1
  if (length(sel) > 1) {
    for (i in 2:length(sel))
      sq[i] <- ifelse(sel[i] == (sel[i-1] + 1), 0 ,1)
  }
  sel2 <- which(sq == 1)
  sel2 <- sel[sel2]
  for (i in 1:length(sel2)) {
    mx <- ifelse(i < length(sel2), sel2[i+1], length(ehour))
    es <- cumsum(ehour[sel2[i]:mx]) * 12 * div / length(ehour)
    it <- which(es<100) + sel2[i] - 1
    gs[it] <- 1
  }
  gs
}

#'mtoraster: Convert a matrix to raster object
#'
#'@description `mtoraster` converts a matrix of global values to raster format.
#'
#'@param m a matrix of values.
#'@param centre if true, centres values to prime meridian (longitude 0°).
#'@param mask option to mask land/sea from plotted values.
#'
#'@import raster
#'
#'@return a raster object.
#'@export
#'
#'@details the function is written for the conversion of matrices of dimensions 73 x 144.
#'
#'@examples
#'
mtoraster <- function(m, centre = TRUE, mask = NA) {
  if (centre) {
    m <- cbind(m[,73:144], m[,1:72])
    if (class(mask) == "matrix" | class(mask) == "array") {
      mask <- cbind(mask[,73:144], mask[,1:72])
      mask[which(mask == 0)] <- NA
      m <- m * mask
    }
    e <- extent(c(-181.25, 178.75, -91.25, 91.25))
    r <- raster(m)
    extent(r) <- e
  }
  else {
    if (class(mask) == "matrix" | class(mask) == "array") {
      mask <- cbind(mask[,73:144], mask[,1:72])
      mask[which(mask == 0)] <- NA
      m <- m * mask
    }
    e <- extent(c(-1.25, 358.75, -91.25, 91.25))
    r <- raster(m)
    extent(r) <- e
  }
  r
}
#'resamplenc: Reformat an nc file to new dimensions.
#'
#'@description `resamplenc` is used to reformat data an nc file to the same dimensions as the data in `filein1`.
#'
#'@param filein1 master nc file.
#'@param filein2 nc file to undergo reformatting.
#'@param var time dimension to be left unchanged.
#'@param varname the name of the variable to be created (character string).
#'@param units the variable's units (character string). Or, pass a zero length string (") to have no units attributed.
#'@param fileout full path name of the new nc file to be created.
#'
#'@import raster
#'@import ncdf4
#'
#'@return a .nc file with the same dimensions as filein1.
#'@export
#'
#'@details recording period for filein2 is left unchanged.
#'
#'@seealso `nctoarray`
#'
#'@examples to follow once there are two nc files with different latitude and longitude dimensions in the data.
#'
resamplenc<-function(filein1, filein2, varname, units, fileout) {
  r1 <-raster(filein1, band=1)
  r2 <- raster(filein2, band = 1)
  nc <- nc_open(filein2)
  tme <- ncvar_get(nc, var = "time")
  nc_close(nc)
  a <- nctoarray(filein2)
  ao <- array(NA, dim = c(dim(r1)[1:2], dim(a)[3]))
  for (i in 1:dim(a)[3]) {
    m <- a[,,i]
    rin <- raster(m, template = r2)
    rou <- resample(rin, r1)
    ao[,,i] <- getValues(rou, format = "matrix")
  }
  ao2 <- aperm(ao,c(2,1,3))
  e <- extent(r1)
  stplat <- (e@ymax - e@ymin)/ dim(r1)[1]
  stplon <- (e@xmax - e@xmin)/ dim(r1)[2]
  latseq <- seq((e@ymin + 0.5 * stplat), (e@ymax - 0.5 * stplat), stplat)
  latseq <- rev(latseq)
  lonseq <- seq((e@xmin + 0.5 * stplon), (e@xmax - 0.5 * stplon), stplon)
  long <- ncdim_def(name="lon",units="degrees_east",vals=lonseq,longname="longitude")
  lat <- ncdim_def(name="lat",units="degrees_north",vals=latseq,longname="latitude")
  ntme <- ncdim_def(name="time",units="hours since 1900-01-01 00:00:00",vals=tme)
  varinfo <- ncvar_def(varname, units, dim=list(long,lat,ntme))
  ncnew <- nc_create(filename="fileout",varinfo)
  ncvar_put(ncnew,varinfo,vals=ao2)
  nc_close(ncnew)
}
#'warna
#'@description Prints a warning message if data span more than one year and indicates method followed.
#'
#'@return If true, returns warning message, "Data spans more than one year. Data aggregated by unique month irrespective of year and one value returned".
#'@export
#'
warna <- function() {
  warning ("Data spans more than one year. Data aggregated by unique month
           irrespective of year and one value returned")
}

#'warnb
#'@description Prints a warning message if data span more than one year and indicates method followed.
#'
#'@return If true, returns warning message, "Data spans more than one year. Calculations performed on all data and single value returned".
#'@export
#'
warnb <- function() {
  warning ("Data spans more than one year. Calculations performed on all data
           and single value returned")
}
