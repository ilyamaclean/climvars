#' bio1: Calculates mean annual temperature
#'
#' @description `bio1` is used to calculated mean annual temperature.
#'
#' @param temps a vector of temperatures, normally for one year (see details).
#' @param tme a `POSIXlt` object representing the date and time of each `temps` value.
#' Ignored if method unspecified.
#' @param method an optional character string describing the method used to
#' calculate mean annual temperature. Options are "anuclim", "dailymaxmin" or unspecified
#' (see details).
#' @return a single numeric value of mean annual temperature.
#' @export
#'
#' @details If `method` is "anuclim" temperatures are aggregated by month, and
#' spline intepolated to weekly before mean annual temperature is calculated,
#' replicating the method used by http://www.worldclim.org/. If `method` is
#' "dailymaxmin", daily mean temperatures are calculated as the mean of daily
#' maximum and minimum temperatures and annual mean calculated from daily
#' means. Otherwise the mean of `temps` is returned. If using `anuclim` method
#' and data span more than one year, data are aggregated by unique month
#' irrespective of year and one value returned. If using `dailymaxmin` method
#' and data span more than one year, calculations will be performed on all data
#' and a single value returned.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' temps <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' tme <- tmecreate(2010, 6)
#' plot(temps~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Temperature")
#' bio1(temps, tme)
#' bio1(temps, tme, method = "anuclim")
#' bio1(temps, tme, method = "dailymaxmin")
#'
bio1 <- function(temps, tme, method = "") {
  if (is.na(sd(temps, na.rm = TRUE)))
    tmean <- NA
  else {
    if (method == "anuclim") {
      if (length(unique(tme$year)) > 1) warna()
      tmth <- aggregate(temps, by = list(tme$mon), mean, na.rm = TRUE)$x
      twk <- spline(tmth, n = length(tmth) / 12 * 52)$y
      tmean <- mean(twk, na.rm = TRUE)
    } else if (method ==  "dailymaxmin") {
      if (length(unique(tme$year)) > 1) warnb()
      tmx <- aggregate(temps, by = list(tme$yday), max, na.rm = TRUE)$x
      tmn <- aggregate(temps, by = list(tme$yday), min, na.rm = TRUE)$x
      tme <- (tmx + tmn) /2
      tmean <- mean(tme, na.rm = TRUE)
    }
    else {
      if (length(unique(tme$year)) > 1) warna()
      tmean <- mean(temps, na.rm = TRUE)
    }
  }
  tmean
}
#' bio2: Calculates mean annual diurnal temperature range
#'
#' @description `bio2` is used to calculate the mean annual diurnal range in
#' temperature (range mean of the maximum-minimum).
#'
#' @param temps a vector of temperatures, normally for one year (see details).
#' @param tme a `POSIXlt` object representing the date and time of each `temps` value.
#' @param method an optional character string describing the method used to
#' mean annual diurnal temperature range. Options are "anuclim" or unspecified (see
#' details).
#'
#' @return a single numeric value of mean diurnal temperature range.
#' @export
#'
#' @details If using `anuclim` method and data span more than one year, data are
#' aggregated by unique month irrespective of year and one value returned. If
#' `method` is "anuclim" temperatures are aggregated by month and spline
#' intepolated to weekly before mean diurnal temperature range is calculated,
#' replicating the method used by http://www.worldclim.org/. If left
#' unspecified and time interval is <= daily, the mean difference between the
#' daily maximum and minimum values is calculated.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' temps <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' tme <- tmecreate(2010, 6)
#' plot(temps~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Temperature")
#' bio2(temps, tme)
#' bio2(temps, tme, method = "anuclim")
#'
bio2 <- function(temps, tme, method = "") {
  if (is.na(sd(temps, na.rm = TRUE)))
    tdtr <- NA
  else {
    if (as.numeric(tme[2]) - as.numeric(tme[1]) > 86400 & method != "anuclim") {
      warning ("time interval > daily. Using anuclim method")
      method <- "anuclim"
    }
    if (method == "anuclim") {
      if (length(unique(tme$year)) > 1) warna()
      tmthmx <- aggregate(temps, by = list(tme$mon), max, na.rm = TRUE)$x
      tmthmn <- aggregate(temps, by = list(tme$mon), min, na.rm = TRUE)$x
      twkmx <- spline(tmthmx, n = length(tmthmx) / 12 * 52)$y
      twkmn <- spline(tmthmn, n = length(tmthmn) / 12 * 52)$y
      tdtr <- mean((twkmx - twkmn), na.rm = TRUE)
    }
    else {
      if (length(unique(tme$year)) > 1) warnb()
      tmx <- aggregate(temps, by = list(tme$yday), max, na.rm = TRUE)$x
      tmn <- aggregate(temps, by = list(tme$yday), min, na.rm = TRUE)$x
      dtr <- tmx - tmn
      tdtr <- mean(dtr, na.rm = TRUE)
    }
  }
  tdtr
}
#' bio4: Calculates temperature seasonality
#'
#' @description `bio4` calculates the variation in temperature over a given year
#' as the coeeficient of variation in the mean temperatures
#'
#' @param temps a vector of temperatures, normally for one year (see details).
#' @param tme a `POSIXlt` object representing the date and time of each `temps` value.
#' @param method an optional character string describing the method used to
#' calculate temperature seasonality. Options are "anuclim" or unspecified (see
#' details).
#'
#' @return a single numeric value representng annual temperature seasonality.
#' @export
#'
#' @details If method is "anuclim" temperatures are aggregated by month and
#' monthly averages are calculated and spline intepolated to weekly values.
#' Temperature seasonality is calculated as the standard deviation of weekly
#' mean temperatures as a percentage of the mean of those temperatures. If
#' using "anuclim" method and data span more than one year, data are aggregated
#' by unique month irrespective of year and one value returned. If method is
#' not specified, calculation is based on the standard deviation of the mean of
#' all temperatures as a percentage of the mean of all temperatures. If method
#' is not specified and data span more than one year, calculations will be
#' performed on all data and a single value returned. For all calculations, the
#' mean in degrees Kelvin is used.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' temps <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' tme <- tmecreate(2010, 6)
#' plot(temps~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Temperature")
#' bio4(temps, tme)
#' bio4(temps, tme, method = "anuclim")
bio4 <- function(temps, tme, method = "") {
  if (is.na(sd(temps, na.rm = TRUE)))
    tcv <- NA
  else {
    if (method == "anuclim") {
      if (length(unique(tme$year)) > 1) warna()
      tmth <- aggregate(temps, by = list(tme$mon), mean, na.rm = TRUE)$x
      twk <- spline(tmth, n = length(tmth) / 12 * 52)$y
      tcv <- sd(twk, na.rm = TRUE) / mean(twk, na.rm = TRUE)*100
    }
    else {
      if (length(unique(tme$year)) > 1) warnb()
      tcv <- sd(temps, na.rm = TRUE) / mean(temps, na.rm = TRUE)*100
    }
  }
  tcv
}
#' bio5: Calculates maximum temperature of the warmest period of the year
#'
#' @description `bio5` is used to calculate the maximum weekly or monthly temperature in each year
#'
#' @param temps a vector of temperatures, normally for one year (see details).
#' @param tme a `POSIXlt` object representing the date and time of each `temps` value.
#' @param method an optional character string describing the method to calculate
#' maximum temperature in the warmest period. Options are "anuclimmean",
#' "anuclimmax" or unspecified (see details).
#'
#' @return a single numeric value of maximum temperature (weekly or monthly) for each year.
#' @export
#'
#' @details If method is "anuclimmean", mean monthly temperatures are spline
#' interpolated to a weekly time period and the maximum weekly value returned.
#' If method is "anuclimmax", monthly maximum temperatures are spline
#' interpolated to a weekly time period and the maximum weekly value returned.
#' If method is unspecified, the maximum temperature across all values for each year is returned.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' temps <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' tme <- tmecreate(2010, 6)
#' plot(temps~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Temperature")
#' bio5(temps, tme)
#' bio5(temps, tme, method = "anuclimmean")
#' bio5(temps, tme, method = "anuclimmax")
#'
bio5 <- function(temps, tme, method = "") {
  if (is.na(sd(temps, na.rm = TRUE)))
    tmx <- NA
  else {
    if (method == "anuclimmean") {
      if (length(unique(tme$year)) > 1) warna()
      tmth <- aggregate(temps, by = list(tme$mon), mean, na.rm = TRUE)$x
      twk <- spline(tmth, n = length(tmth) / 12 * 52)$y
      tmx <- max(twk, na.rm = TRUE)
    } else if (method == "anuclimmax") {
      if (length(unique(tme$year)) > 1) warna()
      tmth <- aggregate(temps, by = list(tme$mon), max, na.rm = TRUE)$x
      twk <- spline(tmth, n = length(tmth) / 12 * 52)$y
      tmx <- max(twk, na.rm = TRUE)
    }
    else {
      if (length(unique(tme$year)) > 1) warnb()
      tmx <- max(temps, na.rm = TRUE)
    }
  }
  tmx
}
#' bio6: Calculates minimum temperature of the coldest period of the year
#'
#' @description `bio6` is used to calculate the minimum temperature value across
#' all months of the year.
#'
#' @param temps a vector of temperatures, normally for one year (see details).
#' @param tme a `POSIXlt` object representing the date and time of each `temps` value.
#' @param method an optional character string describing the method to calculate
#' minimum temperature. Options are "anuclimmean", "anuclimmin" or unspecified (see
#' details).
#'
#' @return a single numeric value of minimum temperature for the defiend period.
#' @export
#'
#' @details "anuclimmean" splines mean monthly temperature to weekly time period
#' across all months and returns the mean minimum weekly temperature. Anuclimmax
#' splines minimum monthly temperatures to weekly time period and returns the
#' minimum weekly temperature. If left unspecified, the minimum temperature
#' across all values of each year is returned.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' temps <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' tme <- tmecreate(2010, 6)
#' plot(temps~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Temperature")
#' bio6(temps, tme)
#' bio6(temps, tme, method = "anuclimmean")
#' bio6(temps, tme, method = "anuclimmin")
bio6 <- function(temps, tme, method = "") {
  if (is.na(sd(temps, na.rm = TRUE)))
    tmn <- NA
  else {
    if (method == "anuclimmean") {
      if (length(unique(tme$year)) > 1) warna()
      tmth <- aggregate(temps, by = list(tme$mon), mean, na.rm = TRUE)$x
      twk <- spline(tmth, n = length(tmth) / 12 * 52)$y
      tmn <- min(twk, na.rm = TRUE)
    } else if (method == "anuclimmin") {
      if (length(unique(tme$year)) > 1) warna()
      tmth <- aggregate(temps, by = list(tme$mon), min, na.rm = TRUE)$x
      twk <- spline(tmth, n = length(tmth) / 12 * 52)$y
      tmn <- min(twk, na.rm = TRUE)
    }
    else {
      if (length(unique(tme$year)) > 1) warnb()
      tmn <- min(temps, na.rm = TRUE)
    }
  }
  tmn
}
#' bio7: Annual temperature range
#'
#' @description `bio7` is used to calculate the annual range in temperature
#'
#' @param temps a vector of temperatures, normally for one year (see details).
#' @param tme a `POSIXlt` object representing the date and time of each `temps` value.
#' @param method An optional character string describing the method used to
#' calculate maximum and minimum temperatures. Options are "anuclimmean",
#' "anuclimmax" or unspecified (see details).
#'
#' @return a single numeric value of annnual temperature range (maximum-minimum
#' temperature values).
#' @export
#'
#' @details If method is "anuclimmean", mean monthly temperatures are spline
#' interpolated to a weekly time period and range calculated from the maximum
#' and minimum weeekly temperature values. If method is "anuclimmaxmin",
#' maximum and minimum monthly temperatues are spline interpolated to a weekly
#' time period and range is calculated from the maximum and minimum mean weekly
#' temperature values. If left unspecified, the range is calculated by
#' subtracting the maximum and minimum temperature alues for each year. To
#' satisfy the requirements for `tme`, a POSIXlt object can be created using
#' the `tmecreate` wrapper function. This calculation should return the value of
#' bio5(temps, tme)-bio6(temps, tme) when methods remain the same.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' temps <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' tme <- tmecreate(2010, 6)
#' plot(temps~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Temperature")
#' bio7(temps, tme)
#' bio7(temps, tme, method = "anuclimmean")
#' bio7(temps, tme, method = "anuclimmaxmin")
#'
bio7 <- function(temps, tme, method = "") {
  if (is.na(sd(temps, na.rm = TRUE)))
    tanr <- NA
  else {
    if (method == "anuclimmean") {
      if (length(unique(tme$year)) > 1) warna()
      tmth <- aggregate(temps, by = list(tme$mon), mean, na.rm = TRUE)$x
      twk <- spline(tmth, n = length(tmth) / 12 * 52)$y
      tanr <- max(twk, na.rm = TRUE) - min(twk, na.rm = TRUE)
    } else if (method == "anuclimmaxmin") {
      if (length(unique(tme$year)) > 1) warna()
      tmthmx <- aggregate(temps, by = list(tme$mon), max, na.rm = TRUE)$x
      tmthmn <- aggregate(temps, by = list(tme$mon), min, na.rm = TRUE)$x
      twkmx <- spline(tmthmx, n = length(tmth) / 12 * 52)$y
      twkmn <- spline(tmthmn, n = length(tmth) / 12 * 52)$y
      tanr <- max(twkmx, na.rm = TRUE) - min(twkmn, na.rm = TRUE)
    }
    else {
      if (length(unique(tme$year)) > 1) warna()
      tanr <- max(temps, na.rm = TRUE) - min(temps, na.rm = TRUE)
    }
  }
  tanr
}
#' bio3: Calculates isothermality
#'
#' @description `bio3` is used to calculate isothermality (day-to-night
#' temperature oscillations relative to annual oscillations).
#'
#' @param temps a vector of temperatures, normally for one year (see details).
#' @param tme a `POSIXlt` object representing the date and time of each `temps` value.
#' @param method An optional character string describing the method used to
#' calculate isothermality. Options include "anuclimmean", "anuclimmaxmin" or
#' unspecified (see details).
#'
#' @return a single numeric value representing isothermality for each year.
#' @export
#'
#' @details
#' If method is "anuclimmean", bio3 is calculated using "anuclim" method.
#' Temperatures are aggregated by month and spline intepolated to weekly before
#' mean diurnal temperature range is calculated, replicating the method used by
#' http://www.worldclim.org/.
#'
#' If method is "anuclimaxmin", bio3 is calculated using "anuclimmaxmin" method.
#' Maximum and minimum monthly temperatues are spline interpolated to a weekly
#' time period and range is calculated from the maximum and minimum mean weekly
#' temperature values.
#'
#' If using method "anuclimmean" or "anuclimmaxmin" and data spans more than one
#' year, data are aggregated by unique month irrespective of year and one value
#' returned.
#'
#' If method is left unspecified, bio3 is calculated using the mean of daily
#' maximum and minimum temperature. If data spans more than one year, data are
#' aggregated by unique month irrespective of year and one value returned. If
#' method is unspecified, bio7 is calculated using the maximum temperature value
#' for the year.If data spans more than one year, bio7 calculations are performed
#' on all data and single value returned.
#'
#' @seealso [bio2()] and [bio7()] for calculating diurnal and annual temperature ranges.
#' [tmecreate()] for creating a 'POSIXlt' object.
#'
#' @examples
#' temps <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' tme <- tmecreate(2010, 24)
#' plot(temps~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Temperature")
#' bio3(temps, tme)
#' bio3(temps, tme, method = "anuclimmean")
#' bio3(temps, tme, method = "anuclimmaxmin")
#'
bio3 <- function(temps, tme, method = "") {
  if (is.na(sd(temps, na.rm = TRUE)))
    tiso <- NA
  else {
    if (method == "anuclimmean") {
      if (length(unique(tme$year)) > 1) warna()
      tiso <- bio2(temps, tme, "anuclim") / bio7(temps, tme, "anuclimean")
    } else if (method == "anuclimmmaxmin") {
      if (length(unique(tme$year)) > 1) warna()
      tiso <- bio2(temps, tme, "anuclim") / bio7(temps, tme, "anuclimmaxmin")
    }
    else {
      if (length(unique(tme$year)) > 1) warnb()
      tiso <- bio2(temps, tme) / bio7(temps, tme)
    }
  }
  tiso
}
#' bio8: Calculates mean temperature of Wettest quarter
#'
#' @description `bio8` is used to calculate the mean temperature in the wettest
#' quarter of the year
#'
#' @param temps a vector of temperatures, normally for one year (see details).
#' @param prec a vector of precipitation values, normally for one year (see
#' details).
#' @param tme1 a `POSIXlt` object representing the date and time of each `temps` value.
#' @param tme2 a `POSIXlt` object representing the date and time of each `prec` value.
#' @param method An optional character string describing the methods used to
#' calculate mean temperature in the wettest quarter. Options are "anuclim" or unspecified
#' (see details).
#'
#' @return a single numeric value of mean temperature of the wettest quarter of
#' the year.
#' @export
#'
#' @details If method is "anuclim", mean monthly temperature and total monthly
#' precipitation is calculated and then spline interpolated to a weekly time
#' period. Precipitation is calculated for each 13-week period and the mean
#' temperature for the wettest period returned. If data spans more than one
#' year, data are aggregated by unique month irrespective of year and one value
#' returned. If method is unspecified, the mean temperature of the wettest
#' quarter is calculated using all `temps` values and precipitation per quarter
#' is calculated using the time interval for measurements. If data span more
#' than one year, calculations are performed on all data and a single value
#' returned.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' temps <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' tme1 <- tmecreate(2010, 6)
#' plot(temps~as.POSIXct(tme1), type = "l", xlab = "Month", ylab = "Temperature")
#' prec <- (10 * sin(c(0:364) * (pi / -360)) + rnorm(365) + 12)
#' tme2 <- tmecreate(2010, 24)
#' plot(prec~as.POSIXct(tme2), type = "l", xlab = "Month", ylab = "Precipitation")
#' bio8(temps, prec, tme1, tme2)
#' bio8(temps, prec, tme1, tme2, method = "anuclim")
#'
bio8 <- function(temps, prec, tme1, tme2, method = "") {
  if (is.na(sd(prec, na.rm = TRUE)) | is.na(sd(temps, na.rm = TRUE)))
    twet <- NA
  else {
    if (method == "anuclim") {
      if (length(unique(tme1$year)) > 1) warna()
      if (length(unique(tme2$year)) > 1) warna()
      tmth <- aggregate(temps, by = list(tme1$mon), mean, na.rm = TRUE)$x
      twk <- spline(tmth, n = length(tmth) / 12 * 52)$y
      pmth <- aggregate(prec, by = list(tme2$mon), sum, na.rm = TRUE)$x
      pwk <- spline(pmth, n = length(pmth) / 12 * 52)$y * 12 / 52
      pwk[which(pwk < 0)] <- 0
      qtr <- function(i) {
        ptw <- c(pwk, pwk)
        psu <- sum(ptw[i: (i + 12)], na.rm = TRUE)
        psu
      }
      wq <- sapply(c(1:length(pwk)), qtr)
      i <- which(wq == max(wq, na.rm = TRUE))[1]
      twk2 <- c(twk, twk)
      twet <- mean(twk2[i:(i + 12)], na.rm = TRUE)
    }
    else {
      if (length(unique(tme1$year)) > 1) warnb()
      qtr <- function(i, int) {
        prec2 <- c(prec, prec)
        psu <- sum(prec2[i: (i + int)], na.rm = TRUE)
        psu
      }
      id <- (as.numeric(tme2[2]) - as.numeric(tme2[1])) / 86400
      dd1 <- 24/(24/(1/id))
      int <- 91 / id
      wq <- sapply(c(1:length(prec)), qtr, int)
      i <- which(wq == max(wq, na.rm = TRUE))[1]
      tid <-(as.numeric(tme1[2]) - as.numeric(tme1[1])) / 86400
      dd2 <- 24/(24/(1/tid))
      tint <- 91 / tid
      ti <- i*(dd1*dd2)
      tte <- c(temps, temps)
      twet <- mean(tte[ti:(ti + tint)], na.rm = TRUE)
    }
  }
  twet
}
#' bio9: Calculates mean temperature of the driest quarter
#'
#' @description `bio9` is used to calculate the mean temperature of the driest
#' quarter of the year
#'
#' @param temps a vector of temperatures, normally for one year (see details).
#' @param prec a vector of precipitation values, normally for one year (see
#' details).
#' @param tme1 a `POSIXlt` object representing the date and time of each `temps` value.
#' @param tme2 a `POSIXlt` object representing the date and time of each `prec` value.
#' @param method character string describing the method used to calculate mean temperature
#' of the driest quarter. Options are "anuclim" or unspecified (see details).
#'
#' @return a single numeric value of mean temperature of the wettest quarter of
#' the year.
#' @export
#'
#' @details If method is "anuclim", mean monthly temperature values are
#' calculated and spline interpolated to a weekly time period. Precipitation
#' values are summed for all months and then spline interpolated to a weekly
#' time period. Mean temeprature of the driest 13-week period is returned.
#' Otherwise, annual precipitation values are used to calculate precipitation
#' in the driest three-month period and mean temperature in this period
#' returned.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' temps <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' tme1 <- tmecreate(2010, 6)
#' plot(temps~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Temperature")
#' prec <- (10 * sin(c(0:364) * (pi / -360)) + rnorm(365) + 12)
#' tme2 <- tmecreate(2010, 24)
#' plot(prec~as.POSIXct(tme2), type = "l", xlab = "Month", ylab = "Precipitation")
#' bio9(temps, prec, tme1, tme2)
#' bio9(temps, prec, tme1, tme2, method = "anuclim")
#'
bio9 <- function(temps, prec, tme1, tme2, method = "") {
  if (is.na(sd(prec, na.rm = TRUE)) | is.na(sd(temps, na.rm = TRUE)))
    tdry <- NA
  else {
    if (method == "anuclim") {
      if (length(unique(tme1$year)) > 1) warna()
      if (length(unique(tme2$year)) > 1) warna()
      tmth <- aggregate(temps, by = list(tme1$mon), mean, na.rm = TRUE)$x
      twk <- spline(tmth, n = length(tmth) / 12 * 52)$y
      pmth <- aggregate(prec, by = list(tme2$mon), sum, na.rm = TRUE)$x
      pwk <- spline(pmth, n = length(pmth) / 12 * 52)$y * 12 / 52
      pwk[which(pwk < 0)] <- 0
      qtr <- function(i) {
        ptw <- c(pwk, pwk)
        psu <- sum(ptw[i: (i + 12)], na.rm = TRUE)
        psu
      }
      dq <- sapply(c(1:length(pwk)), qtr)
      i <- which(dq == min(dq, na.rm = TRUE))[1]
      twk2 <- c(twk, twk)
      tdry <- mean(twk2[i:(i + 12)], na.rm = TRUE)
    }
    else {
      if (length(unique(tme1$year)) > 1) warnb()
      qtr <- function(i, int) {
        prec2 <- c(prec, prec)
        psu <- sum(prec2[i: (i + int)], na.rm = TRUE)
        psu
      }
      id <- (as.numeric(tme2[2]) - as.numeric(tme2[1])) / 86400
      dd1 <- 24/(24/1/id)
      int <- 91 / id
      dq <- sapply(c(1:length(prec)), qtr, int)
      i <- which(dq == min(dq, na.rm = TRUE))[1]
      tid <-(as.numeric(tme1[2]) - as.numeric(tme1[1])) / 86400
      dd2 <- 24/(24/1/tid)
      tint <- 91 / tid
      ti <- i*(dd1*dd2)
      tte <- c(temps, temps)
      tdry <- mean(tte[ti:(ti + tint)], na.rm = TRUE)
    }
  }
  tdry
}
#' bio10: Calculates mean temperature of the Warmest quarter
#'
#' @description `bio10` is used to calculate the mean temperature of the warmest
#' quarter (three months) of the year
#'
#' @param temps a vector of temperatures, normally for one year (see details).
#' @param tme a `POSIXlt` object representing the date and time of each `temps` value.
#' @param method An optional character string describing the method used to
#' calculate mean temperature of the warmest quarter. Options are "anuclim" or
#' unspecified (see details).
#'
#' @return a single numeric value of mean temperature in the warmest quarter of
#' the year.
#' @export
#'
#' @details If method is "anuclim", warmest quarter is determined to the nearest
#' week. Mean monthly temperature values are calculated and spline interpolated
#' to a weekly time period. Precipitation values are summed for all months and
#' then spline interpolated to a weekly time period. Otherwise, the mean
#' temperature of the warmest 3-month period is calculated from annual values.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' temps <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' tme <- tmecreate(2010, 6)
#' plot(temps~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Temperature")
#' bio10(temps, tme)
#' bio10(temps, tme, method = "anuclim")
#'
bio10 <- function(temps, tme, method = "") {
  if (is.na(sd(temps, na.rm = TRUE)))
    thot <- NA
  else  {
    if (method == "anuclim") {
      if (length(unique(tme$year)) > 1) warna()
      tmth <- aggregate(temps, by = list(tme$mon), mean, na.rm = TRUE)$x
      twk <- spline(tmth, n = length(tmth) / 12 * 52)$y
      qtr <- function(i) {
        tw <- c(twk, twk)
        me <- mean(tw[i: (i + 12)], na.rm = TRUE)
        me
      }
      hq <- sapply(c(1:length(twk)), qtr)
      i <- which(hq == max(hq, na.rm = TRUE))[1]
      twk2 <- c(twk, twk)
      thot <- mean(twk2[i:(i + 12)], na.rm = TRUE)
    }
    else {
      if (length(unique(tme$year)) > 1) warnb()
      qtr <- function(i, int) {
        tw <- c(temps, temps)
        me <- mean(tw[i: (i + int)], na.rm = TRUE)
        me
      }
      id <- (as.numeric(tme[2]) - as.numeric(tme[1])) / 86400
      int <- 91 / id
      hq <- sapply(c(1:length(temps)), qtr, int)
      i <- which(hq == max(hq, na.rm = TRUE))[1]
      tte <- c(temps, temps)
      thot <- mean(tte[i:(i + int)], na.rm = TRUE)
    }
  }
  thot
}
#' bio11: Calculates mean temperature of the coldest quarter
#'
#' @description `bio11` is used to calculate the mean temperature of the coldest
#' quarter (three months) of the year
#'
#' @param temps a vector of temperatures, normally for one year (see details).
#' @param prec a vector of precipitation values, normally for one year (see
#' details).
#' @param tme a `POSIXlt` object representing the date and time of each `temps` value.
#' @param method An optional character vector describing the method used to calculate
#' the mean temperature of the coldest quarter. Options  are "anuclim" or unpsecified (see
#' details).
#'
#' @return a single numeric value of mean temperature of the warmest quarter of
#' the year.
#' @export
#'
#' @details If method is "anuclim", mean monthly temperature values are
#' calculated and spline interpolated to a weekly time period. Mean temperature of the coldest 13-week period is
#' determined. If method is left unspecified, mean temperature of the coldest 3-month (91-day) period is
#' calculated from annual temperature values.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' temps <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' tme <- tmecreate(2010, 6)
#' plot(temps~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Temperature")
#' bio11(temps, tme)
#' bio11(temps, tme, method = "anuclim")
#'
bio11 <- function(temps, tme, method = "") {
  if (is.na(sd(temps, na.rm = TRUE)))
    tcold <- NA
  else {
    if (method == "anuclim") {
      if (length(unique(tme$year)) > 1) warna()
      tmth <- aggregate(temps, by = list(tme$mon), mean, na.rm = TRUE)$x
      twk <- spline(tmth, n = length(tmth) / 12 * 52)$y
      qtr <- function(i) {
        tw <- c(twk, twk)
        me <- mean(tw[i: (i + 12)], na.rm = TRUE)
        me
      }
      cq <- sapply(c(1:length(twk)), qtr)
      i <- which(cq == min(cq, na.rm = TRUE))[1]
      twk2 <- c(twk, twk)
      tcold <- mean(twk2[i:(i + 12)], na.rm = TRUE)
    }
    else {
      if (length(unique(tme$year)) > 1) warnb()
      qtr <- function(i, int) {
        tw <- c(temps, temps)
        me <- mean(tw[i: (i + int)], na.rm = TRUE)
        me
      }
      id <- (as.numeric(tme[2]) - as.numeric(tme[1])) / 86400
      int <- 91 / id
      cq <- sapply(c(1:length(temps)), qtr, int)
      i <- which(cq == min(cq, na.rm = TRUE))[1]
      tte <- c(temps, temps)
      tcold <- mean(tte[i:(i + int)], na.rm = TRUE)
    }
  }
  tcold
}

#' bio12: Calculates total annual precipitation
#' @description `bio12` is used to calculate total precipitation in the year
#' @param prec a vector of precipitation values, normally for one year (see
#'  details).
#' @param tme a `POSIXlt` object representing the date and time of each `temps` value.
#' @param method An optional character string describing the method used to
#' calculate total annual precipitation. Options are "anuclim" or unspecified (see
#' details).
#'
#' @return a single numeric value of total annual precipitation.
#' @export
#'
#' @details If method is "anuclim", monthly precipitation values are spline
#' interpolated to a weekly time period and the total for each year returned.
#' Otherwise, all precipitation values for each year are summed.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' prec <- (10 * sin(c(0:364) * (pi / -360)) + rnorm(365) + 12)
#' tme <- tmecreate(2010, 24)
#' plot(prec~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Precipitation")
#' bio12(prec, tme)
#' bio12(prec, tme, method="anuclim")
#'
bio12 <- function(prec, tme, method = "") {
  if (is.na(sd(prec, na.rm = TRUE)))
    map <- NA
  else {
    if (method == "anuclim") {
      if (length(unique(tme$year)) > 1) warna()
      pmth <- aggregate(prec, by = list(tme$mon), sum, na.rm = TRUE)$x
      pwk <- spline(pmth, n = length(pmth) / 12 * 52)$y * 12 / 52
      pwk[which(pwk < 0)] <- 0
      map <- sum(pwk, na.rm = TRUE) / length(unique(tme$year))
    }
    else {
      if (length(unique(tme$year)) > 1) warnb()
      map <- sum(prec, na.rm = TRUE) / length(unique(tme$year))
    }
  }
  map
}
#' bio13: Calculates precipitation of the wettest period
#'
#' @description `bio13` is used to calculate the precipitation of the wettest
#' week or month of the year, depending on the time step.
#'
#' @param prec a vector of precipitation values, normally for one year (see
#' details).
#' @param tme a `POSIXlt` object representing the date and time of each `prec` value.
#' @param method An optional character string describing how the maximum weekly or
#' monthly precipitation is calculated. Options are "week" and "month" (see
#' details).
#'
#' @return a single numeric value of total precipitation in the wettest week or month of the year.
#' @export
#'
#' @details
#' If method is "week", monthly precipitation values are spline interpolated to a
#' weekly time period and the maximum weekly precipitation is returned. If method
#' is "month", monthly precipitation values are summed and the maximum monthly
#' precipitation is returned.
#'
#' If data span more than one year, data are aggregated by unique month
#' irrespective of year and one value returned.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' prec <- (10 * sin(c(0:364) * (pi / -360)) + rnorm(365) + 12)
#' tme <- tmecreate(2010, 24)
#' plot(prec~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Precipitation")
#' bio13(prec, tme, method="week")
#' bio13(prec, tme, method="month")
#' bio13(prec, tme)
#'
bio13 <- function(prec, tme, method = "week") {
  if (is.na(sd(prec, na.rm = TRUE)))
    wp <- NA
  else {
    if (method == "month"){
      if (length(unique(tme$year)) > 1) warna()
      pmth <- aggregate(prec, by = list(tme$mon), sum, na.rm = TRUE)$x
      pmth[which(pmth < 0)] <- 0
      wp <- max(pmth, na.rm = TRUE)
    }
    else {
      if(method == "week"){
        if (length(unique(tme$year)) > 1) warna()
        pmth <- aggregate(prec, by = list(tme$mon), sum, na.rm = TRUE)$x
        pwk <- spline(pmth, n = length(pmth) / 12 * 52)$y * 12 / 52
        pwk[which(pwk < 0)] <- 0
        wp <- max(pwk, na.rm = TRUE)
      }
      else {
        if (length(unique(tme$year)) > 1) warnb()
        pmth <- aggregate(prec, by = list(tme$mon), sum, na.rm = TRUE)$x
        wp <- max(pmth, na.rm = TRUE)
      }
    }
  }
  wp
}
#' bio14: Calculates precipitation of the driest period
#'
#' @description `bio14` is used to calculate the precipitation in the driest
#' period of the year
#' @param prec a vector of precipitation values, normally for one year (see
#' details).
#' @param tme a `POSIXlt` object representing the date and time of each `prec` value.
#' @param method an optional character string describing how the minimum weekly
#'  or monthly precipitation is calculated. Options include"week", "month" or unspecified
#'  (see details).
#'
#' @return a single numeric value of total precipitation in the driest week or
#' month of the year.
#' @export
#'
#' @details If method is "week" or left unspecified, monthly precipitation values
#' are spline interpolated to a weekly time period and the minimum weekly
#' precipitation is returned. If method is "month", the minimum monthly
#' precipitation is returned.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' prec <- (10 * sin(c(0:364) * (pi / -360)) + rnorm(365) + 12)
#' tme <- tmecreate(2010, 24)
#' plot(prec~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Precipitation")
#' bio14(prec, tme)
#' bio14(prec, tme, method="week")
#' bio14(prec, tme, method="month")
#'
bio14 <- function(prec, tme, method = "week") {
  if (is.na(sd(prec, na.rm = TRUE)))
    dp <- NA
  else {
    if (method == "month"){
      if (length(unique(tme$year)) > 1) warna()
      pmth <- aggregate(prec, by = list(tme$mon), sum, na.rm = TRUE)$x
      pmth[which(pmth < 0)] <- 0
      dp <- min(pmth, na.rm = TRUE)
    }
    else {
      if(method == "week"){
        if (length(unique(tme$year)) > 1) warna()
        pmth <- aggregate(prec, by = list(tme$mon), sum, na.rm = TRUE)$x
        pwk <- spline(pmth, n = length(pmth) / 12 * 52)$y * 12 / 52
        pwk[which(pwk < 0)] <- 0
        dp <- min(pwk, na.rm = TRUE)
      }
      else {
        if (length(unique(tme$year)) > 1) warnb()
        pmth <- aggregate(prec, by = list(tme$mon), sum, na.rm = TRUE)$x
        dp <- min(pmth, na.rm = TRUE)
      }
    }
  }
  dp
}
#' bio15: Calculates precipitation seasonality
#'
#' @description `bio15` is used to calculate precipitation seasonality, which is
#' the standard deviation of weekly or monthly precipitation values as a
#' percentage of the mean of those values.
#'
#' @param prec a vector of precipitation values, normally for one year (see
#' details).
#' @param tme a `POSIXlt` object representing the date and time of each `prec` value.
#' @param method an optional character string describing the method used to
#'  calculate precipitation seasonality. Options include "anuclim" or unspecified (see
#'  details).
#'
#' @return a single numeric value representing precipitation seasonality.
#' @export
#'
#' @details
#' If method is "anuclim", monthly precipitation is spline interpolated to a
#' weekly time period and precipitation seasonality calculated using these
#' values, replicating the method used by http://www.worldclim.org/. Otherwise,
#' precipitation seasonality is calculated using yearly values.
#'
#' If using `anuclim` method and data span more than one year, data are
#' aggregated by unique month irrespective of year and one value returned. If
#' method is left unspecified and data span more than one year, calculations
#' will be performed on all data and a single value returned.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' prec <- (10 * sin(c(0:364) * (pi / -360)) + rnorm(365) + 12)
#' tme <- tmecreate(2010, 24)
#' plot(prec~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Precipitation")
#' bio15(prec, tme, method="week")
#' bio15(prec, tme, method="month")
#'
#'
bio15 <- function(prec, tme, method = "anuclim") {
  if (is.na(sd(prec, na.rm = TRUE)))
    cvp <- NA
  else {
    if (method == "anuclim"){
      if (length(unique(tme$year)) > 1) warna()
      pmth <- aggregate(prec, by = list(tme$mon), sum, na.rm = TRUE)$x
      pwk <- spline(pmth, n = length(pmth) / 12 * 52)$y * 12 / 52
      pwk[which(pwk <= 0)] <- 1
      cvp <- (sd(pwk, na.rm = TRUE)  / mean(pwk, na.rm = TRUE)) * 100
    }
    else {
      if (length(unique(tme$year)) > 1) warnb()
      cvp <- (sd(prec, na.rm = TRUE)  / mean (prec, na.rm = TRUE)) * 100
    }
  }
  cvp
}
#' bio16: Calculates precipitation of the wettest quarter
#'
#' @description `bio16` is used to calculate the total precipitation of the
#' wettest quarter of the year
#'
#' @param prec a vector of precipitation values, normally for one year (see
#'  details).
#' @param tme a `POSIXlt` object representing the date and time of each `prec` value.
#' @param method an optional character string describing how precipitation of the wettest
#' quarter is calculated. Options include "anuclim" or unspecified (see details).
#'
#' @return a single numeric value for precipitation in the wettest quarter of the year.
#' @export
#'
#' @details If method is "anuclim", monthly precipitation is spline interpolated
#' to a weekly time period. Precipitation for each 13-week period is calculated
#' and total precipitation in the wettest quarter returned. If data span more
#' than one year, data are aggregated by unique month irrespective of year and
#' one value returned. Otherwise, precipitation for each three-month (91-day)
#' period is calculated and total precipitation in the wettest quarter
#' returned. If data span more than one year, calculations are performed on all
#' data and a single value returned.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' prec <- (10 * sin(c(0:364) * (pi / -360)) + rnorm(365) + 12)
#' tme <- tmecreate(2010, 24)
#' plot(prec~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Precipitation")
#' bio16(prec, tme)
#' bio16(prec, tme, method="anuclim")
#'
bio16 <- function(prec, tme, method = "") {
  if (is.na(sd(prec, na.rm = TRUE)))
    pwet <- NA
  else {
    if (method == "anuclim") {
      if (length(unique(tme$year)) > 1) warna()
      pmth <- aggregate(prec, by = list(tme$mon), sum, na.rm = TRUE)$x
      pwk <- spline(pmth, n = length(pmth) / 12 * 52)$y * 12 / 52
      pwk[which(pwk < 0)] <- 0
      qtr <- function(i) {
        pw <- c(pwk, pwk)
        su <- sum(pw[i: (i + 12)], na.rm = TRUE)
        su
      }
      wq <- sapply(c(1:length(pwk)), qtr)
      i <- which(wq == max(wq, na.rm = TRUE))[1]
      pwk2 <- c(pwk, pwk)
      pwet <- sum(pwk2[i:(i + 12)], na.rm = TRUE)
    }
    else {
      if (length(unique(tme$year)) > 1) warnb()
      qtr <- function(i, int) {
        pw <- c(temps, temps)
        su <- sum(pw[i: (i + int)], na.rm = TRUE)
        su
      }
      id <- (as.numeric(tme[2]) - as.numeric(tme[1])) / 86400
      int <- 91 / id
      wq <- sapply(c(1:length(prec)), qtr, int)
      i <- which(wq == max(wq, na.rm = TRUE))[1]
      pre2 <- c(prec, prec)
      pwet <- sum(pre2[i:(i + int)], na.rm = TRUE)
    }
  }
  pwet
}
#' bio17: Calculates precipitation of the driest quarter
#'
#' @description `bio17` is used to calculate the precipitation in the driest
#' quarter of the year
#'
#' @param prec a vector of precipitation values, normally for one year (see
#' details).
#' @param tme a `POSIXlt` object representing the date and time of each `prec` value.
#' @param method an optional character string describing how quarterly
#' precipitation is calculated. Options include "anuclim" or unspecified (see details).
#'
#' @return a single numeric value of precipitation of the driest quarter.
#' @export
#'
#' @details If method is "anuclim", monthly precipitation is spline interpolated
#' to a weekly time period and precipitation of each 13-week period is
#' calculated. The precipitation in the driest quarter is then found. If data
#' spans more than one year, data are aggregated by unique month irrespective
#' of year and one value returned. Otherwise, precipitation in each three-month
#' period is calculated and total precipitation in the driest quarter
#' returned. If data span more than one year, calculations are performed on all
#' data and single value returned.
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @examples
#' prec <- (10 * sin(c(0:364) * (pi / -360)) + rnorm(365) + 12)
#' tme <- tmecreate(2010, 24)
#' plot(prec~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Precipitation")
#' bio17(prec, tme)
#' bio17(prec, tme, method="anuclim")
#'
bio17 <- function(prec, tme, method = "") {
  if (is.na(sd(prec, na.rm = TRUE)))
    pdry <- NA
  else {
    if (method == "anuclim") {
      if (length(unique(tme$year)) > 1) warna()
      pmth <- aggregate(prec, by = list(tme$mon), sum, na.rm = TRUE)$x
      pwk <- spline(pmth, n = length(pmth) / 12 * 52)$y * 12 / 52
      pwk[which(pwk < 0)] <- 0
      qtr <- function(i) {
        pw <- c(pwk, pwk)
        su <- sum(pw[i: (i + 12)], na.rm = TRUE)
        su
      }
      dq <- sapply(c(1:length(pwk)), qtr)
      i <- which(dq == min(dq, na.rm = TRUE))[1]
      pwk2 <- c(pwk, pwk)
      pdry <- sum(pwk2[i:(i + 12)], na.rm = TRUE)
    }
    else {
      if (length(unique(tme$year)) > 1) warnb()
      qtr <- function(i, int) {
        pw <- c(prec, prec)
        su <- sum(pw[i: (i + int)], na.rm = TRUE)
        su
      }
      id <- (as.numeric(tme[2]) - as.numeric(tme[1])) / 86400
      int <- 91 / id
      dq <- sapply(c(1:length(prec)), qtr, int)
      i <- which(dq == min(dq, na.rm = TRUE))[1]
      pre2 <- c(prec, prec)
      pdry <- sum(pre2[i:(i + int)], na.rm = TRUE)
    }
  }
  pdry
}
#' bio18: Precipitation of the warmest quarter
#'
#' @description `bio18` is used to calculate the precipitation in the warmest
#' quarter of the year
#'
#' @param temps a vector of temperature values, normally for one year (see
#' details).
#' @param prec a vector of precipitation values, normally for one year (see
#' details).
#' @param tme1 a `POSIXlt` object representing the date and time of each `temps` value.
#' @param tme2 a `POSIXlt` object representing the date and time of each `prec` value.
#' @param method an optional character string describing how quarterly
#' temperature and precipitation are calculated. Options are "anuclim" or unspecified
#' (see details).
#'
#' @return a single numeric value of precipitation in the warmest quarter of the
#' year.
#' @export
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @details If method is "anuclim", monthly mean temperature and total monthly
#' precipitation are interpolated to a weekly time period before calculating
#' mean temperature of each 13-week period. The precipitation in the warmest
#' quarter is then found. If data span more than one year, data are aggregated
#' by unique month irrespective of year and one value returned. If method is
#' unspecified, the mean temperature in each three-month period is calculated
#' and precipitation in the coldest quarter returned. If data span more than
#' one year, calculations are performed on all data and single value returned.
#'
#' @examples
#' prec <- (10 * sin(c(0:364) * (pi / -360)) + rnorm(365) + 12)
#' temps <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' tme1 <- tmecreate(2010, 6)
#' tme2<- tmecreate(2010, 24)
#' plot(temps~as.POSIXct(tme1), type = "l", xlab = "Month", ylab = "Temperature")
#' bio18(temps, prec, tme1, tme2)
#' bio18(temps, prec, tme1, tme2, method="anuclim")
bio18 <- function(temps, prec, tme1, tme2, method = "") {
  if (is.na(sd(prec, na.rm = TRUE)) | is.na(sd(temps, na.rm = TRUE)))
    pwarm <- NA
  else {
    if (method == "anuclim") {
      if (length(unique(tme1$year)) > 1) warna()
      if (length(unique(tme2$year)) > 1) warna()
      tmth <- aggregate(temps, by = list(tme1$mon), mean, na.rm = TRUE)$x
      twk <- spline(tmth, n = length(tmth) / 12 * 52)$y
      pmth <- aggregate(prec, by = list(tme2$mon), sum, na.rm = TRUE)$x
      pwk <- spline(pmth, n = length(pmth) / 12 * 52)$y * 12 / 52
      pwk[which(pwk < 0)] <- 0
      qtr <- function(i) {
        tw <- c(twk, twk)
        me <- mean(tw[i: (i + 12)], na.rm = TRUE)
        me
      }
      wq <- sapply(c(1:length(twk)), qtr)
      i <- which(wq == max(wq, na.rm = TRUE))[1]
      pwk2 <- c(pwk, pwk)
      pwarm <- sum(pwk2[i:(i + 12)], na.rm = TRUE)
    }
    else {
      if (length(unique(tme1$year)) > 1) warnb()
      if (length(unique(tme2$year)) > 1) warnb()
      qtr <- function(i, int) {
        tw <- c(temps, temps)
        me <- mean(tw[i: (i + int)], na.rm = TRUE)
        me
      }
      id <- (as.numeric(tme1[2]) - as.numeric(tme1[1])) / 86400
      dd1 <- 24/(24/(1/id))
      int <- 91 / id
      wq <- sapply(c(1:length(temps)), qtr, int)
      i <- which(wq == max(wq, na.rm = TRUE))[1]
      pid <-(as.numeric(tme2[2]) - as.numeric(tme2[1])) / 86400
      dd2 <- 24/(24/(1/pid))
      pint <- 91 / pid
      pi <- i*(dd2/dd1)
      pre2 <- c(prec, prec)
      pwarm <- sum(pre2[pi:(pi + pint)], na.rm = TRUE)
    }
  }
  pwarm
}
#' bio19: Precipitation of the coldest quarter
#'
#' @description `bio19` is used to calculate the precipitation in the coldest
#' quarter of the year.
#'
#' @param temps a vector of temperature values, normally for one year (see
#' details)
#' @param prec a vector of precipitation values, normally for one year (see
#' details).
#' @param tme1 a `POSIXlt` object representing the date and time of each `temps` value.
#' @param tme2 a `POSIXlt` object representing the date and time of each `prec` value.
#' @param method an optional character string describing how quarterly mean
#' temperature and precipitation are calculated. Options are "anuclim" or unspecified
#' (see details).
#'
#' @return a single numeric value of precipitation in the coldest quarter of the
#' year.
#' @export
#'
#' @seealso the [tmecreate()] function can be used to create a POSIXlt object.
#'
#' @details If method is "anuclim", monthly mean temeprature and total monthly
#' precipitation are interpolated to weekly time period before calculating mean
#' temperature for each 13-week period. The precipitation in the coldest
#' quarter is then calculated. If data spans more than one year, data are
#' aggregated by unique month irrespective of year and one value returned If
#' method is left unspecified, the mean temperature in each three-month period
#' is calculated and precipitation in the coldest quarter returned. If data
#' spans more than one year, calculations are performed on all data and single
#' value returned.
#'
#' @examples
#'
#' prec <- 10 * sin(c(0:364) * (pi / -360)) + (rnorm(365) + 12)
#' temps <- 10 * sin(c(0:1459) / (pi * 150)) + rnorm(1460)
#' tme1 <- tmecreate(2010, 6)
#' tme2<- tmecreate(2010, 24)
#' plot(prec~as.POSIXct(tme), type = "l", xlab = "Month", ylab = "Precipitation")
#' bio19(temps, prec, tme1, tme2)
#' bio19(temps, prec, tme1, tme2, method="anuclim")
bio19 <- function(temps, prec, tme1, tme2, method = "") {
  if (is.na(sd(prec, na.rm = TRUE)) | is.na(sd(temps, na.rm = TRUE)))
    pcld <- NA
  else {
    if (method == "anuclim") {
      if (length(unique(tme1$year)) > 1) warna()
      if (length(unique(tme2$year)) > 1) warna()
      tmth <- aggregate(temps, by = list(tme1$mon), mean, na.rm = TRUE)$x
      twk <- spline(tmth, n = length(tmth) / 12 * 52)$y
      pmth <- aggregate(prec, by = list(tme2$mon), sum, na.rm = TRUE)$x
      pwk <- spline(pmth, n = length(pmth) / 12 * 52)$y * 12 / 52
      pwk[which(pwk < 0)] <- 0
      qtr <- function(i) {
        tw <- c(twk, twk)
        me <- mean(tw[i: (i + 12)], na.rm = TRUE)
        me
      }
      cq <- sapply(c(1:length(twk)), qtr)
      i <- which(cq == min(cq, na.rm = TRUE))[1]
      pwk2 <- c(pwk, pwk)
      pcld <- sum(pwk2[i:(i + 12)], na.rm = TRUE)
    }
    else {
      if (length(unique(tme1$year)) > 1) warnb()
      if (length(unique(tme2$year)) > 1) warnb()
      qtr <- function(i, int) {
        tw <- c(temps, temps)
        me <- mean(tw[i: (i + int)], na.rm = TRUE)
        me
      }
      id <- (as.numeric(tme1[2]) - as.numeric(tme1[1])) / 86400
      dd1 <- 24/(24/(1/id))
      int <- 91 / id
      cq <- sapply(c(1:length(temps)), qtr, int)
      i <- which(cq == min(cq, na.rm = TRUE))[1]
      pid <-(as.numeric(tme2[2]) - as.numeric(tme2[1])) / 86400
      dd2 <- 24/(24/(1/pid))
      pint <- 91 / pid
      pi <- i*(dd2/dd1)
      pre2 <- c(prec, prec)
      pcld <- sum(pre2[pi:(pi + pint)], na.rm = TRUE)
    }
  }
  pcld
}

