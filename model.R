# Full TempEst 2 implementation.
# Author: Daniel Philippus
# Last update: March 15, 2024 (code); Feb. 11, 2025 (documentation).
#
# This file contains the full implementation of TempEst 2/SCHEMA.  Most functions are used internally;
# the user is concerned with `full.schema`, and possibly `krig.ssn` and `krig.anom`, which implement
# the main model logic.
#
# `full.schema` contains the "glue code" for the SCHEMA logic, combinging the seasonal and anomaly terms.
# The two component functions handle data processing and prediction.
# `krig.ssn` is the default implementation of seasonality estimation.
# `krig.anom` is the default implementation of anomaly estimation.
#
# The functions `sigmoider`, `sigsmooth`, and those starting with `physfit` are used in anomaly
# estimation for fitting sensitivity coefficients and smoothing inputs.
#
# The `fit.sins` and `make.sints` functions implement three-sine seasonality
# (Philippus, Corona, and Hogue 2024, https://doi.org/10.1111/1752-1688.13228).
# `fit.sins` is used to fit seasonality coefficients to a timeseries (for training),
# and `make.sints` reproduces a seasonality timeseries from the coefficients.


library(tidyverse)
library(fields)

sigmoider <- function(lst, sigrng = 3/10) {
  # Implements a sigmoid function with a mean of 0.
  # lst is the input number.
  # sigrng is the "range" of the sigmoid, scaling the width.  A default of 0.3
  # works well for most sites, but this can be tuned.
  # Returns: sigmoid-scaled LST anomaly.
  1 / (1 + exp(-lst * sigrng)) - 0.5
}

sigsmooth <- function(lst, k = 3/10) {
  # Implements sigmoidal smoothing for LST anomaly - that is,
  # smoothing of recent variation and application of the sigmoid scaling.
  # The smoother (1, 1/2, 1/4, 1/6, 1/8, 1/10, 1/12) performs well for most sites.
  # lst is the vector of LST anomalies.
  # k is the sigmoid "range", described in `sigmoider`.
  # Returns: smoothed and transformed LST anomaly.
  lstsig <- sigmoider(lst, k)
  filt <- 1/c(1, 1:6*2)
  len <- length(filt)
  # Pad with first value, then cut the padding off again
  stats::filter(c(rep(lstsig[1], len), lstsig), filt, sides=1)[(len+1):(len + length(lstsig))]
}

humsmooth <- function(hum) {
  # Implements the simple smoothing function for recent humidity anomaly.
  # The smoother (1, 1/10) works well for most sites.
  # hum is the vector of humidity anomalies.
  # Returns: smoothed humidity anomaly.
  hum + lag(hum) / 10
}

physfit.obj <- function(gd) {
  # Fit anomaly sensitivity, returning an `lm` object.
  # gd is the data frame containing *smoothed anomaly* timeseries: temperature, lst, and humidity.
  # Returns: lm object for temp. anomaly as a linear function of smoothed LST and humidity anomaly.
  temp <- gd$temperature
  lst <- gd$lst
  hum <- gd$humidity
  
  lm(temp ~ 0 + lst + hum)
}

physfit <- function(gd) {
  # Fit anomaly sensitivity, returning the sensitivity coefficients.
  # See description of physfit.obj.
  # Returns: single-row tibble of coefficients.
  obj <- physfit.obj(gd)
  
  as_tibble_row(obj$coefficients)
}

physfit.prd <- function(lst, hum, temp) {
  # Compute the temperature anomaly that would be predicted by the sensitivity fit.
  # Sensitivity fit is described in `physfit.obj`.
  # lst: vector of smoothed LST anomaly
  # hum: vector of smoothed humidity anomaly
  # temp: vector of temperature anomaly
  # Returns: vector of predicted temperature anomalies
  physfit.obj(tibble(lst=lst, humidity=hum, temperature=temp))$fitted.values
}

physfit.vec <- function(lst, hum, temp) {
  # Same as `physfit`, but renamed.
  rename(physfit(tibble(lst=lst, humidity=hum, temperature=temp)),
         Coef.LSTSigmoid = lst,
         Coef.Humidity = hum)
}

physfit.vecplus <- function(lst, hum, temp) {
  # Same as `physfit.vec`, but adds an intercept term for daily *max* prediction.
  # Daily mean anomaly has an intercept of 0; daily max does not.
  as_tibble_row(lm(temp ~ lst + hum)$coefficients) %>%
    rename(InterceptP = `(Intercept)`,
           Coef.LSTSigmoidP = lst,
           Coef.HumidityP = hum)
}

full.schema <- function(sche=krig.ssn, ma=krig.anom, rtn.model=FALSE,
                        use.max=FALSE) {
  # Implements SCHEMA logic.  This is a higher-order function that returns a model training function.
  # Arguments are used to specify SCHEMA behavior.
  # sche: Seasonal component.  Accepts training data as an input and returns a seasonality prediction
  # function (data with id, day, inputs --> data with id, day, temp.doy).
  # ma: Anomaly component.  Accepts training data as an input and returns an anomaly prediction
  # function (data with id, date, day, inputs --> data with id, date, day, temp.anom; also temp.plus if use.max).
  # rtn.model: Instead of returning prediction function, returns the actual models for inspection.
  #   sche and ma functions also take rtn.model as an input and behave likewise.
  # use.max: Also predict daily max temperature.  Affects anomaly, but not seasonality.
  function(indat) {
    # Returns a function that takes data --> prediction
    # Required columns in indat and data:
    # id, lon [deg], lat [deg],
    # day [Julian day], date,
    # elevation [m, mean in 500 m radius],
    # water, shrubland, grassland, barren [abundance in 500 m radius],
    # lst [C, tested with daytime MODIS in 500 m radius],
    # humidity [specific humidity, kg/kg]
    #
    # indat additionally requires temperature [C, daily mean observed] for training
    #
    # if rtn.model, return the model components instead of a predictor function

    ssn <- sche(indat, rtn.model=rtn.model)  # ssn(data) -> tibble(id, day, temp.doy)
    anom <- ma(indat, rtn.model=rtn.model,
               use.max=use.max)  # anom(data) -> tibble(id, date, temp.anom)
    
    if (rtn.model) {
      # Just combine the model vectors returned by the components.
      c(ssn, anom)
    } else {
      function(prdat) {
        # For prediction, the SCHEMA function just needs to join the component results
        # appropriately, then compute the actual prediction.
        doypr <- ssn(prdat)
        anpr <- anom(prdat)
        
        left_join(prdat, anpr, by=c("id", "date")) %>%
          left_join(doypr, by=c("id", "day")) %>%
          (\(x) {if (use.max) {
            mutate(x, temp.mod = temp.doy + temp.anom,
                 temp.max = temp.doy + temp.plus)
          } else {
            mutate(x, temp.mod = temp.doy + temp.anom)
          }})
      }
    }
  }
}

krig.anom <- function(indat, rtn.model=FALSE, use.max=FALSE) {
  # Returns function: data --> predictions
  # Required indat and prediction columns:
  # id, lon [deg], lat [deg],
  # day [Julian day], date,
  # elevation [m, mean in 500 m radius],
  # water [abundance in 500 m radius],
  # lst [C, tested with daytime MODIS in 500 m radius],
  # humidity [specific humidity, kg/kg]
  #
  # indat additionally requires temperature [C, daily mean observed] for training
  # and temperature.max if use.max is TRUE
  # Smoother: computes temperature anomaly and smoothed input anomaly.
  # Or just input, if this is prediction not training.
  smoother <- \(x) x %>%
    drop_na() %>%
    group_by(id, day) %>%
    (
      \(y) {
        # Add temperature max anomaly if applicable
        if ("temperature.max" %in% names(y) & "temperature" %in% names(y)) {
          mutate(y, temperature.plus = temperature.max - mean(temperature,
                                                              na.rm=T))
        } else y
      }
    ) %>%
    mutate(across(c(any_of("temperature"), lst, humidity), # Compute anomalies
                  ~(.x - mean(.x, na.rm=T)))) %>%
    group_by(id) %>%
    arrange(date) %>%
    mutate(lst = sigsmooth(lst),  # apply smoothers
           humidity = humsmooth(humidity)) %>%
    drop_na()

  # Compute anomaly sensitivities as well as the required predictor variables
  fitted <- smoother(indat) %>%
    group_by(id) %>%
    summarize(
      physfit.vec(lst, humidity, temperature),
      if(use.max) physfit.vecplus(lst, humidity, temperature.plus),
      across(c(lon, lat, water, elevation), mean),
      humidity_sd = sd(humidity)
    ) %>%
    ungroup() %>%
    drop_na()

  # Convenience functions to select X and Z inputs for kriging (Z = non-spatial covariates).
  xer <- \(x) as.matrix(select(x, lon, lat))
  zer <- \(x) as.matrix(select(x, elevation, water, humidity_sd))

  # This setup allows pre-fitted spatial covariance coefficients (the ...) to be used for prediction,
  # but not for returning model components.  That's because fitting them can be a bit crashy and slow, so
  # we want actually building a model to be reliable (by using pre-fitted ones), but we also want to
  # be able to re-evaluate them when we're inspecting components.
  fn <- function(y, ...) {
    if (rtn.model) {
      spatialProcess(xer(fitted), y, Z = zer(fitted), Distance = "rdist.earth")
    } else {
      spatialProcess(xer(fitted), y, Z = zer(fitted), Distance = "rdist.earth", ...)
    }
  }

  # Coefficient predictor fits.  aRange and lambda have been pre-selected by maximum likelihood estimation
  # during model development.  Covariance function type (e.g., Matern, exponential) was also tested,
  # but the `fields` default (Matern, number=1) works well.
  coef.lst <- fn(fitted$Coef.LSTSigmoid
                             , aRange = 23, lambda = 1.1
                 )
  coef.lstp <- if (use.max) fn(fitted$Coef.LSTSigmoidP
                  , aRange = 28, lambda = 1.1
  ) else NULL
  
  coef.humidity <- fn(fitted$Coef.Humidity
                                  , aRange = 23, lambda = 1.4
                      )
  coef.humidityp <- if (use.max) fn(fitted$Coef.HumidityP
                      , aRange = 34, lambda = 1.6
  ) else NULL
  
  coef.intp <- if (use.max) fn(fitted$InterceptP
                                    , aRange = 32, lambda = 0.9
  ) else NULL
  
  if (rtn.model) {
    # Just return the actual spatial fit objects.
    if (use.max) {
      list(
        "LST" = coef.lst,
        "Humidity" = coef.humidity,
        "LSTMax" = coef.lstp,
        "HumidityMax" = coef.humidityp,
        "InterceptMax" = coef.intp
      )
    } else {
      list(
        "LST" = coef.lst,
        "Humidity" = coef.humidity
      )
    }
  } else {
    # Actually predict the coefficients.
    function(prdat) {
      # Prepare covariates
      smdat <- prdat %>%
        smoother
      
      codat <- smdat %>%
        group_by(id) %>%
        summarize(
          across(c(lon, lat, water, elevation), mean),
          humidity_sd = sd(humidity)
        ) %>%
        drop_na %>%
        ungroup()

      # Now, predict all of the coefficients.
      codat$lstco <- predict(coef.lst, xer(codat), Z=zer(codat))[,1]
      codat$humco <- predict(coef.humidity, xer(codat), Z=zer(codat))[,1]
      codat$lstpco <- if (use.max) predict(coef.lstp, xer(codat), Z=zer(codat))[,1] else NULL
      codat$humpco <- if (use.max) predict(coef.humidityp, xer(codat), Z=zer(codat))[,1] else NULL
      codat$intpco <- if (use.max) predict(coef.intp, xer(codat), Z=zer(codat))[,1] else NULL

      # Finally, compute the estimated anomaly.
      smdat <- left_join(smdat, codat, by="id") %>%
        (\(y) { if (use.max) {
          mutate(y,
          temp.anom = lst * lstco + humidity * humco,
          temp.plus = intpco + lst * lstpco + humidity * humpco
          )
        } else {
          mutate(y,
                 temp.anom = lst * lstco + humidity * humco)
        }
        })
      
      select(smdat, id, date, temp.anom, any_of("temp.plus"))
    }
  }
}

krig.ssn <- function(indat, rtn.model=FALSE, ...) {
  # Returns a function that takes data --> prediction
  # Required columns in indat and data:
  # id, lon [deg], lat [deg],
  # day [Julian day], date,
  # elevation [m, mean in 500 m radius],
  # water, shrubland, grassland, barren [abundance in 500 m radius],
  # lst [C, tested with daytime MODIS in 500 m radius],
  # humidity [specific humidity, kg/kg]
  #
  # indat additionally requires temperature [C, daily mean observed] for training
  # 
  # This implementation is based on the "three-sine" annual temperature cycle function
  # from Philippus, Corona and Hogue (2024): https://doi.org/10.1111/1752-1688.13228

  # Baseline cosine timeseries (main sine in three-sine).
  costs <- \(day) cos((day - 210) * 2 * pi / 365)

  # We need coverage for some months, but not all, for use in coefficient estimation.
  # Make sure the training data has them.
  req.mos <- c("01", "03", "05", "07", "08", "09", "10", "11")
  mos <- \(x) unique(format(as.Date(x$day-1, "1970-01-01"), "%m"))
  has.mos <- \(x) mean(req.mos %in% mos(x)) == 1
  
  if (!has.mos(indat))
    stop("Missing months in training data; some data from each month are required")

  # Preprocessing.
  preproc <- function(data, fit3s) {
    # fit3s: Are we fitting three-sine coefficients (training) or only estimating them (prediction)?
    data %>%
      group_by(id, day) %>%
      summarize(
        # Static variables - first
        across(c(lon, lat, elevation, date,
                 water, grassland, shrubland, barren), first),
        # Dynamic variables - mean
        across(c(lst, humidity, any_of("temperature")), ~mean(.x, na.rm=TRUE))
      ) %>%
      group_by(id) %>%
      arrange(day) %>%
      mutate(
        # Lots of computed covariates here.  Maxima/minima, annual amplitudes, etc.
        maxT = max(lst, na.rm=T),
        minT = min(lst, na.rm=T),
        amphum = safely(\() lm(humidity ~ costs(day))$coefficients[[2]],
                        otherwise=NA)()$result,
        meanT = mean(lst, na.rm=T),
        mean_hum = mean(humidity, na.rm=T),
        # Days below freezing
        freeze_days = sum(lst < 0, na.rm=T),
        # Cumulative degree-days
        # cdd = cumsum(lst, na.rm=T)
        # Weird thing accounts for missing days:
        # 1. Cumulative sum, ignoring NAs (-->0)
        # 2. Mean value based on non-NAs only
        # 3. Multiply cumulative mean by day of year
        cdd = cumsum(case_match(lst, NA ~ 0, .default = lst)) / cumsum(!is.na(lst)) * day
      ) %>%
      group_modify( ~ {
        if (fit3s) {
          # Add fitted three-sine coefficients for training.
          view <- .x  # This is included for debugging purposes.
          cbind(.x, fit.sins(.x))
        } else
          .x
      }) %>%
      group_by(id, month = format(date, "%m")) %>%
      summarize(
        # Monthly means - or just the first values for static quantities
        across(-c(lst, humidity, cdd, day, any_of("temperature")), first),
        across(c(lst, humidity, cdd), ~mean(.x, na.rm=TRUE))
      ) %>%
      ungroup() %>%
      pivot_wider(names_from = "month", values_from = c("lst", "humidity", "cdd"),
                  id_cols = id, unused_fn = first) %>%
      drop_na(
        # All required columns
        id, lat, lon, elevation,
        water, grassland, barren, shrubland,
        minT, maxT, amphum, freeze_days,
        lst_05, lst_09, lst_10,
        humidity_01, humidity_03, humidity_05, humidity_07, humidity_08,
        humidity_09, humidity_11,
        cdd_09, cdd_10
      ) %>%
      (\(x) {if (fit3s) drop_na(x,
        Intercept, Amplitude, FallWinter, SpringSummer, WinterDay
      ) else x})
  }

  # Prepare training data
  train <- preproc(indat, TRUE)
  
  xer <- \(x) as.matrix(select(x, lon, lat))
  # Each component gets its own set of Z-variables (non-spatial covariates).
  # These have been pre-tuned.
  zer.itx <- \(x) select(x, elevation, water, minT, maxT,
                               freeze_days, lst_09, lst_10, humidity_01,
                               humidity_03, humidity_09, humidity_11)
  zer.amp <- \(x) as.matrix(select(x, elevation, water, minT, amphum,
                                   cdd_09, cdd_10))
  zer.fw <- \(x) as.matrix(select(x, elevation, lst_05, humidity_03, humidity_05,
                                  grassland, barren, water, amphum, freeze_days))
  zer.ssu <- \(x) as.matrix(select(x, lst_05, humidity_07, elevation, water))
  zer.wid <- \(x) as.matrix(select(x, elevation, shrubland, grassland, barren, water, amphum,
                                   freeze_days, lst_05, humidity_03, humidity_05, humidity_08))
  # See explanation on krig.anom for this thing.
  fn <- function(y, z, ...) {
    if (rtn.model) {
      spatialProcess(xer(train), y, Z = z, Distance = "rdist.earth")
    } else {
      spatialProcess(xer(train), y, Z = z, Distance = "rdist.earth", ...)
    }
  }

  # Prepare spatial fits.
  kix <- fn(train$Intercept, zer.itx(train)
            ,aRange = 5.7, lambda = 0.14
            )
  amp <- fn(train$Amplitude, zer.amp(train)
            ,aRange = 31, lambda = 0.6
            )
  fw <- fn(train$FallWinter, zer.fw(train)
           ,aRange = 25, lambda = 0.78
           )
  # Turns out we need the FallWinter coefficient as a predictor of SpringSummer.
  ssu <- fn(train$SpringSummer, cbind(zer.ssu(train), train$FallWinter)
            ,aRange = 15, lambda = 1.2
            )
  wid <- fn(train$WinterDay, zer.wid(train)
            ,aRange = 28, lambda = 2
            )
  
  if (rtn.model) {
    # Just return components.
    list(
      "Intercept" = kix,
      "Amplitude" = amp,
      "AutumnWinter" = fw,
      "SpringSummer" = ssu,
      "WinterDay" = wid
    )
  } else {
    # Actually make predictions
    function(data) {
      if (!has.mos(data)) {
        warning("Missing months in prediction data; some data from each month are required")
        tibble(id=NA, day=NA, temp.doy=NA)
      } else {
        ppd <- preproc(data, FALSE)
        # Generate predictions.  We have to do FallWinter first so we can use it as an
        # input to SpringSummer.
        ppd %>%
          group_by(id) %>%
          group_modify(~{
            fallwint <- predict(fw, xer(.x), Z=zer.fw(.x))[,1]
            sinfit <- tibble(
              Intercept = predict(kix, xer(.x), Z=zer.itx(.x))[,1],
              Amplitude = predict(amp, xer(.x), Z=zer.amp(.x))[,1],
              FallWinter = fallwint,
              SpringSummer = predict(ssu, xer(.x), Z=cbind(zer.ssu(.x), fallwint))[,1],
              WinterDay = predict(wid, xer(.x), Z=zer.wid(.x))[,1],
              # These have so much uncertainty that predicting them doesn't improve performance.
              # Just use baseline values.
              SpringDay = 160,
              FallDay = 330,
              SummerDay = 220
            )
            # Build the timeseries
            make.sints(sinfit)
          }) %>%
          ungroup() %>%
          rename(temp.doy = actemp) %>%
          select(id, day, temp.doy)
      }
    }
  }
}

#### 3-sin functions
# See documentation here: https://github.com/quantum-dan/seasonality

applysin <- function(days, start, domain) {
  # Apply one of the truncated sine functions, within the domain and offset
  # from start
  sin(((days - start) %% 365) * (2*pi)/length(domain)) * (days %in% domain)
}

fit.sins <- function(data, mod=FALSE, aic=FALSE,
                     erv=tibble(
                       Intercept = NA,
                       Amplitude = NA,
                       FallDay = NA,
                       WinterDay = NA,
                       FallWinter = NA,
                       SpringDay = NA,
                       SummerDay = NA,
                       SpringSummer = NA,
                       R2 = NA
                     )) {
  # Fit across all data.  Use three-sinusoid approach.  Data must have day
  # and temperature.
  # mod = return model instead of tibble, unless aic
  # is true, in which case return tibble(sin, anomaly) with AICs
  inp <- data %>%
    group_by(day) %>%
    summarize(temperature = mean(temperature, na.rm=T)) %>%
    drop_na() %>%
    arrange(day)
  
  if (nrow(inp) < 180) {
    warning("Insufficient data coverage for 3-sine fit")
    erv
  } else {
    day <- inp$day
    Ts <- inp$temperature
    
    # First, fit the raw cosine.
    index <- cos((day - 210) * 2 * pi / 365)
    cosfit <- lm(Ts ~ index)$coefficients
    meant <- cosfit[[1]]
    amplitude <- cosfit[[2]] / meant
    cospr <- 1 + index * amplitude
    anomaly <- Ts/meant - cospr
    
    # Now identify coordinates for sines.
    # Use rolling mean to find zeros and peaks.
    rolled <- stats::filter(anomaly, rep(1/31, 31), circular=T)
    
    falld <- day[day >= 300]
    fally <- rolled[day %in% falld]
    fallt <- if (length(falld) > 0) mean(falld[fally == min(fally)][1]) else 330
    
    wind <- day[day <= 110]
    winy <- rolled[day %in% wind]
    wint <- if (length(wind) > 0) mean(wind[winy == max(winy)]) else 80
    
    spd <- day[day >= 120 & day <= 180]
    spy <- rolled[day %in% spd]
    spt <- if (length(spd) > 0) mean(spd[spy == min(spy)]) else 150
    
    sumd <- day[day >= 200 & day <= 240]
    sumy <- rolled[day %in% sumd]
    sumt <- if (length(sumd) > 0) mean(sumd[sumy == max(sumy)]) else 220
    
    # Define sine functions
    # Domain = half of the gap between peaks before and after
    sin1width <- round(((wint - fallt) %% 365)/2)
    sin1dom <- c((fallt - sin1width):366, 1:(wint + sin1width))
    
    sin2width <- round((sumt - spt)/2)
    sin2dom <- (spt - sin2width):(sumt + sin2width)
    
    sin1 <- -applysin(day, fallt - sin1width, sin1dom)
    sin2 <- -applysin(day, spt - sin2width, sin2dom)
    
    fullfit <- lm((Ts - meant * cospr) ~ 0 + sin1 + sin2)
    co <- fullfit$coefficients
    
    if (mod) {
      # Return lm for fit comparison (AIC)
      if (aic) {
        tibble(
          Sin = 2*(2 - as.numeric(logLik(lm(Ts ~ index)))),
          Anom = 2*(8 - as.numeric(logLik(fullfit)))
        )
      } else fullfit
    } else {
      tibble(
        Intercept = meant,
        Amplitude = amplitude * meant,
        FallDay = fallt,
        WinterDay = wint,
        FallWinter = if (!is.na(co[[1]])) co[[1]] else 0,
        SpringDay = spt,
        SummerDay = sumt,
        SpringSummer = co[[2]],
        R2 = cor(Ts, (fullfit$fitted.values + meant * cospr))^2
      )
    }
  }
}

make.sints <- function(sinfit) {
  # Sin fit --> tibble(day, actemp)
  # sinfit: tibble(Intercept, Amplitude, FallDay, WinterDay, FallWinter,
  # SpringDay, SummerDay, SpringSummer, R2)
  # Coefficients fit the actual temperature, not the normalized temperature.
  
  day <- 1:366
  
  index <- cos((day - 210) * 2 * pi / 365)
  
  # Prepare sin functions
  wint <- sinfit$WinterDay
  fallt <- sinfit$FallDay
  sumt <- sinfit$SummerDay
  spt <- sinfit$SpringDay
  
  sin1 <- if (sinfit$FallWinter != 0 & !is.na(wint + fallt)) {
    sin1width <- round(((wint - fallt) %% 365)/2)
    sin1dom <- c((fallt - sin1width):366, 1:(wint + sin1width))
    -applysin(day, fallt - sin1width, sin1dom)
  } else 0
  
  sin2 <- if (sinfit$SpringSummer != 0 & !is.na(sumt + spt)) {
    sin2width <- round((sumt - spt)/2)
    sin2dom <- (spt - sin2width):(sumt + sin2width)
    -applysin(day, spt - sin2width, sin2dom)
  } else 0
  
  actemp <- sinfit$Intercept + sinfit$Amplitude * index +
    sinfit$FallWinter * sin1 + sinfit$SpringSummer * sin2
  actemp[actemp < 0] <- 0
  
  # Compute fit
  tibble(
    day = day,
    actemp = actemp
  )
}

