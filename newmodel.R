library(tidyverse)
library(fields)

# Define functions for smoothing and modeling
sigmoider <- function(lst, sigrng = 3/10) {
  1 / (1 + exp(-lst * sigrng)) - 0.5
}

sigsmooth <- function(lst, k = 3/10) {
  lstsig <- sigmoider(lst, k)
  filt <- 1/c(1, 1:6*2)
  len <- length(filt)
  # Pad with first value, then cut the padding off again
  stats::filter(c(rep(lstsig[1], len), lstsig), filt, sides = 1)[(len + 1):(len + length(lstsig))]
}

humsmooth <- function(hum) {
  hum + lag(hum) / 10
}

physfit.obj <- function(gd) {
  temp <- gd$temperature
  lst <- gd$lst
  hum <- gd$humidity
  
  lm(temp ~ 0 + lst + hum)
}

physfit <- function(gd) {
  obj <- physfit.obj(gd)
  as_tibble_row(obj$coefficients)
}

physfit.prd <- function(lst, hum, temp) {
  physfit.obj(tibble(lst = lst, humidity = hum, temperature = temp))$fitted.values
}

physfit.vec <- function(lst, hum, temp) {
  rename(physfit(tibble(lst = lst, humidity = hum, temperature = temp)),
         Coef.LSTSigmoid = lst,
         Coef.Humidity = hum)
}

physfit.vecplus <- function(lst, hum, temp) {
  as_tibble_row(lm(temp ~ lst + hum)$coefficients) %>%
    rename(InterceptP = `(Intercept)`,
           Coef.LSTSigmoidP = lst,
           Coef.HumidityP = hum)
}

# Model functions
full.schema <- function(sche = krig.ssn, ma = krig.anom, rtn.model = FALSE,
                        use.max = FALSE) {
  function(indat) {
    ssn <- sche(indat, rtn.model = rtn.model)  # ssn(data) -> tibble(id, day, temp.doy)
    anom <- ma(indat, rtn.model = rtn.model, use.max = use.max)  # anom(data) -> tibble(id, date, temp.anom)
    
    if (rtn.model) {
      c(ssn, anom)
    } else {
      function(prdat) {
        doypr <- ssn(prdat)
        anpr <- anom(prdat)
        
        left_join(prdat, anpr, by = c("id", "date")) %>%
          left_join(doypr, by = c("id", "day")) %>%
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

krig.anom <- function(indat, rtn.model = FALSE, use.max = FALSE) {
  smoother <- \(x) x %>%
    drop_na() %>%
    group_by(id, day) %>%
    (
      \(y) {
        if ("temperature.max" %in% names(y) & "temperature" %in% names(y)) {
          mutate(y, temperature.plus = temperature.max - mean(temperature, na.rm = TRUE))
        } else y
      }
    ) %>%
    mutate(across(c(any_of("temperature"), lst, humidity),
                  ~(.x - mean(.x, na.rm = TRUE)))) %>%
    group_by(id) %>%
    arrange(date) %>%
    mutate(lst = sigsmooth(lst),
           humidity = humsmooth(humidity)) %>%
    drop_na()
  
  fitted <- smoother(indat) %>%
    group_by(id) %>%
    summarize(
      physfit.vec(lst, humidity, temperature),
      if (use.max) physfit.vecplus(lst, humidity, temperature.plus),
      across(c(lon, lat, water, elevation), mean),
      humidity_sd = sd(humidity)
    ) %>%
    ungroup() %>%
    drop_na()
  
  xer <- \(x) as.matrix(select(x, lon, lat))
  zer <- \(x) as.matrix(select(x, elevation, water, humidity_sd))
  
  fn <- function(y, ...) {
    if (rtn.model) {
      spatialProcess(xer(fitted), y, Z = zer(fitted), Distance = "rdist.earth")
    } else {
      spatialProcess(xer(fitted), y, Z = zer(fitted), Distance = "rdist.earth", ...)
    }
  }
  
  coef.lst <- fn(fitted$Coef.LSTSigmoid, aRange = 23, lambda = 1.1)
  coef.lstp <- if (use.max) fn(fitted$Coef.LSTSigmoidP, aRange = 28, lambda = 1.1) else NULL
  coef.humidity <- fn(fitted$Coef.Humidity, aRange = 23, lambda = 1.4)
  coef.humidityp <- if (use.max) fn(fitted$Coef.HumidityP, aRange = 34, lambda = 1.6) else NULL
  coef.intp <- if (use.max) fn(fitted$InterceptP, aRange = 32, lambda = 0.9) else NULL
  
  if (rtn.model) {
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
    function(prdat) {
      smdat <- prdat %>%
        smoother
      
      codat <- smdat %>%
        group_by(id) %>%
        summarize(
          across(c(lon, lat, water, elevation), mean),
          humidity_sd = sd(humidity)
        ) %>%
        drop_na() %>%
        ungroup()
      
      codat$lstco <- predict(coef.lst, xer(codat), Z = zer(codat))[,1]
      codat$humco <- predict(coef.humidity, xer(codat), Z = zer(codat))[,1]
      codat$lstpco <- if (use.max) predict(coef.lstp, xer(codat), Z = zer(codat))[,1] else NULL
      codat$humpco <- if (use.max) predict(coef.humidityp, xer(codat), Z = zer(codat))[,1] else NULL
      codat$intpco <- if (use.max) predict(coef.intp, xer(codat), Z = zer(codat))[,1] else NULL
      
      smdat <- left_join(smdat, codat, by = "id") %>%
        (\(y) { if (use.max) {
          mutate(y,
                 temp.anom = lst * lstco + humidity * humco,
                 temp.plus = intpco + lst * lstpco + humidity * humpco)
        } else {
          mutate(y,
                 temp.anom = lst * lstco + humidity * humco)
        }
        })
      
      select(smdat, id, date, temp.anom, any_of("temp.plus"))
    }
  }
}

krig.ssn <- function(indat, rtn.model = FALSE, ...) {
  costs <- \(day) cos((day - 210) * 2 * pi / 365)
  
  req.mos <- c("01", "03", "05", "07", "08", "09", "10", "11")
  mos <- \(x) unique(format(as.Date(x$day - 1, "1970-01-01"), "%m"))
  has.mos <- \(x) mean(req.mos %in% mos(x)) == 1
  
  if (!has.mos(indat))
    stop("Missing months in training data; some data from each month are required")
  
  preproc <- function(data, fit3s) {
    data %>%
      group_by(id, day) %>%
      summarize(
        across(c(lon, lat, elevation, date,
                 water, grassland, shrubland, barren), first),
        across(c(lst, humidity, any_of("temperature")), ~mean(.x, na.rm = TRUE))
      ) %>%
      group_by(id) %>%
      arrange(day) %>%
      mutate(
        maxT = max(lst, na.rm = TRUE),
        minT = min(lst, na.rm = TRUE),
        amphum = safely(\() lm(humidity ~ costs(day))$coefficients[[2]], otherwise = NA)()$result,
        meanT = mean(lst, na.rm = TRUE),
        mean_hum = mean(humidity, na.rm = TRUE),
        freeze_days = sum(lst < 0, na.rm = TRUE),
        cdd = cumsum(case_match(lst, NA ~ 0, .default = lst)) / cumsum(!is.na(lst)) * day
      ) %>%
      group_modify( ~ {
        if (fit3s) {
          view <- .x
          cbind(.x, fit.sins(.x))
        } else
          .x
      }) %>%
      group_by(id, month = format(date, "%m")) %>%
      summarize(
        across(-c(lst, humidity, cdd, day, any_of("temperature")), first),
        across(c(lst, humidity, cdd), ~mean(.x, na.rm = TRUE))
      ) %>%
      ungroup() %>%
      pivot_wider(names_from = "month", values_from = c("lst", "humidity", "cdd"),
                  id_cols = id, unused_fn = first) %>%
      drop_na(
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
  
  train <- preproc(indat, TRUE)
  
  xer <- \(x) as.matrix(select(x, lon, lat))
  zer.itx <- \(x) select(x, elevation, water, minT, maxT, freeze_days, lst_09, lst_10, humidity_01, humidity_03, humidity_09, humidity_11)
  zer.amp <- \(x) as.matrix(select(x, elevation, water, minT, amphum, cdd_09, cdd_10))
  zer.fw <- \(x) as.matrix(select(x, elevation, lst_05, humidity_03, humidity_05, grassland, barren, water, amphum, freeze_days))
  zer.ssu <- \(x) as.matrix(select(x, lst_05, humidity_07, elevation, water))
  zer.wid <- \(x) as.matrix(select(x, elevation, shrubland, grassland, barren, water, amphum, freeze_days, lst_05, humidity_03, humidity_05, humidity_08))
  
  fn <- function(y, z, ...) {
    if (rtn.model) {
      spatialProcess(xer(train), y, Z = z, Distance = "rdist.earth")
    } else {
      spatialProcess(xer(train), y, Z = z, Distance = "rdist.earth", ...)
    }
  }
  
  kix <- fn(train$Intercept, zer.itx(train), aRange = 5.7, lambda = 0.14)
  amp <- fn(train$Amplitude, zer.amp(train), aRange = 31, lambda = 0.6)
  fw <- fn(train$FallWinter, zer.fw(train), aRange = 25, lambda = 0.78)
  ssu <- fn(train$SpringSummer, cbind(zer.ssu(train), train$FallWinter), aRange = 15, lambda = 1.2)
  wid <- fn(train$WinterDay, zer.wid(train), aRange = 28, lambda = 2)
  
  if (rtn.model) {
    list(
      "Intercept" = kix,
      "Amplitude" = amp,
      "AutumnWinter" = fw,
      "SpringSummer" = ssu,
      "WinterDay" = wid
    )
  } else {
    function(data) {
      if (!has.mos(data)) {
        warning("Missing months in prediction data; some data from each month are required")
        tibble(id = NA, day = NA, temp.doy = NA)
      } else {
        ppd <- preproc(data, FALSE)
        ppd %>%
          group_by(id) %>%
          group_modify(~{
            fallwint <- predict(fw, xer(.x), Z = zer.fw(.x))[,1]
            sinfit <- tibble(
              Intercept = predict(kix, xer(.x), Z = zer.itx(.x))[,1],
              Amplitude = predict(amp, xer(.x), Z = zer.amp(.x))[,1],
              FallWinter = fallwint,
              SpringSummer = predict(ssu, xer(.x), Z = cbind(zer.ssu(.x), fallwint))[,1],
              WinterDay = predict(wid, xer(.x), Z = zer.wid(.x))[,1],
              SpringDay = 160,
              FallDay = 330,
              SummerDay = 220
            )
            make.sints(sinfit)
          }) %>%
          ungroup() %>%
          rename(temp.doy = actemp) %>%
          select(id, day, temp.doy)
      }
    }
  }
}

# Sinusoidal fit functions
applysin <- function(days, start, domain) {
  sin(((days - start) %% 365) * (2 * pi) / length(domain)) * (days %in% domain)
}

fit.sins <- function(data, mod = FALSE, aic = FALSE,
                     erv = tibble(
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
  inp <- data %>%
    group_by(day) %>%
    summarize(temperature = mean(temperature, na.rm = TRUE)) %>%
    drop_na() %>%
    arrange(day)
  
  if (nrow(inp) < 180) {
    warning("Insufficient data coverage for 3-sine fit")
    erv
  } else {
    day <- inp$day
    Ts <- inp$temperature
    
    index <- cos((day - 210) * 2 * pi / 365)
    cosfit <- lm(Ts ~ index)$coefficients
    meant <- cosfit[[1]]
    amplitude <- cosfit[[2]] / meant
    cospr <- 1 + index * amplitude
    anomaly <- Ts / meant - cospr
    
    rolled <- stats::filter(anomaly, rep(1/31, 31), circular = TRUE)
    
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
    
    sin1width <- round(((wint - fallt) %% 365)/2)
    sin1dom <- c((fallt - sin1width):366, 1:(wint + sin1width))
    
    sin2width <- round((sumt - spt)/2)
    sin2dom <- (spt - sin2width):(sumt + sin2width)
    
    sin1 <- -applysin(day, fallt - sin1width, sin1dom)
    sin2 <- -applysin(day, spt - sin2width, sin2dom)
    
    fullfit <- lm((Ts - meant * cospr) ~ 0 + sin1 + sin2)
    co <- fullfit$coefficients
    
    if (mod) {
      if (aic) {
        tibble(
          Sin = 2 * (2 - as.numeric(logLik(lm(Ts ~ index)))),
          Anom = 2 * (8 - as.numeric(logLik(fullfit)))
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
  day <- 1:366
  
  index <- cos((day - 210) * 2 * pi / 365)
  
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
  
  tibble(
    day = day,
    actemp = actemp
  )
}

# New function to graph temperature predictions along the river
graph_temp_predictions <- function(data, model_func, day_of_year) {
  data_with_predictions <- data %>%
    mutate(predicted_temp = model_func(data)$temp.doy)
  
  ggplot(data_with_predictions, aes(x = lon, y = lat, color = predicted_temp)) +
    geom_path() +
    scale_color_viridis_c(option = "magma") +
    labs(title = paste("Predicted Stream Temperature on Day", day_of_year),
         x = "Longitude",
         y = "Latitude",
         color = "Temperature (°C)") +
    theme_minimal()
}

# Example usage with your data
# Assuming 'data' is your dataset and 'model_func' is the trained model function
# day_of_year <- 150  # Example day of year
# graph_temp_predictions(data, model_func, day_of_year)
