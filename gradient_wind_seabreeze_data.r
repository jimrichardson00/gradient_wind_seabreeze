# gradient_wind_seabreeze.r

options(scipen=3)
require(doParallel)
registerDoParallel(cores = 3)

# -------------------------------------------------------------------

setwd("/media/jim/FAT323/Projects/gradient_wind_seabreeze")
require(data.table)
metar <- fread("HM01X_Data_040842_999999999445572.txt"
  # , nrow = 0
  # , header = FALSE
  , stringsAsFactors = FALSE
  , data.table = FALSE)
# metar < - as.data.frame(metar, optional = TRUE)
head(metar)

setnames(metar, "hm", "hm")
setnames(metar, "Station Number", "stnnumb")
setnames(metar, "Day/Month/Year in DD/MM/YYYY format", "dateloc")
setnames(metar, "Hour24:Minutes  in HH24:MI format in Local standard time", "timeloc")
setnames(metar, "Precipitation in last 10 minutes in mm", "preci10")
setnames(metar, "Quality of precipitation in last 10 minutes", "prec10q")
setnames(metar, "Precipitation since 9am local time in mm", "prec9am")
setnames(metar, "Quality of precipitation since 9am local time", "prec9aq")
setnames(metar, "Air Temperature in degrees C", "tempera")
setnames(metar, "Quality of air temperature", "tempqua")
setnames(metar, "Wet bulb temperature in degrees C", "wetbulb")
setnames(metar, "Quality of Wet bulb temperature", "wetbulq")
setnames(metar, "Dew point temperature in degrees C", "dewpntt")
setnames(metar, "Quality of dew point temperature", "dewpntq")
setnames(metar, "Relative humidity in percentage %", "relhumi")
setnames(metar, "Quality of relative humidity", "relhumq")
setnames(metar, "Vapour pressure in hPa", "vappres")
setnames(metar, "Quality of vapour pressure", "vappreq")
setnames(metar, "Saturated vapour pressure in hPa", "satvapp")
setnames(metar, "Quality of saturated vapour pressure", "satvapq")
setnames(metar, "Wind speed in km/h", "windkmh")
setnames(metar, "Wind speed quality", "windkmq")
setnames(metar, "Wind direction in degrees true", "winddir")
setnames(metar, "Wind direction quality", "winddrq")
setnames(metar, "Speed of maximum windgust in last 10 minutes in  km/h", "windgkm")
setnames(metar, "Quality of speed of maximum windgust in last 10 minutes", "wingkmq")
setnames(metar, "Cloud amount(of first group) in eighths", "cldoct1")
setnames(metar, "Quality of first group of cloud amount", "cldamq1")
setnames(metar, "Cloud type(of first group) in in_words", "cldtyp1")
setnames(metar, "Quality of first group of cloud type", "cldtyq1")
setnames(metar, "Cloud height (of first group) in feet", "cldhgt1")
setnames(metar, "Quality of first group of cloud height", "cldhtq1")
setnames(metar, "Cloud amount(of second group) in eighths", "cldoct2")
setnames(metar, "Quality of second group of cloud amount", "cldamq2")
setnames(metar, "Cloud type(of second group) in in_words", "cldtyp2")
setnames(metar, "Quality of second group of cloud type", "cldtyq2")
setnames(metar, "Cloud height (of second group) in feet", "cldhgt2")
setnames(metar, "Quality of second group of cloud height", "cldhtq2")
setnames(metar, "Cloud amount(of third group) in eighths", "cldoct3")
setnames(metar, "Quality of third group of cloud amount", "cldamq3")
setnames(metar, "Cloud type(of third group) in in_words", "cldtyp3")
setnames(metar, "Quality of third group of cloud type", "cldtyq3")
setnames(metar, "Cloud height (of third group) in feet", "cldhgt3")
setnames(metar, "Quality of third group of cloud height", "cldhtq3")
setnames(metar, "Cloud amount(of fourth group) in eighths", "cldoct4")
setnames(metar, "Quality of fourth group of cloud amount", "cldamq4")
setnames(metar, "Cloud type(of fourth group) in in_words", "cldtyp4")
setnames(metar, "Quality of fourth group of cloud type", "cldtyq4")
setnames(metar, "Cloud height (of fourth group) in feet", "cldhgt4")
setnames(metar, "Quality of fourth group of cloud height", "cldhtq4")
setnames(metar, "Ceilometer cloud amount(of first group)", "cldocc1")
setnames(metar, "Quality of first group of ceilometer cloud amount", "cldacq1")
setnames(metar, "Ceilometer cloud height (of first group) in feet", "cldhtc1")
setnames(metar, "Quality of first group of ceilometer cloud height", "cldhcq1")
setnames(metar, "Ceilometer cloud amount(of second group)", "cldocc2")
setnames(metar, "Quality of second group of ceilometer cloud amount", "cldacq2")
setnames(metar, "Ceilometer cloud height (of second group) in feet", "cldhtc2")
setnames(metar, "Quality of second group of ceilometer cloud height", "cldhcq2")
setnames(metar, "Ceilometer cloud amount(of third group)", "cldocc3")
setnames(metar, "Quality of third group of ceilometer cloud amount", "cldacq3")
setnames(metar, "Ceilometer cloud height (of third group) in feet", "cldhtc3")
setnames(metar, "Quality of third group of ceilometer cloud height", "cldhcq3")
setnames(metar, "Horizontal visibility in km", "visobsk")
setnames(metar, "Quality of horizontal visibility", "visobsq")
setnames(metar, "Direction of minimum visibility in degrees", "visobmn")
setnames(metar, "Quality of direction of minimum visibility", "visobmq")
setnames(metar, "AWS visibility in km", "visawsk")
setnames(metar, "Quality of AWS(Automatic Weather Station) visibility", "visawsq")
setnames(metar, "Present weather in text", "presewx")
setnames(metar, "Quality of present weather", "preswxq")
setnames(metar, "Intensity of first present weather in text", "prwxin1")
setnames(metar, "Quality of intensity of first present weather", "prwxiq1")
setnames(metar, "Descriptor of first present weather in text", "prwxde1")
setnames(metar, "Quality of descriptor of first present weather", "prwxdq1")
setnames(metar, "Type of first present weather in text", "prwxty1")
setnames(metar, "Quality of type of first present weather", "prwxtq1")
setnames(metar, "Intensity of second present weather in text", "prwxin2")
setnames(metar, "Quality of intensity of second present weather", "prwxiq2")
setnames(metar, "Descriptor of second present weather in text", "prwxde2")
setnames(metar, "Quality of descriptor of second present weather", "prwxdq2")
setnames(metar, "Type of second present weather in text", "prwxty2")
setnames(metar, "Quality of type of second present weather", "prwxtq2")
setnames(metar, "Intensity of third present weather in text", "prwxin3")
setnames(metar, "Quality of intensity of third present weather", "prwxiq3")
setnames(metar, "Descriptor of third present weather in text", "prwxde3")
setnames(metar, "Quality of descriptor of third present weather", "prwxdq3")
setnames(metar, "Type of third present weather in text", "prwxty3")
setnames(metar, "Quality of type of third present weather", "prwxtq3")
setnames(metar, "Descriptor of first recent weather in text", "rewxde1")
setnames(metar, "Quality of descriptor of first recent weather", "rewxdq1")
setnames(metar, "Type of first recent weather in text", "rewxty1")
setnames(metar, "Quality of type of first recent weather", "rewxtq1")
setnames(metar, "Descriptor of second recent weather in text", "rewxde2")
setnames(metar, "Quality of descriptor of second recent weather", "rewxdq2")
setnames(metar, "Type of second recent weather in text", "rewxty2")
setnames(metar, "Quality of type of second recent weather", "rewxtq2")
setnames(metar, "Descriptor of third recent weather in text", "rewxde3")
setnames(metar, "Quality of descriptor of third recent weather", "rewxdq3")
setnames(metar, "Type of third recent weather in text", "rewxty3")
setnames(metar, "Quality of type of third recent weather", "rewxtq3")
setnames(metar, "AWS present weather in text", "prwxaws")
setnames(metar, "Quality of AWS present weather", "prwxawq")
setnames(metar, "AWS weather for last 15 minutes in text", "awwx15m")
setnames(metar, "Quality of AWS weather for last 15 minutes", "awwx15q")
setnames(metar, "AWS weather for last 60 minutes in text", "awwx60m")
setnames(metar, "Quality of AWS weather for last 60 minutes", "awwx60q")
setnames(metar, "Mean sea level pressure in hPa", "mslphPa")
setnames(metar, "Quality of mean sea level pressure", "mslphPq")
setnames(metar, "Station level pressure in hPa", "stprehP")
setnames(metar, "Quality of station level pressure", "stprehq")
setnames(metar, "QNH pressure in hPa", "QNHphPa")
setnames(metar, "Quality of QNH pressure", "QNHhPaq")
setnames(metar, "AWS Flag", "awsflag")
setnames(metar, "#", "x")

metar$visobsm <- metar$visobsk*1000
metar$windkno <- as.numeric(metar$windkmh)*0.539957

# -------------------------------------------------------------------

# creates year, month, day columns from date columns
require(stringr)
metar$year_loc <- str_match(as.character(metar$dateloc), pattern = "(\\d*)/(\\d*)/(\\d*)")[, 4]
metar$month_loc <- str_match(as.character(metar$dateloc), pattern = "(\\d*)/(\\d*)/(\\d*)")[, 3]
metar$day_loc <- str_match(as.character(metar$dateloc), pattern = "(\\d*)/(\\d*)/(\\d*)")[, 2]
metar$hour_loc <- str_match(as.character(metar$timeloc), pattern = "(\\d*)\\:(\\d*)")[, 2]
metar$minute_loc <- str_match(as.character(metar$timeloc), pattern = "(\\d*)\\:(\\d*)")[, 3]
head(metar)

# creates new date column dateloc as date type 
metar$dateloc <- paste(metar$year_loc, "-", metar$month_loc, "-", metar$day_loc, sep = "")
head(metar)

# creates new datetime column dateloc as datetime type 
head(paste(metar$dateloc, " ", metar$timeloc, sep = ""))
metar$datetime_loc <- as.POSIXct(paste(metar$dateloc, " ", metar$timeloc, sep = ""), tz = "Australia/Brisbane", format = "%Y-%m-%d %H:%M")

# pulls out year, month, day etc from loc datetime
require(lubridate)
metar$year_loc <- as.integer(metar$year_loc)
metar$month_loc <- as.integer(metar$month_loc)
metar$day_loc <- as.integer(metar$day_loc)
metar$hour_loc <- as.integer(metar$hour_loc)
metar$minute_loc <- as.integer(metar$minute_loc)

# creates datetime utc by converting from datetime loc
metar$datetime_utc <- format(metar$datetime_loc, usetz = TRUE, tz = "UTC", format = "%Y-%m-%d %H:%M")
metar$datetime_utc <- as.POSIXct(metar$datetime_utc, tz = "UTC", format = "%Y-%m-%d %H:%M")
head(metar$datetime_utc)

# # pulls out year, month, day etc from utc datetime
# require(lubridate)
# metar$year_utc <- year(metar$year_utc)
# metar$month_utc <- month(metar$month_utc)
# metar$day_utc <- day(metar$day_utc)
# metar$hour_utc <- hour(metar$hour_utc)
# metar$minute_utc <- minute(metar$minute_utc)

head(metar)

# -------------------------------------------------------------------

setwd("/media/jim/FAT323/Projects/gradient_wind_seabreeze")
require(data.table)
upperair <- fread("UA01D_Data_040842_999999999446605.txt"
  # , nrow = 0
  # , header = FALSE
  , stringsAsFactors = FALSE
  , data.table = FALSE)
head(upperair)

setnames(upperair, "ua", "ua" )
setnames(upperair, "Station Number", "stnnumb" )
setnames(upperair, "Station Name", "stnname" )
setnames(upperair, "Locality", "localit" )
setnames(upperair, "State", "stateau" )
setnames(upperair, "Latitude", "latitud" )
setnames(upperair, "Longitude", "longitu" )
setnames(upperair, "Month/Year site opened (MM/YYYY)", "monyeao" )
setnames(upperair, "Month/Year site closed (MM/YYYY)", "monyeac" )
setnames(upperair, "WMO Index Number", "wmoindn" )
setnames(upperair, "Rainfall district code", "raindst" )
setnames(upperair, "River station ID", "rivstid" )
setnames(upperair, "Aviation ID", "aviatid" )
setnames(upperair, "Height above MSL", "habvmsl" )
setnames(upperair, "Day/Month/Year Hour24:Minutes in DD/MM/YYYY HH24:MI format in Local Time", "datetime_loc" )
setnames(upperair, "Day/Month/Year Hour24:Minutes in DD/MM/YYYY HH24:MI format in Local Standard Time", "datetime_std" )
setnames(upperair, "Day/Month/Year Hour24:Minutes in DD/MM/YYYY HH24:MI format in UTC - Coordinated Universal Time", "datetime_utc" )
setnames(upperair, "Air temperature in Degrees C", "tempera" )
setnames(upperair, "Quality of air temperature", "tempqua" )
setnames(upperair, "Dew point temperature in Degrees C", "dewpntt" )
setnames(upperair, "Quality of dew point temperature", "dewpntq" )
setnames(upperair, "Relative humidity in percentage %", "relhumi" )
setnames(upperair, "Quality of relative humidity", "relhumq" )
setnames(upperair, "Wind speed measured in knots", "windkno" )
setnames(upperair, "Quality of wind speed", "windknq" )
setnames(upperair, "Wind direction measured in degrees", "winddir" )
setnames(upperair, "Quality of wind direction", "winddiq" )
setnames(upperair, "Pressure in hPa", "preshPa" )
setnames(upperair, "Quality of pressure", "preshPq" )
setnames(upperair, "Geopotential height in gpm to nearest 0.1m", "geopgpm" )
setnames(upperair, "Quality of geopotential height", "geopgpq" )
setnames(upperair, "Level type", "levelt" )
setnames(upperair, "#", "x" )

head(upperair)

upperair$geopgpf <- 3.28084*as.numeric(upperair$geopgpm)

# for(name in names(upperair)) {
#   print(name)
# }

# -------------------------------------------------------------------

require(stringr)
# creates datetime utc by converting from datetime_loc, DD/MM/YYYY HH24:MI
upperair$year_loc <- str_match(as.character(upperair$datetime_loc), pattern = "(\\d*)/(\\d*)/(\\d*) (\\d*)\\:(\\d*)")[, 4]
upperair$month_loc <- str_match(as.character(upperair$datetime_loc), pattern = "(\\d*)/(\\d*)/(\\d*) (\\d*)\\:(\\d*)")[, 3]
upperair$day_loc <- str_match(as.character(upperair$datetime_loc), pattern = "(\\d*)/(\\d*)/(\\d*) (\\d*)\\:(\\d*)")[, 2]
upperair$hour_loc <- str_match(as.character(upperair$datetime_loc), pattern = "(\\d*)/(\\d*)/(\\d*) (\\d*)\\:(\\d*)")[, 5]
upperair$minute_loc <- str_match(as.character(upperair$datetime_loc), pattern = "(\\d*)/(\\d*)/(\\d*) (\\d*)\\:(\\d*)")[, 6]
head(upperair)

paste(upperair$year_loc, "-", upperair$month_loc, "-", upperair$day_loc, " ", upperair$hour_loc, ":", upperair$minute_loc, sep = "")
upperair$datetime_loc <- format(paste(upperair$year_loc, "-", upperair$month_loc, "-", upperair$day_loc, " ", upperair$hour_loc, ":", upperair$minute_loc, sep = ""), usetz = TRUE, tz = "Australia/Brisbane", format = "%Y-%m-%d %H:%M")
upperair$datetime_loc <- as.POSIXct(upperair$datetime_loc, tz = "Australia/Brisbane", format = "%Y-%m-%d %H:%M")
upperair$datetime_loc <- format(upperair$datetime_loc, usetz = TRUE, tz = "Australia/Brisbane", format = "%Y-%m-%d %H:%M")

# creates new date column dateloc as date type 
upperair$dateloc <- paste(upperair$year_loc, "-", upperair$month_loc, "-", upperair$day_loc, sep = "")
head(upperair$dateloc)

# creates new date column dateloc as date type 
upperair$timeloc <- paste(formatC(x = upperair$hour_loc, flag = "0", width = 2), ":", formatC(x = upperair$minute_loc, flag = "0", width = 2), sep = "")
head(upperair$timeloc)

upperair$year_loc <- as.integer(upperair$year_loc)
upperair$month_loc <- as.integer(upperair$month_loc)
upperair$day_loc <- as.integer(upperair$day_loc)
upperair$hour_loc <- as.integer(upperair$hour_loc)
upperair$minute_loc <- as.integer(upperair$minute_loc)

# -------------------------------------------------------------------

require(stringr)
# creates datetime utc by converting from datetime_loc, DD/MM/YYYY HH24:MI
upperair$year_utc <- str_match(as.character(upperair$datetime_utc), pattern = "(\\d*)/(\\d*)/(\\d*) (\\d*)\\:(\\d*)")[, 4]
upperair$month_utc <- str_match(as.character(upperair$datetime_utc), pattern = "(\\d*)/(\\d*)/(\\d*) (\\d*)\\:(\\d*)")[, 3]
upperair$day_utc <- str_match(as.character(upperair$datetime_utc), pattern = "(\\d*)/(\\d*)/(\\d*) (\\d*)\\:(\\d*)")[, 2]
upperair$hour_utc <- str_match(as.character(upperair$datetime_utc), pattern = "(\\d*)/(\\d*)/(\\d*) (\\d*)\\:(\\d*)")[, 5]
upperair$minute_utc <- str_match(as.character(upperair$datetime_utc), pattern = "(\\d*)/(\\d*)/(\\d*) (\\d*)\\:(\\d*)")[, 6]
head(upperair)

paste(upperair$year_utc, "-", upperair$month_utc, "-", upperair$day_utc, " ", upperair$hour_utc, ":", upperair$minute_utc, sep = "")
upperair$datetime_utc <- format(paste(upperair$year_utc, "-", upperair$month_utc, "-", upperair$day_utc, " ", upperair$hour_utc, ":", upperair$minute_utc, sep = ""), usetz = TRUE, tz = "Australia/Brisbane", format = "%Y-%m-%d %H:%M")
upperair$datetime_utc <- as.POSIXct(upperair$datetime_utc, tz = "Australia/Brisbane", format = "%Y-%m-%d %H:%M")
upperair$datetime_utc <- format(upperair$datetime_utc, usetz = TRUE, tz = "Australia/Brisbane", format = "%Y-%m-%d %H:%M")

# creates new date column dateloc as date type 
upperair$dateloc <- paste(upperair$year_utc, "-", upperair$month_utc, "-", upperair$day_utc, sep = "")
head(upperair$dateloc)

# creates new date column dateloc as date type 
upperair$timeloc <- paste(formatC(x = upperair$hour_loc, flag = "0", width = 2), ":", formatC(x = upperair$minute_loc, flag = "0", width = 2), sep = "")
head(upperair$timeloc)

upperair$year_utc <- as.integer(upperair$year_utc)
upperair$month_utc <- as.integer(upperair$month_utc)
upperair$day_utc <- as.integer(upperair$day_utc)
upperair$hour_utc <- as.integer(upperair$hour_utc)
upperair$minute_utc <- as.integer(upperair$minute_utc)

# ------------------------------------------

head(upperair[upperair$timeloc == "09:00", ])
head(upperair)

upperair$dateloc <- as.character(upperair$dateloc)
unique(upperair$dateloc)
unique(upperair[upperair$timeloc == "09:00", ]$dateloc)

require(plyr)
upperair_900hPa <- ddply(.data = upperair[upperair$timeloc == "09:00", ]
	, .variables = c("dateloc")
	, .fun = summarize
	, preshPa_upp = preshPa[which.min(abs(preshPa - 900))]
	, winddir_upp = winddir[which.min(abs(preshPa - 900))]
	, windkno_upp = windkno[which.min(abs(preshPa - 900))]
	)
upperair_900hPa$windkno_upp_rnd <- round(upperair_900hPa$windkno_upp)
head(upperair_900hPa)

# ----------------------------------------------------------------

head(metar)
metar_summerday <- metar[is.na(metar$hour_loc) == FALSE, ]
metar_summerday <- metar[is.na(metar$month_loc) == FALSE, ]
metar_summerday <- metar_summerday[(metar_summerday$hour_loc >= 9) & (metar_summerday$hour_loc <= 18), ]
metar_summerday <- metar_summerday[(metar_summerday$month_loc %in% c(8, 9, 10, 11, 12, 1, 2)), ]
head(metar_summerday)

metar_summerday_df <- ddply(.data = metar_summerday
	, .variables = c("dateloc")
	, .fun = summarize
	, windkno_met = windkno[which.max(as.numeric(na.omit(windkno)))][1]
	, winddir_met = winddir[which.max(as.numeric(na.omit(windkno)))][1]
	, windkno_utc = datetime_utc[which.max(as.numeric(na.omit(windkno)))][1]
	, windkno_avg = mean(as.numeric(na.omit(windkno)))
	, winddir_avg = sum(windkno*winddir)/sum(windkno)
	, tempera_met = tempera[which.max(as.numeric(na.omit(tempera)))][1]
	, tempera_utc = datetime_utc[which.max(as.numeric(na.omit(tempera)))][1]
	, dewpntt_met = dewpntt[which.min(as.numeric(na.omit(dewpntt)))][1]
	, dewpntt_utc = datetime_utc[which.min(as.numeric(na.omit(dewpntt)))][1]
	, windwes_met = ifelse(any(winddir >= 230 & winddir <= 330 & windkno >= 5), 1, 0)
	, n_records_met = length(unique(datetime_utc))
	)
# metar_summerday_df$windwes <- ifelse(metar_summerday_df$winddir_avg >= 230 & metar_summerday_df$winddir_avg >= 330, 1, 0)
head(metar_summerday_df)

metar_summerday_df$dateloc <- as.character(metar_summerday_df$dateloc)
upperair_900hPa$dateloc <- as.character(upperair_900hPa$dateloc)

metar_summerday_mer <- merge(x = metar_summerday_df, y = upperair_900hPa, by = "dateloc")
head(metar_summerday_mer)

metar_summerday_mer$winddir_dif <- metar_summerday_mer$winddir_upp - metar_summerday_mer$winddir_avg
metar_summerday_mer$winddir_dif <- unlist(lapply(metar_summerday_mer$winddir_dif, FUN = deg))
head(metar_summerday_mer)

head(metar_summerday_mer)
setwd("/home/jim/Dropbox/Projects/gradient_wind_seabreeze")
write.csv(metar_summerday_mer, "metar_summerday_mer.csv")

metar_summerday_sum <- ddply(.data = metar_summerday_mer
	, .variables = c("windkno_upp_rnd")
	, .fun = summarize
	, windkno_met = mean(na.omit(windkno_met))
	, winddir_met = mean(na.omit(winddir_met))
	, windkno_utc = mean(na.omit(windkno_utc))
	, windkno_avg = mean(na.omit(windkno_avg))
	, winddir_avg = mean(na.omit(winddir_avg))
	, tempera_met = mean(na.omit(tempera_met))
	, dewpntt_met = mean(na.omit(dewpntt_met))
	, winddir_dif = mean(na.omit(winddir_dif))
	, windwes_met = mean(na.omit(windwes_met))
	, n_recor_sum = length(dateloc)
	)
head(metar_summerday_sum)

setwd("/home/jim/Dropbox/Projects/gradient_wind_seabreeze")
write.csv(metar_summerday_sum, "metar_summerday_sum.csv")

png("winddir_avg.png", width = 1000, height = 1000, pointsize = 30)
plot(metar_summerday_sum$windkno_upp_rnd, metar_summerday_sum$winddir_avg)
lines(metar_summerday_sum$windkno_upp_rnd, metar_summerday_sum$n_recor_sum)
dev.off()

png("winddir_met.png", width = 1000, height = 1000, pointsize = 30)
plot(metar_summerday_sum$windkno_upp_rnd, metar_summerday_sum$winddir_met)
lines(metar_summerday_sum$windkno_upp_rnd, metar_summerday_sum$n_recor_sum)
dev.off()

png("windwes_met.png", width = 1000, height = 1000, pointsize = 30)
plot(metar_summerday_sum$windkno_upp_rnd, metar_summerday_sum$windwes_met)
lines(metar_summerday_sum$windkno_upp_rnd, metar_summerday_sum$n_recor_sum)
dev.off()

png("winddir_dif.png", width = 1000, height = 1000, pointsize = 30)
plot(metar_summerday_sum$windkno_upp_rnd, metar_summerday_sum$winddir_dif)
lines(metar_summerday_sum$windkno_upp_rnd, metar_summerday_sum$n_recor_sum)
dev.off()

png("tempera_met.png", width = 1000, height = 1000, pointsize = 30)
plot(metar_summerday_sum$windkno_upp_rnd, metar_summerday_sum$tempera_met)
lines(metar_summerday_sum$windkno_upp_rnd, metar_summerday_sum$n_recor_sum)
dev.off()

png("dewpntt_met.png", width = 1000, height = 1000, pointsize = 30)
plot(metar_summerday_sum$windkno_upp_rnd, metar_summerday_sum$dewpntt_met)
lines(metar_summerday_sum$windkno_upp_rnd, metar_summerday_sum$n_recor_sum)
dev.off()

png("winddir_avg.png", width = 1000, height = 1000, pointsize = 30)
plot(metar_summerday_sum$windkno_upp_rnd, metar_summerday_sum$windwes)
lines(metar_summerday_sum$windkno_upp_rnd, metar_summerday_sum$n_recor_sum)
dev.off()

# ----------------------------------------------------------------

metar_westerly <- metar[is.na(metar$winddir) == FALSE, ]
metar_westerly <- metar_westerly[(metar_westerly$winddir >= 230 & metar_westerly$winddir <= 330), ]
metar_westerly <- metar_westerly[(metar_westerly$windkno >= 5), ]
metar_westerly <- metar_westerly[(as.numeric(metar_westerly$hour_loc) >= 9) & (as.numeric(metar_westerly$hour_loc) <= 18), ]
metar_westerly <- metar_westerly[(as.integer(metar_westerly$month_loc) %in% c(8, 9, 10, 11, 12, 1, 2)), ]
head(metar_westerly)

metar_westerly_df <- ddply(.data = metar_westerly
	, .variables = c("dateloc")
	, .fun = summarize
	, winddir_met = winddir[which.max(as.numeric(na.omit(windkno)))][1]
	, windkno_met = windkno[which.max(as.numeric(na.omit(windkno)))][1]
	, timeloc_met = timeloc[which.max(as.numeric(na.omit(windkno)))][1]
	, n_records_met = length(unique(datetime_loc))
	)
head(metar_westerly_df)

metar_westerly_df <- metar_westerly_df[metar_westerly_df$n_records >= 5, ]
head(metar_westerly_df)

metar_westerly_mer <- merge(x = metar_westerly_df, y = upperair_900hPa, by = "dateloc")
head(metar_westerly_mer)

setwd("/home/jim/Dropbox/Projects/gradient_wind_seabreeze")
write.csv(metar_westerly_mer, "metar_westerly_mer.csv")

ecdf_westerly <- ecdf(metar_westerly_mer$windkno_upp)
plot(seq(1, 30, 1), ecdf_westerly(seq(1, 30, 1)), type = "l")

# ----------------------------------------------------------------

metar_temperature <- metar[is.na(metar$tempera) == FALSE, ]
metar_temperature <- metar_temperature[(metar_temperature$tempera >= 30), ]
metar_temperature <- metar_temperature[(metar_temperature$hour_loc >= 9) & (metar_temperature$hour_loc <= 18), ]
metar_temperature <- metar_temperature[(metar_temperature$month_loc %in% c(8, 9, 10, 11, 12, 1, 2)), ]
head(metar_temperature)

metar_temperature_df <- ddply(.data = metar_temperature
	, .variables = c("dateloc")
	, .fun = summarize
	, winddir_met = winddir[which.max(as.numeric(na.omit(tempera)))][1]
	, windkno_met = windkno[which.max(as.numeric(na.omit(tempera)))][1]
	, dewpntt_met = windkno[which.max(as.numeric(na.omit(tempera)))][1]
	, tempera_met = windkno[which.max(as.numeric(na.omit(tempera)))][1]
	, timeloc_met = timeloc[which.max(as.numeric(na.omit(tempera)))][1]
	, n_records_met = length(unique(datetime_loc))
	)
head(metar_temperature_df)

metar_temperature_df <- metar_temperature_df[metar_temperature_df$n_records >= 5, ]
head(metar_temperature_df)

metar_temperature_mer <- merge(x = metar_temperature_df, y = upperair_900hPa, by = "dateloc")
head(metar_temperature_mer)

setwd("/home/jim/Dropbox/Projects/gradient_wind_seabreeze")
write.csv(metar_temperature_mer, "metar_temperature_mer.csv")
head(metar_temperature_mer)

ecdf_temperature <- ecdf(metar_temperature_mer$windkno_upp)
plot(seq(1, 30, 1), ecdf_temperature(seq(1, 30, 1)), type = "l")

# ----------------------------------------------------------------

metar_dewpoint <- metar[is.na(metar$dewpntt) == FALSE, ]
metar_dewpoint <- metar_dewpoint[metar_dewpoint$dewpntt <= 10, ]
metar_dewpoint <- metar_dewpoint[(metar_dewpoint$hour_loc >= 9) & (metar_dewpoint$hour_loc <= 18), ]
metar_dewpoint <- metar_dewpoint[(metar_dewpoint$month_loc %in% c(8, 9, 10, 11, 12, 1, 2)), ]
head(metar_dewpoint)

metar_dewpoint_df <- ddply(.data = metar_dewpoint
	, .variables = c("dateloc")
	, .fun = summarize
	, winddir_met = winddir[which.min(as.numeric(na.omit(dewpntt)))][1]
	, windkno_met = windkno[which.min(as.numeric(na.omit(dewpntt)))][1]
	, dewpntt_met = windkno[which.min(as.numeric(na.omit(dewpntt)))][1]
	, tempera_met = windkno[which.min(as.numeric(na.omit(dewpntt)))][1]
	, timeloc_met = timeloc[which.min(as.numeric(na.omit(dewpntt)))][1]
	, n_records_met = length(unique(datetime_loc))
	)
head(metar_dewpoint_df)

metar_dewpoint_df <- metar_dewpoint_df[metar_dewpoint_df$n_records >= 5, ]
head(metar_dewpoint_df)

metar_dewpoint_mer <- merge(x = metar_dewpoint_df, y = upperair_900hPa, by = "dateloc")
head(metar_dewpoint_mer)

setwd("/home/jim/Dropbox/Projects/gradient_wind_seabreeze")
write.csv(metar_dewpoint_mer, "metar_dewpoint_mer.csv")
head(metar_dewpoint_mer)

ecdf_dewpoint <- ecdf(metar_dewpoint_mer$windkno_upp)
plot(seq(1, 30, 1), ecdf_dewpoint(seq(1, 30, 1)), type = "l")


