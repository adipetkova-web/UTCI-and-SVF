#' Data analysis for UA
#' Route Park
#' Author: Adelina Petkova
#' Date: 09.2026

#### load libraries ####
library(ggplot2)
library(dplyr)
library(lubridate)

#### load data ####
data <- read.csv("data/Kestrel_park.csv", skip = 9, header = T)
data <- data[-1, ]

protocol <- read.csv("data/measurement_protocol_route_park.csv",  skip = 6, header = T)

protocol$start..UTC.<- as.POSIXct(paste(protocol$date, protocol$start..UTC.), 
                                  format = "%Y-%m-%d %H:%M:%S",
                                  tz = "UTC")

protocol$end..UTC. <- as.POSIXct(paste(protocol$date, protocol$end..UTC.), 
                                 format = "%Y-%m-%d %H:%M:%S",
                                 tz = "UTC")

data$Time <- as.POSIXct(data$Time, 
                                       format = "%Y-%m-%d %H:%M:%S",
                                       tz = "UTC")

svf <- read.csv("data/svf_fisheye_all_sites.csv", skip = 1)
svf <- svf[-1 , -1]
svf <- svf[ ,-2]

#### Filtering ####
filtered_data <- data.frame()

for(i in seq_len(nrow(protocol))) {
  
  res <- data[
    data$Time >= protocol$start..UTC.[i] &
      data$Time <= protocol$end..UTC.[i],
  ]
  
  if(nrow(res) > 0) {
    res$site.name <- protocol$site.name[i]
    res$site.no <- protocol$site.no[i]
    
    filtered_data <- rbind(filtered_data, res)
  }
}

# write.csv(filtered_data, "data/filetered_data_park.csv", row.names = F)


filtered_data <- filtered_data %>%
  arrange(Time) %>%
  filter({
    dt         <- as.POSIXct(Time, format = "%Y-%m-%d %H:%M:%S")
    gap        <- c(0, diff(as.numeric(dt)))
    block_id   <- cumsum(gap > 30)
    block_start <- ave(as.numeric(dt), block_id, FUN = min)
    as.numeric(dt) >= block_start + 5 * 60
  })

filtered_data<- filtered_data %>% filter(site.name != "courtyard tree")

filtered_data <- filtered_data %>%
  mutate(across(
    -c(Time, site.name, site.no),
    ~ as.numeric(.x)))

#### aggregate ####

aggregated_data <- filtered_data %>%
  mutate(hour = lubridate::floor_date(Time, "hour")) %>%
  group_by(site.name, hour) %>%
  summarise(across(where(is.numeric), ~ round(mean(.x, na.rm = TRUE), 2)),
            .groups = "drop")


#### Tmrt ####

calc_tmrt <- function(Tg, Ta, Va, D = 0.15, epsilon = 0.95) {
  Tmrt <- ((Tg + 273.15)^4 + 
             (1.1e8 * Va^0.6 / (epsilon * D^0.4)) * (Tg - Ta))^0.25 - 273.15
  return(Tmrt)
  }

aggregated_data <- aggregated_data %>%
  mutate(Tmrt = calc_tmrt(Tg = Globe.Temp, Ta = Temp, Va = Wind.Speed))

#### UTCI ####

source("calc_UTCI.R")

aggregated_data <- aggregated_data %>%
  mutate(UTCI = calc_utci(ta = aggregated_data$Temp, tmrt = aggregated_data$Tmrt,
          hur = aggregated_data$Rel..Hum., ws = aggregated_data$Wind.Speed))

#### SVF ####
model_data1 <- aggregated_data %>%
  left_join(svf, by = "site.name")

#### plots ####

#ggplot(data = filtered_data, aes(x = start_time, y = Celsius, color = site.name)) +
#  geom_line()

#ggplot(filtered_data, aes(x = site.name, y = Celsius)) +
 # geom_boxplot() +
  #labs(
   # x = "Site",
    #y = "Air temperature (°C)",
    #title = "Temperature distribution by site") +
  #theme_minimal()


