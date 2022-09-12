library(hkweather)
enq_eqx = function(year){
  
}
enq_rjd = function(day, tz = "HongKong"){
  day = as.POSIXct(day, tz = tz)
  year = year(day)
  days = floor(day - ISOdatetime(year = year, 1, 1, 0, 0, 0))+1
  message(days)
}
draw_sele = function(lon, lat, RJD, d_year = 365, tz = 0, sum_solstice = 172){
  hkweather::hkw_lib()
  phi_r = 23.44 * pi /180
  dr    = sum_solstice
  del_S =  phi_r * cos(2 * pi * (RJD - dr) / d_year)
  lon_rad = lon * pi/180
  lat_rad = lat * pi/180
  
  
  data = data.frame(Time_l = seq(0, 23, 0.5)) %>%
    mutate(UTC = Time_l - tz) %>%
    mutate(elevation = asin(sin(lat_rad)*sin(del_S)-cos(lat_rad)*cos(del_S)*cos(2*pi*UTC/24+lon_rad))) %>%
    mutate(elevation = elevation / pi *180)
  
  data %>%
    ggplot(aes(x = Time_l, y = elevation)) +
    geom_path() +
    theme_bw() +
    scale_x_continuous(breaks = seq(0, 23, 1))+
    scale_y_continuous(breaks = seq(-90, +90, 10))+
    coord_cartesian(expand = F)+
    labs(x = paste0("Local time, at ", lon, "°, ", lat, "°", ", RJD = ", RJD), y = "Solar elevation (degrees)")
}
draw_planck = function(temp, wave_range = c(0, 50)){
  C1 = 3.74*10^8
  C2 = 1.44*10^4
  
  data_wavelength = data.frame(wavelength = seq(min(wave_range), max(wave_range), 0.1))
  data = data.frame(temp = temp) %>%
    expand_grid(data_wavelength) %>%
    mutate(watt = C1 / (wavelength ^ 5 * (exp(C2 / (wavelength * temp)) - 1)))

  data %>%
    ggplot(aes(x = wavelength, y = watt, group = temp, color = as.factor(temp)))+
    geom_path()+
    theme_bw()+
    labs(y = "Planck radiant exitance (W/m2 um-1, log-scale)",
         x = "Wavelength (um, log-scale)",
         color = "Temperature\n(K)")+
    scale_y_continuous(trans = "log10")+
    scale_x_continuous(trans = "log10")
}
draw_aira_noon = function(lon, lat, year = 365, perihelion = 2, sum_solstice = 172){
  #additional variable
  C = 6.283185
  P = 365.256363
  dP = perihelion
  dR = sum_solstice
  e = 0.0167
  a = 149.457
  R = 149.6
  phi_r = (23.44) * pi / 180
  S0 = 1361
  
  loc_lon = 0
  loc_lon_rad = loc_lon * pi / 180
  
  data_lat = data.frame(lat = lat)
  data = data.frame(RJD = 1:year) %>%
    mutate(M = C * (RJD - dP) / P) %>%
    mutate(v = M+0.0333988*sin(M)+0.0003486*sin(2*M)+0.000005*sin(3*M)) %>%
    mutate(RR = a * (1 - e^2) / (1 + e * cos(v))) %>%
    mutate(S = 1361 * (RR/R)^2) %>%
    expand_grid(data_lat) %>%
    mutate(lat_rad = lat * pi / 180) %>%
    mutate(del_s   = phi_r * cos(C*(RJD-dR)/365)) %>%
    mutate(psi     = sin(lat_rad)*sin(del_s)-cos(lat_rad)*cos(del_s)*cos(C/2+loc_lon_rad)) %>%
    mutate(Frad    = S0*sin(psi))
  
  plot1 = data %>%
    ggplot(aes(x = RJD, y = S)) +
    geom_path()+
    theme_bw()+
    scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 365))+
    labs(y = "Solar Irradiance, top of atmosphere (W/m2)",
         x = "RJD")
  
  plot2 = data %>%
    ggplot(aes(x = RJD, y = Frad, group = lat, color = lat))+
    geom_path()+
    theme_bw()+
    scale_color_gradientn(colours = c("Blue", "Red", "Green"), breaks = seq(-90, 90, 45))+
    scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 365))+
    labs(y = "Radiative flux (W/m2)",
         color = "Latitude")
  return(list(list(plot1), list(plot2)))
}


enq_rjd(Sys.time())
leap_year(Sys.time())
draw_sele(lon = 114.1694, lat = 22.3193, tz = 8, RJD = 255)
ggsave("Local Solar Elevation.svg", width = 6, height = 4)
draw_planck(temp = c(6000, 1000, 250), wave_range = c(0, 1000))
ggsave("Planck curve.svg", width = 6, height = 4)
plot = draw_aira_noon(lon = 0, lat = c(-90, -45, 0, +45, +90))
plot[1]
ggsave("Solar irrdaiance aoa.svg", width = 6, height = 4)
plot[2]
ggsave("Radiative flux.svg", width = 6, height = 4)

