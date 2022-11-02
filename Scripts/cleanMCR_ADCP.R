# For editing the large MCR ADCP offshore datasets
library(tidyverse)
library(here)
library(lubridate)

data<-read_table(here("Data","IslandData","MCR_LTER_CTD_ADCP_FOR05.txt"), col_names = c("matlab_datenum", "year", "month", "decimal_day", "U_Vel_02.5m_bin", "U_Vel_03.5m_bin", "U_Vel_04.5m_bin", "U_Vel_05.5m_bin", "U_Vel_06.5m_bin", "U_Vel_07.5m_bin", "U_Vel_08.5m_bin", "U_Vel_09.5m_bin", "U_Vel_10.5m_bin", "U_Vel_11.5m_bin", "U_Vel_12.5m_bin", "U_Vel_13.5m_bin", "U_Vel_14.5m_bin", "U_Vel_15.5m_bin", "U_Vel_16.5m_bin", "U_Vel_17.5m_bin", "U_Vel_18.5m_bin", "U_Vel_19.5m_bin", "U_Vel_20.5m_bin", "U_Vel_21.5m_bin", "V_Vel_02.5m_bin", "V_Vel_03.5m_bin", "V_Vel_04.5m_bin", "V_Vel_05.5m_bin", "V_Vel_06.5m_bin", "V_Vel_07.5m_bin", "V_Vel_08.5m_bin", "V_Vel_09.5m_bin", "V_Vel_10.5m_bin", "V_Vel_11.5m_bin", "V_Vel_12.5m_bin", "V_Vel_13.5m_bin", "V_Vel_14.5m_bin", "V_Vel_15.5m_bin", "V_Vel_16.5m_bin", "V_Vel_17.5m_bin", "V_Vel_18.5m_bin", "V_Vel_19.5m_bin", "V_Vel_20.5m_bin", "V_Vel_21.5m_bin", "Intensity_02.5m_bin", "Intensity_03.5m_bin", "Intensity_04.5m_bin", "Intensity_05.5m_bin", "Intensity_06.5m_bin", "Intensity_07.5m_bin", "Intensity_08.5m_bin", "Intensity_09.5m_bin", "Intensity_10.5m_bin", "Intensity_11.5m_bin", "Intensity_12.5m_bin", "Intensity_13.5m_bin", "Intensity_14.5m_bin", "Intensity_15.5m_bin", "Intensity_16.5m_bin", "Intensity_17.5m_bin", "Intensity_18.5m_bin", "Intensity_19.5m_bin", "Intensity_20.5m_bin", "Intensity_21.5m_bin", "percentGood_02.5m_bin", "percentGood_03.5m_bin", "percentGood_04.5m_bin", "percentGood_05.5m_bin", "percentGood_06.5m_bin", "percentGood_07.5m_bin", "percentGood_08.5m_bin", "percentGood_09.5m_bin", "percentGood_10.5m_bin", "percentGood_11.5m_bin", "percentGood_12.5m_bin", "percentGood_13.5m_bin", "percentGood_14.5m_bin", "percentGood_15.5m_bin", "percentGood_16.5m_bin", "percentGood_17.5m_bin", "percentGood_18.5m_bin", "percentGood_19.5m_bin", "percentGood_20.5m_bin", "percentGood_21.5m_bin", "Temp_adcp", "ADCP_depth", "adcp_wave_measurment87", "adcp_wave_measurment88", "adcp_wave_measurment89", "adcp_wave_measurment90", "adcp_wave_measurment91", "Temp_01m_HeightAboveBottom", "Temp_02m_HeightAboveBottom", "Temp_03m_HeightAboveBottom", "Temp_04m_HeightAboveBottom", "Temp_05m_HeightAboveBottom", "Temp_06m_HeightAboveBottom", "Temp_07m_HeightAboveBottom", "Temp_08m_HeightAboveBottom", "Temp_09m_HeightAboveBottom", "Temp_10m_HeightAboveBottom", "Temp_11m_HeightAboveBottom", "Temp_12m_HeightAboveBottom", "Temp_13m_HeightAboveBottom", "Temp_14m_HeightAboveBottom", "Temp_15m_HeightAboveBottom", "Temp_16m_HeightAboveBottom", "Temp_17m_HeightAboveBottom", "Temp_18m_HeightAboveBottom", "Temp_19m_HeightAboveBottom", "Temp_20m_HeightAboveBottom", "Pressure_01m_HeightAboveBottom", "Pressure_02m_HeightAboveBottom", "Pressure_03m_HeightAboveBottom", "Pressure_04m_HeightAboveBottom", "Pressure_05m_HeightAboveBottom", "Pressure_06m_HeightAboveBottom", "Pressure_07m_HeightAboveBottom", "Pressure_08m_HeightAboveBottom", "Pressure_09m_HeightAboveBottom", "Pressure_10m_HeightAboveBottom", "Pressure_11m_HeightAboveBottom", "Pressure_12m_HeightAboveBottom", "Pressure_13m_HeightAboveBottom", "Pressure_14m_HeightAboveBottom", "Pressure_15m_HeightAboveBottom", "Pressure_16m_HeightAboveBottom", "Pressure_17m_HeightAboveBottom", "Pressure_18m_HeightAboveBottom", "Pressure_19m_HeightAboveBottom", "Pressure_20m_HeightAboveBottom", "Temperature_BTM_20m", "Temperature_BTM_10m", "Available_column_for_unforseen_measurement", "Pressure_deep", "Temperature_deep", "Conductivity_deep", "Salinity_deep", "sigmatheta_deep", "Pressure_shallow", "Temperature_shallow", "Conductivity_shallow", "Salinity_shallow", "sigmatheta_shallow", "Pressure", "Temperature", "Conductivity", "Salinity", "sigmatheta", "Significant_wave_height", "Dominant_wave_period")
                 
)

# function to convert matlab date to datetime
matlab2POS = function(x, timez = "UTC") {
  days = x - 719529 	# 719529 = days from 1-1-0000 to 1-1-1970
  secs = days * 86400 # 86400 seconds in a day
  # This next string of functions is a complete disaster, but it works.
  # It tries to outsmart R by converting the secs value to a POSIXct value
  # in the UTC time zone, then converts that to a time/date string that 
  # should lose the time zone, and then it performs a second as.POSIXct()
  # conversion on the time/date string to get a POSIXct value in the user's 
  # specified timezone. Time zones are a goddamned nightmare.
  return(as.POSIXct(strftime(as.POSIXct(secs, origin = '1970-1-1', 
                                        tz = 'UTC'), format = '%Y-%m-%d %H:%M', 
                             tz = 'UTC', usetz = FALSE), tz = timez))
}

data<-data %>%
  filter(year %in% c(2021,2022)) %>%
  select(matlab_datenum, year, month, decimal_day, Temperature, Conductivity, Salinity, Significant_wave_height, Dominant_wave_period)%>%
  mutate_all(~na_if(., 9999))%>% # replace all 9999 with NA
  mutate(datetime = matlab2POS(matlab_datenum)) %>%
  drop_na(Significant_wave_height) %>%
  mutate(datetime = round_date(datetime)) %>%
  select(datetime, Significant_wave_height, Dominant_wave_period)

write_csv(data, here("Data","IslandData","WestSideADCPMCR.csv"))
