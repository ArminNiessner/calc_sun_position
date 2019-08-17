#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 20:57:50 2019

@author: armin


 Copyright (C) 2019  Armin Niessner
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <https://www.gnu.org/licenses/>.

"""

import pandas as pd
import numpy as np

def calc_sun(latitude, longitude, Datetime, timezone, output):          # https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
    jul_date = pd.Timestamp(Datetime).to_julian_date() - (timezone / 24)
    time = Datetime.split(" ")[1]
    (h, m, s) = time.split(':')
    jul_time = (int(h) * 3600 + int(m) * 60 + int(s)) /  (24 * 60 * 60)
    jul_century = (jul_date - 2451545) / 36525
    Geom_mean_long_sun = (280.46646 + jul_century * (36000.76983 + jul_century * 0.0003032)) % 360
    Geom_mean_anom_sun = 357.52911 + jul_century * (35999.05029 - 0.0001537 * jul_century)
    Eccent_earth_orbit = 0.016708634 - jul_century * (0.000042037 + 0.0000001267 * jul_century)
    Sun_eq_of_ctr = np.sin(np.radians(Geom_mean_anom_sun)) * (1.914602 - jul_century * (0.004817 + 0.000014 * jul_century)) + np.sin(np.radians(2 * Geom_mean_anom_sun)) * (0.019993 - 0.000101 * jul_century) + np.sin(np.radians(3 * Geom_mean_anom_sun)) * 0.000289
    Sun_true_long = Geom_mean_long_sun + Sun_eq_of_ctr
    Sun_true_anom = Geom_mean_anom_sun + Sun_eq_of_ctr
    Sun_rad_vector = (1.000001018 * (1 - Eccent_earth_orbit**2)) / (1 + Eccent_earth_orbit * np.cos(np.radians(Sun_true_anom)))
    Sun_app_long = Sun_true_long - 0.00569 - 0.00478 * np.sin(np.radians(125.04 - 1934.136 * jul_century))
    Mean_obliq_eclip = 23 + (26 + ((21.448 - jul_century * (46.815 + jul_century * (0.00059 - jul_century * 0.001813)))) / 60) / 60
    Obliq_corr = Mean_obliq_eclip + 0.00256 * np.cos(np.radians(125.04 - 1934.136 * jul_century))
    Sun_at_ascen = np.degrees(np.arctan2(np.cos(np.radians(Obliq_corr)) * np.sin(np.radians(Sun_app_long)), np.cos(np.radians(Sun_app_long))))
    Sun_declin = np.degrees(np.arcsin(np.sin(np.radians(Obliq_corr)) * np.sin(np.radians(Sun_app_long))))
    var_y = np.tan(np.radians(Obliq_corr / 2)) * np.tan(np.radians(Obliq_corr / 2))
    Eq_of_time = 4 * np.degrees(var_y * np.sin(2 * np.radians(Geom_mean_long_sun)) - 2 * Eccent_earth_orbit * np.sin(np.radians(Geom_mean_anom_sun)) + 4 * Eccent_earth_orbit * var_y * np.sin(np.radians(Geom_mean_anom_sun)) * np.cos(2 * np.radians(Geom_mean_long_sun)) - 0.5 * var_y**2 * np.sin(4 * np.radians(Geom_mean_long_sun)) - 1.25 * Eccent_earth_orbit**2 * np.sin(2 * np.radians(Geom_mean_anom_sun)))
    HA_sunrise = np.degrees(np.arccos(np.cos(np.radians(90.833)) / (np.cos(np.radians(latitude)) * np.cos(np.radians(Sun_declin))) - np.tan(np.radians(latitude)) * np.tan(np.radians(Sun_declin))))
    Solar_noon = (720 - 4 * longitude - Eq_of_time + timezone * 60) / 1440
    Solar_noon_t = pd.to_datetime(Solar_noon)
    Sunrise = Solar_noon - HA_sunrise * 4/1440
    Sunrise_t = pd.to_datetime(Sunrise)
    Sunset = Solar_noon + HA_sunrise * 4/1440
    Sunset_t = pd.to_datetime(Sunset)
    Sunlight_dur = 8 * HA_sunrise
    True_solar_time = (jul_time * 1440 + Eq_of_time + 4 * longitude - 60 * timezone) % 1440
    if True_solar_time / 4 < 0: 
        Hour_angle = True_solar_time / 4 + 180
    else:
        Hour_angle = True_solar_time / 4 - 180
    Solar_zen_angle = np.degrees(np.arccos(np.sin(np.radians(latitude)) * np.sin(np.radians(Sun_declin)) + np.cos(np.radians(latitude)) * np.cos(np.radians(Sun_declin)) * np.cos(np.radians(Hour_angle))))
    Solar_elev_angle = 90 - Solar_zen_angle
    if Solar_elev_angle > 85:
        Approx_atm_refr = 0
    else:
        if Solar_elev_angle > 5:
            Approx_atm_refr = 58.1 / np.tan(np.radians(Solar_elev_angle)) - 0.07 / np.tan(np.radians(Solar_elev_angle))**3 + 0.000086 / np.tan(np.radians(Solar_elev_angle))**5
        else:
            if Solar_elev_angle > -0.575:
                Approx_atm_refr = 1735 + Solar_elev_angle * (-518.2 + Solar_elev_angle * (103.4 + Solar_elev_angle * (-12.79 + Solar_elev_angle * 0.711)))
            else:
                Approx_atm_refr = -20.772 / np.tan(np.radians(Solar_elev_angle))
    Approx_atm_refr = Approx_atm_refr / 3600
    Solar_elev_corr = Solar_elev_angle + Approx_atm_refr
    if Hour_angle > 0:
        Solar_azi_angle = (np.degrees(np.arccos(((np.sin(np.radians(latitude)) * np.cos(np.radians(Solar_zen_angle))) - np.sin(np.radians(Sun_declin))) / (np.cos(np.radians(latitude)) * np.sin(np.radians(Solar_zen_angle))))) + 180) % 360
    else:
        Solar_azi_angle = (540 - np.degrees(np.arccos(((np.sin(np.radians(latitude)) * np.cos(np.radians(Solar_zen_angle))) - np.sin(np.radians(Sun_declin))) / (np.cos(np.radians(latitude)) * np.sin(np.radians(Solar_zen_angle)))))) % 360
    
    dic = {"jul_date": jul_date,
           "jul_time": jul_time,
           "jul_century": jul_century,
           "Geom_mean_long_sun": Geom_mean_long_sun,
           "Geom_mean_anom_sun": Geom_mean_anom_sun,
           "Eccent_earth_orbit": Eccent_earth_orbit,
           "Sun_eq_of_ctr": Sun_eq_of_ctr,
           "Sun_true_long": Sun_true_long,
           "Sun_true_anom": Sun_true_anom,
           "Sun_rad_vector": Sun_rad_vector,
           "Sun_app_long": Sun_app_long,
           "Mean_obliq_eclip": Mean_obliq_eclip,
           "Obliq_corr": Obliq_corr,
           "Sun_at_ascen": Sun_at_ascen,
           "Sun_declin": Sun_declin,
           "var_y": var_y,
           "Eq_of_time": Eq_of_time,
           "HA_sunrise": HA_sunrise,
           "Solar_noon": Solar_noon,
           "Solar_noon_t": Solar_noon_t, #
           "Sunrise": Sunrise,
           "Sunrise_t": Sunrise_t,  #
           "Sunset": Sunset,
           "Sunset_t": Sunset_t, #
           "Sunlight_dur": Sunlight_dur,
           "True_solar_time": True_solar_time,
           "Hour_angle": Hour_angle,
           "Solar_zen_angle": Solar_zen_angle,
           "Solar_elev_angle": Solar_elev_angle,
           "Approx_atm_refr": Approx_atm_refr,
           "Solar_elev_corr": Solar_elev_corr,
           "Solar_azi_angle": Solar_azi_angle}
    
    return dic[output]

# example:
# Angle of the sun over the Horizon
# Solar_elev = calc_sun(48.5, 9.3, "2019-08-17 16:00:00", 2, "Solar_elev_corr")