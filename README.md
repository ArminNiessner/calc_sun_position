# calc_sun_position
Calculates the position of the sun for a given location and time

## Usage
```
Solar_elev = calc_sun(48.5, 9.3, "2019-08-17 16:00:00", 2, "Solar_elev_corr")
```
Latitude and Longitude in decimal degrees (can be negative e.g. for southern hemisphere). Timezone as integer between -12 and 12.

"Solar_elev_corr": angle of sun above horizon corrected for atmospheric refraction
"Solar_azi_angle": azimuth angle of sun, 0 = norht, 180 = south

## Referene
[NOAA Solar Calculation](https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html)
