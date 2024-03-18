library(raster)
library(sebkc)
library(rgdal)

setwd("E:\\B\\sebal_2019")

NDVI=raster("ndvi2019-12.tif") # Ndvi
Ts=raster("lst2019-12.tif") # Land surface temperature
dem=raster("E:\\B\\sebal\\dem.tif") # Digital elevation model
albedo=raster("alb2019-12.tif") # Albedo
SAVI=raster("savi2019-12.tif") # Soil-Adjusted Vegetation Index 
tmin=raster("tmin2019-12.tif") # Minimum atmospheric temperature 
tmax=raster("tmax2019-12.tif") # Maximum atmospheric temperature 
ETo=raster("eto2019-12.tif") # FAO Penman-monteith reference evapotranspiration 
u=raster("u_component_of_wind_10m.tif") # v component of wind
v=raster("v_component_of_wind_10m.tif") # u component of wind
uw=raster("u22019-12.tif") # wind speed 
latitude = raster("rad2019-1.tif") # raster image containing latitude values
shapefile = readOGR("E:\\B\\spi\\study_area_updated_utm.shp") 
date="2019-01-01" # Date 
Gsc=0.0820##solar constant MJ/m2
zx=10 # hieght at which wind speed measured

###Regression constant, intercept, expressing the fraction of extra-terrestrial radiation reaching the earth on overcast days
as=0.25
bs=0.50
model = "SEBAL" # SEBAL\METRIC
iter.max =2   # maximum iterations of sensible heat calculation.
t1 = 8 # The length of the calculation period in hour; 1 for hour, 0.5 for 30 minutes 0.25 for 15 minutes
i = 1  # Loop counter
rah_hot_change = -10000
coldDT = 0

# Date of year (DOY)
j=yday(date)


ETr=ETo*11.574/24
tmean = (tmax + tmin)/2
etmax = 0.6108 * exp(17.27 * tmax/(tmax + 237.3))
etmin = 0.6108 * exp(17.27 * tmin/(tmin +237.3))
es = (etmax + etmin)/2
ea = 0.6108 * exp(17.27 * tmin/(tmin +237.3))
p = 101.3 * ((293 - 0.0065 * dem)/293)^5.26
sd =  0.409 * sin(2 * pi/365 * J - 1.39)
dr = 1 + 0.033 * cos(2 * pi/365 * J)
psi = p*0.000665
ws = acos(-tan(lat_rad) * tan(sd))
Ra = (1440/pi) * dr * Gsc *((ws * sin(lat_rad) * sin(sd))+ (cos(lat_rad) * cos(sd) * sin(ws)))
Rso = (0.75 + (2 * 10^-5) * dem) * Ra
N = 24/pi * ws
n = (tmean*0.232)+4.352
Rs = (as + bs * (n/N)) * Ra
Rns = (1 - albedo) * Rs
Rnl = 0.000000004903 * (0.34 - 0.14 * sqrt(ea))*((tmax + 273.2)^4 +(tmin + 273.2)^4)/2 * (1.35 * Rs/Rso - 0.35)
Rn = (Rns - Rnl)/0.0864
LAI = (5.434 * SAVI) + 1.5416
zom=0.018*LAI#momentum roughnss length


extent=extent(NDVI)# A extent of the raster file that include c(xmin, xmax; serastercond row: ymin, ymax)

modhot=hotTs(Ts, NDVI, albedo, dem, cluster = 8,
             extent = extent, upper = 0.8, lower = 0.2, plot = TRUE,
             layout = "landscape", draw = "rect", folder = NULL, welev = NULL,
             clip = NULL)


modcold=coldTs(Ts, NDVI, albedo, sunangle = 1.418128, dem,
               cluster = 8, extent = extent, upper = 0.95, lower = 0.2,
               plot = TRUE, layout = "portrait", draw = "rect", folder = NULL,
               welev = NULL, clip = NULL, iter.max = 100)

# x and y coordinates of a hot pixel in the form c(x,y) -------------------


xhot=modhot$x
yhot=modhot$y

# x and y coordinates of a cold pixel in the form c(x,y) -------------------

xcold=modcold$x
ycold=modcold$y

######soil heat flux by G/Rn ratio

G_by_Rn=(Ts)*(0.0038+0.007*albedo)*(1-0.98*NDVI^4)
G_by_Rn[NDVI<0,]=0.5
G_by_Rn[Ts<4&albedo>0.45,]=0.5
G=G_by_Rn*Rn
G=(((1.80*(Ts-273.15))/Rn)+0.084)*84

##############aerodynamic resistance to heat transport#########

k=0.41 #von Karman's constant
## z1 and z2 are heights in meters above the zero plane displacement (d) of the vegetation
z1=0.1
z2=2

u200=((uw*log(200))/zom)/log(zx/zom)
u_star=(k*u200)/(log(200/zom))# friction velocity
rah=log(z2/z1)/(u_star*k)# aerodynamic resistance to heat tranpsort


df = data.frame(a= numeric(), b = numeric(), u_star = numeric(),
                L = numeric(),pairhot=numeric(),paircold=numeric(),
                rah_hot=numeric(),
                rah_cold=numeric(),dThot=numeric(),dTcold=numeric(),
                Hhot=numeric(),Hcold=numeric(),Ts_hot=numeric(),
                Ts_cold=numeric(),
                Rn_hot=numeric(),Rn_cold=numeric(),G_hot=numeric(),
                G_cold=numeric(),
                LEcold=numeric(),LEhot=numeric(),ETint_hot=numeric(),
                ETint_cold=numeric())
Ts_cold=getValues(Ts)[cellFromXY(Ts,c(xcold,ycold))]
Ts_hot=getValues(Ts)[cellFromXY(Ts,c(xhot,yhot))]
Rn_cold=getValues(Rn)[cellFromXY(Rn,c(xcold,ycold))]
Rn_hot=getValues(Rn)[cellFromXY(Rn,c(xhot,yhot))]
G_cold=getValues(G)[cellFromXY(G,c(xcold,ycold))]
G_hot=getValues(G)[cellFromXY(G,c(xhot,yhot))]
print("Computing sensible heat flux........................")


p = 1.03547033 # air pressure in kilopascals (kPa)
cp = 1004.14# air specific heat
pressure = 101.325 * ((293 - (0.0065 * dem))/293)^5.26
g = 9.81 # gravity 
vaporisation = (2.501 - (0.00236 * (Ts - 273.15))) * 10^6

###sensible heat flux#######

while (i <= iter.max) {
  print(paste("Monin-Obukhov length iteration", i, "of",
              iter.max))
  rah_hot = getValues(rah)[cellFromXY(rah, c(xhot, yhot))]
  rah_cold = getValues(rah)[cellFromXY(rah, c(xcold, ycold))]
  if (i > 1) {
    rah_hot_change = round((((rah_hot - rah2)/rah2) * 100), 2)
  }
  hotDT = ((Rn_hot - G_hot) * rah_hot)/(p * cp)
  a = (hotDT - coldDT)/(Ts_hot - Ts_cold)
  b = (hotDT - a)/Ts_hot
  b =(hotDT - (a * Ts_hot))
  dT = (a * Ts) + b
  pair = (1000 * pressure)/(1.01 * ((Ts) * 287))
  p = getValues(pair)[cellFromXY(pair, c(xhot, yhot))]
  H = (pair * cp * dT)/rah
  H_hot = getValues(H)[cellFromXY(H, c(xhot, yhot))]
  H_cold = getValues(H)[cellFromXY(H, c(xcold, ycold))]
  if (model == "METRIC") {
    if (is.null(ETr)) {
      return(print("ETr is needed"))
    }
    LE_cold = 1.05 *getValues(ETr)[cellFromXY(pair,
              c(xcold, ycold))]  * getValues(pair)[cellFromXY(pair,c(xcold, ycold))]
    H_cold = (Rn_cold - LE_cold - G_cold)/3600
    coldDT = (H_cold * rah_cold)/(getValues(pair)[cellFromXY(pair,c(xcold, ycold))] * cp)
  }
  else {
    coldDT = (H_cold * rah_cold)/(getValues(pair)[cellFromXY(pair,c(xcold, ycold))] * cp)
  }
  u_star_hot = getValues(u_star)[cellFromXY(u_star, c(xhot,yhot))]
  L = -((p * cp * (u_star_hot^3) * Ts_hot)/(k * g * H_hot))
  x_200m = (1 - (16 * (200/L)))^0.25
  x_2m = (1 - (16 * (2/L)))^0.25
  x_01m = (1 - (16 * (0.1/L)))^0.25
  w_200m = ifelse(L < 0, 2 * log((1 + x_200m)/2) + log((1 +x_200m^2)/2) - 2 * atan(x_200m) + 0.5 * pi, (ifelse(L > 0, -5 * (200/L), 0)))
  w_2m = ifelse(L < 0, 2 * log((1 + x_2m^2)/2), ifelse(L > 0, -5 * (2/L), 0))
  w_01m = ifelse(L < 0, 2 * log((1 + x_01m^2)/2), ifelse(L > 0, -5 * (2/L), 0))
  u_star = (u200 * k)/(log(200/zom) - w_200m)
  rah = (log(z2/z1) - w_2m + w_01m)/(u_star * k)
  rah2 = getValues(rah)[cellFromXY(rah, c(xhot, yhot))]
  
  df = rbind(df, data.frame(a = a, b = b, u_star = u_star_hot, L = L, paircold = getValues(pair)[cellFromXY(pair,
                            c(xcold, ycold))], pairhot = p, rah_hot = rah_hot, rah_cold = rah_cold, dThot = hotDT, dTcold = coldDT,
                            Hhot = H_hot, Hcold = H_cold, Ts_hot = Ts_hot, Ts_cold = Ts_cold, Rn_hot = Rn_hot, Rn_cold = Rn_cold,
                            G_hot = G_hot, G_cold = G_cold, LEcold = Rn_cold - G_cold - H_cold, LEhot = Rn_hot - G_hot - H_hot,
                            ETint_hot = t1 * 3600 * ((Rn_hot - G_hot - H_hot)/(getValues(vaporisation)[cellFromXY(vaporisation,
                            c(xhot, yhot))])), ETint_cold = t1 * 3600 * ((Rn_cold - G_cold - H_cold)/(getValues(vaporisation)
                                                                                                                                                                                            [cellFromXY(vaporisation, c(xcold, ycold))]))))
  
  i = i + 1
}

LE = Rn - G - H
ETint = (t1 * 3600 * LE)/vaporisation


if (model == "METRIC") {
  EF = ETint / ETr
  EF[EF > 1.1, ] = 1.1
  EF[EF < 0, ] = 0
  if (!is.null(ETr)) {
    ET24 = EF * ETr
  }
} else {
  EF = 1.1 * (LE / (Rn - G))
  EF[EF < 0, ] = 0
  EF[EF > 1.1, ] = 1.1
  if (!is.null(Rn)) {
    ET24 = EF * Rn * 0.035
  }
}



