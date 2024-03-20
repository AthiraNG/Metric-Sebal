library(raster)
library(sebkc)
library(lubridate)

setwd("D:\\sebal") # Set Directory


# Raster Layers (Inputs) --------------------------------------------------

NDVI = raster("ndvi2019-1.tif") # Ndvi
Ts = raster("lst2019-1.tif") # Land surface temperature
Z = raster("dem.tif") # Digital elevation model
albedo = raster("alb2019-1.tif") # Albedo
LAI = raster("LAI2019-1.tif") # Leaf area index
Tmin = raster("tmin2019-1.tif") # Minimum atmospheric temperature 
Tmax = raster("tmax2019-1.tif") # Maximum atmospheric temperature 
u = raster("u_component_of_wind_10m.tif") # v component of wind
v = raster("v_component_of_wind_10m.tif") # u component of wind
uw = raster("uh2019-1.tif") # wind speed 
latitude = raster("rad2019-1.tif") # raster image containing latitude values
Tdew = raster("tdew2019-1.tif") # Dewpoint temperature
RHmax = raster("RH_max.tif") # Maximum relative humidity
RHmin = raster("RH_min.tif") # Minimun relative humidity
ETo = raster("ETo2019-1.tif") # Reference evapotranspiration

# Constants Used ----------------------------------------------------------

gsc = 0.0820 # solar constant MJ/m2
h = 10 # hieght at which wind speed measured
k=0.41 #von Karman's constant
# z1 and z2 are heights in meters above the zero plane displacement (d) of the vegetation
z1=0.1
z2=2
#Regression constant, intercept, expressing the fraction of extra-terrestrial radiation reaching the earth on overcast days
as=0.25
bs=0.50
p = 1.03547033 # air pressure in kilopascals (kPa)
cp = 1004.14# air specific heat
g = 9.81 # gravity 


# Specify the model (SEBAL/METRIC) ----------------------------------------

model = "METRIC" 

# Other values used in the model ------------------------------------------

iter.max =4   # maximum iterations of sensible heat calculation.
t1 = 8 # The length of the calculation period in hour; 1 for hour, 0.5 for 30 minutes 0.25 for 15 minutes
i = 1  # Loop counter
rah_hot_change = -10000
coldDT = 0

# Date Of Year ------------------------------------------------------------

date = "2019-01-01" # Date 
j=yday(date) # Date of year (DOY)

# Functions to calculate all the variables --------------------------------

# function to convert DMS to degree decimal

dms_to_decimal = function(degrees, minutes, seconds) {
  decimal = degrees + (minutes / 60) + (seconds / 3600)
  return(decimal)
}

# function to calculate saturation vapour pressure curve

calc_satvpr_slope <- function(Tmean) {
  d=(4098*(0.6108*(2.718282*(((17.27*Tmean)/(Tmean+237.3))))))/((Tmean+237.2)^2)
  return(d)
}

# function to calculate saturation vapor pressure at the air temperature

calc_eTmax <- function(Tmax) {
  eTmax=0.6108*(2.718282*((17.27*Tmax)/(Tmax+237.3)))
  return(eTmax)
}
calc_eTmin <- function(Tmin) {
  eTmin=0.6108*(2.718282*((17.27*Tmin)/(Tmin+237.3)))
  return(eTmin)
}

# function to calculate actual vapor pressure

calc_ea <- function(Tdew = NULL, RHmax = NULL, RHmin = NULL, eTmin = NULL, eTmax = NULL) {
  if (!is.null(Tdew)) {
    # Calculate ea using Tdew
    ea = 0.6108 * (2.718282 * ((17.27 * Tdew) / (Tdew + 237.3)))
  } else if (!is.null(RHmax) && !is.null(RHmin) && !is.null(eTmin) && !is.null(eTmax)) {
    # Calculate ea using RHmax, RHmin, eTmin, and eTmax
    ea = ((eTmin * (RHmax / 100)) + (eTmax * (RHmin / 100))) / 2
  }
  
  return(ea)
}

# function to calculate inverse relative distance Earth-Sun (dr)

calc_dr <- function(j) {
  dr=1+(0.033*(cos(((2*pi)/365)*j)))
  return(dr)
}

# function to calculate solar declination 

calc_sd <- function(j) {
  sd=0.409*sin((((2*pi)/365)*j)-1.39)
  return(sd)
}

# function to calculate solar hour angle
calc_ωs <- function(radians,sd) {
  ωs=acos((-tan(radians))*tan(sd))
  return(ωs)
}

# function to calculate atmospheric pressure 

calc_atm_pressure <- function(z) {
  atm_pressure= (((293-(0.0065*z))/293)^5.26)*101.3
  return(atm_pressure)
}

# Function to calculate wind speed at 2 meters height (u2)

calc_u2 =  function(ωs, h) {
  u2 = ωs * (4.87 / log((67.8 * h) - 5.42))
  return(u2)
}

# Function to calculate extraterrestrial radiation (Ra)

calc_Ra = function(gsc, dr, ωs, radians, sd) {
  Ra = (1440/pi)*(gsc)*(dr)*((ωs*(sin(radians))*
                                (sin(sd)))+((cos(radians))*(cos(sd))*(sin(ωs))))
  return(Ra)
}

# Function to calculate Mean daily solar radiation (Rs)

calc_Rs <- function(Tmax,Tmin,Ra) {
  Rs=(Tmax-Tmin)^0.5*Ra*0.16
  return(Rs)
}
# Function to calculate Net solar or net shortwave radiation

calc_Rns <- function(alb,Rs) {
  Rns=(1-alb)*Rs
  return(Rns)
}

# Function to calculate extraterrestrial radiation (Rnl)

calc_Rnl = function(Tmax,Tmin,ωs,ea,Rs,Rso) {
  Rnl=(((Tmax+273.16)^4+(Tmin+273.16)^4)/2)*4.903*10^(-9)*
    (0.34-(0.14*sqrt(ea)))*((1.35*(Rs/Rso))-0.35)
  return(Rnl)
}

# Function to calculate FAO-56 Reference evapotranspiration (ETo)

calc_ETo = function(Rn,d,r,u2,Tmean,Rso,es,ea) {
  ETo=((0.409*Rn*d)/(d+(r*(1+(0.34*u2)))))+
    ((900*(u2/(Tmean+273)))*((es-ea)*r)*(d+
                                           (r*(1+(0.34*u2)))))
  return(ETo)
}

# Radians 

if (file.exists("rad2019-1.tif")) {
  latitude <- raster("rad2019-1.tif")
  
  latitude_radians <- (pi / 180) * latitude
} else {
  # Assuming you have defined 'degrees', 'minutes', and 'seconds'
  decimal_degrees <- dms_to_decimal(degrees, minutes, seconds)
  
  latitude_radians <- (pi / 180) * decimal_degrees
}

# Mean daily temperature

stack_tminmax=stack(Tmax,Tmin)# Stack of min and max temperature
Tmean=calc(stack_tminmax,mean)

# Slope of saturation vapor pressure curve

d=calc_satvpr_slope(Tmean)

#Atmospheric Pressure

atm_pressure= calc_atm_pressure(Z)

# Psychrometric constant

r=0.000665*atm_pressure

# Saturation vapor pressure at the air temperature

eTmax=calc_eTmax(Tmax)
eTmin=calc_eTmin(Tmin)

# The mean saturation vapor pressure

eTmax_Tmin_stcak=stack(eTmax,eTmin)
es=calc(eTmax_Tmin_stcak,mean)

# Actual vapor pressure 

if (file.exists("tdew2019-1.tif")) {
  Tdew = raster("tdew2019-1.tif")
  ea = calc_ea(Tdew)
} else if (file.exists("RH_max.tif") && file.exists("RH_min.tif") && file.exists("tmin2019-1.tif") 
           && file.exists("tmax2019-1.tif")) {
  RHmax = raster("RH_max.tif")
  RHmin = raster("RH_min.tif")
  Tmax = raster("tmax2019-1.tif")
  Tmin = raster("tmin2019-1.tif")
  # Calculate actual vapor pressure from RHmax and RHmin, if available
  ea <- calc_u2(RHmax,Rmin,Tmax,Tmin)
}

# The inverse relative distance Earth-Sun (dr)

dr=calc_dr(j)

# Solar declination

sd=calc_sd(j)

# Sunset hour angle

ωs=calc_ωs(latitude_radians,sd)

# Extraterrestrial radiation

Ra=calc_Ra(gsc, dr, ωs, latitude_radians, sd)

# Mean daily solar radiation

Rs=calc_Rs(Tmax,Tmin,Ra)

# Clear sky solar radiation

Rso=0.759*Ra

# Net solar or net shortwave radiation

Rns=calc_Rns(albedo,Rs)

# Incoming net long wave radiation

Rnl=calc_Rnl(Tmax,Tmin,ωs,ea,Rs,Rso)

# Net radiation

Rn=Rns-Rnl

# Average wind speed at 2 meter

if (file.exists("uh2019-1.tif")) {
  WS <- raster("uh2019-1.tif")
  u2 <- calc_u2(WS, h)
} else if (file.exists(u) && file.exists(v)) {
  u <- raster(u)
  v <- raster(v)
  # Calculate wind speed from u and v components, if available
  WS <- sqrt(u^2 + v^2)
  u2 <- calc_u2(WS, h)
}

ETr = ETo*11.574/24 # alfalfa reference crop

zom=0.018*LAI # momentum roughnss length

vaporisation = (2.501 - (0.00236 * (Ts - 273.15))) * 10^6




extent=extent(NDVI) # A extent of any raster file that include c(xmin, xmax; serastercond row: ymin, ymax)

modhot=hotTs(Ts, albedo, Z, cluster = 8,
             extent = extent, upper = 0.8, lower = 0.2, plot = TRUE,
             layout = "landscape", draw = "rect", folder = NULL, welev = NULL,
             clip = NULL)


modcold=coldTs(Ts,albedo, sunangle = 1.418128, Z,
               cluster = 8, extent = extent, upper = 0.95, lower = 0.2,
               plot = TRUE, layout = "portrait", draw = "rect", folder = NULL,
               welev = NULL, clip = NULL, iter.max = 100)

# x and y coordinates of a hot pixel in the form c(x,y) -------------------


xhot=modhot$x
yhot=modhot$y

# x and y coordinates of a cold pixel in the form c(x,y) -------------------

xcold=modcold$x
ycold=modcold$y


# Soil heat flux by G/Rn ratio --------------------------------------------

G_by_Rn=(Ts)*(0.0038+0.007*albedo)*(1-0.98*NDVI^4)
G_by_Rn[NDVI<0,]=0.5
G_by_Rn[Ts<4&albedo>0.45,]=0.5
G=G_by_Rn*Rn
G=(((1.80*(Ts-273.15))/Rn)+0.084)*84


# Aerodynamic resistance to heat transport --------------------------------

u200=((uw*log(200))/zom)/log(h/zom) # Wind speed at a blending height assumed to be 200 m

u_star=(k*u200)/(log(200/zom)) # friction velocity

rah=log(z2/z1)/(u_star*k) # aerodynamic resistance to heat tranpsort


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



# Sensible Heat Flux ------------------------------------------------------


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
  pair = (1000 * atm_pressure)/(1.01 * ((Ts) * 287))
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

LE = Rn - G - H # Latent heat flux
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
ET24


