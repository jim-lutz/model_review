# ExtTempFunc.R
# functions called by ExtTempCalc.Rmd

# Properties of Air
# linear fit to cover ~ 20 °C to 75 °C
# from https://www.engineeringtoolbox.com/
# see ../Peer Review/Properties_of_Air.ods

beta <- function(
  # coefficient of thermal expansion of air, 1/K  
  Tair = 20 # air temperature, °C 
  # default beta = 3.43e-3  
){
  beta = (-0.01048275862069 * Tair + 3.63793103448276) * 1e-3
  return(beta)
}

nu <- function(
  # kinematic viscosity of air, (m2/s), using air at 22.5 °C
  Tair = 22.5 # air temperature, °C 
  # default nu = 15.28e-6 
){
  nu = (0.094*Tair + 13.17) * 1e-6
  return(nu)
  }

alpha <- function(
  # Air thermal diffusivity at atmospheric pressure, (m2/s), using air at 22.5 °C
  Tair = 22.5 # air temperature, °C 
  # default alpha =  22.045e-6 
){
  alpha = (0.143853658536585*Tair + 18.7818536585366) * 1e-6
  return(alpha)
}


k_air <- function(
    # Air thermal conductivity at atmospheric pressure, (W/m-K)
    Tair = 22.5 # air temperature, °C 
    # default k_air = 26.06e-3
  ){
  k_air = 0.071967176004528 * Tair + 24.4579117147708
  return(k_air)
  }



