# ExtTempFunc.R
# functions called by ExtTempCalc.Rmd

# Functions of Properties of Air
# linear fit to cover ~ 20 °C to 75 °C
# from https://www.engineeringtoolbox.com/
# see ../Peer Review/Properties_of_Air.ods

beta <- function(
  # coefficient of thermal expansion of air, 1/K 
  # used in Grashof number
  Tair = 20 # air temperature, °C 
  # default beta = 3.43e-3  
){
  beta <- (-0.01048275862069 * Tair + 3.63793103448276) * 1e-3
  return(beta)
}

nu <- function(
  # kinematic viscosity of air, (m2/s), using air at 22.5 °C
  # used in Grashof number and Prandtl number
  Tair = 22.5 # air temperature, °C 
  # default nu = 15.28e-6 
){
  nu <- (0.094*Tair + 13.17) * 1e-6
  return(nu)
  }

alpha <- function(
  # Air thermal diffusivity at atmospheric pressure, (m2/s), using air at 22.5 °C
  # used in Prandtl number
  Tair = 22.5 # air temperature, °C 
  # default alpha =  22.045e-6 
){
  alpha <- (0.143853658536585*Tair + 18.7818536585366) * 1e-6
  return(alpha)
}


k_air <- function(
  # Air thermal conductivity at atmospheric pressure, (W/m-K)
  # used with Nusselt number to calculate convective heat transfer coefficients
  Tair = 22.5 # air temperature, °C 
  # default k_air = 26.06e-3
  ){
  k_air <- 0.071967176004528 * Tair + 24.4579117147708
  return(k_air)
  }


# Functions for conductive heat transfer
# from 2021 ASHRAE® HANDBOOK FUNDAMENTALS, SI Edition
# CHAPTER 4 HEAT TRANSFER
# Table 2 One-Dimensional Conduction Shape Factors

Qcond.cyl <- function( 
  # function to calculate the heat transfer through 
  # hollow cylinder of length L with negligible heat transfer from end surfaces
  length,      # length of cylinder, m
  r_outer,     # outer radius of cylinder, m
  r_inner,     # inner radius of cylinder, m
  T_inner = 75, # temperature of inner surface, °C
  T_outer,      # temperature of outer surface, °C
  k = 0.025     # thermal conductivity of cylinder,  W/m-K, default to polyurethane
)
{
  Qcond.cyl <- ( 2 * pi * k * length *(T_inner - T_outer)) / 
    ( log(r_outer/r_inner) )
  
  return(Qcond.cyl)
}

Qcond.disc <- function( 
  # function to calculate conductive heat transfer through a disc
  radius,       # radius of disc, m
  thick,        # thickness of cylinder, m
  Thot = 75,    # hot surface temperature, °C
  Tcold,        # cold surface temperature, °C
  k = k_ins     # thermal conductivity of disc, W/m-K
)
{
  Qcond.disc <- ( (k_ins * pi * radius^2) / (thick)) * (Thot - Tcold)
  
  return(Qcond.disc)
}
  

# Functions for natural convective heat transfer
# from 2021 ASHRAE® HANDBOOK FUNDAMENTALS, SI Edition
# CHAPTER 4 HEAT TRANSFER
# Table 9 Natural Convection Correlations

Gr <- function(
  # calculate the Grashof number
  g     = 9.81,     # acceleration due to Earth's gravity (m/s^2)
  beta  = 3.43e-3,  # the coefficient of thermal expansion of air at 20°C (1/K)
  Tsurf = 25,       # the surface temperature (°C )
  Tamb  = 20,       # the bulk temperature (°C)
  L ,               # characteristic length (m)
  mu    = 15.28e-6  # the kinematic viscosity (m2/s), using air at 22.5°C
) 
{
  Gr <- ( g * beta * (Tsurf - Tamb) * L ^ 3 )/ mu^2
  return(Gr)
}

Pr <- function(
  # calculate the Prandtl number
  mu    = 15.28e-6,  # the kinematic viscosity (m2/s), using air at 22.5°C
  alpha = 22.045e-6 # the thermal diffusivity (m2/s) for air at 22.5°C)
) 
{
  Pr <- mu / alpha
  return(Pr)
}

Ra <- function(
  # calculate the Rayleigh number
  Gr,  # Grashof number
  Pr   # Prandtl number
) 
{
  Ra <- Gr * Pr
  return(Ra)
}

Nu.vert.nat <- function(
  # Natural Convection Correlation for Vertical plate, 
  # Tsurf = constant
  # Characteristic dimension: L = height
  # Properties at (ts + t∞)/2 except β at t∞  
  Ra,  # Rayleigh number
  Pr   # Prandtl number
) 
{
  if (Ra < 10^9) {
    Nu.vert.nat <- 0.68 + 
      (0.67 * Ra^(1/4)) /
      ( 1 + (0.492/Pr)^(9/16)  )^(4/9)
    
  } else {
    if (Ra < 10^12) {
      Nu.vert.nat <- (0.825 + 
                   (0.387 * Ra^(1/6)) /
                   ( 1 + (0.437/Pr)^(9/16)  )^(8/27)
      )^2
    } else {
      Nu.vert.nat <- NA
    }
  }
  
  return(Nu.vert.nat)
  
}



# Functions for forced convective heat transfer
# from 2021 ASHRAE® HANDBOOK FUNDAMENTALS, SI Edition
# CHAPTER 4 HEAT TRANSFER
# Table 8 Forced-Convection Correlations

Re.cyl.frc <- function( 
  # Reynolds number for cross flow over cylinder
  V, # air velocity, m/s
  D, # diameter, m    this is the characteristic length
  nu   # kinematic viscosity, (m2/s)
){
  Re.cyl.frc <- V * D / nu
  
  return(Re.cyl.frc)
}



Nu.cyl.frc <- function(
  # Nusselt number for forced flow over cylinder
  # IV. External Flows for Cross Flow over Cylinder: 
  # Characteristic length = D = diameter. Re = VD/ν.
  # All properties at arithmetic mean of surface and fluid temperatures.
  # Average value of h
  Re,  # Reynolds number
  Pr   # Prandtl number
){
  Nu.cyl.frc <- 0.3 + (
    ( 0.62 * Re^(1/2) * Pr^(1/3) ) / 
    ( 1 + (0.4/Pr)^(2/3) ) ^(1/4)
                      ) *
    ( 1 + (Re/282000)^(5/8) ) ^(4/5)

  return(Nu.cyl.frc)

} 
  

