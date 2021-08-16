function [Q, Q_convect, Q_rad] = Q_ext_top_losses(k, d_o, dT, T_amb, emissivity, stefan_boltzmann)
% FUNCTION NAME:
%   Q_ext_top_losses
%
% DESCRIPTION:
% Compute the heat losses from the top of the hot water cylinder.
%
% The heat losses are assumed to be due to a mixture of convection and radiation
% Conductive losses are assumed not to play a significant contribution to the 
% thermal losses.
%
% INPUT:
%   k - Thermal conductivity of surrounding air (W.m-1.K-1).
%   d_o - Outer diameter of hot water cylinder (m).
%   dT - Temperature rise of the hot water cylinder exterior to ambient (K).
%   T_amb - Ambient temperature of the surroundings (K).
%   emissivity - Radiative emissivity of the hot water cylinder's exterior
%      surface (dimensionless).
%   stefan_boltzmann - Stefan-Boltzmann constant to be used for computing
%      the radiation losses. Used to allow for different power units
%      (W.m-2.K-4).
%
% OUTPUT:
%   Q - Total losses from top of hot water cylinder (W).
%   Q_convect - Convective losses (W).
%   Q_rad - Radiative losses (W).
%
%   The function may be made to give losses in kWh.day-1 by using the following
%   units for k and stefan_boltzmann parameters:
%   k - kwh.day-1.K-1
%   stefan_boltzmann - kWh.day-1.m-2.K-4
%
% ASSUMPTIONS AND LIMITATIONS:
%   This function has not been vectorised. Only one argument may be a vector.
%
%   The hot water cylinder is assumed to be in the vertical position.
%
%   Penerations to the exterior surface of the hot water cylinder are not taken
%   into account.
%
%   This function is an implementation of scenario 2.1.22 in Conduction Heat 
%   Transfer Solutions, James H. VanSant, August 1983.
%
% REVISION HISTORY:
%   20210709 - bjem
%       * Initial implementation


% Assume laminar convection on flat horizontal plate using the small curvature approximation.
    c = 0.27;
    m = 0.25;
    
    % Characteristic dimension is taken as the ratio of area to perimeter.
    % For a circle, this quantity is simply d_o.
    
    %stefan_boltzmann = 5.670374419e-8 * 0.024;
    %emissivity = 0.05;
    Pr = 0.71;
    nu = 1.516e-5;
    A = (pi * d_o^2) / 4;
    x = A / (pi * d_o);
    
    T_rad_eff = ((dT + T_amb + 273.15)^4 - (T_amb + 273.15)^4)^0.25;
    T_film = T_amb + 0.5 * dT;
    Gr_d = (9.81 * dT * x^3) / ((T_film + 273.15) * nu^2);
    
    
    Nu = c * (Gr_d * Pr)^m;
    h = (k / x) * Nu;
    
    Q_convect = h * A * dT;
    Q_rad = stefan_boltzmann * emissivity * A * T_rad_eff^4;
    Q = Q_convect + Q_rad;
    
    1;

end