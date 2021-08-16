function [Q, Q_convect, Q_rad] = Q_ext_bottom_losses(k, d_o, dT)
% FUNCTION NAME:
%   Q_ext_bottom_losses
%
% DESCRIPTION:
% Compute the heat losses from the base of the hot water cylinder.
%
% The losses model assumes the bottom of the hot water cylinder is
% in contact with a semi-infinite solid with thermal
% conductivity 'k'.
%
% The heat losses are assumed entirely due to thermal conduction. Convective and
% radiative losses are always assumed to be zero.
%
% INPUT:
%   k - Thermal conductivity of the base (W.m-1.K-1).
%   d_o - Outer diameter of hot water cylinder (m).
%   dT - Temperature rise of the hot water cylinder exterior to ambient (K).
%
% OUTPUT:
%   Q - Total losses through base of hot water cylinder (W).
%   Q_convect - Convective losses (W). This value is always zero.
%   Q_rad - Radiative losses (W). This value is always zero.
%
%   The function may be made to give losses in kWh.day-1 by using the following
%   units for k parameter:
%   k - kWh.day-1.K-1
%
% ASSUMPTIONS AND LIMITATIONS:
%   This function has not been vectorised. Only one argument may be a vector.
%
%   This function is an implementation of scenario 2.1.22 in Conduction Heat 
%   Transfer Solutions, James H. VanSant, August 1983.
%
% REVISION HISTORY:
%   20210709 - bjem
%       * Initial implementation
    
    Q = 2 * d_o * k * dT;
    Q_convect = 0;
    Q_rad = 0;
    
    1;

end