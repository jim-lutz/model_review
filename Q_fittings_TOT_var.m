function [Q_var, Q_var_fittings, Q_var_TPR, Q_var_thermopocket] = Q_fittings_TOT_var(n_fittings, n_TPR, n_thermopocket, dT, k_per_fitting, k_per_TPR, h_thermopocket)
% FUNCTION NAME:
%   Q_fittings_TOT_var
%
% DESCRIPTION:
% Compute the sub-variance of the heat losses from the fittings attached to the
% hot water cylinder per Table E.9 of AS/NZS 4234:2021.
%
% This function is used for the situation where the variance is assumed to be
% related to the square of the square of the total losses.
%
% This function assumes the ratio is one. The true variance must have this ratio
% supplied. One basis is the penalty factor and an assumption of what percentile
% the penalty factor represents. The default penalty factor in AS/NZS 4234:2021
% is 1.35. If the range of 1 - 1.35 - 1.7 times the insulation calculated via
% Q_fittings_TOT is assumed to represent 5 - 50 - 95 percentiles, then the ratio
% is (1.7 - 1.35) / normcdfinv(0.95) ~= 0.2103. The output of this function must
% be multiplied by the square of this ratio to get the true variance.
%
% The variance in this function is not assumed to be the simple square of the
% losses because the losses between the different fittings is uncorrelated.
%
% INPUT:
%   n_fittings - The effective number of pipe fittings (dimensionless).
%   n_TPR - The effective number of TPR valves (dimensionless).
%   n_thermopocket - The effective number of thermostat pockets (dimensionless).
%   dT - Temperature rise of the hot water cylinder exterior to ambient (K).
%   k_per_fitting - Thermal conductivity of an individual fitting or pipe
%      (W.K-1).
%   k_per_TPR - Thermal conductivity of an individual TPR valve (W.K-1).
%   h_per_thermopocket - Heat transfer coefficient of a single thermostat
%      pocket (W.K-1).
%
% OUTPUT:
%   Q - Total losses through base of hot water cylinder unscaled variance (W2).
%   Q_fittings - Thermal losses through pipe fittings unscaled variance (W2).
%   Q_TPR - Thermal losses through TPR valves unscaled variance (W2).
%   Q_thermopocket - Thermal losses through thermostat pockets unscaled variance
%      (W2).
%
% The total losses variance is the sum of the three loss components variances.
%
% The function may be made to give losses in kWh.day-1 by using the following
% units for k parameter:
% k_per_fitting - kWh.day-1.K-1
% k_per_TPR - kWh.day-1.K-1
% h_per_thermopocket - kWh.day-1.K-1
%
% ASSUMPTIONS AND LIMITATIONS:
%   This function has not been vectorised. Only one argument may be a vector.
%
% REVISION HISTORY:
%   20210709 - bjem
%       * Initial implementation
    
    Q_var_fittings = k_per_fitting^2 * n_fittings * dT^2;
    Q_var_TPR = k_per_TPR^2 * n_TPR * dT^2;
    Q_var_thermopocket = n_thermopocket * h_thermopocket^2 * dT^2;
    
    Q_var = Q_var_fittings + Q_var_TPR + Q_var_thermopocket;
    
end