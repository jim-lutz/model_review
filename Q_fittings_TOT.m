function [Q, Q_fittings, Q_TPR, Q_thermopocket] = Q_fittings_TOT(n_fittings, n_TPR, n_thermopocket, dT, k_per_fitting, k_per_TPR, h_thermopocket)
% FUNCTION NAME:
%   Q_fittings_TOT
%
% DESCRIPTION:
% Compute the heat losses from the fittings attached to the hot water cylinder
% per Table E.9 of AS/NZS 4234:2021.
%
% The losses in Table E.9 can be broadly summarised into three categories:
% * Losses that depend on the surface area of a tank opening
%   (categorised as thermostat pockets).
% * Losses that depend on the diameter of the pipe or fitting.
%
% This function does not distinguish between insulated and non-insulated
% fittings. This is done by selecting the values of k_per_fitting and k_per_TPR
% as appropriate. Table E.9 allows losses from pipes or fittings to be deducted
% by 30% if they are insulated.
%
% The temperature rise may be to ambient, or to the exterior surface of the
% hot water cylinder.
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
%   Table E.9 is for the following conditions:
%   * Temperature rise (dT) of 55 K.
%   * The effective number of a fitting or TPR valve of a different effective
%     diameter than the 'base case' used to establish k_per_fitting or k_per_TPR
%     is scaled by the ratio of the base case diameter to the actual diameter.
%   * The effective number of a thermopocket a different effective area than the
%    'base case' used to establish h_per_thermopocket is scaled by the ratio of
%     the base case area to the actual area.
%
% It is possible to interchange the fittings to different categories e.g. a 
% capped pipe outlet or inlet may be modelled as an area-dependent loss similar
% to a thermopocket.
%
% OUTPUT:
%   Q - Total losses through base of hot water cylinder (W).
%   Q_fittings - Thermal losses through pipe fittings (W).
%   Q_TPR - Thermal losses through TPR valves (W).
%   Q_thermopocket - Thermal losses through thermostat pockets (W).
%
% The total losses is the sum of the three loss components.
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
    
    T_fittings = k_per_fitting * n_fittings;
    T_TPR = k_per_TPR * n_TPR;
    T_thermopocket = h_thermopocket * n_thermopocket;
    
    Q_fittings = T_fittings * dT;
    Q_TPR = T_TPR * dT;
    Q_thermopocket = T_thermopocket * dT;
    
    Q = Q_fittings + Q_TPR + Q_thermopocket;
    
end