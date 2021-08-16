function [Q_total, d_i, A_i, h] = Q_ins_end_losses(d_o, Dx1, Dx_end, k, DT)
% FUNCTION NAME:
%   Q_ins_end_losses
%
% DESCRIPTION:
% Compute the internal heat losses for the end of the tank depicted in Figure
% E.4 of AS/NZS 4234:2021.
%
%
% INPUT:
%   d_o - Outer diameter of hot water cylinder (m).
%   h_o - Outer height of hot water cylinder (m).
%   Dx1 - Thickness of the top of radial insulation of the hot water
%         cylinder (m).
%   Dx_end - Thickness of the end insulation (either top or bottom) (m).
%   k - Thermal conductivity of the hot water insulation (W.m-1.K-1).
%   dT - Temperature rise of the hot water cylinder exterior to ambient (K).
%
% OUTPUT:
%   Q_total - Total internal losses of the hot water cylinder through the
%     insulation.through base of hot water cylinder (W).
%   d_i - Internal diameter of the hot water cylinder (m).
%   A_i - Internal circular area of the hot water cylinder (m2).
%   h - Heat transfer coefficient through the end insulation (W.m-2.K-1).
%
%   The function may be made to give losses in kWh.day-1 by using the following
%   units for k:
%   k - kwh.day-1.K-1
%
% ASSUMPTIONS AND LIMITATIONS:
%   This function has not been vectorised. Only one argument may be a vector.
%
%   The insulation material has uniform and homogenous thermal conductivity
%   all over the tank, including at the base (compare note just above Table
%   E.6).
%
% REVISION HISTORY:
%   20210709 - bjem
%       * Initial implementation

d_i = (d_o - 2 * Dx1);

A_i = pi * (d_i / 2)^2;

h = k / Dx_end;

Q_total = h * A_i * DT;

end