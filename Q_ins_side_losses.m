function [Q_total, d_i, h_i, A_rad, h] = Q_ins_side_losses(d_o, h_o, Dx1, Dx2, Dx3, k, dT)
% FUNCTION NAME:
%   Q_ins_side_losses
%
% DESCRIPTION:
% Compute the internal radial heat losses for the side the tank depicted in
% Figure E.4 of AS/NZS 4234:2021.
%
%
% INPUT:
%   d_o - Outer diameter of hot water cylinder (m).
%   h_o - Outer height of hot water cylinder (m).
%   Dx1 - Thickness of the top of radial insulation of the hot water
%         cylinder (m).
%   Dx2 - Thickness of the top insulation (m).
%   Dx3 - Thickness of the bottom insulation (m).
%   k - Thermal conductivity of the hot water insulation (W.m-1.K-1).
%   dT - Temperature rise of the hot water cylinder exterior to ambient (K).
%
% OUTPUT:
%   Q_total - Total internal losses of the hot water cylinder through the
%     insulation.through base of hot water cylinder (W).
%   d_i - Internal diameter of the hot water cylinder (m).
%   h_i - Internal height of the hot water cylinder (m).
%   A_rad - Effective internal area used for radial heat transfer (m2).
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

h_i = h_o - (Dx2 + Dx3);
d_i = d_o - 2 * Dx1;
A_rad = pi * d_i * h_i;

h = (2 * k) / (log(d_o / d_i) * d_i);

Q_total = h * A_rad * dT;

end