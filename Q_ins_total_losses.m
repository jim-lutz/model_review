function Q_total = Q_ins_total_losses(d_o, h_o, Dx1, Dx2, Dx3, k, DT)
% FUNCTION NAME:
%   Q_ins_total_losses
%
% DESCRIPTION:
% Compute the total internal heat losses for the tank depicted in Figure E.4 of
% AS/NZS 4234:2021.
%
%
% INPUT:
%   d_o - Outer diameter of hot water cylinder (m).
%   h_o - Outer height of hot water cylinder (m).
%   Dx1 - Thickness of the radial insulation of the hot water cylinder (m).
%   Dx2 - Thickness of the top insulation of the hot water cylinder (m).
%   Dx3 - Thickness of the bottom insulation of the hot water cylinder (m).
%   k - Thermal conductivity of the hot water insulation (W.m-1.K-1).
%   dT - Temperature rise of the hot water cylinder exterior to ambient (K).
%
% OUTPUT:
%   Q_total - Total internal losses of the hot water cylinder through the
%     insulation.through base of hot water cylinder (W).
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

[Q_top, ~, ~] = Q_ins_end_losses(d_o, Dx1, Dx2, k, DT);
[Q_bottom, ~, ~] = Q_ins_end_losses(d_o, Dx1, Dx3, k, DT);
[Q_side, ~, ~, ~, ~] = Q_ins_side_losses(d_o, h_o, Dx1, Dx2, Dx3, k, DT);

Q_total = Q_top + Q_bottom + Q_side;

end