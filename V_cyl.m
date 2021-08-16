function V = V_cyl(d_o, h_o, Dx1, Dx2, Dx3)
% FUNCTION NAME:
%   V_cyl
%
% DESCRIPTION:
% Compute the internal volume of a hot water cylinder as depicted in Figure E.4
% of AS/NZS 4234:2021.
%
% INPUT:
%   d_o - Outer diameter of hot water cylinder (m).
%   h_o - Outer height of the hot water cylinder (m).
%   Dx1 - Thickness of the radial insulation of the hot water cylinder (m).
%   Dx2 - Thickness of the top insulation of the hot water cylinder (m).
%   Dx3 - Thickness of the bottom insulation of the hot water cylinder (m).
%
% OUTPUT:
%   V - Hot water cylinder storage capacity (m3).
%
% ASSUMPTIONS AND LIMITATIONS:
%   This function has not been vectorised. Only one argument may be a vector.
%
% REVISION HISTORY:
%   20210709 - bjem
%       * Initial implementation
V = (pi * (((d_o - 2 * Dx1) / 2)^2) * (h_o - (Dx2 + Dx3)));

1;

end