function Q_MEPS = Q_MEPS_4606_Table_1(V)
% FUNCTION NAME:
%   Q_MEPS_4606_Table_1
%
% DESCRIPTION:
% Compute the maximum allowable hot water cylinder heat losses (MEPS losses)
% from standard NZS 4606:1989 according to clause 2.4.1(c) and Table 1.
% This also covers the equivalent losses in NZS 4602:1988.
%
% The temperature rise is 55.6 K.
%
% The heat loss calculation is performed according to 2.4.1(c) and Table 1,
% with the following extra notes:
% * despite anything suggested (or not made clear) by 2.4.1(c), the heat loss
%   is computed so that it never decreases as a function of the hot water
%   cylinder capacity.
%
% To implement this condition, the heat loss model is treated in the following
% ways:
% 1. The allowable heat loss from the equations is pre-computed and put into a
%    lookup table for linear interpolation. The resulting interpolated value is
%    then rounded to the nearest 0.1 kWh.day-1.
% 2. The allowable losses in Table 1 proceed in a "slab" fashion. The Table 1
%    losses are assumed to be the same from one capacity listed to just below
%    another e.g. the maximum allowed loss for a hot water cylinder with
%    capacity of 90 L to just less than 135 L (135 * (1 - eps)) is assumed to be
%    1.2 kWh.day-1 by this method.
% 3. The two canditate loss values are compared, and the maximum value is taken
%    as the maximum allowed losses.
%
% The loss function is partially optimised to remove unneccessary comparisons
% where one type of loss will clearly always be chosen.
%
% INPUT:
%   V - Hot water cylinder storage capacity (m3).
%
% OUTPUT:
%   Q_MEPS - Maximum allowed losses that are allowed for the hot water cylinder
%      with storage volume V (kWh.day-1). The losses are rounded to the nearest
%      0.1 kWh.day-1.
%
% ASSUMPTIONS AND LIMITATIONS:
%   This function has not been vectorised. Only one argument may be a vector.
%
%
% REVISION HISTORY:
%   20210709 - bjem
%       * Initial implementation
Q_table =[...
0.0000, 0.4; ...
% 0.0065, 0.5; ...
% 0.0130, 0.5; ...
% 0.0220, 0.6; ...
% 0.0450, 0.8; ...
0.0900, 1.2; ...
% 0.1350, 1.4; ...
% 0.1800, 1.6; ...
% 0.2250, 1.8; ...
% 0.2700, 2.0; ...
0.3600, 2.5; ...
%0.4500, 2.9; ...
0.5400, 3.4; ...
0.6300, 3.8; ...
];

Q_table2 =[...
0.0000, 0.4;
0.0900 - eps(0.09), 1.156; ...
0.09, 1.152; ...
0.63, 3.744; ...
0.71, 4.128; ...
];


Q_MEPS_spot = interp1(Q_table(:,1), Q_table(:,2), V, 'previous');

Q_MEPS_interp = ceil(10 * interp1(Q_table2(:,1), Q_table2(:,2), V)) / 10;

Q_MEPS = max(Q_MEPS_spot, Q_MEPS_interp);

end