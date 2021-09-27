%
%   Description: Optimisation file for hot water cylinder losses evaluation.
%   The computations of losses are based on those in AS/NZS 4234:2021
%   Appendix E.
%   
%   Author: Bart Milne (bart.milne@eeca.govt.nz)
%
%   Date: 14 June 2021
%
%   Comment: Optimiser for various models of hot water cylinders.
%   
%   In addition to any assumptions in Appendix E; the following assumptions are
%   made about the cylinder construction:
%   * inner cylinder walls are of negligible thickness compared to the
%   * insulation;
%   * inner and outer cylinder wall thermal resistance is negligible
%     compared to the insulation;
%   * cylinder insulation is the same in all parts of the cylinder;
%   * cylinder bottom is insulated with the same material as the top and side;
%   * insulation thickness at the top and bottom of the cylinder is equal;
%   * one element only in the cylinder;
%   * two pipe connections are present that are not insulated;
%   * heat losses occur from cold connection the same as for the hot connection;
%   * TPR valve is present and not insulated;
%   * thermostat pocket has an effective diameter of 35 mm;
%   * temperature rise above ambient is 55 K;
%   * cylinder exterior walls are held at a constant temperature of 20?C;
%   * cylinder insulation is polyurethane with a bulk thermal conductivity of
%     0.025 W.m-1.K-1.

clear; close; clc
  # not sure what this does, but didn't squawk
pkg load statistics
  # ran with no output
pkg load optim
  # ran with no output

filename = 'cylinder_data_202021_Check_Testing.xlsx'
  #{ 
  this is the name of the input file. 
  see also data_maps.ods for summary of data sets
  this spreadsheet doesn't have test results in it?
  #}
%filename = 'cylinder_data.xlsx'
filename_output = 'out.xlsx'
  # name of results file

[M_opt, T_data, R_data] = xlsread(filename);
  #{  
'xlsread' is a function from the file /home/jiml/snap/octave/78/octave/io-2.6.3/xlsread.m

Return argument NUMARR (M_opt) contains the numeric data, optional return
arguments TXTARR (T_data) and RAWARR (R_data) contain text strings and the raw
spreadsheet cell data, respectively.
  T_data has text in rows 1 & 2   
  #}  
% Get rid of the first three columns in the numerical data array.
% They contain a serial counter and name and model information which is not
% used for optimisation.
M_opt(:,1:3) = [];

% The remaining columns in M_opt are: 
% Capacity (m3),
% Diameter (m),
% Height (m),
% Height/Diameter Ratio,
% Number of fittings,
% Number of TPR valves,
% Number of thermostat pockets,
% Thermal conductivity of insulation (W.m-1.K-1),
% Standard (4692 means AS/NZS 4692, 4606 means NZS 4602)
  # seems to have done that

% Some physical constants used in the modelling.
% The quantities in which watts appear in the SI units
% are converted to kWh.day-1 by multiplying by 0.024.
%
stefan_boltzmann = 5.670374419e-8 * 0.024
emissivity = 0.05
  #{
  This doesn't look right!!
    from ASHRAE Fundamentals, 
    CHAPTER 4. HEAT TRANSFER
    3. THERMAL RADIATION
    Table 5 Emissivities and Absorptivities of Some Surfaces
    Emissivity of white paint is 0.90 - 0.85  
  from report 'Optimisation and Statistical Estimation of Hot water cylinder Losses According to ASNZS 4234 Appendix E'  
    "� is the surface emissivity (dimensionless, 
    assumed to be 0.4 for metal hot water cylinder exteriors)."
  #}
k_wood = 0.130 * 0.024
k_air  = 0.024 * 0.024
  #{
  from ASHRAE Fundamentals, CHAPTER 33. Physical PROPERTIES OF MATERIALS
  Table 3 Properties of Solids
  Thermal Conductivity, W/(m-K)
    Softwoods 0.11/0.16
    0.130 looks plausible 
  Table 1 Properties of Vapor
  Thermal Conductivity, W/(m-K)
    Nitrogen    0.0240
    Oxygen      0.0244
    Water vapor 0.0247
    0.024 looks fine
  #}
% The temperature rise and thermal conductivity. The thermal conductivities are
% taken from the model losses in Appendix E. The multiplication by 0.024 is to
% convert W (watts) into kWh.day-1 so the thermal conductivity becomes in units
% of kWh.day-1.m-1.K-1.
dT_4692 = 55;
dT_4606 = 55.6;
k_per_fitting = (5 * 0.0269) / 55;
  # pipe diameter = 0.0269 = 1.059 inches
  # assumed insulated, OK?
k_per_TPR = 0.12 / 55;
  # 0.12 = d_pipe * 5 ?
  # d_pipe = 0.024 = 0.945 inches, Probably OK
h_thermopocket = 27 * (pi * 0.035^2 / 4) / 55;
  # 0.035 m = 1.38 inches radius?
T_amb = 20;
  #{
    from 4234, Table E.9 ? Heat loss through tank fitttings 
        Tank opening (e.g. themostat pocket)            A_o * 27
          A_o denotes the area of the opening in the insulation 
          in square metres (e.g. for a thermostat opening).
        Upinsulated pipe gr fitting (e.g PTR valve)     d_pipe * 5
        Insulated pipe or fitting                       d_pipe * 3.5
  #}
#
% Volume ratio. Assume target volume is actually a multiple of the volume.
% This is based on the provision AS/NZS 4692 that the volume must be from 1
% to 1.05 times the contained volume.
V_ratio_base = 1;
  #{
    This is  from 4692.1, 5.2.2 Nominal capacity?New Zealand
    In standard this is the allowed ratio of measured to nominal volume
    to keep mfrs from claiming more volume than they provide
  #}
% Get MEPS level for every cylinder from Table A5 of AS/NZS 4692.1. This MEPS
% level does not include the TPR valve allowance. The offset is added as
% required.
Q_MEPS_4692 = Q_MEPS_4692_Table_A5(M_opt(:,1)) + (M_opt(:,5) > 0) * 0.12;
  # M_opt(:,1) is Capacity (m3)
  # M_opt(:,5) is number of TPR Valves
  # Q_MEPS_4692_Table_A5() is in Q_MEPS_4692_Table_A5.m
Q_MEPS_4606 = Q_MEPS_4606_Table_1(M_opt(:,1));
  # the calculations in Table 1 and nearby equations are a rather strange
  # way of doing MEPS, value in table if that volume, equation if not
  # seems like this is implemented OK
Q_MEPS = Q_MEPS_4692;
  # looks like it matchs TABLE A5 from 4692.2
Q_MEPS(M_opt(:, 10) == 4606) = Q_MEPS_4606(M_opt(:, 10) == 4606);


% Compute the z-score that 95% of cylinders would fall below. The value is based
% on having a mean insulation losses of 1.35 times the "basic" (unpenalised)
% insulation losses. Therefore, at least 95% of cylinders are expected to have
% losses below 1.7 times the basic losses. Use the difference between the
% minimum losses and the mean losses to get the standard deviation.
z_target = norminv(0.95);
r_penalty_ = 1.35
r_penalty_mult = 1
r_penalty = 1 + (r_penalty_ - 1) * r_penalty_mult
r95_penalty = 1 + 2 * (r_penalty_ - 1) * r_penalty_mult
Q_sd_ratio = round(1e6 * ((r_penalty - 1) / z_target)) / 1e6;
# haven't checked this, he's using 1.35 as a variable.

% Iterate through the optimisation for every cylinder model (row in M_opt
% matrix).
N = size(M_opt, 1)
for ii=1:N
    % Extract the cylinder volume, outer diameter, outer height, number of
    % fittings, number of TPR valves and thermostat pocket effective diameter.
    V_cyl_i = M_opt(ii, 1);
    d_o = M_opt(ii, 2);
    h_o = M_opt(ii, 3);
    n_fittings = M_opt(ii, 4);
    n_TPR = M_opt(ii,5);
    n_thermopocket = M_opt(ii,6);
    k = M_opt(ii,7) * 0.024;
    emissivity = M_opt(ii,8);
    v_forced = M_opt(ii,9);
    std_ = M_opt(ii, 10);
    
    % Compute the dimensions of the inner cylinder as if it were a uniformly
    % scaled version of the outer cylinder with the design volume. This 
    % computation is used to get initial values for the optimal insulation
    % thickness optimisation guess.
    V_cyl_o = pi * (d_o^2 * h_o) / 4;
    tank_ratio = cbrt(V_cyl_i / V_cyl_o);
    # Compute the real cube root of each element of X.
    
    % Use the tank volume ratio dimensions to provide initial estimates for the
    % optimum insulation thickness. Make sure the upper bound is never greater
    % than half the outer diameter or half the outer height.
    t0_ratio = 5;
    t0 = (1 - tank_ratio) * [d_o; h_o] / 2;
    lb = t0 * 1 / t0_ratio;
    ub = min(t0 * t0_ratio, 0.5 * [d_o; h_o]);
    
    % Optimisation and equality functions. The first function is the insulation
    % losses to optimise. The second implements the constraint that the volume
    % of the inner tank must match the volume as designed.
    if (std_ == 4606)
        dT = dT_4606; # 55.6
    else
        dT = dT_4692; # 55
    end
    
    % Define the optimisation functions and equality constraint functions.
    % The actual temperature rise used does not affect the optimisation,
    % as the temperature rise is a simple multiplier on the quantity to be
    % optimised.
    opt_func = @(t) Q_ins_total_losses(d_o, h_o, t(1), t(2), t(2), k, dT);
    eq_func  = @(t) V_cyl(d_o, h_o, t(1), t(2), t(2)) - V_ratio_base * V_cyl_i;
    
    % Perform a sequential quadratic programme to optimise the insulation
    % thickness. This seems quite well-behaved numerically.
    [x, obj, info, iter, nf, lambda] = sqp(t0, opt_func, [], eq_func, lb, ub);
    
    % Extract the radial end end insulation thicknesses that yield the minimum
    % losses. The thickness is in metres. Round down to nearest mm to ensure
    % volume is preserved.
    t_rad_min  = floor(x(1) * 1000) / 1000;
    t_ends_min = floor(x(2) * 1000) / 1000;
    
    % Compute a "check volume" for comparison to the nominal volume. This is
    % intended to never be less than the nominal volume.
    V_cyl_check = floor(V_cyl(d_o, h_o, t_rad_min, t_ends_min, t_ends_min) * 1000) / 1000;

    % Get the inner cylinder dimensions that correspond to the minimum losses.
    % I also calculate the aspect ratio, as there seems to be a correlation
    % between the aspect ratio and chance of passing.
    d_i = d_o - 2 * t_rad_min;
    h_i = h_o - 2 * t_ends_min;
    o_ratio = round(10000 * h_o / d_o) / 10000;
    i_ratio = round(10000 * h_i / d_i) / 10000;
    
    % Q_ins_min is the "basic" insulation losses. Q_ins_total is Q_ins_min with an
    % added penalty to take into account insulation imperfections. This penalty
    % is set at 35% in the standard, but the actual penalty level can be
    % adjusted using the r_penalty_mult variable.
    Q_ins_min = Q_ins_total_losses(d_o, h_o, t_rad_min, t_ends_min, t_ends_min, k, dT);
    Q_ins_tot = r_penalty * Q_ins_min;

    % Fittings losses. These are computed using a function for convenience,
    % but use the parameters in Table E.9.
    Q_fittings = Q_fittings_TOT(n_fittings, n_TPR, n_thermopocket, dT, k_per_fitting, k_per_TPR, h_thermopocket);

    % The variance is assumed not to be a simple multiple of the fittings
    % value, as it is assumed that the individual differences between the 
    % fittings are random.
    Q_fittings_var = Q_fittings_TOT_var(n_fittings, n_TPR, n_thermopocket, dT, k_per_fitting, k_per_TPR, h_thermopocket) * Q_sd_ratio;
      # this function is from Q_fittings_TOT_var.m
      # dT - Temperature rise of the hot water cylinder exterior to ambient (K).
      # shouldn't this be from water temp to ambient???
      # in 4342 E.9 Fittings heat loss, these are from water temp to ambient
      # these are losses through these feed throughs. the thermal path
      # is through the metal of the pipe, TPR, or reduced insulation around thermopocket

    % The total losses are the sum of the insulation losses (with penalty)
    % and the fittings losses.
    Q_TOT = round(1000 * (Q_ins_tot + Q_fittings)) / 1000;
    
    % The standard deviation is computed by adding variances.
    Q_ins_tot_var = (Q_sd_ratio * Q_ins_tot)^2;
    Q_TOT_var = Q_ins_tot_var + Q_fittings_var; 
    Q_TOT_sd = round(1000 * sqrt(Q_TOT_var)) / 1000;

    % The 95th percentile heat losses are given to compare these to the MEPS
    % levels.
    Q_95_min = round(1000 * ((2 * r_penalty - 1) * Q_ins_min + Q_fittings)) / 1000;
    Q_95 = round(1000 * (Q_TOT + z_target * Q_TOT_sd)) / 1000;
    Q_95_max = round(1000 * ((2 * r_penalty - 1) * (Q_ins_min + Q_fittings))) / 1000;

    % Compute the z score and the probability of passing the MEPS threshold.
    Q_TOT_z_score = round(10000 * ((Q_MEPS(ii) - Q_TOT) / Q_TOT_sd)) / 10000;
    
    % Q_probability is the chance of passing Stage 1 check testing.
    % Q_probability2 is the probability of passing Stage 2 check testing.
    % The probabilities are different because the hot water cylinder must
    % fail Stage 1 testing before Stage 2 testing is performed. The chance
    % of passing Stage 2 testing is slightly higher because there is a finite
    % chance that the two extra cylinders purchased may both pass.
    Q_probability = round(10000 * normcdf(Q_TOT_z_score)) / 10000;
    Q_probability2 = round(10000 * (Q_probability * (1 + Q_probability - Q_probability^2))) / 10000;
    
    % Compute the extra insulation thickness required on the top and side to
    % achieve MEPS compliance.
    % Find the insulation losses that would be required. for the 95th percentile
    % to achieve MEPS. The solving function is simplified in light of
    % r_penalty - 1 representing the 45th percentile above the mean implicitly.
    % This problem can be solved directly for obtaining the thickness
    % of the insulation, but it is easier to split the problem into two steps.
    %
    % Step 1: Find the insulation losses that would be required.
    %
    % Step 2: Find the insulation thickness that is required to achieve the
    % target losses.
    %
    % The equation is more complicated because the standard deviations of the
    % insulation losses and the fittings losses are assumed not to add together
    % directly.
    fsolve_fcn_losses95 = @(Q__) r_penalty * Q__ + Q_fittings + z_target * sqrt(Q_sd_ratio^2 * r_penalty^2 * Q__^2 + Q_fittings_var) - Q_MEPS(ii);

    % Step 1: Get the maximum insulation losses via numerical solution.
    Q_MEPS_ins95 = fsolve(fsolve_fcn_losses95, Q_TOT);
    
    % Step 2a: Find the required insulation thickness, assuming the extra
    % insulation thickness required is uniform all around the hot water
    % cylinder. The insulation thickness is always rounded up to the
    % nearest mm.
    %
    % The initial guess is no extra insulation at all.
    %
    fsolve_fcn = @(t) (Q_ins_total_losses(d_o + 2 * t, h_o + 2 * t, t_rad_min + t, t_ends_min + t, t_ends_min + t, k, dT)) - Q_MEPS_ins95;
    t_meps_uniform = fsolve(fsolve_fcn, 0);
    t_meps_uniform = ceil(1000 * t_meps_uniform) / 1000;
    
    t_rad_uniform_meps = ceil(1000 * (t_rad_min  + t_meps_uniform)) / 1000;
    t_ends_uniform_meps = ceil(1000 * (t_ends_min + t_meps_uniform)) / 1000;
    d_o_meps_uniform = ceil(1000 * (d_o + 2 * t_meps_uniform)) / 1000;
    h_o_meps_uniform = ceil(1000 * (h_o + 2 * t_meps_uniform)) / 1000;
    aspect_ratio_meps_uniform = h_o_meps_uniform / d_o_meps_uniform;
    
    % Step 2b: Find the insulation thickness required if the radial and end
    % thicknesses of insulation are allowed to vary independently of each
    % other. This means that the top and bottom insulation thicknesses are
    % assumed to be identical. The solver uses the individual variations
    % of heat losses with respect to each thickness (the Jacobian) to
    % find the most optimal path to proceed on. There is no guarantee that
    % the solution found is optimal in any sense (e.g. minimal insulation
    % volume) except that the sum of squares of the insulation thicknesses is
    % minimal.
    %
    % The initial guess is no extra insulation at all.
    %
    fsolve_fcn2 = @(t) (Q_ins_total_losses(d_o + 2 * t(1), h_o + 2 * t(2), t_rad_min + t(1), t_ends_min + t(2), t_ends_min + t(2), k, dT)) - Q_MEPS_ins95;
    
    t_meps = fsolve(fsolve_fcn2, [0 0]');
    t_meps = ceil(1000 * t_meps) / 1000;
    t_rad_extra_meps = t_meps(1);
    t_top_bottom_extra_meps = t_meps(2);
    
    t_rad_meps    = ceil(1000 * (t_rad_min  + t_rad_extra_meps)) / 1000;
    t_top_meps    = ceil(1000 * (t_ends_min + t_top_bottom_extra_meps)) / 1000;
    t_bottom_meps = ceil(1000 * (t_top_meps)) / 1000;
    d_o_meps = ceil(1000 * (d_o + 2 * t_rad_extra_meps)) / 1000;
    h_o_meps = ceil(1000 * (h_o + 2 * t_top_bottom_extra_meps)) / 1000;
    aspect_meps = h_o_meps / d_o_meps;
    
    Q_ins_compliance = Q_ins_total_losses(d_o_meps, h_o_meps, t_rad_meps, t_top_meps, t_bottom_meps, k, dT);
    
    Q_ins_tot_compliance = r_penalty * Q_ins_compliance;
    Q_ins_tot_compliance_var = (Q_sd_ratio * Q_ins_tot_compliance)^2;
    Q_ins_tot_compliance_sd = round((1000 * sqrt(Q_ins_tot_compliance_var))) / 1000;
    
    Q_meps_compliance = round(1000 * (Q_ins_tot_compliance + Q_fittings)) / 1000;
    Q_meps_compliance_var = sum([Q_ins_tot_compliance_sd^2, Q_fittings_var]);
    Q_meps_compliance_sd = round(1000 * sqrt(Q_meps_compliance_var)) / 1000;
    
    Q_meps_compliance_check = Q_meps_compliance + z_target * Q_meps_compliance_sd;

    % Hot water cylinder Exterior Temperature rise
    %
    % In this case, the hot water cylinder temperature rise is assumed to be
    % variable.
    %
    % Heat dissipation is modelled in the following ways:
    % * Top: Convection & Radiation
    % * Side: Convection & Radiation
    % * Bottom: Conduction
    %
    % Assume the outer cylinder has a constant temperature rise above ambient.
    % Find an exterior tank temperature that balances the interior heat
    % transfer and the exterior heat transfer.
    %
    % This is the basic function for balancing the heat transfer. Q__ is the
    % heat dissipated when the hot water cylinder exterior has zero temperature
    % rise above ambient. It is assumed this heat transfer scles linearly with
    % the temperature rise from the hot water cylinder exterior to the exterior.
    %
    % The nett effect is that the losses during testing will decrease if the
    % exterior surface of the hot water cylinder is allowed to rise above
    % ambient. This is because the temperature difference from the hot water
    % cylinder interior tank to the exterior is reduced.
    %
    T_ext_fcn_basic = @(Q__, r, T) ...
    Q__ * ((dT - T) / dT) - ...
    (Q_ext_side_losses(k_air, d_o, h_o, T, T_amb, emissivity, stefan_boltzmann, v_forced) + ...
     Q_ext_top_losses (k_air, d_o, T, T_amb, emissivity, stefan_boltzmann) + ...
     Q_ext_bottom_losses(k_wood, d_o, T));
    
    % Create the zero finding function. The total losses (including fittings
    % and TPR valve losses) are included because it is assumed these will
    % dissipate into the same air temperature as the hot water cylinder
    % exterior.
    T_ext_fcn_TOT = @(T) T_ext_fcn_basic(Q_TOT, 1, T);
    
    % Find the value of the exterior hot water cylinder temperature rise
    % that gives an energy balance. I use lsqnonlin here rather than
    % fsolve to prevent numerical problems if fsolve tries a negative
    % value of dT.
    %
    dT_ext = lsqnonlin(T_ext_fcn_TOT, dT/2, 0, dT);
      # Solve nonlinear least-squares (nonlinear data-fitting) 
      # problems min [EuclidianNorm(f(x))] .^ 2 x
      # part of optim package
    dT_ext_(ii,1) = dT_ext;
    
    dT_mod(ii,1) = dT - dT_ext;
    
    Q_TOT_tr(ii,1) = round(1000 * Q_TOT * (dT - dT_ext) / dT) / 1000;
    Q_TOT_tr_sd(ii,1) = round(1000 * Q_TOT_sd * (dT - dT_ext) / dT) / 1000;
    Q_TOT_tr_z_score(ii,1) = round(1000 * (Q_MEPS(ii,1) - Q_TOT_tr(ii,1)) / Q_TOT_tr_sd(ii,1)) / 1000;
    Q_TOT_tr_probability(ii,1) = round(10000 * normcdf(Q_TOT_tr_z_score(ii,1))) / 10000;
    Q_TOT_tr_probability2(ii,1) = round(10000 * (Q_TOT_tr_probability(ii,1) * (1 + Q_TOT_tr_probability(ii,1) - Q_TOT_tr_probability(ii,1)^2))) / 10000;
        
    [Q_side_tr  , ~, ~] = Q_ext_side_losses(k_air, d_o, h_o, dT_ext, T_amb, emissivity, stefan_boltzmann, v_forced);
    [Q_top_tr   , ~, ~] = Q_ext_top_losses (k_air, d_o, dT_ext, T_amb, emissivity, stefan_boltzmann);
    [Q_bottom_tr, ~, ~] = Q_ext_bottom_losses(k_wood, d_o, dT_ext);
    
    1;
    
    % Create the row of the results matrix.
    V(ii,:) = [...
    dT, ...
    t_rad_min, ...
    t_ends_min, ...
    d_i, ...
    h_i, ...
    V_cyl_check, ...
    o_ratio, ...
    i_ratio, ...
    Q_ins_min, ...
    Q_ins_tot, ...
    Q_fittings, ...
    Q_TOT, ...
    Q_TOT_sd, ...
%    Q_ins_tot_var, ...
%    Q_fittings_var, ...
%    Q_TOT_var, ...
%    Q_95, ...
    Q_MEPS(ii), ...
    Q_TOT_z_score, ...
    Q_probability, ...
    Q_probability2, ...
    dT_mod(ii,1), ...
    dT_ext, ...
    Q_TOT_tr(ii,1), ...
    Q_TOT_tr_sd(ii,1), ...
    Q_TOT_tr_z_score(ii,1), ...
    Q_TOT_tr_probability(ii,1), ...
    Q_TOT_tr_probability2(ii,1), ...
%    Q_ins_compliance, ...
%    Q_ins_tot_compliance,...
%    Q_fittings, ...
%    Q_meps_compliance, ...
%    Q_meps_compliance_sd, ...
%    Q_fittings_var, ...
%    Q_ins_tot_compliance_var, ...
%    Q_meps_compliance_var, ...
%    t_meps_uniform, ...
%    d_o_meps_uniform, ...
%    h_o_meps_uniform, ...
%    aspect_ratio_meps_uniform, ...
%    t_rad_extra_meps, ...
%    t_top_bottom_extra_meps, ...
%    t_rad_meps, ...
%    t_top_meps, ...
%    d_o_meps, ...
%    h_o_meps, ...
%    round(1000 * aspect_meps) / 1000, ...
%    round(1000 * Q_meps_compliance_check) / 1000, ...
    ];
   
    1;
end


% The columns in M_opt are:
cols = {};

cols(1,:) = { ...
'Temperature Rise used for computing losses statistics (K)', ...
'Optimal Radial Insulation Thickness (m)', ...
'Optimal End Insulation Thickness (m)', ...
'Inner Diameter (m)', ...
'Inner Height (m)', ...
'Inner Volume Check (m3)', ...
'External Aspect (height/diameter) ratio', ...
'Internal Aspect (height/diameter) ratio', ...
'Minimum insulation losses (without penalty) (kWh.day-1)', ...
'Minimum insulation losses (with penalty) (kWh.day-1)', ...
'Losses from fittings (including TPR valve) (kWh.day-1)', ...
'Total losses (kWh.day-1)', ...
'Total losses standard deviation (kWh.day-1)', ...
'MEPS losses including TPR valve where applicable (kWh.day-1)', ...
'Calculated z score to pass MEPS', ...
'Probability of any cylinder complying with MEPS given the z score calculated', ...
'Probability of any cylinder passing stage 2 check testing', ...
'Modified temperature rise to hot water cylinder exterior (K)', ...
'Exterior to ambient temperature rise (K)', ...
'Thermal losses at modified temperature rise (kWh.day-1)', ...
'Thermal losses at modified temperature rise standard deviation (kWh.day-1)', ...
'Calculated z score to pass MEPS with modified T rise', ...
'Probability of any cylinder with modified T rise complying with MEPS given the z score calculated', ...
'Probability of any cylinder with modified T rise passing stage 2 check testing', ...
};




cols(2,:) = { ...
'dT', ...
't_rad_min', ...
't_ends_min', ...
'd_i', ...
'h_i', ...
'V_cyl_check', ...
'o_ratio', ...
'i_ratio', ...
'Q_ins_min', ...
'Q_ins_tot', ...
'Q_fittings', ...
'Q_TOT', ...
'Q_TOT_sd', ...
'Q_MEPS', ...
'Q_TOT_z_score', ...
'Q_probability', ...
'Q_probability2', ...
'dT_mod', ...
'dT_ext', ...
'Q_TOT_tr', ...
'Q_TOT_tr_sd', ...
'Q_TOT_tr_z_score', ...
'Q_TOT_tr_probability', ...
'Q_TOT_tr_probability2', ...
};

% Append the extra columns to a new output spreadsheet. The existing data
% in the input spreadsheet is copied across first, then the new data is
% appended to that.
xlswrite(filename_output, R_data)
xlswrite(filename_output, cols, 1, 'N1')
xlswrite(filename_output, V, 1, 'N3')

%%%%% END OF FILE %%%%%
