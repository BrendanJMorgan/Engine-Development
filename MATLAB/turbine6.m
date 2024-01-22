load('impulse_turbine_efficiency.mat');

cea_gg = CEA('problem','rocket', ...
     'equilibrium', ...  % equilibrium flow, infinite chamber,
     'o/f',OF_gg, ...          % Define OF ratio
     'case','CEAM-rocket1', ...
     'p(psi)',p_gg/6894.76,'pi/p',p_gg/p_amb, ... % Define nozzle with pressure ratio
     'reactants', ... 
     'fuel','C2H5OH(L)','wt%',proof,'t(k)',T_amb,'H2O(L)','wt%',1-proof,'t(k)',T_amb, ... % 190-Proof Ethanol for Fuel, at ambient
     'oxid','O2(L)','wt%',100, ... % Liquid Oxygen for Oxidizer, at boiling point
     'output','short','massf','transport','mks', ...
     'end');

turbine_shaft_speed = shaft_speed*gear_ratio; % rad/s

mach_spouting = sqrt(2/(gamma_gg-1)*((p_amb/p_gg)^((1-gamma_gg)/gamma_gg)-1)); % unitless
T_spouting = Tt_gg*(p_amb/p_gg)^((gamma_gg-1)/gamma_gg); % K - static temperature for perfectly expanded flow
mw_gg = 0.001*cea_gg.output.eql.mw(3); % kg/mol - molecular weight (fully combusted)
R_gg = 8.3145 / mw_gg; % J/mol*K - specific gas constant
sonic_spouting = sqrt(gamma_gg*R_gg*T_spouting); % m/s - speed of sound at nozzle exit plane

v_spouting_norm = mach_spouting*sonic_spouting; % m/s - speed at which gas leaves the nozzles

%% Rotor Dimensions
% One stage impulse rotor - does not work for multistage and/or velocity driven turbines
nozzle_angle = 70*pi/180; % empirical recommendation by Farokhi, page 691
v_pitchline = sin(nozzle_angle)*v_spouting_norm/2; % m/s
v_spouting = [2*v_pitchline; sqrt(v_spouting_norm^2 - 4*v_pitchline^2)]; % [m/s, m/s]
r_pitchline = v_pitchline / turbine_shaft_speed; 

isentropic_v_ratio = v_pitchline/v_spouting_norm; % unitless
stage_efficiency = interp1(impulse_turbine_efficiency(:,1), impulse_turbine_efficiency(:,2), isentropic_v_ratio, 'spline'); % unitless
turbine_power = shaft_power/stage_efficiency; % W

mdot_gg = turbine_power/(2*v_pitchline^2) % kg/s

%% Rotor Blades
blade_number_rotor = 64; % unitless - total number of blades on the rotor
blade_gap_rotor = 1/4*0.0254; % m
dens_gas_spouted = cea_gg.output.eql.density(3); % kg/m3
blade_length_rotor = sqrt(4*r_pitchline*mdot_gg/(blade_number_rotor*blade_gap_rotor*dens_gas_spouted*v_spouting(2))); % m - blades are unshrouded

v_blade_inlet_rotor = [v_pitchline; v_spouting(2)]; % [m/s, m/s]
v_blade_outlet_rotor = [-v_pitchline; v_spouting(2)]; % [m/s, m/s]
blade_angle_inlet_rotor = atan2(v_blade_inlet_rotor(1), v_blade_inlet_rotor(2)); % rad
blade_angle_outlet_rotor = atan2(v_blade_outlet_rotor(1), v_blade_outlet_rotor(2)); % rad

r_rotor_tip = r_pitchline + blade_length_rotor/2; % m
r_rotor_base = r_pitchline - blade_length_rotor/2; % m



%% Mass Flow Adjustments
mdot_total = mdot_cc + mdot_gg; 
gg_fraction = mdot_gg / mdot_total; % not totally sure if denominator should be cc instead of total
mdot_fuel_gg = mdot_gg*(1/(1+OF_gg)); % kg/s - Fuel Mass Flow Rate
mdot_ox_gg = mdot_gg*(OF_gg/(1+OF_gg)); % kg/s - Oxidizer Mass Flow Rate
mdot_fuel_total = mdot_fuel_cc+mdot_fuel_gg; % kg/s
mdot_ox_total = mdot_ox_cc+mdot_ox_gg; % kg/s




