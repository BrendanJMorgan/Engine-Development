load('impulse_turbine_efficiency.mat');

turbine_shaft_speed = shaft_speed*gear_ratio; % rad/s

mach_spouting = sqrt(2/(gamma_gg-1)*((p_amb/p_gg)^((1-gamma_gg)/gamma_gg)-1)); % unitless
T_spouting = Tt_gg*(p_amb/p_gg)^((gamma_gg-1)/gamma_gg); % K - static temperature for perfectly expanded flow
mw_gg = 0.001*cea_gg.output.eql.mw(3); % kg/mol - molecular weight (fully combusted)
R_gg = 8.3145 / mw_gg; % J/mol*K - specific gas constant
sonic_spouting = sqrt(gamma_gg*R_gg*T_spouting); % m/s - speed of sound at nozzle exit plane

v_spouting_norm = mach_spouting*sonic_spouting; % m/s - speed at which gas leaves the nozzles

%% Rotor Dimensions
% One stage impulse rotor - does not work for multistage and/or velocity driven turbines
v_pitchline = sin(nozzle_angle)*v_spouting_norm/2; % m/s
v_spouting = [2*v_pitchline; sqrt(v_spouting_norm^2 - 4*v_pitchline^2)]; % [m/s, m/s]
r_pitchline = v_pitchline / turbine_shaft_speed; 

isentropic_v_ratio = v_pitchline/v_spouting_norm; % unitless
stage_efficiency = interp1(impulse_turbine_efficiency(:,1), impulse_turbine_efficiency(:,2), isentropic_v_ratio, 'spline'); % unitless
turbine_power = turbine_shaft_power/stage_efficiency; % W

mdot_gg = turbine_power/(2*v_pitchline^2); % kg/s

%% Rotor Blades

v_blade_inlet_rotor = [v_pitchline; v_spouting(2)]; % [m/s, m/s] - relative velocity 
v_blade_outlet_rotor = [-v_pitchline; v_spouting(2)]; % [m/s, m/s] - relative_velocity
blade_angle_inlet_rotor = atan2(v_blade_inlet_rotor(1), v_blade_inlet_rotor(2)); % rad
blade_angle_outlet_rotor = atan2(v_blade_outlet_rotor(1), v_blade_outlet_rotor(2)); % rad

dens_gas_spouted = cea_gg.output.eql.density(3); % kg/m3
blade_length_rotor = sqrt(4*r_pitchline*mdot_gg/(blade_number_rotor*blade_gap_rotor*sec(blade_angle_inlet_rotor)*dens_gas_spouted*v_spouting(2))); % m - blades are unshrouded

r_rotor_tip = r_pitchline + blade_length_rotor/2; % m
r_rotor_base = r_pitchline - blade_length_rotor/2; % m

%% Mass Flow Adjustments
mdot_total = mdot_cc + mdot_gg; 
gg_fraction = mdot_gg / mdot_total; % not totally sure if denominator should be cc instead of total
mdot_fuel_gg = mdot_gg*(1/(1+OF_gg)); % kg/s - Fuel Mass Flow Rate
mdot_ox_gg = mdot_gg*(OF_gg/(1+OF_gg)); % kg/s - Oxidizer Mass Flow Rate
mdot_fuel_total = mdot_fuel_cc+mdot_fuel_gg; % kg/s
mdot_ox_total = mdot_ox_cc+mdot_ox_gg; % kg/s

%% Nozzle Plate
A_throat_np = mdot_gg*c_star_gg/p_gg/nozzle_number; % m2
d_throat_np = 2*sqrt(A_throat_np/pi); % m
d_throat_np/0.0254




