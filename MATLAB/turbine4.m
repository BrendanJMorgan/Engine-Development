load('impulse_turbine_efficiency.mat');
% gear_teeth = [30; 35; 36; 40; 48; 50; 60; 70; 90]; % unitless - module 1 keyed shaft from KHK
% gear_teeth = [20,22,24,25,26,28,30,32,36,40,48,50,52,60,80,90,100]; % unitless - module 1.5 keyed shaft from KHK

%% Rotor Dimensions
% One stage impulse rotor - does not work for multistage and/or velocity driven turbines

turbine_params_initial = [3, 3*0.0254, 400*6894.76];

calculate_cost(turbine_params_initial, impulse_turbine_efficiency, shaft_speed, shaft_power, OF_gg, p_amb, T_amb, proof)

lb = [1, 0.5*0.0254, 10*6894.76]; % Lower bounds for gear_ratio, r_pitchline, p_gg
ub = [100/18, 12*0.0254, 500*6894.76]; % Upper bounds for gear_ratio, r_pitchline, p_gg

options = optimset('Display', 'iter');
cost_function = @(turbine_params) calculate_cost(turbine_params, impulse_turbine_efficiency, shaft_speed, shaft_power, OF_gg, p_amb, T_amb, proof)
[turbine_params_optimal] = fmincon(cost_function, turbine_params_initial, [], [], [], [], lb, ub, [], options);

% Final Optimized Parameters
gear_ratio = turbine_params_optimal(1) % unitless
r_pitchline = turbine_params_optimal(2) % m
r_pitchline/0.0254
p_gg = turbine_params_optimal(3) % Pa
p_gg/6894.76
[mdot_gg, v_spouting, v_pitchline, Tt_gg] = calculate_rotor(gear_ratio, r_pitchline, p_gg, impulse_turbine_efficiency, shaft_speed, shaft_power, OF_gg, p_amb, T_amb, proof)
mdot_gg/mdot_cc
nozzle_angle = asin(v_pitchline/v_spouting); % rad
nozzle_angle*180/pi

% turbine_thrust = mdot_gg*v_spouting(2)

%% Rotor Blades
% blade_number_rotor = 64; % unitless - total number of blades on the rotor
% blade_gap_rotor = 1/4*0.0254; % m
% blade_length_rotor = sqrt(4*pitchline_radius*mdot_gg/(blade_number_rotor*blade_gap_rotor*dens_gas_spouted*v_spouting(2))); % m - blades are unshrouded

%% Helper Functions

function cost = calculate_cost(params, impulse_turbine_efficiency, shaft_speed, shaft_power, OF_gg, p_amb, T_amb, proof)
    gear_ratio = params(1);
    r_pitchline = params(2);
    p_gg = params(3);

    mdot_gg = calculate_rotor(gear_ratio, r_pitchline, p_gg, impulse_turbine_efficiency, shaft_speed, shaft_power, OF_gg, p_amb, T_amb, proof);

    cost = 1E12*mdot_gg + 1E12*r_pitchline;
end

function [mdot_gg, v_spouting, v_pitchline, Tt_gg] = calculate_rotor(gear_ratio, r_pitchline, p_gg, impulse_turbine_efficiency, shaft_speed, shaft_power, OF_gg, p_amb, T_amb, proof)
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
    cea_gg.output.eql.mach(3); % unitless - mach number for perfectly expanded flow
    Tt_gg = cea_gg.output.eql.temperature(1); % K - stagnation / gg chamber pressure
    gamma_gg = cea_gg.output.eql.gamma(3); % unitless - ratio of specific heats
    mach_spouting = sqrt(2/(gamma_gg-1)*((p_amb/p_gg)^((1-gamma_gg)/gamma_gg)-1)); % unitless
    T_spouting = Tt_gg*(p_amb/p_gg)^((gamma_gg-1)/gamma_gg); % K - static temperature for perfectly expanded flow
    mw_gg = 0.001*cea_gg.output.eql.mw(3); % kg/mol - molecular weight (fully combusted)
    R_gg = 8.3145 / mw_gg; % J/mol*K - specific gas constant
    sonic_spouting = sqrt(gamma_gg*R_gg*T_spouting); % m/s - speed of sound at nozzle exit plane
    
    v_spouting = mach_spouting*sonic_spouting; % m/s    
    v_pitchline = turbine_shaft_speed*r_pitchline; % m/s
    
    isentropic_v_ratio = v_pitchline/norm(v_spouting); % unitless
    stage_efficiency = interp1(impulse_turbine_efficiency(:,1), impulse_turbine_efficiency(:,2), isentropic_v_ratio); % unitless
    turbine_power = shaft_power/stage_efficiency; % W
    
    mdot_gg = turbine_power/(2*v_pitchline^2); % kg/s
end


%% Estimations
% degR = 0; % turbine degree of reaction; 0 for impulse design
% 
% turbine_efficiency = 0.8; % Guess
% power_turbine = shaft_power / turbine_efficiency; % W
% 
% power_rotor = power_turbine / turbine_stage_number; % W
% specific_power_rotor = power_rotor / mdot_gg; % J/kg
% flow_coeff_rotor = 0.7; % Guess
% velocity_rotor = sqrt(specific_power_rotor / (1-flow_coeff_rotor^2)); % m/s
% diameter_rotor = 2*velocity_rotor/shaft_speed; % m
% vdot_gg = mdot_gg / dens_gas(end); % m3/s - using gas density from CC not GG, just for estimation purposes
% blade_height_rotor = vdot_gg / (pi*diameter_rotor*velocity_rotor*flow_coeff_rotor); % m

