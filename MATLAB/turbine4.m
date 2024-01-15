turbine_shaft_speed = shaft_speed*gear_ratio % rad/s


%% Rotor
rotor_torque = shaft_power/turbine_shaft_speed % N*m
flow_coeff_rotor = 0.7; % Guess
v_pitchline1 = sqrt(shaft_power / (mdot_gg*(1-flow_coeff_rotor^2))) % m/s
v_pitchline2 = sqrt(shaft_power/(2*mdot_gg*(1-reaction_degree)))
% Flow Coefficient?
r_pitchline = v_pitchline/turbine_shaft_speed % m
r_pitchline/0.0254
nozzle_angle = 30*pi/180 % rad - empirical recommendation by Farokhi is 70 deg, page 691

v_nozzle_abs = turbine_shaft_speed*r_pitchline * [1, acot(nozzle_angle)] % [m/s, m/s]
mdot_gg/mdot_cc

%%
% gas_generator


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

