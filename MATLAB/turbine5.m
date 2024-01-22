load('impulse_turbine_efficiency.mat');

gear_ratio = 100/18; % unitless - higher makes for a smaller, faster turbine
p_gg = 10000*6894.76; % Pa - gas generator chamber pressure

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
Tt_gg = cea_gg.output.eql.temperature(1); % K - stagnation, gg chamber temperature
gamma_gg = cea_gg.output.eql.gamma(3); % unitless - ratio of specific heats
mach_spouting = sqrt(2/(gamma_gg-1)*((p_amb/p_gg)^((1-gamma_gg)/gamma_gg)-1)); % unitless
T_spouting = Tt_gg*(p_amb/p_gg)^((gamma_gg-1)/gamma_gg); % K - static temperature for perfectly expanded flow
mw_gg = 0.001*cea_gg.output.eql.mw(3); % kg/mol - molecular weight (fully combusted)
R_gg = 8.3145 / mw_gg; % J/mol*K - specific gas constant
sonic_spouting = sqrt(gamma_gg*R_gg*T_spouting); % m/s - speed of sound at nozzle exit plane

v_spouting_norm = mach_spouting*sonic_spouting; % m/s - speed at which gas leaves the nozzles

%% Rotor Dimensions
% One stage impulse rotor - does not work for multistage and/or velocity driven turbines
objective_function = @(r_pitchline) calculate_rotor(r_pitchline, turbine_shaft_speed,impulse_turbine_efficiency, shaft_power, v_spouting_norm);
[r_pitchline, mdot_gg] = fminunc(objective_function, 2*0.0254);

% Final Optimized Parameters
[mdot_gg, v_pitchline] = calculate_rotor(r_pitchline, turbine_shaft_speed,impulse_turbine_efficiency, shaft_power, v_spouting_norm);
nozzle_angle = asin(2*v_pitchline/v_spouting_norm); % rad

if ~(nozzle_angle < 70*pi/180)
    v_pitchline
    v_spouting_norm
    error("Insufficient pressure ratio across nozzle plate")
end
 
% turbine_thrust = mdot_gg*v_spouting(2)

%% Rotor Blades
% blade_number_rotor = 64; % unitless - total number of blades on the rotor
% blade_gap_rotor = 1/4*0.0254; % m
% blade_length_rotor = sqrt(4*pitchline_radius*mdot_gg/(blade_number_rotor*blade_gap_rotor*dens_gas_spouted*v_spouting(2))); % m - blades are unshrouded

%% Helper Functions

function [mdot_gg, v_pitchline] = calculate_rotor(r_pitchline, turbine_shaft_speed,impulse_turbine_efficiency, shaft_power, v_spouting_norm)
    v_pitchline = turbine_shaft_speed*r_pitchline; % m/s
    
    isentropic_v_ratio = v_pitchline/v_spouting_norm; % unitless
    stage_efficiency = interp1(impulse_turbine_efficiency(:,1), impulse_turbine_efficiency(:,2), isentropic_v_ratio, 'spline'); % unitless
    if stage_efficiency < 0.0001
        stage_efficiency = 0.0001; % helps the solver converge
    end
    turbine_power = shaft_power/stage_efficiency; % W
    
    mdot_gg = turbine_power/(2*v_pitchline^2) % kg/s

    if 2*v_pitchline < v_spouting_norm
        mdot_gg = 100000000*mdot_gg
    end
end

