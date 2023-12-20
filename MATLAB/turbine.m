%% Notes
% Blade Shape is generally dictated by stresses, not flow

%% Requirements

% specific_power_turbine = power_turbine / mdot_gg; % W/kg
% torque_turbine = power_turbine / v_ang; % N-m
% % v_inlet = M_exit_gg*c_exit_gg; % m/s
% turbine_p_ratio = p_amb/p_gg; % Pa - turbine pressure ratio

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
% 
% 
% diameter_rotor / 0.0254;
% blade_height_rotor / 0.0254;

%% Nozzle Plate (np)
nozzle_number = 16;
A_throat_np = mdot_gg*c_star_gg/p_gg/nozzle_number; % m2
d_throat_np = 2*sqrt(A_throat_np/pi); % m

p_back_turbine = p_amb; % Pa - assuming the backpressure is ambient, I have no idea if this is valid
M_np = sqrt( 2/(gamma_gg-1) * ((p_back_turbine/p_gg)^((1-gamma_gg)/gamma_gg) - 1) ); % Absolute mach post nozzle, pre rotor

A_exit_np = A_throat_np/M_np * ( (2+(gamma_gg-1)*M_np^2)/(gamma_gg+1) ) ^ (0.5*(gamma_gg+1)/(gamma_gg-1)); % m2
d_exit_np = 2*sqrt(A_throat_np/pi); % m - more of a virtual diameter because the nozzle is so tilted

T_exit_np = Tt_gg/(1+(gamma-1)/2*M_np^2); % K temperature at the virtual exit; reality may be more complicated

alpha_nozzle = 70*pi/180; % empirical recommendation by Farokhi, page 691
v_np = M_np*sqrt(gamma_gg*R_gg*T_exit_np)*[cos(alpha_nozzle), sin(alpha_nozzle)]

%% Rotor

% 1) Selection of rotational speed and annulus dimensions
% Rotational speed requirements are driven by the pumps


%% 2) Selection of the number of stages

% flow_coeff = v_nozzle / (shaft_speed*r_rotor);night
% 
% Tt_rotor = Tt_nozzle - power_rotor / (mdot_gg*cp_gg)
% stage_loading_coeff = g*cp_gg*(Tt_nozzle-Tt_rotor)/(shaft_speed*r_rotor)^2;

%% 3) Calculation of airflow angles for each stage at the mean radius 


%% 4) Calculation of airflow angle variations from the hub to tip for each stage 


%% 5) Selection of blade material

%% 6) Selection of blading using experimental cascade data 


%% 7) Selection of turbine cooling, if needed 


%% 8) Verification of turbine efficiency based on cascade loss data 


%% 9) Prediction of off-design performance 


%% 10) Rig testing of design


