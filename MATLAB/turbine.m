%% Requirements

% specific_power_turbine = power_turbine / mdot_gg; % W/kg
% torque_turbine = power_turbine / v_ang; % N-m
% % v_inlet = M_exit_gg*c_exit_gg; % m/s
% turbine_p_ratio = p_amb/p_gg; % Pa - turbine pressure ratio

%% Estimations
turbine_efficiency = 0.8; % Guess
power_turbine = shaft_power / turbine_efficiency % W

power_rotor = power_turbine / turbine_stage_number; % W
specific_power_rotor = power_rotor / mdot_gg % J/kg
flow_coeff_rotor = 0.7; % Guess
velocity_rotor = sqrt(specific_power_rotor / (1-flow_coeff_rotor^2)) % m/s
diameter_rotor = 2*velocity_rotor/shaft_speed % m
vdot_gg = mdot_gg / dens_gas(end); % m3/s - using gas density from CC not GG, just for estimation purposes
blade_height_rotor = vdot_gg / (pi*diameter_rotor*velocity_rotor*flow_coeff_rotor) % m


diameter_rotor / 0.0254
blade_height_rotor / 0.0254
%% 1) Selection of rotational speed and annulus dimensions
% Rotational speed requirements are driven by the pumps


%% 2) Selection of the number of stages


%% 3) Calculation of airflow angles for each stage at the mean radius 


%% 4) Calculation of airflow angle variations from the hub to tip for each stage 


%% 5) Selection of blade material


%% 6) Selection of blading using experimental cascade data 


%% 7) Selection of turbine cooling, if needed 


%% 8) Verification of turbine efficiency based on cascade loss data 


%% 9) Prediction of off-design performance 


%% 10) Rig testing of design


