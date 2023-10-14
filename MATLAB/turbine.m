%% Requirements

power_turbine = 8000; % placeholder power_fuel + power_ox; % W
specific_power_turbine = power_turbine / mdot_gg; % W/kg
torque_turbine = power_turbine / v_ang; % N-m
% v_inlet = M_exit_gg*c_exit_gg; % m/s
turbine_p_ratio = p_amb/p_gg; % Pa - turbine pressure ratio

% turbine_eff = specific_power_turbine / ( T_gg * cp_gg * (1-turbine_p_ratio^((gamma_gg-1)/gamma_gg)) ); % Required turbine efficiency
% if turbine_eff >= 1
%     error("Gas generator insufficient to drive turbopump");
% end

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


