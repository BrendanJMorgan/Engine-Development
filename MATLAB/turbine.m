%% Requirements

power_turbine = power_fuel + power_ox; % W
specific_power_turbine = power_turbine / mdot_gg; % W/kg
torque_turbine = power_turbine / v_ang; % N-m
% v_inlet = M_exit_gg*c_exit_gg; % m/s
turbine_p_ratio = p_amb/p_gg; % Pa - turbine pressure ratio

%%

turbine_eff = specific_power_turbine / ( T_gg * cp_gg * (1-turbine_p_ratio^((gamma_gg-1)/gamma_gg)) ); % Required turbine efficiency

format short g
T_gg
turbine_eff
if turbine_eff >= 1
    error("Gas generator insufficient to drive turbopump");
end