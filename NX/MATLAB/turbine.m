%% Requirements

power_turbine = power_fuel + power_ox; % W

torque_turbine = power_turbine / ang_vel; % W

v_inlet = M_exit_gg*c_exit_gg; % m/s


%%

turbine_eff = 0.2;

pt_out = p_amb