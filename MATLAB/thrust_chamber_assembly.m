combustion_chamber
cc_geometry
cc_gas_flow
coolant_flow

%% Often commented out because of long runtime. thermal_balance.m does not affect any inputs into powerhead.m
% T_wall_hot = 800*ones(1,length(x)); % K - initial GUESS for the hot wall temperatures
% thermal_balance
% thermal_balance % intentionally run twice (greater convergence)