%% Turbopump

% LOx Pump
density_ox = PropsSI('D','Q',0,'P',p_amb,'Oxygen'); % kg/m3 - lox at boiling
p_out = p_cc + 25*6894.76; % Pa - GUESS RIGHT NOW
dp_inducer = 25*6894.76; % Pa - PLACEHOLDER
p_in = p_amb + dp_inducer; % Pa - inlet pressure
T_in = PropsSI("T", "P", p_amb, "Q", 0, "Oxygen");
head = (p_out - p_in) / (density_ox*g); % m
vdot = mdot_ox_total/density_ox; % m3/s - Oxidizer Mass Flow Rate
density = density_ox; % kg/m3
vapor_pressure = PropsSI("P", "T", T_in, "Q", 0, "Oxygen");
clock = 1; % Counterclockwise
front_drive = false;

pump
% inducer

shaft_power = pump_power;

writematrix([impeller_curve/0.0254, zeros(length(impeller_curve(:,1)),1)], 'Curves/ox_impeller_curve_inches.txt', 'Delimiter', ',')  
writematrix([shroud_curve/0.0254, zeros(length(shroud_curve(:,1)),1)], 'Curves/ox_shroud_curve_inches.txt', 'Delimiter', ',')
writematrix([blade_curve(:,1)/0.0254, zeros(length(impeller_curve(:,1)),1), -blade_curve(:,2)/0.0254], 'Curves/ox_blade_curve_inches.txt', 'Delimiter', ',')
writematrix([r_volute.*cos(theta_volute)/0.0254, r_volute.*sin(theta_volute)/0.0254, zeros(length(theta_volute(:,1)),1)], 'Curves/ox_volute_curve_inches.txt', 'Delimiter', ',')

% Fuel Pump
density_fuel = PropsSI('D','T',T_amb,'P',p_amb,['Ethanol[',num2str(proof),']&Water[',num2str(1-proof),']']); % kg/m3 - ethanol at STP
p_out = p_cool(1) + 50*6894.76; % Pa - outlet pressure plus 50 psi of margin for plumbing losses and various inefficiencies, TOTAL GUESS RIGHT NOW
dp_inducer = 25*6894.76; % Pa - PLACEHOLDER
p_in = p_amb + dp_inducer; % Pa - inlet pressure
T_in = T_amb;
head = (p_out - p_in) / (density_fuel*g); % m
vdot = mdot_ox_total/density_ox; % m3/s - Oxidizer Mass Flow Rate
density = density_fuel; % kg/m3
vapor_pressure = PropsSI("P", "T", T_in, "Q", 0, ['Ethanol[',num2str(proof),']&Water[',num2str(1-proof),']']);
clock = -1;
front_drive = true;

pump

shaft_power = shaft_power + pump_power;

writematrix([impeller_curve/0.0254, zeros(length(impeller_curve(:,1)),1)], 'Curves/fuel_impeller_curve_inches.txt', 'Delimiter', ',')  
writematrix([shroud_curve/0.0254, zeros(length(shroud_curve(:,1)),1)], 'Curves/fuel_shroud_curve_inches.txt', 'Delimiter', ',')
writematrix([blade_curve(:,1)/0.0254, zeros(length(impeller_curve(:,1)),1), -blade_curve(:,2)/0.0254], 'Curves/fuel_blade_curve_inches.txt', 'Delimiter', ',')
writematrix([r_volute.*cos(theta_volute)/0.0254, r_volute.*sin(theta_volute)/0.0254, zeros(length(theta_volute(:,1)),1)], 'Curves/fuel_volute_curve_inches.txt', 'Delimiter', ',')

% Gas Generator
% gas_generator

% Turbine
turbine4

% Shaft
shaft_torque = shaft_power/shaft_speed; % N*m
shaft_polar_moment = pi/2*r_shaft^4; % m4
shaft_stress = shaft_torque*r_shaft/shaft_polar_moment; % Pa
shaft_torsion_FSu = shaft_shear_strength / shaft_stress
shaft_torsion_MSu = shaft_torsion_FSu / FSu_design - 1 





