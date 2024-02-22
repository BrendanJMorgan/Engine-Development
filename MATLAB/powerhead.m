%% Total Mass Flow Rate GUESSES
mdot_gg = mdot_cc*gg_fraction_guess/(1-gg_fraction_guess); % kg/s - Propellant mass flow rate into gas generator
mdot_fuel_gg = mdot_gg*(1/(1+OF_gg)); % kg/s - Fuel Mass Flow Rate
mdot_ox_gg = mdot_gg*(OF_gg/(1+OF_gg)); % kg/s - Oxidizer Mass Flow Rate
mdot_fuel_total = mdot_fuel_cc+mdot_fuel_gg; % kg/s
mdot_ox_total = mdot_ox_cc+mdot_ox_gg; % kg/s
% actual mdot_gg is not solved until turbine is called

%% LOx Pump
density_pump = PropsSI('D','Q',0,'P',p_amb,'Oxygen'); % kg/m3 - lox at boiling
p_out = p_cc + 25*6894.76; % Pa - GUESS RIGHT NOW
dp_inducer = 25*6894.76; % Pa - PLACEHOLDER
p_in = p_amb + dp_inducer; % Pa - inlet pressure
T_in = PropsSI("T", "P", p_amb, "Q", 0, "Oxygen");
mdot_pump = mdot_ox_total; % kg/s - Oxidizer Mass Flow Rate
vapor_pressure_pump = PropsSI("P", "T", T_in, "Q", 0, "Oxygen");
clock = 1; % Counterclockwise

pump
inducer

shaft_power = pump_power; % W
shaft_power_theory = pump_power_theory; % W

% Write Results
shroud_curve_ox = shroud_curve;
impeller_curve_ox = impeller_curve; 
blade_curve_ox = blade_curve;
blade_number_ox = blade_number;
volute_curve_ox = volute_curve;
t_tongue_ox = t_tongue;
h_volute_ox = h_volute;

writematrix([impeller_curve/0.0254, zeros(length(impeller_curve(:,1)),1)], 'Curves/ox_impeller_curve_inches.txt', 'Delimiter', ',')  
writematrix([shroud_curve/0.0254, zeros(length(shroud_curve(:,1)),1)], 'Curves/ox_shroud_curve_inches.txt', 'Delimiter', ',')
writematrix([blade_curve(:,1)/0.0254, zeros(length(impeller_curve(:,1)),1), -blade_curve(:,2)/0.0254], 'Curves/ox_blade_curve_inches.txt', 'Delimiter', ',')
writematrix([volute_curve_ox(:,1)/0.0254, volute_curve_ox(:,2)/0.0254, zeros(length(volute_curve_ox(:,1)),1)], 'Curves/ox_volute_curve_inches.txt', 'Delimiter', ',')

%% Fuel Pump
density_pump = PropsSI('D','T',T_amb,'P',p_amb,['Ethanol[',num2str(proof),']&Water[',num2str(1-proof),']']); % kg/m3 - ethanol at STP
p_out = p_cool(1) + 50*6894.76; % Pa - outlet pressure plus 50 psi of margin for plumbing losses and various inefficiencies, TOTAL GUESS RIGHT NOW
dp_inducer = 25*6894.76; % Pa - PLACEHOLDER
p_in = p_amb + dp_inducer; % Pa - inlet pressure
T_in = T_amb; % K
mdot_pump = mdot_fuel_total; % kg/s - Fuel Mass Flow Rate
vapor_pressure_pump = PropsSI("P", "T", T_in, "Q", 0, ['Ethanol[',num2str(proof),']&Water[',num2str(1-proof),']']);
clock = -1;

pump

shaft_power = shaft_power + pump_power; % W
shaft_power_theory = shaft_power_theory + pump_power_theory; % W

% Write Results
shroud_curve_fuel = shroud_curve;
impeller_curve_fuel = impeller_curve; 
blade_curve_fuel = blade_curve;
blade_number_fuel = blade_number; 
volute_curve_fuel = volute_curve;
t_tongue_fuel = t_tongue;
h_volute_fuel = h_volute; 

writematrix([impeller_curve/0.0254, zeros(length(impeller_curve(:,1)),1)], 'Curves/fuel_impeller_curve_inches.txt', 'Delimiter', ',')  
writematrix([shroud_curve/0.0254, zeros(length(shroud_curve(:,1)),1)], 'Curves/fuel_shroud_curve_inches.txt', 'Delimiter', ',')
writematrix([blade_curve(:,1)/0.0254, zeros(length(impeller_curve(:,1)),1), -blade_curve(:,2)/0.0254], 'Curves/fuel_blade_curve_inches.txt', 'Delimiter', ',')
writematrix([volute_curve_fuel(:,1)/0.0254, volute_curve_fuel(:,2)/0.0254, zeros(length(volute_curve_fuel(:,1)),1)], 'Curves/fuel_volute_curve_inches.txt', 'Delimiter', ',')

% Gas Generator
gas_generator2

%% Turbine
turbine_shaft_power = shaft_power/gear_efficiency; % unitless
turbine6

% Efficiency
transmission_efficiency = shaft_power_theory / turbine_shaft_power; % unitless
isp_real = isp_ideal*mdot_cc/mdot_total;

% Shaft
shaft_torque = shaft_power/shaft_speed; % N*m
shaft_polar_moment = pi/2*r_shaft^4; % m4
shaft_stress = shaft_torque*r_shaft/shaft_polar_moment; % Pa
shaft_torsion_FSu = shaft_shear_strength / shaft_stress; % unitless
shaft_torsion_MSu = shaft_torsion_FSu / FSu_design - 1; % unitless





