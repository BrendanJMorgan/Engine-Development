% Lox Pump, Fuel Pump (each with impeller/blades/volute/inducer), Turbine, Gas Generator

%% Total Mass Flow Rate GUESSES
mdot_gg = mdot_cc*gg_fraction_guess/(1-gg_fraction_guess); % kg/s - Propellant mass flow rate into gas generator
mdot_fuel_gg = mdot_gg*(1/(1+OF_gg)); % kg/s - Fuel Mass Flow Rate
mdot_ox_gg = mdot_gg*(OF_gg/(1+OF_gg)); % kg/s - Oxidizer Mass Flow Rate
mdot_fuel_total = mdot_fuel_cc+mdot_fuel_gg; % kg/s
mdot_ox_total = mdot_ox_cc+mdot_ox_gg; % kg/s
% actual mdot_gg is not solved until turbine is called

%% Liquid Oxygen Pump
density_pump = PropsSI('D','Q',0,'P',p_amb,'Oxygen'); % kg/m3 - lox at boiling
p_out = p_cc*(1+cc_stiffness) + 10*6894.76; % Pa - injector inlet pressure plus 10 psi of margin for plumbing losses
p_in = p_amb; % Pa - inlet pressure; assumes zero dp across inducer for conservatism
T_in = PropsSI("T", "P", p_amb, "Q", 0, "Oxygen");
mdot_pump = mdot_ox_total; % kg/s - Oxidizer Mass Flow Rate
vapor_pressure_pump = PropsSI("P", "T", T_in, "Q", 0, "Oxygen");
clock = 1; % Counterclockwise
hub_tip_ratio_inducer = 0.4; % unitless - between 0.2 and 0.4 for rear drive and between 0.5 and 0.6 for front drive

pump

shaft_power_pump = power_pump; % W
shaft_power_pump_theory = power_pump_theory; % W

    % Write results to ox specific variables
    p_in_ox = p_in;
    p_out_ox = p_out;

    shroud_curve_ox = shroud_curve;
    impeller_curve_ox = impeller_curve; 
    blade_curve_ox = blade_curve;
    blade_number_ox = blade_number;
    
    volute_curve_ox = volute_curve;
    h_volute_ox = h_volute;
    t_tongue_ox = t_tongue;
    r_volute_outlet_ox = r_volute_outlet;
    OD_volute_outlet_ox = OD_volute_outlet;
    
    h_min_ox = h_min_inducer;
    r_tip_ox = r_tip_inducer;
    r_hub_base_ox = r_hub_base;
    r_hub_top_ox = r_hub_top;
    blade_lead_ox = blade_lead_inducer;
    clearance_ox = clearance_axial_inducer;

    writematrix([impeller_curve/0.0254, zeros(length(impeller_curve(:,1)),1)], 'Curves/ox_impeller_curve_inches.txt', 'Delimiter', ',')  
    writematrix([shroud_curve/0.0254, zeros(length(shroud_curve(:,1)),1)], 'Curves/ox_shroud_curve_inches.txt', 'Delimiter', ',')
    writematrix([blade_curve(:,1)/0.0254, zeros(length(impeller_curve(:,1)),1), -blade_curve(:,2)/0.0254], 'Curves/ox_blade_curve_inches.txt', 'Delimiter', ',')
    writematrix([volute_curve_ox(:,1)/0.0254, volute_curve_ox(:,2)/0.0254, zeros(length(volute_curve_ox(:,1)),1)], 'Curves/ox_volute_curve_inches.txt', 'Delimiter', ',')
    
    data_labels = {'h_volute_ox';'t_tongue_ox';'r_volute_outlet_ox';'OD_volute_outlet_ox'};
    data = table(data_labels, {h_volute/0.0254;t_tongue/0.0254;r_volute_outlet/0.0254;OD_volute_outlet/0.0254});
    writetable(data, 'Output Parameters/ox_volute_inches.txt', 'WriteVariableNames',0);
    
    data_labels = {'h_min_ox';'r_tip_ox';'r_hub_base_ox';'r_hub_top_ox';'blade_lead_ox';'clearance_ox'};
    data = table(data_labels, {h_min_inducer/0.0254;r_tip_inducer/0.0254;r_hub_base/0.0254;r_hub_top/0.0254;blade_lead_inducer/0.0254;clearance_axial_inducer/0.0254});
    writetable(data, 'Output Parameters/ox_inducer_inches.txt', 'WriteVariableNames',0);

%% Fuel Pump (0.75/0.25 Ethanol/Water)
density_pump = PropsSI('D','T',T_amb,'P',p_amb,['Ethanol[',num2str(proof),']&Water[',num2str(1-proof),']']); % kg/m3 - ethanol at STP
p_out = p_cool(1) + 10*6894.76; % Pa - regen inlet pressure plus 10 psi of margin for plumbing losses
p_in = p_amb; % Pa - inlet pressure; assumes zero dp across inducer for conservatism
T_in = T_amb; % K
mdot_pump = mdot_fuel_total; % kg/s - Fuel Mass Flow Rate
vapor_pressure_pump = PropsSI("P", "T", T_in, "Q", 0, ['Ethanol[',num2str(proof),']&Water[',num2str(1-proof),']']);
clock = -1;
hub_tip_ratio_inducer = 0.39; % unitless - between 0.2 and 0.4 for rear drive and between 0.5 and 0.6 for front drive

pump

shaft_power_pump = shaft_power_pump + power_pump; % W
shaft_power_pump_theory = shaft_power_pump_theory + power_pump_theory; % W

    % Write results to fuel specific variables
    p_in_fuel = p_in;
    p_out_fuel = p_out;

    shroud_curve_fuel = shroud_curve;
    impeller_curve_fuel = impeller_curve; 
    blade_curve_fuel = blade_curve;
    blade_number_fuel = blade_number;
    
    volute_curve_fuel = volute_curve;
    h_volute_fuel = h_volute;
    t_tongue_fuel = t_tongue;
    r_volute_outlet_fuel = r_volute_outlet;
    OD_volute_outlet_fuel = OD_volute_outlet;
    
    h_min_fuel = h_min_inducer;
    r_tip_fuel = r_tip_inducer;
    r_hub_base_fuel = r_hub_base;
    r_hub_top_fuel = r_hub_top;
    blade_lead_fuel = blade_lead_inducer;
    clearance_fuel = clearance_axial_inducer;
    
    
    writematrix([impeller_curve/0.0254, zeros(length(impeller_curve(:,1)),1)], 'Curves/fuel_impeller_curve_inches.txt', 'Delimiter', ',')  
    writematrix([shroud_curve/0.0254, zeros(length(shroud_curve(:,1)),1)], 'Curves/fuel_shroud_curve_inches.txt', 'Delimiter', ',')
    writematrix([blade_curve(:,1)/0.0254, zeros(length(impeller_curve(:,1)),1), -blade_curve(:,2)/0.0254], 'Curves/fuel_blade_curve_inches.txt', 'Delimiter', ',')
    writematrix([volute_curve_fuel(:,1)/0.0254, volute_curve_fuel(:,2)/0.0254, zeros(length(volute_curve_fuel(:,1)),1)], 'Curves/fuel_volute_curve_inches.txt', 'Delimiter', ',')
    
    data_labels = {'h_volute_fuel';'t_tongue_fuel';'r_volute_outlet_fuel';'OD_volute_outlet_fuel'};
    data = table(data_labels, {h_volute/0.0254;t_tongue/0.0254;r_volute_outlet/0.0254;OD_volute_outlet/0.0254});
    writetable(data, 'Output Parameters/fuel_volute_inches.txt', 'WriteVariableNames',0);
    
    data_labels = {'h_min_fuel';'r_tip_fuel';'r_hub_base_fuel';'r_hub_top_fuel';'blade_lead_fuel';'clearance_fuel'};
    data = table(data_labels, {h_min_inducer/0.0254;r_tip_inducer/0.0254;r_hub_base/0.0254;r_hub_top/0.0254;blade_lead_inducer/0.0254;clearance_axial_inducer/0.0254});
    writetable(data, 'Output Parameters/fuel_inducer_inches.txt', 'WriteVariableNames',0);

%% Gas Generator
gas_generator2

%% Turbine
shaft_power_turbine = shaft_power_pump/gear_efficiency; % unitless - the power through the turbine SHAFT - not the turbine itself, which will be subject to a stage inefficiency

turbine2

    % Write results to text file
    data_labels = {'Base Radius (in)'; 'Pitchline Radius (in)'; 'Tip Radius (in)'; 'Blade Width (in)'; 'Blade Depth (in)';
        'Blade Angle (deg)'; 'Blade Opening (in)'; 'Blade Number'; 'Nozzle Angle (deg)'; 'Nozzle Diameter (in)'; 'Nozzle Number'};
    data = table(data_labels, {r_rotor_base/0.0254; r_pitchline_rotor/0.0254; r_rotor_tip/0.0254; blade_width_rotor/0.0254; blade_length_rotor/0.0254; ...
        incidence_angle_blade*180/pi; blade_opening_rotor/0.0254; blade_number_rotor; angle_nozzle*180/pi; d_throat_nozzle/0.0254; nozzle_number});
    writetable(data, 'Output Parameters/turbine_inches.txt', 'WriteVariableNames',0);

% Efficiency
transmission_efficiency = shaft_power_pump_theory / shaft_power_turbine; % unitless
isp_real = isp_ideal*mdot_cc/mdot_total;

% Shaft
shaft_torque_pump = shaft_power_pump/shaft_speed; % N*m
shaft_polar_moment = pi/2*r_shaft^4; % m4
shaft_stress = shaft_torque_pump*r_shaft/shaft_polar_moment; % Pa
shaft_torsion_FSu = shaft_shear_strength / shaft_stress; % unitless
shaft_torsion_MSu = shaft_torsion_FSu / FSu_design - 1; % unitless





