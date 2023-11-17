%% Requirements
p_out_fuel = p_cool(1) + 50*6894.76; % Pa - outlet pressure plus 50 psi of margin for plumbing losses and various inefficiencies
dp_pre_regen = 50*6894.76; % Pa - TOTAL GUESS RIGHT NOW - pressure loss across lines between pump outlet and regen inlet
p_out_fuel = p_cool(1) + dp_pre_regen; % Pa - outlet pressure
p_in_fuel = p_amb; % Pa - inlet pressure
p_out_ox = pc; % Pa - outlet pressure
p_in_ox = p_amb; % Pa - inlet pressure

head_fuel = (p_out_fuel - p_in_fuel) / (density_fuel*g); % m

vdot_fuel = mdot_fuel_cc/density_fuel; % m3/s - Fuel Volumetric Flow Rate
vdot_ox = mdot_ox_cc/density_ox; % m3/s - Oxidizer Mass Flow Rate


%% Establishing the Pump Configuration
vapor_pressure_fuel = py.CoolProp.CoolProp.PropsSI("P", "T", T_amb, "Q", 0, "Ethanol");
NPSH_fuel = p_amb/(density_fuel*g) - vapor_pressure_fuel/(density_fuel*g); % Net positive suction head. Should be updated with plumbing water column and steady state inlet velocity
specific_speed_fuel = shaft_speed*sqrt(vdot_fuel) / (g*head_fuel)^0.75;
suction_specific_speed_fuel = shaft_speed * sqrt(vdot_fuel) / (NPSH_fuel*g)^0.75;
suction_specific_speed_imperial_fuel = suction_specific_speed_fuel*2733.016;


if specific_speed_fuel < 0.1 || specific_speed_fuel > 0.6
    fprintf("Pump Specific Speed is %g. Recommended is 0.1 to 0.6 for centrifugal geometries.", specific_speed_fuel)
end

%% Sizing the Pump
specific_speed_fuel = shaft_speed*sqrt(vdot_fuel)/(head_fuel*g)^0.75;
outlet_flow_coeff = 0.175*sqrt(specific_speed_fuel); % Phi_i in pump handbook

head_coeff_fuel = 0.4 / specific_speed_fuel^0.25; % Omega in pump handbook. Valid for specific speeds < 1 (which you should be anyways; centrifugal pumps should be between 0.1 and 0.6)

blockage = 0.85; % 1 is completely open, pump handbook says 0.85 is typical?


%% Subfunctions

impeller
blades
volute

%% Plotting

% Impeller
figure(1)
line(shroud_curve(:,1), shroud_curve(:,2))
hold on
plot(control_points(:,1),control_points(:,2),'o','color','r')
axis equal
plot(impeller_curve(:,1), impeller_curve(:,2))

% Blades
figure(2);
hold on;
delta_angle = 2 * pi / blade_number; % Calculate the angle to rotate each blade
for i = 0:(blade_number-1)
    rotation_matrix = [cos(i * delta_angle), -sin(i * delta_angle); sin(i * delta_angle), cos(i * delta_angle)];
    rotated_curve = blade_curve * rotation_matrix';
    plot(rotated_curve(:, 1), rotated_curve(:, 2), 'LineWidth', 2); 
    plot(NaN, NaN); % Prevent connection between different blades
end
hold off;
title('Impeller Blades');
axis equal;
grid on;
line([0 0], ylim);  %x-axis
line(xlim, [0 0]);  %y-axis


%% Requirements for Turbine
fuel_pump_efficiency = 0.5; % VERY ARBITRARY PLACEHOLDER
fuel_pump_power = vdot_fuel * head_fuel * density_fuel * g / fuel_pump_efficiency;
ox_pump_power = fuel_pump_power; % PLACEHOLDER

shaft_power = fuel_pump_power + ox_pump_power;






