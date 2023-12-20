%% Requirements
p_out_fuel = p_cool(1) + 50*6894.76; % Pa - outlet pressure plus 50 psi of margin for plumbing losses and various inefficiencies
dp_pre_regen = 50*6894.76; % Pa - TOTAL GUESS RIGHT NOW - pressure loss across lines between pump outlet and regen inlet
p_out_fuel = p_cool(1) + dp_pre_regen; % Pa - outlet pressure
p_in_fuel = p_amb; % Pa - inlet pressure
p_out_ox = p_cc; % Pa - outlet pressure
p_in_ox = p_amb; % Pa - inlet pressure

head_fuel = (p_out_fuel - p_in_fuel) / (density_fuel*g); % m

vdot_fuel = mdot_fuel_cc/density_fuel; % m3/s - Fuel Volumetric Flow Rate
vdot_ox = mdot_ox_cc/density_ox; % m3/s - Oxidizer Mass Flow Rate


%% Establishing the Pump Configuration
dr = 0.0001;

vapor_pressure_fuel = py.CoolProp.CoolProp.PropsSI("P", "T", T_amb, "Q", 0, "Ethanol");
NPSH_fuel = p_amb/(density_fuel*g) - vapor_pressure_fuel/(density_fuel*g); % Available Net positive suction head. Measure hoe close fluid is to flashing/cavitating. Should be updated with plumbing water column and steady state inlet velocity
specific_speed_fuel = shaft_speed*sqrt(vdot_fuel) / (g*head_fuel)^0.75; % Specific speed, used to categorize impeller types
suction_specific_speed_fuel = shaft_speed * sqrt(vdot_fuel) / (NPSH_fuel*g)^0.75;
suction_specific_speed_imperial_fuel = suction_specific_speed_fuel*2733.016;


if specific_speed_fuel < 0.1 || specific_speed_fuel > 0.6
    fprintf("Pump Specific Speed is %g. Recommended is 0.1 to 0.6 for centrifugal geometries.\n", specific_speed_fuel)
end

%% Sizing the Pump
specific_speed_fuel = shaft_speed*sqrt(vdot_fuel)/(head_fuel*g)^0.75;
outlet_flow_coeff = 0.175*sqrt(specific_speed_fuel); % phi_i in pump handbook
head_coeff_fuel = 0.4 / specific_speed_fuel^0.25; % psi in pump handbook. Valid for specific speeds < 1 (which you should be anyways; centrifugal pumps should be between 0.1 and 0.6)
blockage = 0.85; % 1 is completely open, pump handbook says 0.85 is typical?

%% Subfunctions

impeller
blades
volute

%% Centrifugal Integrity - this is a cursory sanity check, a full FEM should be made to verify
% andddddd it's completely fucking wrong. use flywheel stress.
% 
% r = (0:dr:r_exit)'; % m - unified radial coordinate system
% impeller_curve = [ [(0:dr:r_shaft)', impeller_height*ones(ceil(r_shaft/dr),1)] ; impeller_curve];
% 
% bronze_density = 7580; % kg/m3 - for Aluminum Bronze Alloy 63000 from Azom spec sheet
% ss_density = 7970; % kg/m3 - for stainless steel 316 average from Azom spec sheet
% 
% impeller_interpolated = interp1(impeller_curve(:,1), impeller_curve(:,2), r);
% shroud_interpolated = interp1(shroud_curve(:,1), shroud_curve(:,2), r);
% shroud_interpolated(~(shroud_interpolated>=0)) = 0;
% 
% % Shell Integrals of Revolution
% dV_impeller = 2*pi*r*dr.*impeller_interpolated;
% dm_impeller = [ss_density*dV_impeller(1:floor(r_shaft_ss/dr)); bronze_density*dV_impeller(ceil(r_shaft_ss/dr):end)];
% dV_shroud = 2*pi*r*dr.*shroud_interpolated;
% dm_shroud = bronze_density*dV_shroud;
% 
% force_centrifugal = cumsum(shaft_speed^2 .*(dm_impeller+dm_shroud).* r, 'reverse');
% stress_centrifugal = force_centrifugal ./ (2*pi*r.*(impeller_interpolated + shroud_interpolated));

%% Plotting

% Impeller
figure(1)
clf
line(shroud_curve(:,1), shroud_curve(:,2))
hold on
plot(control_points(:,1),control_points(:,2),'o','color','r')
axis equal
plot(impeller_curve(:,1), impeller_curve(:,2))
line([0 0], ylim);  %x-axis
line(xlim, [0 0]);  %y-axis
hold off

% Blades
figure(2);
clf
hold on;

delta_angle = 2 * pi / blade_number; % Calculate the angle to rotate each blade
for i = 0:(blade_number-1)
    rotation_matrix = [cos(i * delta_angle), -sin(i * delta_angle); sin(i * delta_angle), cos(i * delta_angle)];
    rotated_curve = blade_curve * rotation_matrix';
    plot(rotated_curve(:, 1), rotated_curve(:, 2), 'LineWidth', 2); 
    plot(NaN, NaN); % Prevent connection between different blades
end
plot(blade_control_points(:,1),blade_control_points(:,2),'o','color','r')
plot(r_volute_midline.*cos(theta_volute), r_volute_midline.*sin(theta_volute), '--'); % Volute Flow Midline
plot(r_volute_wall.*cos(theta_volute), r_volute_wall.*sin(theta_volute)); % Volute Wall
line([r_volute_wall(1), r_volute_wall(1)+t_tongue], [0, 0]); % Cutwater
fplot(@(t) r_shaft*sin(t), @(t) r_shaft*cos(t)); % Shaft
fplot(@(t) r_eye_inner*sin(t), @(t) r_eye_inner*cos(t)); % Eye Inlet

hold off;
title('Impeller Blades and Volute');
axis equal;
grid on;

% Structural 
% figure(3)
% plot(r, stress_centrifugal/1E6)
% xlabel("Radial Distance (m)")
% ylabel("Centrifugal Tensile Stress (MPa)")


%% Requirements for Turbine
fuel_pump_efficiency = 0.5; % VERY ARBITRARY PLACEHOLDER
fuel_pump_power = vdot_fuel * head_fuel * density_fuel * g / fuel_pump_efficiency;
ox_pump_power = fuel_pump_power; % PLACEHOLDER

shaft_power = fuel_pump_power + ox_pump_power;






