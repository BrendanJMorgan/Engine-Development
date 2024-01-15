%% Establishing the Pump Configuration
dr = 0.0001;

% specific_speed = shaft_speed*sqrt(vdot) / (g*head)^0.75;
% suction_specific_speed = shaft_speed * sqrt(vdot) / (NPSH*g)^0.75;
% suction_specific_speed_imperial = suction_specific_speed*2733.016;

%% Sizing the Pump


specific_speed = shaft_speed*sqrt(vdot)/(head*g)^0.75;
outlet_flow_coeff = 0.175*sqrt(specific_speed); % phi_i in pump handbook
head_coeff = 0.4 / specific_speed^0.25; % psi in pump handbook. Valid for specific speeds < 1 (which you should be anyways; centrifugal pumps should be between 0.1 and 0.6)
blockage = 0.85; % 1 is completely open, pump handbook says 0.85 is typical?
new_blockage = 0; % placeholder

if specific_speed < 0.1 || specific_speed > 0.6
    fprintf("Pump Specific Speed is %g. Recommended is 0.1 to 0.6 for centrifugal geometries.\n", specific_speed)
end

%% Subfunctions

impeller
while (abs(new_blockage-blockage) > 0.01)
    blockage = new_blockage;
    blades
end
volute

%% Requirements for Turbine
pump_efficiency = 0.5; % VERY ARBITRARY PLACEHOLDER
pump_power = vdot * head * density * g / pump_efficiency;








