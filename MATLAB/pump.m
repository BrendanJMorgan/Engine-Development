%% Establishing the Pump Configuration
dr = 0.0001;

%% Sizing the Pump


specific_speed = shaft_speed*sqrt(vdot_pump)/(head_pump*g)^0.75;
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


% Other
d_inlet_pump = 3/8*0.0254; % m
v_inlet_pump = vdot_pump / (pi/4*d_inlet_pump^2) % m/s
cavitation_inlet = (p_amb - vapor_pressure_pump) / (0.5*density_pump*v_inlet_pump^2) % unitless

d_outlet_pump = 3/8*0.0254; % m
v_outlet_pump = vdot_pump / (pi/4*d_outlet_pump^2); % m/s
cavitation_outlet = (p_out - vapor_pressure_pump) / (0.5*density_pump*v_outlet_pump^2); % unitless


%% Requirements for Turbine
pump_efficiency = 0.5; % VERY ARBITRARY PLACEHOLDER
pump_power = vdot_pump * head_pump * density_pump * g / pump_efficiency;








