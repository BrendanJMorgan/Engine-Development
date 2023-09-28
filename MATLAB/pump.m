%% Requirements
p_out_fuel = p_cool(1); % Pa - outlet pressure
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
