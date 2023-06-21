%% Requirements

p_out_fuel = p_cool(1); % Pa - outlet pressure
p_in_fuel = p_amb; % Pa - inlet pressure
p_out_ox = pc; % Pa - outlet pressure
p_in_ox = p_amb; % Pa - inlet pressure

fuel_head = (p_out_fuel - p_in_fuel) / (density_fuel*g); % m

vdot_fuel = mdot_fuel_cc/density_fuel; % m3/s - Fuel Volumetric Flow Rate
vdot_ox = mdot_ox_cc/density_ox; % m3/s - Oxidizer Mass Flow Rate


%% Establishing the Pump Configuration
vapor_pressure_fuel = py.CoolProp.CoolProp.PropsSI("P", "T", T_amb, "Q", 0, "Ethanol");
NPSH_fuel = p_amb/(density_fuel*g) - vapor_pressure_fuel/(density_fuel*g); % Net positive suction head. Should be updated with plumbing water column and steady state inlet velocity
suction_specific_speed_fuel = shaft_speed * sqrt(vdot_fuel) / (NPSH_fuel*g)^0.75;
suction_specific_speed_imperial_fuel = suction_specific_speed_fuel*2733.016;

%% Sizing the Pump
specific_speed_fuel = shaft_speed*sqrt(vdot_fuel)/(fuel_head*g)^0.75;
outlet_flow_coeff = 0.175*sqrt(specific_speed_fuel); % Phi in pump handbook

exit_radius = sqrt(fuel_head*g/outlet_flow_coeff)/shaft_speed; % r2 in pump handbook

head_coeff_fuel = g*fuel_head/(shaft_speed*exit_radius)^2; % Omega in pump handbook

blockage = 0.85; % 1 is completely open, pump handbook says 0.85 is typical?
exit_width = vdot_fuel / (2*pi*shaft_speed*exit_radius^2*outlet_flow_coeff*blockage);

%% Impeller Design

% The eye
% eye_flow_coeff = 


% The exit radius r2 (or diameter D2)


% The exit width b2


% Hub and shroud profiles


% Construction of the blades


% %% Inputs
% 
% r_pump = 2*0.0254; % m - blade radius
% h_exit_pump = 1*0.0254; % m - exit gap height
% 
% pump_eff = 0.55;
% blade_angle = 30*pi/180;
% 
% 
% %%
% 
% head_fuel = (p_out-p_in)/(density_fuel*9.80665); % m - pump pressure head
% head_ox = (p_out-p_in)/(density_ox*9.80665); % m - pump pressure head
% 
% v_rad_fuel = vdot_fuel/(2*pi*r_pump*h_exit_pump); % m/s - radial velocity
% v_rad_ox = vdot_ox/(2*pi*r_pump*h_exit_pump); % m/s - radial velocity
% 
% v_ang_fuel_ideal = ((pump_eff*v_rad_fuel*tan(blade_angle*pi/180) + sqrt(pump_eff*(4*9.80665*head_fuel+pump_eff*v_rad_fuel^2*tan(blade_angle*pi/180))) / (2*pump_eff*r_pump)));
% v_ang_ox_ideal = ((pump_eff*v_rad_ox*tan(blade_angle*pi/180) + sqrt(pump_eff*(4*9.80665*head_ox+pump_eff*v_rad_fuel^2*tan(blade_angle*pi/180))) / (2*pump_eff*r_pump)));
% 
% v_ang = (v_ang_fuel_ideal + v_ang_ox_ideal); % rad/s - angular velocity is a compromise between the pumps
% 
% power_fuel = mdot_fuel_cc*9.80665*head_fuel/pump_eff; % W
% power_ox = mdot_ox_cc*9.80665*head_ox/pump_eff; % W