%% Requirements

p_out = p_cool(1); % Pa - outlet pressure
p_in = p_amb; % Pa - inlet pressure

vdot_fuel = mdot_fuel_cc/density_fuel; %m3/s - Fuel Volumetric Flow Rate
vdot_ox = mdot_ox_cc/density_ox; %m3/s - Oxidizer Mass Flow Rate

%% Inputs

r_pump = 2*0.0254; % m - blade radius
h_exit_pump = 1*0.0254; % m - exit gap height

pump_eff = 0.55;
blade_angle = 30*pi/180;


%%

head_fuel = (p_out-p_in)/(density_fuel*9.80665); % m - pump pressure head
head_ox = (p_out-p_in)/(density_ox*9.80665); % m - pump pressure head

v_rad_fuel = vdot_fuel/(2*pi*r_pump*h_exit_pump); % m/s - radial velocity
v_rad_ox = vdot_ox/(2*pi*r_pump*h_exit_pump); % m/s - radial velocity

v_ang_fuel_ideal = ((pump_eff*v_rad_fuel*tan(blade_angle*pi/180) + sqrt(pump_eff*(4*9.80665*head_fuel+pump_eff*v_rad_fuel^2*tan(blade_angle*pi/180))) / (2*pump_eff*r_pump)));
v_ang_ox_ideal = ((pump_eff*v_rad_ox*tan(blade_angle*pi/180) + sqrt(pump_eff*(4*9.80665*head_ox+pump_eff*v_rad_fuel^2*tan(blade_angle*pi/180))) / (2*pump_eff*r_pump)));

v_ang = (v_ang_fuel_ideal + v_ang_ox_ideal); % rad/s - angular velocity is a compromise between the pumps

power_fuel = mdot_fuel_cc*9.80665*head_fuel/pump_eff; % W
power_ox = mdot_ox_cc*9.80665*head_ox/pump_eff; % W