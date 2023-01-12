%% Inputs
clear all

thrust = 1000*4.44822; % N - Thrust
pc = 250*6894.76; % Pa - Stagnation / Chamber Pressure
p_amb = 13.49*6894.76; % psi - ambient pressure at 2400 feet elevation
Tamb = 293; % K - Ambient Temperature
OF = 1.4; % Oxidizer/Fuel Ratio
proof = 0.95; % How much ethanol in fuel

% Geometry
converge_angle = 45*pi/180; % rad
diverge_angle = 15*pi/180; % rad
l_chamber = 12*0.0254; % m       
r_throat = 0.25*0.0254; % m - radius of curvature around the throat
d2_chamber = 4*0.0254; % m
thickness = 0.125*0.0254; % m
d1_chamber = d2_chamber - 2*thickness;

% Coolant Channels
n_pipe = 8; % number of channels
gap_pipe = 0.125*0.0254; % Gap between channels
h_pipe = 0.25*0.0254; % m - coolant channel height

%% Properties

k_al6061 = 167; % W/m-K - thermal conductivity of aluminum 6061

density_fuel =	795.965; % kg/m3 - ethanol at STP
cp_fuel = 2570; % J/kg-K - specific heat
kin_visc_fuel = 1.34E-06; % m2/s - kinematic viscosity
k_fuel = 0.167; % W/m-K - thermal conjductivity
boiling_fuel = 351.5; % K - boiling point at atmospheric pressure

density_ox = 1141; % kg/m3 - lox at boiling

%% Runs

combustion
exhaust_flow
thermal_balance
pump
gas_generator

%% Results

figure(1)
plot(x,T_wall_cold,x,T_wall_hot,x,T_cool)
legend("Cold Wall","Hot Wall","Coolant",'Location','southeast');
xlabel("Distance from Injector (m)");
ylabel("Temperature (K)");

% figure(2)
% plot(x,p_cool*0.000145038);
% xlabel("Distance from Injector (m)");
% ylabel("Pressure (psi)");
% 
% figure(3)
% plot(x,v_cool);
% xlabel("Distance from Injector (m)");
% ylabel("Velocity (m/s)");
% 
% figure(4)
% plot(x,r1,x,r2,x,-r1,x,-r2, 'color','blue');
% title("Combustion Chamber Contours")
% axis equal
