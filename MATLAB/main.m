%% Inputs
clear all
% clf

thrust = 1000*4.44822; % N - Thrust
pc = 250; % psi - Stagnation / Chamber Pressure
pamb = 13.49; % psi - ambient pressure at 2400 feet elevation
Tamb = 293; % K - Ambient Temperature
OF = 1.4; % Oxidizer/Fuel Ratio
proof = 0.95; % How much ethanol in fuel

% Geometry
converge_angle = 45*pi/180;
diverge_angle = 15*pi/180;
l_chamber = 12*0.0254; % m       
d_chamber = 3.875*0.0254; % m
r_throat = 0.5*0.0254; % m - radius of curvature around the throat
thickness = 0.0025; % m


%% Runs

combustion
exhaust_flow
thermal_balance



%% Results

figure(1)
plot(x,T_wall_cold,x,T_wall_hot,x,T_cool)
legend("Cold Wall","Hot Wall","Coolant",'Location','southeast');
xlabel("Distance from Injector (m)");
ylabel("Temperature (K)");

figure(2)
plot(x,p_cool*0.000145038);
xlabel("Distance from Injector (m)");
ylabel("Pressure (psi)");

figure(3)
plot(x,v_cool);
xlabel("Distance from Injector (m)");
ylabel("Velocity (m/s)");

figure(4)
plot(x,r1,x,r2,x,-r1,x,-r2, 'color','blue');
title("Combustion Chamber Contours")
axis equal
