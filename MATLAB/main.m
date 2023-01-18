%% Inputs 

clear all

thrust = 3600*4.44822; % N - Thrust
pc = 300*6894.76; % Pa - Stagnation / Chamber Pressure
p_amb = 13.49*6894.76; % psi - ambient pressure at 2400 feet elevation
Tamb = 293; % K - Ambient Temperature
OF = 1.2; % Oxidizer/Fuel Ratio
proof = 0.95; % How much ethanol in fuel
c_star_eff = 0.75; % Characteristic Vel Efficiency, experimental
c_tau_eff = 0.96; % Thrust Coefficient Efficiency Factor

p_gg = 1000*6894.76; % Pa - chamber pressure inside gas generator

% Geometry
dx = 0.001; % m - position step
converge_angle = 45*pi/180; % rad
diverge_angle = 15*pi/180; % rad
l_chamber = 12*0.0254; % m
r_throat = 0.5*0.0254; % m - radius of curvature around the throat
d2_chamber = 6*0.0254; % m
thickness = 0.2*0.0254; % m
d1_chamber = d2_chamber - 2*thickness; % m

% Coolant Channels
n_pipe1 = 16; % number of channels along barrel
n_pipe2 = 8; % number of channels near throat
n_pipe3 = 16; % number of channels along lower nozzle section
gap_pipe = 1/16*0.0254; % Gap between channels (fin thickness)
h_pipe = 0.25*0.0254; % m - coolant channel height
merge_radius = 0.45*d1_chamber; % m - when contour is below this radius, transition to n_pipe2

%% Properties

k_al6061 = 167; % W/m-K - thermal conductivity of aluminum 6061
yield_al6061 = 276E6; % Pa - tensile strength of aluminum 6061
FS_design = 2.0; % Minimum design factor of safety

density_fuel =	795.965; % kg/m3 - ethanol at STP
cp_fuel = 2570; % J/kg-K - specific heat
kin_visc_fuel = 1.34E-06; % m2/s - kinematic viscosity
k_fuel = 0.167; % W/m-K - thermal conductivity

density_ox = 1141; % kg/m3 - lox at boiling

%% Runs and Plots

combustion
geometry
exhaust_flow
coolant_flow
thermal_balance
structures
pump
gas_generator
turbine


%% Results

% figure(1)
% clf
% yyaxis left
% colororder('default')
% plot(x,r1,x,r2,x,-1*r1,x,-1*r2, 'color','blue');
% hold on
% yyaxis right
% colororder('default')
% plot(x,n_pipe,'color','red');
% hold off
% xlabel("Distance from Injector (m)");
% title("Combustion Chamber Contours")

figure(2)
clf
plot(x,T_wall_cold,x,T_wall_hot,x,T_cool,x,Tf,x,Tab,x,Tref)
yline(0)
legend("Cold Wall","Hot Wall","Coolant","Free-Stream Gas","Adiabatic (no cooling)","Gas Property Reference",'Location','northwest');
xlabel("Distance from Injector (m)");
ylabel("Temperature (K)");
title("Engine Steady-State Temperatures")

figure(3)
clf
yyaxis left
colororder('default')
plot(x,p_gas/6894.76, 'color','blue');
ylabel("Chamber Pressure (psi)")
hold on
yyaxis right
colororder('default')
plot(x,yield_cc./hoop,'color','red');
hold on
hold off
xlabel("Distance from Injector");
plot(x,FS_design*ones(1,length(x)),'--','color','red')
ylabel("Chamber Wall FSy");
title("Combustion Chamber Structural Integrity")


% isp
% mdot_cc

% figure(2)
% plot(x,p_cool*0.000145038);
% xlabel("Distance from Injector (m)");
% ylabel("Pressure (psi)");
% 
% figure(3)
% plot(x,v_cool);
% xlabel("Distance from Injector (m)");
% ylabel("Velocity (m/s)");
