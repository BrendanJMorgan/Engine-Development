%% Inputs 

clear all

thrust = 3000*4.44822; % N - Thrust
pc = 250*6894.76; % Pa - Stagnation / Chamber Pressure
p_amb = 13.49*6894.76; % psi - ambient pressure at 2400 feet elevation
Tamb = 293; % K - Ambient Temperature
OF = 1.2; % Oxidizer/Fuel Ratio
proof = 0.95; % How much ethanol in fuel
c_star_eff = 0.75; % Characteristic Vel Efficiency, experimental
c_tau_eff = 0.96; % Thrust Coefficient Efficiency Factor

film_fraction = 0.1; % Fraction of the fuel mass flow dedicated to film cooling orifices - typically 3%-10% (Huzel and Huang)
film_width = 0.5*0.0254;

p_gg = 1000*6894.76; % Pa - chamber pressure inside gas generator

% Geometry
dx = 0.001; % m - position step
converge_angle = 45*pi/180; % rad
diverge_angle = 15*pi/180; % rad
l_chamber = 16*0.0254; % m
r_throat = 1*0.0254; % m - radius of curvature around the throat
d2_chamber = 8*0.0254; % m
thickness = 1/4*0.0254; % m
d1_chamber = d2_chamber - 2*thickness; % m

% Coolant Channels 
n_pipe1 = 16; % number of channels along barrel
n_pipe2 = 16; % number of channels near throat
n_pipe3 = 16; % number of channels along lower nozzle section
gap_pipe = 1/8*0.0254; % Gap between channels (fin thickness)
h_pipe = 0.5*0.0254; % m - coolant channel height
merge_radius = 0.45*d1_chamber; % m - when contour is below this radius, transition to n_pipe2
flow_direction = -1; % 1 = forward flow (injector to nozzle), -1 = counter flow (nozzle to injector)

%% Properties

wall = "aluminum";

switch wall
    case "aluminum" % 6061
        k_wall = 253; % W/m-K - thermal conductivity (at 571 C)
    case "steel"
        k_wall = 42.6; % W/m-K - thermal conductivity (at 100 C)
        case "copper"
        k_wall = 398; % W/m-K - thermal conductivity
end


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

% figure(2)
% clf
% colororder('default')
% plot(x,r1,x,r2,x,-1*r1,x,-1*r2, 'color','blue');
% axis equal
% xlabel("Distance from Injector (m)");
% title("Combustion Chamber Contours")

figure(1)
clf
plot(x,T_wall_cold,x,T_wall_hot,x,T_cool,x,Tf,x,Tab,x,Tref)
yline(0)
legend("Cold Wall","Hot Wall","Coolant","Free-Stream Gas","Adiabatic (no cooling)","Gas Property Reference",'Location','northwest');
xlabel("Distance from Injector (m)");
ylabel("Temperature (K)");
title("Engine Steady-State Temperatures")

% figure(3)
% clf
% yyaxis left
% colororder('default')
% plot(x,p_gas/6894.76, 'color','blue');
% ylabel("Chamber Pressure (psi)")
% hold on
% yyaxis right
% colororder('default')
% plot(x,yield_cc./hoop,'color','red');
% hold on
% hold off
% xlabel("Distance from Injector");
% plot(x,FS_design*ones(1,length(x)),'--','color','red')
% ylabel("Chamber Wall FSy");
% title("Combustion Chamber Structural Integrity")

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
