clear all

%% Inputs 

% Constants
g = 9.81;
p_amb = 13.49*6894.76; % psi - ambient pressure at 2400 feet elevation
T_amb = 293; % K - Ambient Temperature

% Overall Engine Performance Targets
thrust = 5000*4.44822; % N - Thrust
pc = 300*6894.76; % Pa - Stagnation / Chamber Pressure

OF = 1.4; % Oxidizer/Fuel Ratio
proof = 0.95; % How much ethanol in fuel
c_star_eff = 0.75; % Characteristic Vel Efficiency, experimental
c_tau_eff = 0.96; % Thrust Coefficient Efficiency Factor

film_fraction = 0.1; % Fraction of the fuel mass flow dedicated to film cooling orifices - typically 3%-10% (Huzel and Huang)

p_gg = 1000*6894.76; % Pa - chamber pressure inside gas generator
gg_fraction = 0.05; % Fraction of total mass flow sent to the gas generator. Context: F1 = 0.030, J2 = 0.014
OF_gg = 0.3; % OF Ratio - "[Most] operate at mixture ratios from 0.2 to 1.0, with hydrocarbons falling in the lower end, about 0.3" (NASA 1972)

% Geometry
dx = 0.002; % m - position step 
converge_angle = 45*pi/180; % rad
diverge_angle = 15*pi/180; % rad
l_star = 2.75; % m
rc_throat = 1*0.0254; % m - radius of curvature around the throat
d2_chamber = 8*0.0254; % m
thickness = 1/16*0.0254; % m
d1_chamber = d2_chamber - 2*thickness; % m
r_chamber = d1_chamber/2;

% Coolant Channels 
n_pipe1 = 16; % number of channels along barrel
n_pipe2 = 8; % number of channels near throat
n_pipe3 = 16; % number of channels along lower nozzle section
gap_pipe = 1/8*0.0254; % Gap between channels (fin thickness)
h_pipe = 0.25*0.0254; % m - coolant channel height
merge_radius = 0.45*d1_chamber; % m - when contour is below this radius, transition to n_pipe2
flow_direction = -1; % 1 = forward flow (injector to nozzle), -1 = counter flow (nozzle to injector)

% Turbomachinery
shaft_speed = 20000*0.1047198; % rad/s - angular velocity of the shaft; thus also the angular velocity of the turbine and both pump impellers+inducers (there is no gearing)
r_shaft = (0.25+0.05)*0.0254; % m
impeller_thickness = 1/8*0.0254; % m - thickness of impeller at the exit point, not including blades
impeller_height = 0.8*0.0254; % m - from base of impeller to eye plane

%% Properties

wall = "steel";

switch wall
    case "aluminum" % 6061
        k_wall = 253; % W/m-K - thermal conductivity (at 571 C)
    case "steel"
        k_wall = 42.6; % W/m-K - thermal conductivity (at 100 C)
    case "stainless"
        k_wall = 14.0;
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
% geometry
% exhaust_flow
% coolant_flow
% thermal_balance
% structures
% pump
% turbine
% gas_generator



%% Results

isp_ideal
mdot_total

% writematrix([impeller_curve/0.0254, zeros(length(impeller_curve(:,1)),1)], 'impeller_curve_inches.txt', 'Delimiter', ',')  
% writematrix([shroud_curve/0.0254, zeros(length(shroud_curve(:,1)),1)], 'shroud_curve_inches.txt', 'Delimiter', ',')

% figure(2)
% clf
% colororder('default')
% plot(x,r1,x,r2,x,-1*r1,x,-1*r2, 'color','blue');
% axis equal
% xlabel("Distance from Injector (m)");
% title("Combustion Chamber Contours")

% figure(1)
% clf
% plot(x,T_wall_cold,x,T_wall_hot,x,T_cool,x,T_free,x,Tab,x,Tref,x,T_film)
% yline(0)
% legend("Cold Wall","Hot Wall","Coolant","Free-Stream Gas","Adiabatic (no cooling)","Gas Property Reference","Fuel Film",'Location','northeast');
% xlabel("Distance from Injector (m)");
% ylabel("Temperature (K)");
% title("Engine Steady-State Temperatures")

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