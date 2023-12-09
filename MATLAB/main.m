clear all
clf

%% TO DO

% Need to adjust all coolprop calls with ethanol to include the water content
% Reevaluate why CEA machs are not even close to isentropic relations
% Decide if c tau should be changed to a bell nozzle
% Compare CEA transport properties to those found from mixture() function

%% Inputs 

% Constants
g = 9.81;
p_amb = 13.49*6894.76;      % psi - ambient pressure at 2400 feet elevation
T_amb = 293;                % K - Ambient Temperature

% Overall Engine Performance Targets
thrust_target = 3000*4.44822;   % N - Thrust
p_cc = 250*6894.76;             % Pa - Stagnation / Chamber Pressure

% Combustion Chamber (CC)
OF = 1.4;                   % Oxidizer/Fuel Ratio (by mass)
proof = 0.95;               % How much ethanol in fuel, remaining part is water (by mass)
c_star_eff = 0.75;          % Characteristic Vel Efficiency, experimental
c_tau_eff = 0.96;           % Thrust Coefficient Efficiency Factor
gamma_guess = 1.22;
c_tau_guess = 0.983*c_tau_eff*sqrt( (2*gamma_guess^2/(gamma_guess-1) * (2/(gamma_guess+1))^((gamma_guess+1)/(gamma_guess-1)) * (1-(p_amb/p_cc)^((gamma_guess-1)/gamma_guess) ) ) ); 
A_throat = thrust_target / (p_cc*c_tau_guess*c_star_eff); % m2 - Throat Area

% Gas Generator (GG)
p_gg = 500*6894.76;         % Pa - chamber pressure inside gas generator
gg_fraction = 0.05;         % Fraction of total mass flow sent to the gas generator. Context: F1 = 0.030, J2 = 0.014
OF_gg = 0.3;                % OF Ratio - "[Most] operate at mixture ratios from 0.2 to 1.0, with hydrocarbons falling in the lower end, about 0.3" (NASA 1972)

% Chamber/Nozzle Geometry
dx = 0.001;                 % m - position step 
converge_angle = 45*pi/180; % rad - half-cone convergene angle of combustion chamber end
diverge_angle = 15*pi/180;  % rad - half-cone divergence angle of nozzle
l_star = 1;                 % m - combustion length
rc_throat = 1*0.0254;       % m - radius of curvature around the throat
d2_chamber = 4.5*0.0254;    % m
thickness = 1/4*0.0254;     % m % SLANT VS VERTICAL THICKNESS
d1_chamber = d2_chamber - 2*thickness; % m
r1_chamber = d1_chamber/2;

% Coolant Channels
n_pipe1 = 16; % number of channels along barrel
n_pipe2 = 16; % number of channels near throat
n_pipe3 = 16; % number of channels along lower nozzle section
gap_pipe = 1/4*0.0254;          % Gap between channels (fin thickness)
h_pipe = 1/16*0.0254;           % m - coolant channel height
merge_radius = 0.45*d1_chamber; % m - when contour is below this radius, transition to n_pipe2
flow_direction = -1;            % 1 = forward flow (injector to nozzle), -1 = counter flow (nozzle to injector)

% Film Cooling
film_fraction = 0.03;   % Fraction of the fuel mass flow dedicated to film cooling orifices - typically 3%-10% (Huzel and Huang)
v_injection = 10;       % m/s - combustion gas must have some initial velocity for injector film cooling to work mathematically
injection_efficiency = 1.0;

% Turbomachinery
shaft_speed = 30000*0.1047198;	    % rad/s - angular velocity of the shaft, rotors, impeller, and inducers (no gearing)
r_shaft = 9/32*0.0254;				% m - a little bit of clearance around a 1/2 inch shaft
r_shaft_ss = 1/4*0.0254;            % m - portion of shaft that is stainless steel
impeller_thickness = 1/8*0.0254;	% m - thickness of impeller at the exit point, not including blades
impeller_height = 0.5*0.0254;		% m - from base of impeller to eye plane
turbine_stage_number = 2;			% number of stages in the turbine (one stage = one rotor + one stator)

%% Properties

wall = "aluminum";

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

% MAKE COOLPROP CALLS FOR PROPERTIES
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

% T_wall_hot = 800*ones(1,length(x)); % K - initial GUESS for the hot wall temperatures
% thermal_balance
% thermal_balance
% structures
pump
% gas_generator
turbine


%% Results

% thrust_lbf = thrust/4.44822;
% thrust_lbf
% isp_ideal
% isp_real
% mdot_total
impeller_diameter_inches = 2*r_exit/0.0254

% writematrix([impeller_curve/0.0254, zeros(length(impeller_curve(:,1)),1)], 'impeller_curve_inches.txt', 'Delimiter', ',')  
% writematrix([shroud_curve/0.0254, zeros(length(shroud_curve(:,1)),1)], 'shroud_curve_inches.txt', 'Delimiter', ',')
% 
% writematrix([blade_curve(:,1)/0.0254, zeros(length(impeller_curve(:,1)),1), blade_curve(:,2)/0.0254], 'blade_curve_inches.txt', 'Delimiter', ',')

% figure(2)
% clf
% colororder('default')
% plot(x,r1,x,r2,x,-1*r1,x,-1*r2, 'color','blue');
% axis equal
% xlabel("Distance from Injector (m)");
% title("Combustion Chamber Contours")

% figure(1)
% clf
% plot(x,T_wall_cold,x,T_wall_hot,x,T_cool,x,T_film,x,T_free,x,T_ab,x,T_recovery,x,T_ref)
% yline(0)
% legend("Cold Wall","Hot Wall","Regen Coolant","Film Coolant","Free-Stream Gas","Adiabatic", "Recovery", "Gas Property Reference",'Location','northeast');
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

toc