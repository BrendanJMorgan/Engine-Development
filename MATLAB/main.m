clear all
format short g
PropsSI = @py.CoolProp.CoolProp.PropsSI;

%% Inputs 

% Constants
g = 9.81;
p_amb = 13.62*6894.76;      % psi - ambient pressure at 2100 feet elevation
T_amb = 293;                % K - Ambient Temperature

% Overall Engine Performance Targets
thrust_target = 3000*4.44822;   % N - Thrust
p_cc = 400*6894.76;             % Pa - Chamber (Stagnation) Pressure

% Combustion Chamber (CC)
OF = 1.4;                   % Oxidizer/Fuel Ratio (by mass)
proof = 0.95;               % How much ethanol in fuel, remaining part is water (by mass)
c_star_eff = 0.75;          % Characteristic Vel Efficiency, experimental
c_tau_eff = 0.96;           % Thrust Coefficient Efficiency Factor
gamma_guess = 1.15;
c_tau_guess = 0.983*c_tau_eff*sqrt( (2*gamma_guess^2/(gamma_guess-1) * (2/(gamma_guess+1))^((gamma_guess+1)/(gamma_guess-1)) * (1-(p_amb/p_cc)^((gamma_guess-1)/gamma_guess) ) ) ); 
A_throat = thrust_target / (p_cc*c_tau_guess*c_star_eff); % m2 - Throat Area

% Gas Generator (GG)
p_gg = 400*6894.76;             % Pa - chamber pressure inside gas generator
gg_fraction_guess = 0.01;       % unitless - Fraction of total mass flow sent to the gas generator. Context: F1 = 0.030, J2 = 0.014
OF_gg = 0.2;                    % unitless - OF Ratio - "[Most] operate at mixture ratios from 0.2 to 1.0, with hydrocarbons falling in the lower end, about 0.3" (NASA 1972)
c_star_eff_gg = 0.75;           % m/s - characteristic velocity efficiency, experimental ANY PAPERS ON THIS?
c_tau_eff_gg = 0.96;            % unitless - Thrust Coefficient Efficiency Factor

% Pumps
shaft_speed = 20000*pi/30;          % rad/s - angular velocity of the pump shaft, impeller, and inducers
gear_ratio = 80/20;                 % unitless - higher makes for a smaller, faster turbine
gear_efficiency = 0.95;             % unitless - arbitrary right now, but typical for spur gears
clock = 1;                          % 1 for counterclockwise and -1 for clockwise (looking down at ox pump)
r_eye_inner = 1/4*0.0254;	        % m - a little bit of clearance around a 1/2 inch shaft
r_shaft = 0.005;                    % m - portion of shaft that is stainless steel
impeller_thickness = 1/8*0.0254;    % m - thickness of impeller at the exit point, not including blades
cavitation_target = 2.0;            % unitless

% Turbine
nozzle_number = 5;              % unitless
nozzle_angle = 45*pi/180;       % rad
diverge_angle_gg = 15*pi/180;   % rad - half-cone divergence angle of nozzle plate
blade_number_rotor = 32;        % unitless - total number of blades on the rotor
blade_gap_rotor = 3/16*0.0254;   % m - perpendicular (i.e. the end mill diameter)

% Chamber/Nozzle Geometry
dx = 0.001;                 % m - position step 
converge_angle = 45*pi/180; % rad - half-cone convergene angle of combustion chamber end
diverge_angle = 15*pi/180;  % rad - half-cone divergence angle of nozzle
l_star = 1;                 % m - combustion length
rc_throat = 1*0.0254;       % m - radius of curvature around the throat
d2_chamber = 5*0.0254;      % m
thickness = 1/4*0.0254;     % m % SLANT VS VERTICAL THICKNESS
d1_chamber = d2_chamber - 2*thickness; % m
r1_chamber = d1_chamber/2;  % m

% Coolant Channels
n_pipe1 = 16;   % number of channels along barrel
n_pipe2 = 16;   % number of channels near throat
n_pipe3 = 16;   % number of channels along lower nozzle section
gap_pipe = 1/4*0.0254;          % m - Gap between channels (fin thickness)
h_pipe = 1/16*0.0254;           % m - coolant channel height
merge_radius = 0.45*d1_chamber; % m - when contour is below this radius, transition to n_pipe2
flow_direction = -1;            % 1 = forward flow (injector to nozzle), -1 = counter flow (nozzle to injector)

% Film Cooling
film_fraction = 0.03;               % unitless - Fraction of the fuel mass flow dedicated to film cooling orifices - typically 3%-10% (Huzel and Huang)
v_injection = 10;                   % m/s - combustion gas must have some initial velocity for injector film cooling to work mathematically
injection_efficiency = 1.0;         % I forgot what this is



%% Properties

k_wall = 253; % W/m-K - thermal conductivity of aluminum (at 571 C)

FSu_design = 2.0; % Minimum design factor of safety, ultimate

shaft_tensile_strength = 510E6; % Pa - stainless steel 304, ultimate
shaft_shear_strength = shaft_tensile_strength/sqrt(3);  % Pa - conservative Tresca criterion


%% Runs

combustion_chamber
cc_geometry
cc_gas_flow
coolant_flow

% T_wall_hot = 800*ones(1,length(x)); % K - initial GUESS for the hot wall temperatures
% thermal_balance
% thermal_balance % intentioanlly run twice (greater convergence)

powerhead

%% Results

% thrust_lbf = thrust/4.44822;
% thrust_lbf
% isp_ideal
% isp_real
% mdot_total
% 
% figure(1); clf
% colororder('default')
% plot(x,r1,x,r2,x,-1*r1,x,-1*r2, 'color','blue');
% axis equal
% xlabel("Distance from Injector (m)");
% title("Combustion Chamber Contours")
% 
% figure(2); clf
% plot(x,T_wall_cold,x,T_wall_hot,x,T_cool,x,T_film,x,T_free,x,T_ab,x,T_recovery,x,T_ref)
% yline(0)
% legend("Cold Wall","Hot Wall","Regen Coolant","Film Coolant","Free-Stream Gas","Adiabatic", "Recovery", "Gas Property Reference",'Location','northeast');
% xlabel("Distance from Injector (m)");
% ylabel("Temperature (K)");
% title("Engine Steady-State Temperatures")
% 
% % Pump Impeller and Blades
% figure(3); clf;
% line(shroud_curve_ox(:,1)/0.0254, shroud_curve_ox(:,2)/0.0254)
% hold on
% axis equal
% plot(impeller_curve_ox(:,1)/0.0254, impeller_curve_ox(:,2)/0.0254)
% line([0 0], ylim);  % x-axis
% line(xlim, [0 0]);  % y-axis
% title("Impeller and Shroud Contours")
% hold off

% Blades
% figure(4); clf
% hold on
% delta_angle = 2 * pi / blade_number_ox; % Calculate the angle to rotate each blade
% for i = 0:(blade_number-1)
%     rotation_matrix = [cos(i * delta_angle), -sin(i * delta_angle); sin(i * delta_angle), cos(i * delta_angle)];
%     rotated_curve = blade_curve_ox * rotation_matrix';
%     plot(rotated_curve(:, 1)/0.0254, rotated_curve(:, 2)/0.0254, 'LineWidth', 2); 
%     plot(NaN, NaN); % Prevent connection between different blades
% end
% plot(volute_curve_ox(:,1)/0.0254, volute_curve_ox(:,2)/0.0254);
% 
% delta_angle = 2 * pi / blade_number_fuel; % Calculate the angle to rotate each blade
% for i = 0:(blade_number-1)
%     rotation_matrix = [cos(i * delta_angle), -sin(i * delta_angle); sin(i * delta_angle), cos(i * delta_angle)];
%     rotated_curve = blade_curve_fuel * rotation_matrix';
%     plot(rotated_curve(:, 1)/0.0254, rotated_curve(:, 2)/0.0254-6, 'LineWidth', 2); 
%     plot(NaN, NaN); % Prevent connection between different blades
% end
% plot(volute_curve_fuel(:,1)/0.0254, volute_curve_fuel(:,2)/0.0254-6);
% 
% plot(r_pitchline/0.0254*cos(theta), r_pitchline/0.0254*sin(theta), 'LineStyle','--')
% plot(r_rotor_tip/0.0254*cos(theta), r_rotor_tip/0.0254*sin(theta), r_rotor_base/0.0254*cos(theta), r_rotor_base/0.0254*sin(theta))
% 
% title('Impellers, Volutes, and Rotor');
% axis equal;
% grid on;
