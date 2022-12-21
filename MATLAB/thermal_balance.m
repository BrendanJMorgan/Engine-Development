tic

%% Inputs
clear all
clf

thrust = 1000*4.44822; % N - Thrust
pc = 250; % psi - Stagnation / Chamber Pressure
pamb = 13.49; % psi - ambient pressure at 2400 feet elevation
Tamb = 293; % K - Ambient Temperature
OF = 1.5; % Oxidizer/Fuel Ratio
proof = 0.95; % How much ethanol in fuel

%% CEA Run Rocket

 cea = CEA('problem','rocket', ...
     'equilibrium', ...  % equilibrium flow, infinite chamber, 
     'o/f',OF, ...          % Define OF ratio
     'case','CEAM-rocket1', ...
     'p(psi)',pc,'pi/p',pc/pamb, ... % Define nozzle with pressure ratio
     'reactants', ... 
     'fuel','C2H5OH(L)','wt%',proof,'t(k)',Tamb,'H2O(L)','wt%',1-proof,'t(k)',Tamb, ... % 190-Proof Ethanol for Fuel, at ambient
     'oxid','O2(L)','wt%',100, ... % Liquid Oxygen for Oxidizer, at boiling point
     'output','short','massf','transport','mks', ...
     'end');

pc = pc*6894.76; % Converting to Pa

Tc = cea.output.eql.temperature(1); % Stagnation / Chamber Temperature
gamma = mean(cea.output.eql.gamma); % Ratio of Specific Heats % MAKE X DEPENDENT

c_star_ideal = cea.output.eql.cstar(3); % Characteristic Velocity
c_star_eff = 0.75; % Char Vel Efficiency, experimental
c_star = c_star_ideal*c_star_eff;

M_inj = cea.output.eql.mach(1); % Mach at injector
M_comb = 0.002; % Mach at start of convergence (combustor?)
M_throat = cea.output.eql.mach(2); % Mach at throat, should be 1 unless something is very wrong
if M_throat < 0.99 || M_throat > 1.01
    error("Sonic point is not at throat");
end
M_exit = cea.output.eql.mach(3); % Mach at exit, should be supersonic

%% Define Geometry

Ae_At = ((gamma+1)/2)^(-0.5*(gamma+1)/(gamma-1))*(1/M_exit)*(1+(gamma-1)/2*M_exit^2)^((gamma+1)/(2*(gamma-1))); % Expansion / Area Ratio
%a_throat = c_star*m_dot/pc;

converge_angle = 45*pi/180;
diverge_angle = 15*pi/180;
l_chamber = 12*0.0254; % m       
d_chamber = 3.875*0.0254; % m

d_throat = 0.0502680; % m - make dependent on CEA run
d_exit = 0.09586; % m - make dependent on CEA run
r_throat = 0.01; % m - radius of curvature around the throat
thickness = 0.0025; % m

c_tau_ideal = gamma * sqrt( (2/(gamma-1) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * (1-(pamb/pc)^((gamma-1)/gamma) ))) + ((pc-pamb)/pc)*Ae_At; % Ideal Thrust Coefficient
lambda_cone = 0.5*(1+cos(diverge_angle)); % Divergence correction factor for conical nozzle
c_tau_eff = 0.96; % Thrust Coefficient Efficiency Factor
c_tau = c_tau_ideal*lambda_cone*c_tau_eff; % Thrust Coefficient

x1_throat = 0.5*(d_chamber-d_throat)/tan(converge_angle) + l_chamber; % Projected point of the converging straight line contour IF throat had no curvature
x2_throat = x1_throat + r_throat*(sin(converge_angle)+cot(converge_angle)*(cos(converge_angle)-1)); % Central point of the throat
x3_throat = x2_throat + r_throat*(sin(diverge_angle)+cot(diverge_angle)*(cos(diverge_angle)-1)); % Projected point of the diverging straight line contour IF throat had no curvature

x_exit = x3_throat + 0.5*(d_exit-d_throat)/tan(diverge_angle); % Position of nozzle exit, aka length of the entire chamber + nozzle

dx = 0.001; % m - position step
x = 0:dx:x_exit; % Position domain

% Define inner contour
r1 = ones(1,length(x))*0.5*d_chamber; % Chamber
    % Convergence
    r1(ceil(l_chamber/dx):end) = 0.5*(d_chamber-d_throat)*(x1_throat-x(ceil(l_chamber/dx):end))/(x1_throat-l_chamber) + 0.5*d_throat;
    % Divergence
    r1(floor(x3_throat/dx):end) = 0.5*(d_exit-d_throat)*(x(floor(x3_throat/dx):end)-x3_throat)/(x_exit-x3_throat) + 0.5*d_throat;
    % Throat Arc
        x_arc = -r_throat*sin(converge_angle):dx:r_throat*sin(diverge_angle);
        throat_arc = 0.5*d_throat - sqrt( r_throat^2 - x_arc.^2 ) + r_throat;
        [~,index] = min(throat_arc);
    r1(round(x2_throat/dx)-index+1:round(x2_throat/dx) + length(x_arc)-index) = throat_arc; 

% Define outer contour
r2 = r1 + thickness;

% figure(1);
% plot(x,r1,x,r2);
% axis equal;
% title("Engine Contour (m)");


%% Exhaust Flow Properties

mdot = thrust/(c_star*c_tau);     % kg/s - Propellant Mass Flow Rate
mdot_fuel = mdot*(1/(1+OF));     % kg/s - Fuel Mass Flow Rate
mdot_ox = mdot*(OF/(1+OF));  % kg/s - Oxidizer Mass Flow Rate

M = zeros(1,length(x)); % Mach Number

% Chamber
M(1:round(l_chamber/dx)) = linspace(M_inj, M_comb, round(l_chamber/dx));
dm = 0.0001;

% Converging
m = 0:dm:1;
for i = round(l_chamber/dx):1:round(x2_throat/dx)
    % Isentropic Area Relation
    [~,index] = min(abs( (1./m)*(0.5*(gamma+1)).^(-0.5*(gamma+1)./(gamma-1)).*(1+0.5*(gamma-1) .* m.^2).^((gamma+1)./(2*(gamma-1))) - (r1(i)*2./d_throat).^2 ) )   ;
    M(i) = (-dm + dm * index);
end
factor1 = (1-M_comb)/(1-M(round(l_chamber/dx)));
M(round(l_chamber/dx):round(x2_throat/dx)) = 1 - ((1-M(round(l_chamber/dx):round(x2_throat/dx))) * factor1); % Stretch to fit CEA data

% Diverging
m = 1:dm:3;
for i = round(x2_throat/dx):1:length(x)
    % Isentropic Area Relation
    [~,index] = min(abs( (1./m)*(0.5*(gamma+1)).^(-0.5*(gamma+1)./(gamma-1)).*(1+0.5*(gamma-1) .* m.^2).^((gamma+1)./(2*(gamma-1))) - (r1(i)*2./d_throat).^2 ) )   ;
    M(i) = (1 -dm + dm * index) * (M_exit);
end
factor2 = M_exit/M(length(x));
M(round(x2_throat/dx):length(x)) = M(round(x2_throat/dx):length(x)) * factor2; % Stretch to fit CEA data


%% Thermal Environment
Twall = 500; % K - wall temperature, independent variable
Tf = Tc ./ (1 + (gamma-1)/2 * M.^2); % Free-Stream Temperature
Tref = Tf .* (1 + 0.032*M.^2 + 0.58*(Twall./Tf-1));

% These all need to be updated with better numbers from CEA
cp = cea.output.eql.cp; % Specific Heat, J/kg-K
cp = interp1([0,x2_throat,x_exit], [cp(1),cp(2),cp(3)], x); % Interpolate
visc = cea.output.eql.viscosity; % Dynamic Viscosity, Pa-s
visc = interp1([0,x2_throat,x_exit], [visc(1),visc(2),visc(3)], x); % Interpolate
Pr = struct2array(cea.output.eql.prandtl(:,1)); % Prandtl Number
Pr = interp1([0,x2_throat,x_exit], [Pr(1),Pr(2),Pr(3)], x); % Interpolate

% Bartz relation
sigma = 1 ./ ( (0.5*(Twall/Tc)*(1+(gamma-1)/2 * M.^2)+0.5).^0.68 .* (1+(gamma-1)/2 * M.^2).^0.12 ); % Meaningless (?) coefficient
h_gas = (0.026./d_throat.^0.2).*(visc.^0.2 .* cp ./ Pr.^0.6 ) .* (pc./c_star).^0.8 .* (d_throat./r_throat)^0.1 .* (0.5*d_throat./r1).^1.8 .* sigma; % Convective Heat Transfer Coefficient

Tab = Tc.*(1+Pr.^0.33.*(gamma-1)/2.*M.^2)./(1+(gamma-1)/2.*M.^2); % Adiabatic Wall Temperature
q_gas = h_gas.*pi.*r1.*(Tab-Twall); % W/m - Linear Heat Flux
Qtotal = sum(q_gas*dx);



%% Coolant Properties
% Ethanol is the coolant
density_cool =	795.965; % kg/m3 - ethanol at STP
cp_cool = 2570; % J/kg-K - specific heat
kin_visc_cool = 1.34E-06; % m2/s - kinematic viscosity
cond_cool = 0.167; % W/m2-K - thermal conductivity

%% Coolant Loop Flow

d_pipe = 0.125*0.0254; % m - coolant pipe diameter
n_pipe = floor(floor(pi*(d_chamber+d_pipe) / d_pipe) / 2); % number of pipes
v_cool = mdot_fuel / (density_cool*n_pipe*pi*d_pipe^2 / 4); % m/s - fluid velocity
Re_cool = v_cool*d_pipe/kin_visc_cool; % Reynold's number
Pr_cool = cp_cool*density_cool*kin_visc_cool/cond_cool; % Prandtl Number
f_cool = (0.79 * log(Re_cool)-1.64)^(-2); % friction factor, smooth pipe approximation

if (Pr_cool < 0.5 || Pr_cool > 2000 || Re_cool < 3000 || Re_cool > 5E6)
    %error("Gnielenski correlation invalid");
end

Nu_cool = (f_cool/8)*(Re_cool-1000)*Pr_cool/(1+12.7*(f_cool/8)^0.5*(Pr_cool^(2/3)-1)); % Nusselt Number, Gnielenski correlation
h_cool = Nu_cool*cond_cool/d_pipe; % Convective heat transfer coefficient

%% Wall Conduction

k_al6061 =	167; % W/m2-K - thermal conductivity of aluminum 6061
fin_thickness = 0.125*0.0254; % m
fin_height = 0.25*0.0254; % m
fin_length  = x_exit; % m
m_coeff = sqrt(2*h_gas/(k_al6061*fin_thickness)); % coefficient for fin efficency equation
fin_eff = tanh(m_coeff*dx)/m_coeff*dx; % fin efficiency for each length step
area_cool = (2*fin_height*fin_eff+fin_thickness)*dx	; % adjusted contact area of coolant on channel walls

%% Thermal Balance

T_cool = zeros(1,length(x)+1);
T_cool(1) = Tamb; % K - coolant inlet temperature
T_wall_cold = zeros(1,length(x));
T_wall_hot = zeros(1,length(x));

toc
tic 

Q_gas = q_gas*pi*2*r1(i)*dx;

for i=1:1:length(x)-1
    T_wall_cold(i) = Q_gas(i) / (h_cool*area_cool*n_pipe) + T_cool(i);
    T_wall_hot(i) = Q_gas(i)*(r2(i)-r1(i))/(k_al6061*pi*2*r1(i)*dx) + T_wall_cold(i);
    
    dT = q_gas(i)*cp_cool*n_pipe / mdot_fuel; % Coolant Temperature change
    T_cool(i+1) = T_cool(i) + dT;
end

toc



















