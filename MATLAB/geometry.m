%% Nozzle Parameters

Ae_At = cea.output.eql.aeat(end); % Expansion / Area Ratio

c_star_ideal = cea.output.eql.cstar(3); % Characteristic Velocity
c_star = c_star_ideal*c_star_eff;

c_tau_ideal = gamma * sqrt( (2/(gamma-1) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * (1-(p_amb/pc)^((gamma-1)/gamma) ))); % Ideal Thrust Coefficient
lambda_cone = 0.5*(1+cos(diverge_angle)); % Divergence correction factor for conical nozzle
c_tau = c_tau_ideal*lambda_cone*c_tau_eff; % Thrust Coefficient

v_exhaust_ideal = cea.output.eql.sonvel(3)*cea.output.eql.mach(3); % m/s - ideal exhaust velocity
v_exhaust = c_star*c_tau; % m/s - actual exhaust velocity
isp_ideal = cea.output.eql.isp(end); % s - specific impulse

%% Mass Flow Rates
mdot_cc = thrust/v_exhaust;     % kg/s - Propellant mass flow rate into combustion chamber
mdot_fuel_cc = mdot_cc*(1/(1+OF)); % kg/s - Fuel Mass Flow Rate
mdot_ox_cc = mdot_cc*(OF/(1+OF)); % kg/s - Oxidizer Mass Flow Rate

mdot_gg = mdot_cc*gg_fraction/(1-gg_fraction); % kg/s - Propellant mass flow rate into gas generator
mdot_fuel_gg = mdot_gg*(1/(1+OF_gg)); % kg/s - Fuel Mass Flow Rate
mdot_ox_gg = mdot_gg*(OF_gg/(1+OF_gg)); % kg/s - Oxidizer Mass Flow Rate

mdot_total = mdot_cc+mdot_gg; % kg/s - Propellant mass flow rate through pumps
mdot_fuel_total = mdot_fuel_cc+mdot_fuel_gg; % kg/s
mdot_ox_total = mdot_ox_cc+mdot_ox_gg; % kg/s

%% Cross Section Areas
R_gas = 8.3145 / (0.001*cea.output.eql.mw(2)); % J/kg-K - Specific Gas Constant (throat)

A_throat = mdot_cc*sqrt(Tc)/pc * sqrt(R_gas/gamma) * ((gamma+1)/2)^((gamma+1)/(2*(gamma-1))); % m2 - throat area
A2_throat = thrust / pc; % m2 - throat area
A3_throat = mdot_cc*c_star/pc; % m2 - throat area
d_throat = sqrt(4*A_throat/pi); % m - throat diameter
r_throat = d_throat/2;
A_exit = A_throat*Ae_At; % m2 - exit area
d_exit = sqrt(4*A_exit/pi); % m - exit diameter

Ac_At = pi*d1_chamber^2 / (4*A_throat);

if Ac_At < 2
    fprintf ("Contraction ratio is %g, below recommended minimum of 2\n", (d1_chamber/d_throat)^2);
elseif Ac_At > 5
    fprintf("Contraction ratio is %g, above recommended maximum of 5\n", (d1_chamber/d_throat)^2);
end

%% Chamber Length
l_chamber = l_star*(r_throat/r_chamber)^2 - r_chamber*(r_throat^2+pi*sqrt(r_chamber*r_throat)/(tan(converge_angle)*3*r_chamber^2) - r_chamber/(3*tan(converge_angle)));

%% Contours
x1_throat = 0.5*(d1_chamber-d_throat)/tan(converge_angle) + l_chamber; % Projected point, onto central axis, of the converging straight line contour IF throat had no curvature
x2_throat = x1_throat + r_throat*(sin(converge_angle)+cot(converge_angle)*(cos(converge_angle)-1)); % Central point of the throat
x3_throat = x2_throat + r_throat*(sin(diverge_angle)+cot(diverge_angle)*(cos(diverge_angle)-1)); % Projected point, onto central axis, of the diverging straight line contour IF throat had no curvature

x_exit = x3_throat + 0.5*(d_exit-d_throat)/tan(diverge_angle); % Position of nozzle exit, aka length of the entire chamber + nozzle
x = 0:dx:x_exit; % Position domain

% Define inner contour
r1 = ones(1,length(x))*0.5*d1_chamber; % Chamber

    % Convergence
    r1(floor(l_chamber/dx) : end) = 0.5 * ( d1_chamber - d_throat ) * ( x1_throat - x(floor(l_chamber/dx):end) ) / ( x1_throat - x(floor(l_chamber/dx)) ) + 0.5*d_throat;

    % Divergence
    r1(floor(x3_throat/dx) : end) = 0.5 * ( d_exit-d_throat ) * ( x(floor(x3_throat/dx):end) - x3_throat ) / ( x(end) - x3_throat ) + 0.5*d_throat;

    % Throat Arc
        x_arc = -r_throat*sin(converge_angle):dx:r_throat*sin(diverge_angle);
        throat_arc = 0.5*d_throat - sqrt( r_throat^2 - x_arc.^2 ) + r_throat;
        [~,index] = min(throat_arc);
    r1(ceil(x2_throat/dx)-index+1 : ceil(x2_throat/dx) + length(x_arc) - index) = throat_arc;

% Define outer contour
r2 = r1 + thickness;

pipe_bound1 = find(r1 <= merge_radius, 1, 'first'); % First index where inner contour is at a given radius
pipe_bound2 = find(r1 <= merge_radius, 1, 'last'); % Last index where innner contour is at a given radius

n_pipe = zeros(1,length(x));
n_pipe( 1 : pipe_bound1 ) = n_pipe1;
n_pipe( pipe_bound1 : pipe_bound2) = n_pipe2;
n_pipe( pipe_bound2 : end) = n_pipe3;

if min(n_pipe) <= 0
    error("Geometric conditions for coolant channels are invalid");
end

if pipe_bound2 < length(x)
    n_pipe( pipe_bound2 : end ) = n_pipe3;
end

w_pipe = 2*pi*r2./n_pipe - gap_pipe; % m - coolant channel width
if max(w_pipe<0) == 1
    error("Sections of coolant channels have zero thickness");
end
