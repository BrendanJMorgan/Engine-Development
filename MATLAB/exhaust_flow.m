%% Define Geometry

Ae_At = cea.output.eql.aeat(end); % Expansion / Area Ratio

c_tau_ideal = gamma * sqrt( (2/(gamma-1) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * (1-(pamb/pc)^((gamma-1)/gamma) ))); % + ((pc-pamb)/pc)*Ae_At; % Ideal Thrust Coefficient
lambda_cone = 0.5*(1+cos(diverge_angle)); % Divergence correction factor for conical nozzle
c_tau_eff = 0.96; % Thrust Coefficient Efficiency Factor
c_tau = c_tau_ideal*lambda_cone*c_tau_eff; % Thrust Coefficient

v_exhaust_ideal = cea.output.eql.sonvel(3)*cea.output.eql.mach(3); % m/s - ideal exhaust velocity
v_exhaust = c_star*c_tau; % m/s - actual exhaust velocity
isp = cea.output.eql.isp(end); % s - specific impulse

mdot = thrust/v_exhaust;     % kg/s - Propellant Mass Flow Rate
mdot_fuel = mdot*(1/(1+OF));     % kg/s - Fuel Mass Flow Rate
mdot_ox = mdot*(OF/(1+OF));  % kg/s - Oxidizer Mass Flow Rate

A_throat = thrust / (pc*c_tau); % m2 - throat area
d_throat = sqrt( A_throat/pi);
A_exit = A_throat*Ae_At;
d_exit = sqrt(A_exit/pi);

% A_exit = (thrust - mdot*v_exhaust) / (pc - pamb);
% if A_exit < 0
%     error("Maximum possible thrust: %g lbf, Isp: %g s", mdot*v_exhaust*0.224809, isp);
% end
% d_exit = sqrt( A_exit/pi );
% A_throat = A_exit / Ae_At;
% d_throat = sqrt( A_throat/pi );

% d_throat = 0.0502680; % m - make dependent on CEA run
% d_exit = 0.09586; % m - make dependent on CEA run



x1_throat = 0.5*(d_chamber-d_throat)/tan(converge_angle) + l_chamber; % Projected point, onto central axis, of the converging straight line contour IF throat had no curvature
x2_throat = x1_throat + r_throat*(sin(converge_angle)+cot(converge_angle)*(cos(converge_angle)-1)); % Central point of the throat
x3_throat = x2_throat + r_throat*(sin(diverge_angle)+cot(diverge_angle)*(cos(diverge_angle)-1)); % Projected point, onto central axis, of the diverging straight line contour IF throat had no curvature

x_exit = x3_throat + 0.5*(d_exit-d_throat)/tan(diverge_angle); % Position of nozzle exit, aka length of the entire chamber + nozzle

dx = 0.001; % m - position step
x = 0:dx:x_exit; % Position domain

% Define inner contour
r1 = ones(1,length(x))*0.5*d_chamber; % Chamber
    % Convergence
    r1(floor(l_chamber/dx) : end) = 0.5 * ( d_chamber - d_throat ) * ( x1_throat - x(floor(l_chamber/dx):end) ) / ( x1_throat - x(floor(l_chamber/dx)) ) + 0.5*d_throat;
    % Divergence
    r1(floor(x3_throat/dx) : end) = 0.5 * ( d_exit-d_throat ) * ( x(floor(x3_throat/dx):end) - x3_throat ) / ( x(end) - x3_throat ) + 0.5*d_throat;
    % Throat Arc
        x_arc = -r_throat*sin(converge_angle):dx:r_throat*sin(diverge_angle);
        throat_arc = 0.5*d_throat - sqrt( r_throat^2 - x_arc.^2 ) + r_throat;
        [~,index] = min(throat_arc);
    r1(ceil(x2_throat/dx)-index+1 : ceil(x2_throat/dx) + length(x_arc) - index) = throat_arc; 

% Define outer contour
r2 = r1 + thickness;

% Coolant Channels
n_pipe = 32; % number of channels
gap_pipe = 0.125*0.0254; % Gap between channels
w_pipe = 2*pi*r2/n_pipe - gap_pipe; % m - coolant channel width
h_pipe = 0.25*0.0254/2; % m - coolant channel height

% figure(1);
% plot(x,r1,x,r2);
% axis equal;
% title("Engine Contour (m)");

%% Exhaust Flow Properties



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
