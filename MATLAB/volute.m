%% Pump Handbook page 2.68

r_tongue = 1.06*r_exit_impeller; % m - from pump handbook but somewhat arbitrary
t_tongue = 0.70*blade_thickness; % m - from pump handbook but somewhat arbitrary
A_throat = ( vdot_pump/(3.8*r_exit_impeller*v_exit_tang) * ( 1 + sqrt(1+15.2*r_exit_impeller*v_exit_tang*(r_tongue+t_tongue)/vdot_pump) ) )^2; % m2
r_throat = r_exit_impeller + 0.5*sqrt(A_throat) + t_tongue; % m

theta = linspace(0, 2*pi, 1000)'; % rad
A_volute = ( vdot_pump.*theta./(7.6*pi*r_exit_impeller*v_exit_tang) .* ( 1 + sqrt(1+30.4*pi*r_exit_impeller*r_tongue*v_exit_tang./(vdot_pump*theta)) ) ).^2; % m2
A_volute(isnan(A_volute)) = 0;

h_volute = sqrt(0.25*t_tongue^2+A_volute(end)) - 0.5*t_tongue; % m - height of volute gap - set so the outlet is a square
r_volute = r_tongue + A_volute/h_volute; % m - radius of volute wall
OD_volute_outlet = 2*sqrt(A_volute(end)/pi); % m - diameter of volute outlet hole
r_volute_outlet = r_volute(end) - h_volute(end)/2;

v_throat = vdot_pump./A_throat; % m/s - fluid exit velocity into plumbing

volute_curve = r_volute.*[cos(clock*theta), sin(clock*theta)]; % [m, m]

