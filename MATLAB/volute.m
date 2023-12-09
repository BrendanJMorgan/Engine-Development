%% Pump Handbook page 2.68

r_tongue = 1.06*r_exit; % m - from pump handbook but somewhat arbitrary
t_tongue = 0*0.70*blade_thickness; % m - from pump handbook but somewhat arbitrary
A_throat = ( vdot_fuel/(3.8*r_exit*v_exit_tang) * ( 1 + sqrt(1+15.2*r_exit*v_exit_tang*(r_tongue+t_tongue)/vdot_fuel) ) )^2; % m2
r_throat = r_exit + 0.5*sqrt(A_throat) + t_tongue; % m

theta_volute = linspace(0, 2*pi, 1000); % rad
A_volute = ( vdot_fuel.*theta_volute./(7.6*pi*r_exit*v_exit_tang) .* ( 1 + sqrt(1+30.4*pi*r_exit*r_tongue*v_exit_tang./(vdot_fuel*theta_volute)) ) ).^2; % m2
A_volute(isnan(A_volute)) = 0;
r_volute = r_tongue + 0.5*sqrt(A_volute); % m

v_throat = vdot_fuel./A_throat
figure(1)
polarplot(theta_volute, r_volute)
