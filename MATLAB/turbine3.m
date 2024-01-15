mdot_gg = 0.1 % kg/s - placeholder
r_pitchline = 5/0.0254; % m - constant radius for rotors/stators/nozzles
reaction_degree = 0; 

p5 = p_amb; % Pa
C5 = [0, Cz]; % [m/s, m/s] - zero swirl at the exit is ideal

% Rotor 2 (station 4 to 5)
rotor2_power = 0.5*shaft_power; % W - 50% to share torque evenly between stages
C4 = [C5(1)-rotor2_power/(mdot_gg*r_pitchline*shaft_speed),Cz]; % [m/s, m/s]

% Stator (station 3 to 4)
C3 = [(1-reaction_degree)*2*r_pitchline_shaft_speed-C4(1),Cz]; % [m/s, m/s]

%% Rotor 1 (station 2 to 3)
rotor1_power = 0.5*shaft_power; % W - 50% to share torque evenly between stages
C2 = [C3(1)-rotor1_power/(mdot_gg*r_pitchline*shaft_speed),Cz]; % [m/s, m/s]


%% Nozzle Plate (station 1 to 2)

nozzle_number = 64;
A_throat_np = mdot_gg*c_star_gg/p_gg/nozzle_number; % m2
d_throat_np = 2*sqrt(A_throat_np/pi); % m

p_back_turbine = p_amb; % Pa - assuming the backpressure is ambient, I have no idea if this is valid
M_np = sqrt( 2/(gamma_gg-1) * ((p_back_turbine/p_gg)^((1-gamma_gg)/gamma_gg) - 1) ); % Absolute mach post nozzle, pre rotor

A_exit_np = A_throat_np/M_np * ( (2+(gamma_gg-1)*M_np^2)/(gamma_gg+1) ) ^ (0.5*(gamma_gg+1)/(gamma_gg-1)); % m2
d_exit_np = 2*sqrt(A_throat_np/pi); % m - more of a virtual diameter because the nozzle is so tilted
	
alpha_np = 55*pi/180; % rad - empirical recommendation by Farokhi is 70 deg, page 691
T_exit_np = Tt_gg/(1+(gamma_gg-1)/2*M_np^2); % K - temperature at the virtual exit; reality may be more complicated
v_flow_np = M_np*sqrt(gamma_gg*R_gg*T_exit_np)*[sin(alpha_np), cos(alpha_np)]; % [m/s, m/s] - tangential and axial velocity components

dens_np = (p_back_turbine/p_gg)^(1/gamma_gg); % kg/m3



% 1-R = (Ct2+Ct1)/(2U)    