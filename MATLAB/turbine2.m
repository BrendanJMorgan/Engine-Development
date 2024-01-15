% Nozzle Plate (np)
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

% Rotor
velocity_ratio = 0.25; % NASA SP-8110 page 83
if velocity_ratio >= 0.45
    n_stages = 1; % Reaction staging
elseif (0.20 <= velocity_ratio) && (velocity_ratio <= 0.30)
    n_stages = 2; % velocity-compounded impulse staging+
elseif (0.30 <= velocity_ratio) && (velocity_ratio <= 0.40)
    n_stages = 1; % impulse staging
elseif velocity_ratio <= 0.15
    n_stages = 3;
end
    

r_rotor = velocity_ratio*norm(v_flow_np)/shaft_speed/sqrt(n_stages); % m
midline_speed_rotor = shaft_speed*r_rotor; % m/s - U in Matt & Ohain
	
stage_loading_coeff = shaft_power/(mdot_gg*midline_speed_rotor^2); % psi in Mattingly & Ohain
flow_coeff = v_flow_np(2)/midline_speed_rotor; 
v_flow_rotor = [shaft_power/(mdot_gg*midline_speed_rotor)+v_flow_np(1), v_flow_np(2)]; % [m/s, m/s] - assuming axial velocity remains constant

blade_length_rotor = mdot_gg/(2*pi*r_rotor*dens_np*v_flow_np(2)); % m - this is just a rearranged continuity equation

thrust_turbine = mdot_gg*v_flow_np(2); % N

thrust_turbine/4.44822
2*(r_rotor+blade_length_rotor/2)/0.0254
blade_length_rotor/0.0254

