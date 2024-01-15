
initial_guess = [0.5, 0, 0.8]; % [mdot_gg, degree_reaction, flow_coeff]

% Define the bounds for the parameters
lb = [0, 0, 0]; % Lower bounds
ub = [Inf, 1, 2]; % Upper bounds

% Define the options for the optimization
options = optimset('fmincon');
options.Algorithm = 'sqp'; % Sequential Quadratic Programming

% Perform the optimization using fmincon
objective_function = @(params) objective_function_rotor(params, shaft_power, shaft_speed, c_star_gg, p_amb, p_gg, gamma_gg, R_gg, Tt_gg);
constraint_function = @(params) triangle_constraints(params, shaft_power, shaft_speed, c_star_gg, p_amb, p_gg, gamma_gg, R_gg, Tt_gg);
[opt_params, fval] = fmincon(objective_function, initial_guess, [], [], [], [], lb, ub, constraint_function, options);
opt_params

mdot_gg = opt_params(1);
degree_reaction = opt_params(2);
flow_coeff = opt_params(3);
[r_rotor, stage_loading_coeff] = turbine_calculations(mdot_gg, degree_reaction, flow_coeff, shaft_power, shaft_speed, c_star_gg, p_amb, p_gg, gamma_gg, R_gg, Tt_gg);
2*r_rotor/0.0254
stage_loading_coeff


function [r_rotor, stage_loading_coeff] = turbine_calculations(mdot_gg, degree_reaction, flow_coeff, shaft_power, shaft_speed, c_star_gg, p_amb, p_gg, gamma_gg, R_gg, Tt_gg)
	% Nozzle Plate (np)
	nozzle_number = 16;
	A_throat_np = mdot_gg*c_star_gg/p_gg/nozzle_number; % m2
	d_throat_np = 2*sqrt(A_throat_np/pi); % m
	
	p_back_turbine = p_amb; % Pa - assuming the backpressure is ambient, I have no idea if this is valid
	M_np = sqrt( 2/(gamma_gg-1) * ((p_back_turbine/p_gg)^((1-gamma_gg)/gamma_gg) - 1) ); % Absolute mach post nozzle, pre rotor
	
	A_exit_np = A_throat_np/M_np * ( (2+(gamma_gg-1)*M_np^2)/(gamma_gg+1) ) ^ (0.5*(gamma_gg+1)/(gamma_gg-1)); % m2
	d_exit_np = 2*sqrt(A_throat_np/pi); % m - more of a virtual diameter because the nozzle is so tilted
		
	alpha_np = 70*pi/180; % empirical recommendation by Farokhi, page 691
	T_exit_np = Tt_gg/(1+(gamma_gg-1)/2*M_np^2); % K temperature at the virtual exit; reality may be more complicated
	v_flow_np = M_np*sqrt(gamma_gg*R_gg*T_exit_np)*[sin(alpha_np), cos(alpha_np)]; % [m/s, m/s] - tangential and axial velocity components
	
	% Rotor
    isentropic_velocity_ratio = 0.35; %  A single-row impulse stage delivers best performance at velocity ratios between 0.30 and 0.40 (NASA 1974)
    r_rotor = isentropic_velocity_ratio*norm(v_flow_np)/shaft_speed; % m
	% midline_speed_rotor = v_flow_np(2)/flow_coeff; % m/s - U in Matt & Ohain
		
	stage_loading_coeff = shaft_power/(mdot_gg*midline_speed_rotor^2); % psi in Mattingly & Ohain
	v_flow_rotor = [shaft_power/(mdot_gg*midline_speed_rotor)+v_flow_np(1), v_flow_np(2)]; % [m/s, m/s] - assuming axial velocity remains constant
	r_rotor = midline_speed_rotor/shaft_speed; % m
	alpha_rotor = atan( norm(v_flow_np)/norm(v_flow_rotor)*(stage_loading_coeff/flow_coeff - tan(alpha_np)) ); % rad
	
	beta_np = atan(tan(alpha_np) - 1/flow_coeff); % rad
	beta_rotor = atan(tan(beta_np) - 2*degree_reaction/flow_coeff); % rad
end

function r_rotor = objective_function_rotor(params, shaft_power, shaft_speed, c_star_gg, p_amb, p_gg, gamma_gg, R_gg, Tt_gg)
    mdot_gg = params(1);
    degree_reaction = params(2);
    flow_coeff = params(3);

    [r_rotor,~] = turbine_calculations(mdot_gg, degree_reaction, flow_coeff, shaft_power, shaft_speed, c_star_gg, p_amb, p_gg, gamma_gg, R_gg, Tt_gg);
end

function [c, ceq] = triangle_constraints(params, shaft_power, shaft_speed, c_star_gg, p_amb, p_gg, gamma_gg, R_gg, Tt_gg)
    mdot_gg = params(1);
    degree_reaction = params(2);
    flow_coeff = params(3);
    
	[~,stage_loading_coeff] = turbine_calculations(mdot_gg, degree_reaction, flow_coeff, shaft_power, shaft_speed, c_star_gg, p_amb, p_gg, gamma_gg, R_gg, Tt_gg);
    
    % Constraints to keep the flow_coeff and stage_loading_coeff within the triangle
    % c(i) < 0 means the i-th constraint is satisfied
	A = [0.3,1];
	B = [1.4,1];
	C = [0.6,4.5];
    c = [line_eq(flow_coeff, stage_loading_coeff, A, B),
         line_eq(flow_coeff, stage_loading_coeff, B, C),
         line_eq(flow_coeff, stage_loading_coeff, C, A)];
    
    % No equality constraints
    ceq = [];
end

function val = line_eq(phi, psi, p1, p2)
    % Calculate the value for the line equation given two points p1 and p2
    % and the current flow_coeff (phi) and stage_loading_coeff (psi)
    val = (p2(2) - p1(2)) * phi - (p2(1) - p1(1)) * psi + p2(1)*p1(2) - p2(2)*p1(1);
end




















