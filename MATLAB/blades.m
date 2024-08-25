%% Impeller Blades

%% Leading and Trailing Edges
r_inlet = r_eye_impeller; % m - radial location of blade inlets - typically vertically aligned with eye radius - this assumption can be changed

%% Surface of Revolution
rms_curve = 0.5 * (shroud_curve + impeller_curve); % [m,m] - halfway between shroud and impeller

%% Inlet Blade Angles
inlet_index = find(rms_curve(:,1) > r_inlet, 1);
inlet_gap = norm( shroud_curve(inlet_index, :) -  impeller_curve(inlet_index, :) );
u_inlet = shaft_speed*r_inlet;				    % m/s - rotational speed at the inlet (U_1 in pump handbook)
v_inlet = vdot_pump/(2*pi*r_inlet*inlet_gap);   % m/s - fluid velocity at the inlet, assumed to be entirely radial. Can be found with CFD instead for more accuracy (V_1 in pump handbook)
blade_angle_inlet = atan2(u_inlet,v_inlet);     % rad - inlet angle of blades relative to tangential azimuth  (beta_1,f in pump handbook)

%% Outlet Blade Angles
slip_factor = 0.15; % Buseman form "typically between 0.1 and 0.2 ... for frictionless flow ... with log-spiral blades ... and a 2D, radial-flow geometry with parallel hub and shroud"
v_blade_exit = shaft_speed*r_exit_impeller;	% m/s - blade tangential velocity (U_2 in pump handbook)
v_slip = slip_factor*v_blade_exit;			% m/s - slip velocity (V_s in pump handbook)
v_merid_exit = vdot_pump/(2*pi*r_exit_impeller*w_exit);	% m/s - meriodonal velocity (V_m,2 in pump handbook)                                                       

hydraulic_efficiency = 1 - 0.071 / vdot_pump^0.25; % Jekat's Empirical Formula - wraps up all the various efficiencies into one formula, valid for all specific speeds

v_exit_tang = g*head_pump / (hydraulic_efficiency*v_blade_exit);			    % m/s - this is Euler's pump equation rearranged (V_theta,2 in pump handbook)                                                                                                                                                                                                                                                                                                                                                                                                                              
blade_angle_outlet = atan(v_merid_exit / (v_blade_exit-v_exit_tang-v_slip));    % rad - angle between blade tip azimuth and local tangential azimuth (beta_2 in pump handbook)
if blade_angle_outlet <= 0
	blade_angle_outlet*180/pi
	error("Cannot converge to solution with valid blade outlet angle")
end

%% Individual Blade Curve
params_initial = [0.01,0.01,-clock]; % starting values; control1, control2
cost_function = @(params) cost_function_blade(params, r_inlet, r_exit_impeller, blade_angle_inlet, blade_angle_outlet, clock);
blade_params_optimal = fminsearch(cost_function, params_initial);

control1 = blade_params_optimal(1); % m - length between inlet point and inlet control point
control2 = blade_params_optimal(2); % m - length between outlet point and outlet control point
sweep = blade_params_optimal(3);    % rad
[blade_curve, blade_control_points, ~, ~] = compute_blade_curve(control1, control2, sweep, r_inlet, r_exit_impeller, blade_angle_inlet, blade_angle_outlet, clock);

%% Number of Blades
blade_z = interp1(rms_curve(:,1), rms_curve(:,2), blade_curve(:,1));                        % m - vertical position component of blade curve
blade_arc_length = sum( vecnorm( diff([blade_curve(:,1), blade_curve(:,2), blade_z])' ) );  % m - arc length of one individual blade, measure at rms line
solidity = interp1([0, 0.4, 3], [1.8, 1.8, 1], specific_speed);					            % solidity is the optimal ratio of blade chord to blade spacing. Pump handbook page 2.36 (sigma)
blade_number = round(solidity * 2*pi*r_exit_impeller / blade_arc_length);	                % number of blades, rounded to the nearest even number (n)

%% Compute Blockage (pump handbook page 2.63)
blade_thickness = 0.04*r_exit_impeller; % m
boundary_layer_thickness = 0.002*blade_arc_length; % m
meridional_length = sum( vecnorm( diff([blade_curve(:,1), blade_z])' ) ); % m - arc length of blade projected onto the meridional plane
new_blockage = 1 - (2*boundary_layer_thickness+blade_thickness)/(0.0254*2*pi * meridional_length/blade_arc_length); 

%% Helper Functions
% Function to compute slopes and minimum curvature for blade curve in a pump
function [blade_curve, control_points, control_error1, control_error2] = compute_blade_curve(control1, control2, sweep, radius1, radius2, polar1, polar2, clock)
	% clock = 1; % 1 for counterclockwise and -1 for clockwise
	inlet_point = [radius1,0];																% [m,m] - blade inner edge starting location
	outlet_point = radius2*[cos(sweep), sin(sweep)];										% [m,m] - blade outer edge ending location
	inlet_azimuth = atan2(inlet_point(2), inlet_point(1)) - clock*(pi/2 - polar1);			% rad - blade inner edge absolute direction
	outlet_azimuth = atan2(outlet_point(2), outlet_point(1)) - clock*(pi/2 - polar2);		% rad - blade outer edge absolute direction
	inlet_control = control1*[cos(inlet_azimuth),sin(inlet_azimuth)] + inlet_point;			% [m,m] - control point for the inlet edge
	outlet_control = -control2*[cos(outlet_azimuth),sin(outlet_azimuth)] + outlet_point;	% [m,m] - control point for the outlet edge

	% Generate the shroud curve using four bezier control points
    control_points = [inlet_point;inlet_control;outlet_control;outlet_point];
    blade_curve = bezier(control_points); % [m,m]

	control_error1 = atan2(blade_curve(2,2)-inlet_point(2), blade_curve(2,1)-inlet_point(1)) - inlet_azimuth;
	control_error2 = atan2(blade_curve(end-1,2)-outlet_point(2), blade_curve(end-1,1)-outlet_point(1)) - outlet_azimuth;
end

function cost = cost_function_blade(params, radius1, radius2, polar1, polar2, clock)
    control1 = params(1);   % m - length between inlet point and inlet control point
	control2 = params(2);   % m - length between outlet point and outlet control point
    sweep = params(3);      % rad
    [blade_curve, ~, control_error1, control_error2] = compute_blade_curve(control1, control2, sweep, radius1, radius2, polar1, polar2, clock);

	blade_arc_length = sum(sqrt(sum(diff(blade_curve).^2, 2))); % m - arc length of one individual blade
    s = [0; cumsum(sqrt(diff(blade_curve(:, 1)).^2 + diff(blade_curve(:, 2)).^2))];
    blade_curvature = abs(gradient(blade_curve(:, 1), s) .* gradient(gradient(blade_curve(:, 2), s), s) - gradient(blade_curve(:, 2), s) ...
		.* gradient(gradient(blade_curve(:, 1), s), s)) ./ (gradient(blade_curve(:, 1), s).^2 + gradient(blade_curve(:, 2), s).^2).^(3/2); % 1/m
	max_blade_curvature = max(blade_curvature); % 1/m
    
	cost = 1*abs(max_blade_curvature) + 1*abs(sweep) + 0.1*abs(blade_arc_length) + abs(control_error1) + abs(control_error2);
end







