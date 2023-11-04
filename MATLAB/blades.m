%% Impeller Blades

%% Leading and Trailing Edges
r_inlet = r_eye; % m - this assumption can be changed

%% Surface of Revolution
rms_curve =  [shroud_curve(:,1) - 0.5*crosswise_gap.*cos(normal_angles), shroud_curve(:,2) - 0.5*crosswise_gap.*sin(normal_angles)]; % [m,m] - halfway between shroud and impeller

%% Inlet Blade Angles
inlet_gap = crosswise_gap(find(rms_curve > r_inlet, 1));
u_inlet = shaft_speed*r_inlet; % m/s - rotational speed at the inlet
v_inlet = vdot_fuel/(2*pi*r_inlet*inlet_gap); % m/s - fluid velocity at the inlet, assumed to be entirely radial. Can be found with CFD instead for more accuracy
blade_angle_inlet = atan(u_inlet/v_inlet); % rad - inlet angle of blades relative to straight radial direction

%% Outlet Blade Angles
slip_factor = 0.15; % Buseman form "typically between 0.1 and 0.2 ... for frictionless flow through impellers with logarithmic-spiral blades ... and a 2D, radial-flow geometry with parallel hub and shroud"
v_blade_exit = shaft_speed*r_exit; % m/s
v_slip = slip_factor*v_blade_exit; % m/s - slip velocity
v_merid_exit = vdot_fuel/(2*pi*r_exit*w_exit); % m/s                                                                                

hydraulic_efficiency = 1 - 0.071 / vdot_fuel^0.25; % Jekat's Empirical Formula

v_flow_exit_tangential = g*head_fuel / (hydraulic_efficiency*v_blade_exit); % m/s - this is Euler's pump equation rearranged                                                                                                                                                                                                                                                                                                                                                                                                                                     
blade_angle_outlet = atan(v_merid_exit / (v_blade_exit-v_flow_exit_tangential - v_slip)); % rad

%% Individual Blade Curve 
params_initial = [0.01,0.01,pi/4]; % starting values; control1, control2
cost_function = @(params) cost_function_blade(params, r_inlet, r_exit, blade_angle_inlet, blade_angle_outlet);
blade_params_optimal = fminsearch(cost_function, params_initial);

control1 = blade_params_optimal(1);
control2 = blade_params_optimal(2);
sweep = blade_params_optimal(3);
[blade_curve, blade_control_points, ~, ~] = compute_blade_curve(control1, control2, sweep, r_inlet, r_exit, blade_angle_inlet, blade_angle_outlet);

%% Number of Blades
blade_arc_length = sum(sqrt(sum(diff(blade_curve).^2, 2))); % m - arc length of one individual blade
solidity = interp1([0, 0.4, 3], [1.8, 1.8, 1], specific_speed_fuel); % solidity is the optimal ratio of blade chord to blade spacing. Pump handbook page 2.36
blade_number = 2*round(0.5*(solidity * 2*pi*r_exit / blade_arc_length)); % number of blade, rounded to the nearest even number

%% Helper Functions
% Function to compute slopes and minimum curvature for blade curve in a pump
function [blade_curve, control_points, control_error1, control_error2] = compute_blade_curve(control1, control2, sweep, radius1, radius2, polar1, polar2)

	inlet_point = [radius1,0]; % [m,m] - blade inner edge starting location
	outlet_point = radius2*[cos(sweep), sin(sweep)]; % [m,m] - blade outer edge ending location
	inlet_azimuth = atan2(inlet_point(2), inlet_point(1)) + pi/2 - polar1; % rad - blade inner edge absolute direction
	outlet_azimuth = atan2(outlet_point(2), outlet_point(1)) - pi/2 - polar2; % rad - blade outer edge absolute direction
	inlet_control = control1*[cos(inlet_azimuth),sin(inlet_azimuth)] + inlet_point; % [m,m] - control point for the inlet edge
	outlet_control = control2*[cos(outlet_azimuth),sin(outlet_azimuth)] + outlet_point; % [m,m] - control point for the outlet edge

	% Generate the shroud curve using bezier control points
    control_points = [inlet_point;inlet_control;outlet_control;outlet_point];
    blade_curve = bezier(control_points); % [m,m]

	control_error1 = atan2(blade_curve(2,2)-inlet_point(2), blade_curve(2,1)-inlet_point(1)) - inlet_azimuth;
	control_error2 = atan2(blade_curve(end-1,2)-outlet_point(2), blade_curve(end-1,1)-outlet_point(1)) - outlet_azimuth;
end

function cost = cost_function_blade(params, radius1, radius2, polar1, polar2)
    control1 = params(1); % m 
	control2 = params(2); % m
    sweep = params(3); % rad
    [blade_curve, control_points, control_error1, control_error2] = compute_blade_curve(control1, control2, sweep, radius1, radius2, polar1, polar2);

	blade_arc_length = sum(sqrt(sum(diff(blade_curve).^2, 2))); % m
    s = [0; cumsum(sqrt(diff(blade_curve(:, 1)).^2 + diff(blade_curve(:, 2)).^2))];
    blade_curvature = abs(gradient(blade_curve(:, 1), s) .* gradient(gradient(blade_curve(:, 2), s), s) - gradient(blade_curve(:, 2), s) ...
		.* gradient(gradient(blade_curve(:, 1), s), s)) ./ (gradient(blade_curve(:, 1), s).^2 + gradient(blade_curve(:, 2), s).^2).^(3/2); % 1/m
	max_blade_curvature = max(blade_curvature); % 1/m

	cost = 1*abs(max_blade_curvature) + 1*abs(sweep) + 0.1*abs(blade_arc_length) + abs(control_error1) + abs(control_error2);
end






