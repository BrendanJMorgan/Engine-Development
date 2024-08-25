%% Impeller Contours - pg 2.29-2.32 in pump handbook

%% The eye
% Eye radius initial guessing
lower_bound = 1.001*r_eye_inner; % m - eye radius must be larger than its inner radius to have physical meaning
upper_bound = 10*r_eye_inner; % m - if eye radius is somehow above this you have bigger problems
r_guess = (lower_bound + upper_bound) / 2;

% Using fzero with bounded search
options = optimset('TolX',1e-12,'TolFun',1e-12);
[r_eye_impeller, fval, exitflag, output] = fzero(@(r_eye_impeller) eye_solver(r_eye_impeller, vdot_pump, shaft_speed, eye_flow_coeff, r_eye_inner), [lower_bound, upper_bound], options);

if exitflag == 0
    error('Pump eye radius solver failed');
end

A_eye = pi*(r_eye_impeller^2 - r_eye_inner^2); % m2 - inlet area of the impeller eye
eye_flow_coeff = vdot_pump/A_eye / (shaft_speed*r_eye_impeller); % unitless - in case solver didn't exactly match input

%% The exit radius r2 (or diameter D2)
r_exit_impeller = 1 / shaft_speed * sqrt (g*head_pump / head_coeff);

%% The exit width b2
w_exit = vdot_pump / (2*pi*shaft_speed*r_exit_impeller^2*outlet_flow_coeff*blockage);

%% Hub and shroud profiles

r_min = 0.5*r_eye_impeller; % m - minimum allowable radius of curvature
bezier_params_initial = [0.01, 0.01, 0.01]; % starting values
cost_function = @(bezier_params) cost_function_impeller(bezier_params, r_eye_impeller, r_exit_impeller, w_exit, impeller_thickness, r_eye_inner, r_min);
bezier_params_optimal = fminsearch(cost_function, bezier_params_initial);

bezier_vert_optimal = bezier_params_optimal(1);
bezier_horiz_optimal = bezier_params_optimal(2);
impeller_height_optimal = bezier_params_optimal(3);

[~, ~, min_shroud_curvature, ~, shroud_curve, impeller_curve, control_points] = compute_curves(bezier_vert_optimal, bezier_horiz_optimal, r_eye_impeller, r_exit_impeller, w_exit, impeller_thickness, impeller_height_optimal, r_eye_inner);

min_radius = 1/max(min_shroud_curvature);
if min_radius < 0.99*r_min    
	fprintf("\nPump shroud curvature reaches %.2f mm radius. Need above ~%.2f mm (half of eye radius) to avoid cavitation.", min_radius*1000, r_min*1000)
end

%% Functions
function f = eye_solver(r_eye_impeller, vdot_pump, shaft_speed, eye_flow_coeff, r_eye_inner)
    f = r_eye_impeller - (vdot_pump / (pi*shaft_speed*eye_flow_coeff*(1-r_eye_inner^2/r_eye_impeller^2)) ) ^ (1/3);
end

% Function to compute slopes and minimum curvature for shroud and impeller curves in a pump
function [slope_shroud, slope_impeller, min_shroud_curvature, flatness_impeller, shroud_curve, impeller_curve, control_points] = compute_curves(bezier_vert, bezier_horiz, r_eye_impeller, r_exit_impeller, w_exit, impeller_thickness, impeller_height, r_eye_inner)

    % Generate the shroud curve using bezier control points
    control_points = [r_eye_impeller, impeller_height; r_eye_impeller, impeller_height - bezier_vert; r_exit_impeller - bezier_horiz, impeller_thickness + w_exit; r_exit_impeller, impeller_thickness + w_exit];
    shroud_curve = bezier(control_points); % [m,m]

    % Compute the flow areas at eye and exit
    A_pump_eye = pi * (r_eye_impeller^2 - r_eye_inner^2); % m2
    A_pump_exit = 2 * pi * r_exit_impeller * w_exit; % m2

    % Generate area array from eye to exit along the shroud curve
    A_pump_flow = linspace(A_pump_eye, A_pump_exit, length(shroud_curve))'; % m2

    % Compute the angles of the normals at each point on the shroud curve
    normal_angles = atan2(gradient(shroud_curve(:,2)), gradient(shroud_curve(:,1))) + pi/2; % rad
    normal_angles(normal_angles > 2 * pi) = normal_angles(normal_angles > 2 * pi) - 2 * pi; % rad

    % Compute the control wise gap between shroud and impeller curves at each point
    crosswise_gap = sec(normal_angles) .* (shroud_curve(:,1) - sqrt(pi^2 * shroud_curve(:,1).^2 - pi * A_pump_flow .* cos(normal_angles)) / pi); % m

    % Generate impeller curve based on the shroud curve and the crosswise gap
    impeller_curve = [shroud_curve(:,1) - crosswise_gap .* cos(normal_angles), shroud_curve(:,2) - crosswise_gap .* sin(normal_angles)]; % [m,m]

    % Calculate the slope of the shroud and impeller curves at their ends
    slope_shroud = (shroud_curve(end, 2) - shroud_curve(end - 1, 2)) / (shroud_curve(end, 1) - shroud_curve(end - 1, 1));
    slope_impeller = (impeller_curve(end, 2) - impeller_curve(end - 1, 2)) / (impeller_curve(end, 1) - impeller_curve(end - 1, 1));

    % Calculate flatness, i.e. the max height dip below the impeller exit plane
    flatness_impeller = min(impeller_curve(:,2)) - impeller_curve(end, 2); % m

    % Compute the minimum curvature of the shroud curve
    s = [0; cumsum(sqrt(diff(shroud_curve(:, 1)).^2 + diff(shroud_curve(:, 2)).^2))];
    min_shroud_curvature = abs(gradient(shroud_curve(:, 1), s) .* gradient(gradient(shroud_curve(:, 2), s), s) - gradient(shroud_curve(:, 2), s) .* gradient(gradient(shroud_curve(:, 1), s), s)) ./ (gradient(shroud_curve(:, 1), s).^2 + gradient(shroud_curve(:, 2), s).^2).^(3/2);
end

function cost = cost_function_impeller(bezier_params, r_eye_impeller, r_exit_impeller, w_exit, impeller_thickness, r_eye_inner, r_min)
    bezier_vert = bezier_params(1);
    bezier_horiz = bezier_params(2);
    impeller_height = bezier_params(3);
    [slope_shroud, slope_impeller, min_shroud_curvature, flatness_impeller] = compute_curves(bezier_vert, bezier_horiz, r_eye_impeller, r_exit_impeller, w_exit, impeller_thickness, impeller_height, r_eye_inner);
    
    if flatness_impeller > 0
        flatness_impeller = 0;
    end

    % The cost is the sum of the absolute slopes (want them to be 0 for a horizontal line) and the curvature
    cost = 100*abs(slope_shroud) + 100*abs(slope_impeller) + 1000*abs(1/max(min_shroud_curvature)-r_min) + 1000*abs(flatness_impeller);
end
