%% Impeller Contours
%% The eye
converged = false;
eye_flow_coeff = 0.25; %  % phi_e in pump handbook; pg. 2.29; 0.2-0.3 for impellers, ~0.1 or less for inducers UPDATE TO DERIVE FROM CAVITATION COEFFICIENT TAU

% Eye radius initial guessing
lower_bound = 1.01*r_shaft;  % m - eye radius must be larger than shaft radius to have physical meaning
upper_bound = 10*r_shaft;    % m - if eye radius is somehow above this you have bigger problems
r_guess = (lower_bound + upper_bound) / 2;

% Using fzero with bounded search
options = optimset('TolX',1e-12,'TolFun',1e-12);
[r_eye, fval, exitflag, output] = fzero(@(r_eye) equation_to_solve(r_eye, vdot_fuel, shaft_speed, eye_flow_coeff, r_shaft), [lower_bound, upper_bound], options);

if exitflag == 0
    error('Pump eye radius solver failed');
end

%% The exit radius r2 (or diameter D2)
r_exit = 1 / shaft_speed * sqrt (g*head_fuel / head_coeff_fuel);

%% The exit width b2
w_exit = vdot_fuel / (2*pi*shaft_speed*r_exit^2*outlet_flow_coeff*blockage);

%% Hub and shroud profiles

r_min = 0.5*r_eye; % m - minimum allowable radius of curvature
bezier_params_initial = [0.008, 0.008]; % starting values
cost_function = @(bezier_params) cost_function_impeller(bezier_params, r_eye, r_exit, w_exit, impeller_thickness, impeller_height, r_shaft, r_min);
bezier_params_optimal = fminsearch(cost_function, bezier_params_initial);

bezier_vert_optimal = bezier_params_optimal(1);
bezier_horiz_optimal = bezier_params_optimal(2);

shroud_points = [r_eye,impeller_height; r_eye, impeller_height-bezier_vert_optimal; r_exit-bezier_horiz_optimal,impeller_thickness+w_exit; r_exit, impeller_thickness+w_exit]; % [m,m] - bezier control points of the shroud
shroud_curve = bezier(shroud_points); % [m,m] - lower curved edge of revolution of the pump shroud

A_eye = pi*(r_eye^2-r_shaft^2); % m2 - area of the annular eye entrance
A_exit = 2*pi*r_exit*w_exit; % m2 - area of the cylindrical exit THIS VARIABLE IS OVERLOADED

A_pump_flow = linspace(A_eye, A_exit, length(shroud_curve))'; % m2 - cross sectional area perpendicular to flow

normal_angles = atan2(gradient(shroud_curve(:,2)), gradient(shroud_curve(:,1))) + pi/2; % rad
normal_angles(normal_angles > 2*pi) = normal_angles(normal_angles > 2*pi) - 2*pi;

crosswise_gap =  sec(normal_angles) .* ( shroud_curve(:,1) - sqrt( pi^2 * shroud_curve(:,1).^2 - pi*A_pump_flow.*cos(normal_angles) ) / pi); % m - gap between shroud curve and impeller curve at each point

impeller_curve = [shroud_curve(:,1) - crosswise_gap.*cos(normal_angles), shroud_curve(:,2) - crosswise_gap.*sin(normal_angles)]; % [m,m] - upper curved edge of revolution of the pump impeller

% Check that shroud curve does not kink too much; ie that it's radius of curvature is not below r_min
s = [0; cumsum(sqrt(diff(shroud_curve(:,1)).^2 + diff(shroud_curve(:,2)).^2))];
curvature = abs(gradient(shroud_curve(:,1), s).*gradient(gradient(shroud_curve(:,2), s), s) - gradient(shroud_curve(:,2), s).*gradient(gradient(shroud_curve(:,1), s), s)) ./ (gradient(shroud_curve(:,1), s).^2 + gradient(shroud_curve(:,2), s).^2).^(3/2);
[max_curvature, idx] = max(curvature);
min_radius = 1/max_curvature;
if min_radius < 0.99*r_min    
    fprintf("Pump shroud curvature reaches %.2f mm radius. Need above ~%.2f mm (half of eye radius) to avoid cavitation.", min_radius*1000, r_min*1000)
end



%% Functions 
function f = equation_to_solve(r_eye, vdot_fuel, shaft_speed, eye_flow_coeff, r_shaft)
    f = r_eye - (vdot_fuel / (pi*shaft_speed*eye_flow_coeff*(1-r_shaft^2/r_eye^2)) ) ^ (1/3);
end

% Function to compute slopes and minimum curvature for shroud and impeller curves in a pump
function [slope_shroud, slope_impeller, min_shroud_curvature] = compute_curves(bezier_vert, bezier_horiz, r_eye, r_exit, w_exit, impeller_thickness, impeller_height, r_shaft)

    % Generate the shroud curve using bezier control points
    control_points = [r_eye, impeller_height; r_eye, impeller_height - bezier_vert; r_exit - bezier_horiz, impeller_thickness + w_exit; r_exit, impeller_thickness + w_exit];
    shroud_curve = bezier(control_points);

    % Compute the flow areas at eye and exit
    A_eye = pi * (r_eye^2 - r_shaft^2);
    A_exit = 2 * pi * r_exit * w_exit;

    % Generate area array from eye to exit along the shroud curve
    A_pump_flow = linspace(A_eye, A_exit, length(shroud_curve))';

    % Compute the angles of the normals at each point on the shroud curve
    normal_angles = atan2(gradient(shroud_curve(:,2)), gradient(shroud_curve(:,1))) + pi/2;
    normal_angles(normal_angles > 2 * pi) = normal_angles(normal_angles > 2 * pi) - 2 * pi;

    % Compute the crosswise gap between shroud and impeller curves at each point
    crosswise_gap = sec(normal_angles) .* (shroud_curve(:,1) - sqrt(pi^2 * shroud_curve(:,1).^2 - pi * A_pump_flow .* cos(normal_angles)) / pi);

    % Generate impeller curve based on the shroud curve and the crosswise gap
    impeller_curve = [shroud_curve(:,1) - crosswise_gap .* cos(normal_angles), shroud_curve(:,2) - crosswise_gap .* sin(normal_angles)];

    % Calculate the slope of the shroud and impeller curves at their ends
    slope_shroud = (shroud_curve(end, 2) - shroud_curve(end - 1, 2)) / (shroud_curve(end, 1) - shroud_curve(end - 1, 1));
    slope_impeller = (impeller_curve(end, 2) - impeller_curve(end - 1, 2)) / (impeller_curve(end, 1) - impeller_curve(end - 1, 1));

    % Compute the minimum curvature of the shroud curve
    s = [0; cumsum(sqrt(diff(shroud_curve(:, 1)).^2 + diff(shroud_curve(:, 2)).^2))];
    min_shroud_curvature = abs(gradient(shroud_curve(:, 1), s) .* gradient(gradient(shroud_curve(:, 2), s), s) - gradient(shroud_curve(:, 2), s) .* gradient(gradient(shroud_curve(:, 1), s), s)) ./ (gradient(shroud_curve(:, 1), s).^2 + gradient(shroud_curve(:, 2), s).^2).^(3/2);
end

function cost = cost_function_impeller(bezier_params, r_eye, r_exit, w_exit, impeller_thickness, impeller_height, r_shaft, r_min)
    bezier_vert = bezier_params(1);
    bezier_horiz = bezier_params(2);
    [slope_shroud, slope_impeller, min_shroud_curvature] = compute_curves(bezier_vert, bezier_horiz, r_eye, r_exit, w_exit, impeller_thickness, impeller_height, r_shaft);

    % The cost is the sum of the absolute slopes (we want them to be 0 for a horizontal line) and the curvature
    cost = 100*abs(slope_shroud) + abs(slope_impeller) + 100*abs(1/max(min_shroud_curvature)-r_min);
end
