%% Requirements

p_out_fuel = p_cool(1); % Pa - outlet pressure
p_in_fuel = p_amb; % Pa - inlet pressure
p_out_ox = pc; % Pa - outlet pressure
p_in_ox = p_amb; % Pa - inlet pressure

head_fuel = (p_out_fuel - p_in_fuel) / (density_fuel*g); % m

vdot_fuel = mdot_fuel_cc/density_fuel; % m3/s - Fuel Volumetric Flow Rate
vdot_ox = mdot_ox_cc/density_ox; % m3/s - Oxidizer Mass Flow Rate


%% Establishing the Pump Configuration
vapor_pressure_fuel = py.CoolProp.CoolProp.PropsSI("P", "T", T_amb, "Q", 0, "Ethanol");
NPSH_fuel = p_amb/(density_fuel*g) - vapor_pressure_fuel/(density_fuel*g); % Net positive suction head. Should be updated with plumbing water column and steady state inlet velocity
specific_speed_fuel = shaft_speed*sqrt(vdot_fuel) / (g*head_fuel)^0.75;
suction_specific_speed_fuel = shaft_speed * sqrt(vdot_fuel) / (NPSH_fuel*g)^0.75;
suction_specific_speed_imperial_fuel = suction_specific_speed_fuel*2733.016;

if specific_speed_fuel < 0.1 || specific_speed_fuel > 0.6
    fprintf("Pump Specific Speed is %g. Recommended is 0.1 to 0.6 for centrifugal geometries.", specific_speed_fuel)
end

%% Sizing the Pump
specific_speed_fuel = shaft_speed*sqrt(vdot_fuel)/(head_fuel*g)^0.75;
outlet_flow_coeff = 0.175*sqrt(specific_speed_fuel); % Phi_i in pump handbook

head_coeff_fuel = 0.4 / specific_speed_fuel^0.25; % Omega in pump handbook. Valid for specific speeds < 1 (which you should be anyways; centrifugal pumps should be between 0.1 and 0.6)

blockage = 0.85; % 1 is completely open, pump handbook says 0.85 is typical?

%% Impeller Contours

% The eye
r_eye = 1*0.0254; % m - initial guess of eye inlet radius
while true
    old_r_eye = r_eye;
    eye_flow_coeff = 0.25; %  % phi_e in pump handbook; pg. 2.29; 0.2-0.3 for impellers, ~0.1 or less for inducers UPDATE TO DERIVE FROM CAVITATION COEFFICIENT TAU
    r_eye = (vdot_fuel / (pi*shaft_speed*eye_flow_coeff*(1-r_shaft^2/r_eye^2)) ) ^ (1/3);
    if abs(r_eye - old_r_eye) < 0.01*0.0254
        break;
    end
end

% The exit radius r2 (or diameter D2)
r_exit = 1 / shaft_speed * sqrt (g*head_fuel / head_coeff_fuel);

% The exit width b2
w_exit = vdot_fuel / (2*pi*shaft_speed*r_exit^2*outlet_flow_coeff*blockage);

% Hub and shroud profiles

r_min = 0.5*r_eye; % m - minimum allowable radius of curvature
bezier_params_initial = [0.008, 0.008]; % starting values
options = optimset('Display', 'iter'); % display progress
cost_function_with_fixed_params = @(bezier_params) cost_function(bezier_params, r_eye, r_exit, w_exit, impeller_thickness, impeller_height, r_shaft, r_min);
bezier_params_optimal = fminsearch(cost_function_with_fixed_params, bezier_params_initial, options);

bezier_vert_optimal = bezier_params_optimal(1);
bezier_horiz_optimal = bezier_params_optimal(2);

shroud_points = [r_eye,impeller_height; r_eye, impeller_height-bezier_vert_optimal; r_exit-bezier_horiz_optimal,impeller_thickness+w_exit; r_exit, impeller_thickness+w_exit]; % [m,m] - bezier control points of the shroud
shroud_curve = bezier(shroud_points); % [m,m] - lower curved edge of revolution of the pump shroud

A_eye = pi*(r_eye^2-r_shaft^2); % m2 - area of the annular eye entrance
A_exit = 2*pi*r_exit*w_exit; % m2 - area of the cylindrical exit

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


figure(1)
line(shroud_curve(:,1), shroud_curve(:,2))
hold on
plot(shroud_points(:,1),shroud_points(:,2),'o','color','r')
axis equal

plot(impeller_curve(:,1), impeller_curve(:,2))

line([0 0], ylim);  %x-axis
line(xlim, [0 0]);  %y-axis

%% Blades

% Leading and Trailing Edges
r_inlet = r_eye; % m - this assumption can be changed

% Surface of Revolution
rms_curve =  [shroud_curve(:,1) - 0.5*crosswise_gap.*cos(normal_angles), shroud_curve(:,2) - 0.5*crosswise_gap.*sin(normal_angles)]; % [m,m] - halfway between shroud and impeller

% Inlet Blade Angles
inlet_gap = crosswise_gap(find(rms_curve > r_inlet, 1));
u_inlet = shaft_speed*r_inlet; % m/s - rotational speed at the inlet
v_inlet = vdot_fuel/(2*pi*r_inlet*inlet_gap) % m/s - fluid velocity at the inlet, assumed to be entirely radial. Can be found with CFD instead for more accuracy
angle_inlet = atan(u_inlet/v_inlet); % rad - inlet angle of blades relative to straight radial direction

% Outlet Blade Angles
slip_factor

%%

% Finds the slopes at the ends of curves and the minimum curvature
function [slope_shroud, slope_impeller, shroud_curvature] = compute_curves(bezier_vert, bezier_horiz, r_eye, r_exit, w_exit, impeller_thickness, impeller_height, r_shaft)
    shroud_points = [r_eye,impeller_height; r_eye, impeller_height-bezier_vert; r_exit-bezier_horiz,impeller_thickness+w_exit; r_exit, impeller_thickness+w_exit];
    shroud_curve = bezier(shroud_points);

    A_eye = pi*(r_eye^2-r_shaft^2);
    A_exit = 2*pi*r_exit*w_exit;
    
    A_pump_flow = linspace(A_eye, A_exit, length(shroud_curve))';

    normal_angles = atan2(gradient(shroud_curve(:,2)), gradient(shroud_curve(:,1))) + pi/2;
    normal_angles(normal_angles > 2*pi) = normal_angles(normal_angles > 2*pi) - 2*pi;
  
    crosswise_gap =  sec(normal_angles) .* ( shroud_curve(:,1) - sqrt( pi^2 * shroud_curve(:,1).^2 - pi*A_pump_flow.*cos(normal_angles) ) / pi ); % m - gap between shroud curve and impeller curve at each point
    
    impeller_curve = [shroud_curve(:,1) - crosswise_gap.*cos(normal_angles), shroud_curve(:,2) - crosswise_gap.*sin(normal_angles)];

    % Calculate slope at end of the curves
    slope_shroud = (shroud_curve(end, 2) - shroud_curve(end - 1, 2)) / (shroud_curve(end,1) - shroud_curve(end-1,1));
    slope_impeller = (impeller_curve(end, 2) - impeller_curve(end - 1, 2)) / (impeller_curve(end,1) - impeller_curve(end-1,1));

    % Find minimum curvature
    s = [0; cumsum(sqrt(diff(shroud_curve(:,1)).^2 + diff(shroud_curve(:,2)).^2))];
    shroud_curvature = abs(gradient(shroud_curve(:,1), s).*gradient(gradient(shroud_curve(:,2), s), s) - gradient(shroud_curve(:,2), s).*gradient(gradient(shroud_curve(:,1), s), s)) ./ (gradient(shroud_curve(:,1), s).^2 + gradient(shroud_curve(:,2), s).^2).^(3/2);

end

function cost = cost_function(bezier_params, r_eye, r_exit, w_exit, impeller_thickness, impeller_height, r_shaft, r_min)
    bezier_vert = bezier_params(1);
    bezier_horiz = bezier_params(2);
    [slope_shroud, slope_impeller, shroud_curvature] = compute_curves(bezier_vert, bezier_horiz, r_eye, r_exit, w_exit, impeller_thickness, impeller_height, r_shaft);

    % The cost is the sum of the absolute slopes (we want them to be 0 for a horizontal line) and the curvature
    cost = abs(slope_shroud) + abs(slope_impeller) + 100*abs(1/max(shroud_curvature)-r_min);
end

function curve = bezier(points) % Credit to Lakshmi Narasimhan of NIT Allahabad, India
    for i=0:1:3
        sigma(i+1)=factorial(3)/(factorial(i)*factorial(3-i));
    end  
    projection=[];
    for u=0:0.001:1
        for d=1:4
            UB(d)=sigma(d)*((1-u)^(4-d))*(u^(d-1));
        end
        projection=cat(1,projection,UB);
    end
    curve = projection*points;
end
