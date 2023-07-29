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

%% Impeller Design

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

bezier_params_initial = [0.008, 0.008]; % starting values
options = optimset('Display', 'iter'); % display progress
cost_function_with_fixed_params = @(bezier_params) cost_function(bezier_params, r_eye, r_exit, w_exit, impeller_thickness, impeller_height, r_shaft);
bezier_params_optimal = fminsearch(cost_function_with_fixed_params, bezier_params_initial, options);

bezier_vert_optimal = bezier_params_optimal(1);
bezier_horiz_optimal = bezier_params_optimal(2);

shroud_points = [r_eye,impeller_height; r_eye, impeller_height-bezier_vert_optimal; r_exit-bezier_horiz_optimal,impeller_thickness+w_exit; r_exit, impeller_thickness+w_exit];
shroud_curve = bezier(shroud_points); % [m,m] - lower curved edge of revolution of the pump shroudjj

A_eye = pi*(r_eye^2-r_shaft^2); % m2 - area of the annular eye entrance
A_exit = 2*pi*r_exit*w_exit; % m2 - area of the cylindrical exit

A_pump_flow = linspace(A_eye, A_exit, length(shroud_curve))'; % m2 - cross sectional area perpendicular to flow

normal_angles = atan2(gradient(shroud_curve(:,2)), gradient(shroud_curve(:,1))) + pi/2;
normal_angles(normal_angles > 2*pi) = normal_angles(normal_angles > 2*pi) - 2*pi;
% normal_angles = [0; (mod(atan2(diff(shroud_curve(:,2)), diff(shroud_curve(:,1))) + pi/2, 2*pi)); pi/2]; % rad - angle of vector normal to each point of shroud curve
% normal_angles = (normal_angles(1:end-1) + normal_angles(2:end))/2;

crosswise_gap =  sec(normal_angles) .* ( shroud_curve(:,1) - sqrt( pi^2 * shroud_curve(:,1).^2 - A_pump_flow.*cos(normal_angles) ) / pi); % m - gap between shroud curve and impeller curve at each point

impeller_curve = [shroud_curve(:,1) - crosswise_gap.*cos(normal_angles), shroud_curve(:,2) - crosswise_gap.*sin(normal_angles)]; % [m,m] - upper curved edge of revolution of the pump impeller

r_min = 0.5*r_eye; % m - minimum allowable radius of curvature



% Check that shroud curve does not kink too much; ie that it's radius of curvature is not below r_min
s = [0; cumsum(sqrt(diff(shroud_curve(:,1)).^2 + diff(shroud_curve(:,2)).^2))];
curvature = abs(gradient(shroud_curve(:,1), s).*gradient(gradient(shroud_curve(:,2), s), s) - gradient(shroud_curve(:,2), s).*gradient(gradient(shroud_curve(:,1), s), s)) ./ (gradient(shroud_curve(:,1), s).^2 + gradient(shroud_curve(:,2), s).^2).^(3/2);
if any(curvature > 1/r_min)
    [max_curvature, idx] = max(curvature);
    min_radius = 1/max_curvature;
    fprintf("Pump shroud curvature reaches %.2f mm radius. Need above ~%.2f mm (half of eye radius) to avoid cavitation.", min_radius*1000, r_min*1000)
    %plot(shroud_curve(idx, 1), shroud_curve(idx, 2), 'go');
end


figure(1)
line(shroud_curve(:,1), shroud_curve(:,2))
hold on
plot(shroud_points(:,1),shroud_points(:,2),'o','color','r')
axis equal

plot(impeller_curve(:,1), impeller_curve(:,2))

%%

function [slope_shroud, slope_impeller] = compute_curves(bezier_vert, bezier_horiz, r_eye, r_exit, w_exit, impeller_thickness, impeller_height, r_shaft)
    shroud_points = [r_eye,impeller_height; r_eye, impeller_height-bezier_vert; r_exit-bezier_horiz,impeller_thickness+w_exit; r_exit, impeller_thickness+w_exit];
    shroud_curve = bezier(shroud_points);

    A_eye = pi*(r_eye^2-r_shaft^2);
    A_exit = 2*pi*r_exit*w_exit;
    
    A_pump_flow = linspace(A_eye, A_exit, length(shroud_curve))';

    normal_angles = atan2(gradient(shroud_curve(:,2)), gradient(shroud_curve(:,1))) + pi/2;
    normal_angles(normal_angles > 2*pi) = normal_angles(normal_angles > 2*pi) - 2*pi;
%     normal_angles = [0; (mod(atan2(diff(shroud_curve(:,2)), diff(shroud_curve(:,1))) + pi/2, 2*pi)); pi/2];
%     normal_angles = (normal_angles(1:end-1) + normal_angles(2:end))/2;
  
    crosswise_gap =  sec(normal_angles) .* ( shroud_curve(:,1) - sqrt( pi^2 * shroud_curve(:,1).^2 - A_pump_flow.*cos(normal_angles) ) / pi ); % m - gap between shroud curve and impeller curve at each point
    
    impeller_curve = [shroud_curve(:,1) - crosswise_gap.*cos(normal_angles), shroud_curve(:,2) - crosswise_gap.*sin(normal_angles)];

    % Calculate slope at end of the curves
    s_shroud = [0; cumsum(sqrt(diff(shroud_curve(:,1)).^2 + diff(shroud_curve(:,2)).^2))];
    s_impeller = [0; cumsum(sqrt(diff(impeller_curve(:,1)).^2 + diff(impeller_curve(:,2)).^2))];

    % Calculate slope at end of the curves
    ds_shroud = s_shroud(end) - s_shroud(end-1);
    ds_impeller = s_impeller(end) - s_impeller(end-1);
    slope_shroud = (shroud_curve(end, 2) - shroud_curve(end - 1, 2)) / ds_shroud;
    slope_impeller = (impeller_curve(end, 2) - impeller_curve(end - 1, 2)) / ds_impeller;
end

function cost = cost_function(bezier_params, r_eye, r_exit, w_exit, impeller_thickness, impeller_height, r_shaft)
    bezier_vert = bezier_params(1);
    bezier_horiz = bezier_params(2);
    [slope_shroud, slope_impeller] = compute_curves(bezier_vert, bezier_horiz, r_eye, r_exit, w_exit, impeller_thickness, impeller_height, r_shaft);

    % The cost is the sum of the absolute slopes (we want them to be 0 for a horizontal line)
    cost = abs(slope_shroud) + abs(slope_impeller);
end

function curve = bezier(points) % Credit to Lakshmi Narasimhan of NIT Allahabad, India

    for i=0:1:3
        sigma(i+1)=factorial(3)/(factorial(i)*factorial(3-i));
    end
    
    projection=[];
    for u=0:0.0001:1
        for d=1:4
            UB(d)=sigma(d)*((1-u)^(4-d))*(u^(d-1));
        end
        projection=cat(1,projection,UB);
    end
    curve = projection*points;

end


%%

% figure(1); clf
% plot(shroud_x, shroud_y);
% hold on
% plot(impeller_x, impeller_y);
% axis square
% 
% figure(2); clf
% plot(r_curve)
% hold on
% plot(sin(theta))
% plot(shroud_y)
% figure(2)
% plot(r_curve)
% hold on
% plot(r_min)

% Construction of the blades

% %% Inputs
% 
% r_pump = 2*0.0254; % m - blade radius
% h_exit_pump = 1*0.0254; % m - exit gap height
% 
% pump_eff = 0.55;
% blade_angle = 30*pi/180;
% 
% 
% %%
% 
% head_fuel = (p_out-p_in)/(density_fuel*9.80665); % m - pump pressure head
% head_ox = (p_out-p_in)/(density_ox*9.80665); % m - pump pressure head
% 
% v_rad_fuel = vdot_fuel/(2*pi*r_pump*h_exit_pump); % m/s - radial velocity
% v_rad_ox = vdot_ox/(2*pi*r_pump*h_exit_pump); % m/s - radial velocity
% 
% v_ang_fuel_ideal = ((pump_eff*v_rad_fuel*tan(blade_angle*pi/180) + sqrt(pump_eff*(4*9.80665*head_fuel+pump_eff*v_rad_fuel^2*tan(blade_angle*pi/180))) / (2*pump_eff*r_pump)));
% v_ang_ox_ideal = ((pump_eff*v_rad_ox*tan(blade_angle*pi/180) + sqrt(pump_eff*(4*9.80665*head_ox+pump_eff*v_rad_fuel^2*tan(blade_angle*pi/180))) / (2*pump_eff*r_pump)));
% 
% v_ang = (v_ang_fuel_ideal + v_ang_ox_ideal); % rad/s - angular velocity is a compromise between the pumps
% 
% power_fuel = mdot_fuel_cc*9.80665*head_fuel/pump_eff; % W
% power_ox = mdot_ox_cc*9.80665*head_ox/pump_eff; % W