NPSH_available = p_amb/(density*g) - vapor_pressure/(density*g) % m - Net positive suction head. Assuming zero plumbing water column and inlet velocity
NPSH_required = 10 % m - completely arbitrary
NPSH_inducer = NPSH_required - NPSH_available

%% Tip Diameter

hub_tip_ratio = 0.30; % unitless - nu

ss_speed = shaft_speed*sqrt(vdot)/(g * NPSH_inducer)^0.75 % suction specific speed - 60/(2*pi)*sqrt(15850.32)/(3.28084/g)
ss_speed_imperial = ss_speed * 60/(2*pi)*sqrt(15850.32)/(3.28084/g)^0.75
corrected_ss_speed_imperial = ss_speed_imperial / sqrt(1-hub_tip_ratio) % corrected suction specific speed

inducer_flow_coeff = 2*3574/corrected_ss_speed_imperial / ( 1 + sqrt( 1+6*(3574/corrected_ss_speed_imperial)^2 ) ) % unitless - phi

max_cavitation_number = 2*inducer_flow_coeff^2 / (1-2*inducer_flow_coeff^2)

tip_diameter = 0.3048*0.37843*(vdot*15850.323 / ((1-hub_tip_ratio^2)*shaft_speed*60/(2*pi)*inducer_flow_coeff))^(1/3) % m

%% Tip Contour

%% Fluid Thermodynamic Effects

%% Blade Profile
blade_angle = 

blade_wedge_angle = blade_angle - atan(1.10*inducer_flow_coeff); % rad

%% Blade Sharpness
leading_edge_radius = 0.01*blade_thickness; % m

%% Blade Sweep

%% Blade Cant
blade_cant_angle = 
blade_angle = blade_cant_angle/0.425;

%% Blade Lead
blade_lead = 2*pi*tip_diameter/2*tan(blade_angle); % m - how much the blade advances per turn

%% Blade Thickness



