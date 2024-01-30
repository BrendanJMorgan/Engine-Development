% NASA SP-8052 "LIQUID ROCKET ENGINE TURBOPUMP INDUCER" 1971

NPSH_available = p_amb/(density_pump*g) - vapor_pressure_pump/(density_pump*g); % m - Net positive suction head. Assuming zero plumbing water column and inlet velocity
NPSH_required = 10; % m - completely arbitrary
NPSH_inducer = NPSH_required - NPSH_available;

%% Tip Diameter

hub_tip_ratio = 0.30; % unitless - nu
vdot_corrected = vdot_pump / (1-hub_tip_ratio^2); % m3/s - Q'


ss_speed = shaft_speed*sqrt(vdot_pump)/(g * NPSH_inducer)^0.75; % pseudo unitless - suction specific speed - 60/(2*pi)*sqrt(15850.32)/(3.28084/g)
ss_speed_imperial = ss_speed * 60/(2*pi)*sqrt(15850.32)/(3.28084/g)^0.75; % pseudo unitless
corrected_ss_speed_imperial = ss_speed_imperial / sqrt(1-hub_tip_ratio); % pseudo unitless - corrected suction specific speed

inducer_flow_coeff = 2*3574/corrected_ss_speed_imperial / ( 1 + sqrt( 1+6*(3574/corrected_ss_speed_imperial)^2 ) ); % unitless - phi

max_cavitation_number = 2*inducer_flow_coeff^2 / (1-2*inducer_flow_coeff^2);

tip_diameter_inducer = 0.3048*0.37843*(vdot_corrected*15850.323 / (shaft_speed*60/(2*pi)*inducer_flow_coeff))^(1/3); % m

%% Blade Construction

v_inlet_inducer = 4*vdot_corrected / (pi*tip_diameter_inducer^2); % m/s - c_m
flow_blade_ratio = 0.425; % unitless

r_pitchline_inducer = tip_diameter_inducer; % m
flow_incidence_angle = atan(v_inlet_inducer / (shaft_speed*r_pitchline_inducer)); % rad - gamma
blade_angle = flow_incidence_angle/0.425;

%% Blade Profile
blade_wedge_angle = blade_angle - atan(1.10*inducer_flow_coeff); % rad
leading_edge_radius = 0.01*blade_thickness; % m

%% Blade Sweep

%% Blade Lead
blade_lead = 2*pi*tip_diameter_inducer/2*tan(blade_angle); % m - how much the blade advances per turn

%% Blade Thickness



