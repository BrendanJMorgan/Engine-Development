% NASA SP-8052 "LIQUID ROCKET ENGINE TURBOPUMP INDUCER" 1971

%% Requirement for Inducer
Ca_target = 2.0;
v_impeller_inlet = vdot_pump / (pi*(r_eye^2-r_eye_inner^2)); % m/s
NPSH_required = Ca_target*v_impeller_inlet^2 / (2*g) % m 
NPSH_available = (p_amb - vapor_pressure_pump)/(density_pump*g); % m - Net positive suction head. Assuming zero plumbing water column and inlet velocity
NPSH_inducer = NPSH_required - NPSH_available; % m

%% Tip Diameter

hub_tip_ratio = 0.30; % unitless - nu
vdot_corrected = vdot_pump / (1-hub_tip_ratio^2); % m3/s - Q'

ss_speed = shaft_speed*sqrt(vdot_pump)/(g * NPSH_inducer)^0.75; % pseudo unitless - suction specific speed - 60/(2*pi)*sqrt(15850.32)/(3.28084/g)
ss_speed_imperial = ss_speed * 60/(2*pi)*sqrt(15850.32)/(3.28084/g)^0.75; % pseudo unitless
corrected_ss_speed_imperial = ss_speed_imperial / sqrt(1-hub_tip_ratio); % pseudo unitless - corrected suction specific speed

inducer_flow_coeff = 2*3574/corrected_ss_speed_imperial / ( 1 + sqrt( 1+6*(3574/corrected_ss_speed_imperial)^2 ) ); % unitless - phi

max_cavitation_number = 2*inducer_flow_coeff^2 / (1-2*inducer_flow_coeff^2);

r_tip_inducer = 0.5*0.3048*0.37843*(vdot_corrected*15850.323 / (shaft_speed*60/(2*pi)*inducer_flow_coeff))^(1/3); % m

%% Blade Construction

v_inlet_inducer = [0; 0];
v_inlet_inducer(2) = vdot_corrected / (pi*(1-hub_tip_ratio^2)*r_tip_inducer^2); % m/s - c_m
v_inlet_inducer(1) = shaft_speed*r_tip_inducer; % m/s - w1u ????
flow_blade_ratio = 0.425; % unitless

% cavitation_number_goal = 2; % unitless - K
% v_cavity = v_inlet_inducer*sqrt(1+K); % [m/s, m/s]

flow_incidence_angle = atan(v_inlet_inducer(2) / v_inlet_inducer(1)); % rad - alpha
blade_angle = flow_incidence_angle/0.425 % rad - beta

%% Blade Lead
blade_number_inducer = 4; % 3 or 4 is considered good
blade_lead_inducer = 2*pi*r_tip_inducer*tan(blade_angle); % m - how much the blade advances per turn
blade_pitch_inducer = blade_lead_inducer / blade_number_inducer % m
v_inlet_inducer

%% Blade Sweep
% solidity_inducer = 2.5;
% h_inducer = blade_lead_inducer*solidity_inducer; % m - height of inducer

%% Clearance Losses

%% Axial Clearance
% Between inducer and impeller

%% 
v_tip_inducer = shaft_speed*r_tip_inducer; % m/s
head_coeff_inducer = NPSH_inducer*g / v_tip_inducer^2; % unitless - camber needed if >0.075




