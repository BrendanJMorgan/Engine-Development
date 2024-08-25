% NASA SP-8052 "LIQUID ROCKET ENGINE TURBOPUMP INDUCER" 1971

%% Net Positive Suction Head (NPSH) Requirements

k1 = 1.2; 
k2 = 0.2334 + (shaft_speed*r_eye_impeller/128.3)^4;
cavitation_limit = (shaft_speed*r_eye_impeller)^2/(2*g); % unitless
NPSH_SE = ((k1+k2)*eye_flow_coeff^2 + k2) * cavitation_limit; % m - NPSH, shockless entry. Pump handbook pg. 2.88
NPSH_required = NPSH_SE; % m - assuming operation at best efficiency point (BEP) and shockless

NPSH_available = (p_amb - vapor_pressure_pump)/(density_pump*g); % m - Net positive suction head. Assuming zero plumbing water column and inlet velocity
NPSH_inducer = (NPSH_required - NPSH_available)*NPSH_margin; % m

%% Tip Diameter
vdot_corrected = vdot_pump / (1-hub_tip_ratio_inducer^2); % m3/s - volumetric flow rate through open portion of inducer

ss_speed = shaft_speed*sqrt(vdot_pump)/(g * NPSH_inducer)^0.75; % unitless - suction specific speed
ss_speed_imperial = ss_speed*2733.016; % pseudo unitless
corrected_ss_speed_imperial = ss_speed_imperial / sqrt(1-hub_tip_ratio_inducer^2); % pseudo unitless - corrected suction specific speed

flow_coeff_inducer = 3574/corrected_ss_speed_imperial / (1 + sqrt( 1+6*(3574/corrected_ss_speed_imperial)^2)/2 ); % unitless - phi
flow_coeff_inducer = 1.3077*sqrt(1-hub_tip_ratio_inducer^2)/ss_speed / (1 + 0.5*sqrt( 1+6*(1.3077*sqrt(1-hub_tip_ratio_inducer^2)/ss_speed)^2)); % unitless - phi

r_tip_inducer = 0.5*0.3048*0.37843*(vdot_corrected*15850.323 / (shaft_speed*30/pi*flow_coeff_inducer))^(1/3); % m
r_hub_inducer = r_tip_inducer*hub_tip_ratio_inducer; % m

r_pitchline_inducer = 0.5*(1+hub_tip_ratio_inducer)*r_tip_inducer; % m

%% Blade Construction
v_inlet_inducer = [0; 0];
v_inlet_inducer(1) = shaft_speed*r_tip_inducer; % m/s - tangential velocity at inducer pitchline tip
v_inlet_inducer(2) = vdot_corrected / (pi*r_tip_inducer^2); % m/s - axial velocity at inducer pitchline

flow_blade_ratio = 0.35; % unitless

flow_incidence_angle = atan(v_inlet_inducer(2) / v_inlet_inducer(1)); % rad - alpha
blade_angle_inducer = flow_incidence_angle/(1-flow_blade_ratio); % rad - beta

blade_lead_inducer = 2*pi*r_tip_inducer*tan(blade_angle_inducer); % m - blade lead - how much the blade advances per turn
solidity_inducer = 2.5;
h_min_inducer = blade_lead_inducer*solidity_inducer/blade_number_inducer*sin(blade_angle_inducer); % m - height of inducer

%% Axial Clearance
% Between inducer and impeller
clearance_axial_inducer = blade_lead_inducer / blade_number_inducer; % m

%% Head
v_tip_inducer = shaft_speed*r_tip_inducer; % m/s
head_coeff_inducer = NPSH_inducer*g / v_tip_inducer^2; % unitless - 0.15 is the cutoff between low and high head - camber needed if >0.075

%% Hub Construction
r_hub_base = r_eye_inner; % m
r_hub_top = 2*hub_tip_ratio_inducer*r_tip_inducer - r_hub_base; % m

%% Clearance Losses - not relevant with a proper shroud
r_cavity_inducer = r_tip_inducer + clearance_radial_inducer; % m
ss_speed_clearanced = ss_speed*(1-0.575*sqrt((r_cavity_inducer-r_tip_inducer)/(r_tip_inducer-r_hub_inducer))); % unitless
head_coeff_inducer_clearanced = head_coeff_inducer*(1-1.0*sqrt((r_cavity_inducer-r_tip_inducer)/(r_tip_inducer-r_hub_inducer))); % unitless

%% Other
inducer_torque = mdot_pump*shaft_speed*r_tip_inducer^2; % N*m - this is a conservative maximum
