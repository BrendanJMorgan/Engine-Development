% Theoretical spouting velocities of high-energy working fluids often force operation at
% velocity ratios that come short of delivering maximum efficiency. This condition occurs
% because the pitchline velocities required for higher turbine efficiencies exceed allowable
% stresses in the turbine rotor. The practice in this case requires developing turbine design
% power with the highest allowable pitchline velocity and least number of turbine stages
% needed to accomplish the job. - NASA turbines 1974, page 83

load('expansion_loss_data.mat');
blade_angle_interp = 0:1:140; % deg
blade_width_interp = [0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 2.0]; % in
load('incidence_loss_data.mat');
load('mach_loss_data.mat');
load('impulse_turbine_efficiency.mat');
load('pitch_chord_zweifel.mat');

turbine_shaft_speed = shaft_speed*gear_ratio; % rad/s

%% Nozzle Plate
mach_spouting = sqrt(2/(gamma_gg-1)*((p_amb/p_gg)^((1-gamma_gg)/gamma_gg)-1)); % unitless
T_spouting = Tt_gg*(p_amb/p_gg)^((gamma_gg-1)/gamma_gg); % K - static temperature for perfectly expanded flow
mw_gg = 0.001*cea_gg.output.eql.mw(3); % kg/mol - molecular weight (fully combusted)
R_gg = 8.3145 / mw_gg; % J/mol*K - specific gas constant
sonic_spouting = sqrt(gamma_gg*R_gg*T_spouting); % m/s - speed of sound at nozzle exit plane

A_throat_np = mdot_gg*c_star_gg/p_gg; % kg/s - cumulative throat area on the nozzle plate
nozzle_number = A_throat_np / (pi/4*d_throat_nozzle^2); % unitless

%% Velocities
% One stage impulse rotor - does not work for multistage and/or velocity driven turbines
isentropic_v_ratio = 0.35; % unitless - "A single-row impulse stage delivers best performance at velocity ratios between 0.30 and 0.40"
v_pitchline_rotor = r_pitchline_rotor*turbine_shaft_speed; % m/s - tangential rotor speed at pitchline

% v_spouting_norm = v_pitchline_rotor/isentropic_v_ratio; % unitless
v_spouting_norm = v_exhaust_gg; % m/s
isentropic_v_ratio = v_pitchline_rotor/v_spouting_norm
v_axial_rotor = sqrt(v_spouting_norm^2 - 4*v_pitchline_rotor^2); % m/s
v_spouting = [2*v_pitchline_rotor; v_axial_rotor]; % [m/s, m/s] - [tangential, axial]
angle_nozzle = atan(v_spouting(1)/v_spouting(2)); % rad
v_blade_relative = [v_pitchline_rotor; v_axial_rotor]; % m/s
incidence_angle_blade = atan(v_blade_relative(1)/v_blade_relative(2)); % rad - this is measured FROM THE VERTICAL

    blade_chord_angle = 90*pi/180; % rad - assuming symmetrical blades until proven otherwise
    exit_vector_angle = pi/2 - incidence_angle_blade ; % rad - this is measured FROM THE HORIZONTAL
    entrance_vector_angle = exit_vector_angle; % rad - assuming symmetrical blades, FROM THE HORIZONTAL

M_spouting = norm(v_spouting)*(gamma_gg*R_gg*Tt_gg - 0.5*(gamma_gg-1)*norm(v_spouting).^2) .^ -0.5; % unitless - mach at nozzle plate exit
M_spouting = M_spouting*v_spouting/norm(v_spouting); % unitless
M_blade_relative = norm(v_blade_relative)/norm(v_spouting)*M_spouting; % unitless - mach in corotating frame of blades

%% Blade Spacing
A_choked = mdot_gg*sqrt(Tt_gg)/p_gg * sqrt(R_gg/gamma_gg) * ((gamma_gg+1)/2)^(0.5*(gamma_gg+1)/(gamma_gg-1)); % m2 - minimum throughput area needed to prevent gas backup
A_exit_rotor = A_throat_np * ( (2 + (gamma_gg - 1)*M_spouting.^2 ) / (gamma_gg + 1) ) .^ ((gamma_gg + 1)/(2*(gamma_gg - 1))) / M_spouting; % m2 - cumulative area of all blade gaps in rotor

if exit_vector_angle < 30*pi/180
    fprintf('Rotor blade angle of %g deg is too steep (Zweifel has no correlation)\n', incidence_angle*180/pi);
elseif exit_vector_angle > 80*pi/180
    fprintf('Rotor blade angle of %g deg is too shallow (Zweifel has no correlation)\n', incidence_angle*180/pi);
else
    pitch_chord_ratio = interp1(pitch_chord_zweifel(:,1), pitch_chord_zweifel(:,2), -exit_vector_angle*180/pi, 'spline'); % unitless - NASA Turbines 1974, figure 26, Zweifel impulse blades curve
end
zweifel_coeff = (0.90+1.15)/2; % unitless
blade_chord_rotor = blade_width_rotor/sin(blade_chord_angle); % m
blade_pitch_rotor = blade_chord_rotor * pitch_chord_ratio; % m - the ideal circumferential distance between neighboring blades
blade_opening_rotor = blade_pitch_rotor * sin(incidence_angle_blade); % m

blade_number_rotor =  4*round(2*pi*r_pitchline_rotor/blade_pitch_rotor/4); % unitless - number of blades along the circumference of the rotor

%% Energies
% Adjust turbine gas-path expansion and kinetic energy, using the respective two-dimensional
% loss coefficients for expansion energy, kinetic energy, Mach number, and incidence angle as
% determined from plots like those in figures 16 through 18.

expansion_loss_coeff = interp2(blade_angle_interp, blade_width_interp, expansion_loss_data, (pi-entrance_vector_angle-exit_vector_angle)*180/pi, blade_width_rotor/0.0254); % unitless - NASA Turbines 1976; fig. 16
kinetic_loss_coeff = 2*expansion_loss_coeff^2 - 1; % unitless, equation 3
incidence_loss_coeff = interp1(incidence_loss_data(:,1), incidence_loss_data(:,2), incidence_angle_blade*180/pi); % unitless
mach_loss_coeff = interp1(mach_loss_data(:,1), mach_loss_data(:,2), norm(M_blade_relative)); % unitless
    if norm(M_blade_relative) <= 1; mach_loss_coeff = 1; end % assuming that mach loss does not apply to a subsonic case?
stage_efficiency = interp1(impulse_turbine_efficiency(:,1), impulse_turbine_efficiency(:,2), isentropic_v_ratio, 'spline'); % unitless
clearance_loss_coeff = -1.63*tip_clearance/blade_length_rotor+1;

p_exit_turbine = p_amb; % Pa
T_exit_turbine = Tt_gg*(p_exit_turbine/p_gg)^((gamma_gg-1)/gamma_gg); % K
expansion_specific_work = R_gg*T_exit_turbine*(gamma_gg/(gamma_gg-1))*(1-(p_exit_turbine/p_gg)^((gamma_gg-1)/gamma_gg)); % J/kg - equation 1b, equation 2
kinetic_specific_work = 0.5*norm(v_spouting)^2; % J/kg
specific_work_rotor = expansion_specific_work*expansion_loss_coeff + kinetic_specific_work*kinetic_loss_coeff*incidence_loss_coeff*mach_loss_coeff; % J/kg - equation 5b

turbine_efficiency_full = specific_work_rotor / (expansion_specific_work + kinetic_specific_work); % unitless
nozzle_arc = 2*pi*admission_fraction*r_pitchline_rotor; % m
wheel_ratio = v_axial_rotor / norm(v_spouting); % unitless - ration between outlet/inlet velocity of blades. Must be absolute, because the realtive ratio would always be 1?
turbine_efficiency_partial = (1+wheel_ratio*(1-blade_pitch_rotor/(3*nozzle_arc)))/(1+wheel_ratio) * turbine_efficiency_full - 1.539E-9*isentropic_v_ratio*dens_gas_spouted*turbine_shaft_speed*(1-admission_fraction)/admission_fraction; % unitless - equation 14

turbine_efficiency = turbine_efficiency_partial*clearance_loss_coeff;

power_turbine = shaft_power_turbine/turbine_efficiency; % W - ideal power needed to meet pump requirements
mdot_gg = 0.5 * power_turbine / (turbine_shaft_speed*r_pitchline_rotor)^2; % kg/s - mass flow rate through the gg / manifold / rotor - this is actually equation 10b simplified with impulse rotor assumptions
% Impulse Rotor Assumption: c1*cos(alpha1) = 2v; w2*cos(alpha2) = v; v = omega*r; P = mdot*delta_h
torque_turbine = shaft_power_turbine / turbine_shaft_speed; % N*m

%% Blade Depth

if blade_length_rotor < 1.5*A_choked/(blade_number_rotor*blade_opening_rotor)
    error("Rotor blades are too short and choking the flow")
end

r_rotor_tip = r_pitchline_rotor + blade_length_rotor/2; % m
r_rotor_base = r_pitchline_rotor - blade_length_rotor/2; % m

%% Mass Flow Adjustments
mdot_total = mdot_cc + mdot_gg;
gg_fraction = mdot_gg / mdot_total; % not totally sure if denominator should be cc instead of total
mdot_fuel_gg = mdot_gg*(1/(1+OF_gg)); % kg/s - Fuel Mass Flow Rate
mdot_ox_gg = mdot_gg*(OF_gg/(1+OF_gg)); % kg/s - Oxidizer Mass Flow Rate
mdot_fuel_total = mdot_fuel_cc+mdot_fuel_gg; % kg/s
mdot_ox_total = mdot_ox_cc+mdot_ox_gg; % kg/s

%% Rotor Stress
AN2 = pi*((r_rotor_tip-r_rotor_base)/0.0254)^2 * (turbine_shaft_speed*30/pi)^2;

rotor_density = 7660; % kg/m3 - 316 stainless - https://bssa.org.uk/bssa_articles/elevated-temperature-physical-properties-of-stainless-steels/
rotor_poissons = 0.31; % 316 at 700 C
rotor_strength = 275E6; % Pa

r_rotor_tip/0.0254;
stress_radial = (3+rotor_poissons)/8*rotor_density*turbine_shaft_speed^2*r_rotor_tip^2; % Pa
stress_radial/6.895e+6; % ksi
stress_FSu = rotor_strength/stress_radial;

%% Other
thrust_turbine = mdot_gg * v_spouting(2);
scfm_gg = mdot_gg*30.42*60 % standard cubic feet per minute of nitrogen
spin_time_ideal = 8.64/mdot_gg; % s - there are 8.64 kg of nitrogen in a k bottle. in reality flow will stop when pressure drops too low.
