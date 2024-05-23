load('impulse_turbine_efficiency.mat');
load('pitch_chord_zweifel');

turbine_shaft_speed = shaft_speed*gear_ratio; % rad/s

%% Nozzle Plate
mach_spouting = sqrt(2/(gamma_gg-1)*((p_amb/p_gg)^((1-gamma_gg)/gamma_gg)-1)); % unitless
T_spouting = Tt_gg*(p_amb/p_gg)^((gamma_gg-1)/gamma_gg); % K - static temperature for perfectly expanded flow
mw_gg = 0.001*cea_gg.output.eql.mw(3); % kg/mol - molecular weight (fully combusted)
R_gg = 8.3145 / mw_gg; % J/mol*K - specific gas constant
sonic_spouting = sqrt(gamma_gg*R_gg*T_spouting); % m/s - speed of sound at nozzle exit plane

% v_spouting_norm = mach_spouting*sonic_spouting; % m/s - speed at which gas leaves the nozzles

A_throat_np = pi/4*(d_throat_nozzle)^2; % m
mdot_gg = A_throat_np*nozzle_number*p_gg/c_star_gg; % kg/s
A_choked = mdot_gg*sqrt(Tt_gg)/p_gg * sqrt(R_gg/gamma_gg) * ((gamma_gg+1)/2)^(0.5*(gamma_gg+1)/(gamma_gg-1));


%% Energy Losses
%expansion_energy_loss = % NASA Turbines 1976; fig. 16

%% Rotor Dimensions
% One stage impulse rotor - does not work for multistage and/or velocity driven turbines
stage_efficiency = 1.00;
error = 1;
while abs(error) > 0.01
    stage_efficiency_old = stage_efficiency;

    power_turbine = shaft_power_turbine/stage_efficiency; % W
    v_pitchline_rotor = sqrt(0.5*power_turbine/mdot_gg); % m/s - tangential rotor speed at pitchline
    
    v_spouting = [2*v_pitchline_rotor; sqrt(v_spouting_norm^2 - 4*v_pitchline_rotor^2)]; % [m/s, m/s] - [tangential, axial]
    angle_nozzle = atan(v_spouting(1)/v_spouting(2)); % rad
    r_rotor_pitchline = v_pitchline_rotor / turbine_shaft_speed; % m
    
    isentropic_v_ratio = v_pitchline_rotor/v_spouting_norm; % unitless
    stage_efficiency = interp1(impulse_turbine_efficiency(:,1), impulse_turbine_efficiency(:,2), isentropic_v_ratio, 'spline'); % unitless

    error = stage_efficiency - stage_efficiency_old;
end

torque_turbine = power_turbine / turbine_shaft_speed;

%% Blade Profile
blade_chord_angle = 90*pi/180; % rad - assuming symmetrical blades until proven otherwise
angle_rotor = atan(v_pitchline_rotor/v_spouting(2)); % rad - this is measured FROM THE VERTICAL
exit_vector_angle = pi/2 - angle_rotor; % rad - this is measured FROM THE HORIZONTAL

if exit_vector_angle < 30*pi/180
    % fprintf('Rotor blade angle of %g deg is too steep (Zweifel has no correlation)\n', angle_rotor*180/pi);
elseif exit_vector_angle > 80*pi/180
    % fprintf('Rotor blade angle of %g deg is too shallow (Zweifel has no correlation)\n', angle_rotor*180/pi);
else
    pitch_chord_ratio = interp1(pitch_chord_zweifel(:,1), pitch_chord_zweifel(:,2), -exit_vector_angle*180/pi, 'spline'); % unitless - NASA Turbines 1974, figure 26, Zweifel impulse blades curve
end

blade_inlet_angle = angle_rotor; % rad
blade_outlet_angle = angle_rotor; % rad

zweifel_coeff = (0.90+1.15)/2; % unitless
blade_chord_rotor = blade_width_rotor/sin(blade_chord_angle); % m
blade_pitch_rotor = blade_chord_rotor * pitch_chord_ratio; % m - the ideal circumferential distance between neighboring blades
blade_opening_rotor = blade_pitch_rotor * sin(angle_rotor); % m


%% Blade Depth
blade_number_rotor = round(2*pi*r_rotor_pitchline / blade_pitch_rotor); % unitless - number of blades along the circumference of the rotor

rotor_tab = 10*pi/180; % rad - there are four gaps around the blisk where blades are excluded
blade_number_rotor =  4*floor(r_rotor_pitchline*(2*pi-4*rotor_tab)/blade_pitch_rotor/4);
                                  
blade_length_rotor = mdot_gg/(dens_gas_spouted*v_pitchline_rotor*cos(angle_rotor)*blade_opening_rotor*blade_number_rotor); % m

r_rotor_tip = r_rotor_pitchline + blade_length_rotor/2; % m
r_rotor_base = r_rotor_pitchline - blade_length_rotor/2; % m

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
