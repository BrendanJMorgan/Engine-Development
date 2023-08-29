%% Combustion Thermal Environment
T_free = Tc ./ (1 + (gamma-1)/2 * M.^2); % Free-Stream Temperature

%% Film Flow

mdot_film = film_fraction*mdot_cc;
mdot_gas = (1-film_fraction)*mdot_cc;

cp_film_liq = py.CoolProp.CoolProp.PropsSI('C', 'P', pc, 'Q', 0, 'Ethanol'); % J/kg-s
visc_film_liq = py.CoolProp.CoolProp.PropsSI('V', 'P', pc, 'Q', 0, 'Ethanol'); % Pa-s
dens_film_liq = py.CoolProp.CoolProp.PropsSI('D', 'P', pc, 'Q', 0, 'Ethanol'); % kg/m3
Pr_film_liq = py.CoolProp.CoolProp.PropsSI('Prandtl', 'P', pc, 'Q', 0, 'Ethanol');
T_boil_film = py.CoolProp.CoolProp.PropsSI( 'T', 'P', pc, 'Q', 0, 'Ethanol' ); % K - Boiling point of film coolant at chamber pressure
T_film_initial = T_amb; % K - THIS SHOULD BE THE COOLANT OUTLET TEMPERATURE

heat_vaporization_film = py.CoolProp.CoolProp.PropsSI('H','P',pc,'Q',1,'Ethanol') - py.CoolProp.CoolProp.PropsSI('H','P',pc,'Q',0,'Ethanol');

cp_film_gas = py.CoolProp.CoolProp.PropsSI( 'C', 'P', pc, 'Q', 1, 'Ethanol');
visc_film_gas = py.CoolProp.CoolProp.PropsSI( 'V', 'P', pc, 'Q', 1, 'Ethanol');
dens_film_gas = py.CoolProp.CoolProp.PropsSI( 'D', 'P', pc, 'Q', 1, 'Ethanol');
Pr_film_gas = py.CoolProp.CoolProp.PropsSI('Prandtl', 'P', pc, 'Q', 1, 'Ethanol');

% Approximations
cp_film_avg = cp_film_gas;
visc_film_avg = visc_film_gas;
dens_film_avg = dens_film_gas;
Pr_film_avg = Pr_film_gas;

%%
% check reynold's number, convergence angle affects laminar/turbulent regime boundaries

%stanton_film = stanton_coeff * re_film^-0.2 * pr_film^-0.6; % turbulent regime

%stanton_film = 0.318 * re_film^-0.5 * pr_film^-0.6; % laminar regime

rho_gas = 1;
v_film = 1;
H_aw = 1;
H_film = 1;

%q_wall = density_gas*v_film*stanton_film*(H_aw - H_film); % W - heat flux to wall, NASA 1977 pg. 82
T_film = 500;
surface_tension_film = py.CoolProp.CoolProp.PropsSI('I', 'T', T_film, 'Q', 0, 'Ethanol'); % N/m
mol_weight_film = py.CoolProp.CoolProp.PropsSI('M', 'T', T_film, 'Q', 0, 'Ethanol'); % N/m
mol_weight_gas

P_film_vapor = Pc / (1 + film_mol_weight)
T_interface = 

%lewis_number_gas = gas_cond / (gas_density*cp_gas*diffusivity_gas);

entrainment_augmentation = 1.3; % empirical factor for coolant injection. For injection parallel to core flow, 1.0 to 1.6; for swirl injection, maybe as low as 0.4
film_surface_tension;
entrainment_parameter = 2.10908 * entrainment_augmentation * (rho_gas/g)^0.5 * v_exhaust * (T_free/T_interface)^0.25 / surface_tension; % lbf^-1/2
surface_vaporization = pi*d_chamber
film_liquid_length = log( 1 + entrainment*mdot_film/surface_vaporization)/entrainment;




