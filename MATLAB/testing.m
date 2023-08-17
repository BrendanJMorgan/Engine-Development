
% check reynold's number, convergence angle affects laminar/turbulent regime boundaries

stanton_film = stanton_coeff * re_film^-0.2 * pr_film^-0.6; % turbulent regime

stanton_film = 0.318 * re_film^-0.5 * pr_film^-0.6; % laminar regime


q_wall = density_gas*v_film*stanton_film*(H_aw - H_film); % W - heat flux to wall, NASA 1977 pg. 82


%lewis_number_gas = gas_cond / (gas_density*cp_gas*diffusivity_gas);

entrainment_augmentation = 1.3; % empirical factor for coolant injection. For injection parallel to core flow, 1.0 to 1.6; for swirl injection, maybe as low as 0.4
film_surface_tension;
entrainment_parameter = 2.10908 * entrainment_augmentation * (rho_gas/g)^0.5 * v_gas * (T_gas/T_interface)^0.25 / surface_tension; % lbf^-1/2
film_liquid_length = log( 1 + entrainment*mdot_film/surface_vaporization)/entrainment;