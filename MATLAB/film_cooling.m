function length = film_cooled_length(dens_exhaust, v_exhaust, mdot_film, d_avg, T_film)
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
end