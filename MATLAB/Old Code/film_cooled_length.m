%% Objective: Finds the length of a film cooled region per NASA paper SP-8124

function length = film_cooled_length(mdot_film, T_film, Prandtl, viscosity, d_avg, l_chamber, dens_gas, v_inj_gas, T_gas, H_gas, p_gas, cp_gas, mol_gas)
    surface_tension = py.CoolProp.CoolProp.PropsSI('I', 'T', T_film, 'Q', 0, 'Ethanol'); % N/m
    mol_film = py.CoolProp.CoolProp.PropsSI('M', 'T', T_film, 'Q', 0, 'Ethanol'); % N/m
    
    H_gas_stag = H_gas + v_inj_gas.^2/2; % J - stagnation enthalpy of combustion gas

    T0 = 0;
    T_interface = T_film;
    p_vapor = 101325;
    
    while(abs(T_interface-T0) > 1)
        H_film_saturation = py.CoolProp.CoolProp.PropsSI('H', 'P', p_vapor, 'Q', 0, 'Ethanol');
        H_film =  py.CoolProp.CoolProp.PropsSI('H', 'T', T_film, 'P', p_vapor, 'Ethanol');
        T0 = T_interface;
        T_interface = py.CoolProp.CoolProp.PropsSI('T', 'P', p_vapor, 'Q', 0, 'Ethanol'); % K - temperature at combustion gas / film coolant interface
        H_gas_film = H_gas + cp_gas*(T_interface-T_gas); % J - enthalpy of combustion gas at temperature of film coolant;
        B = (H_gas_stag - H_gas_film) / (H_film_saturation - H_film);
        p_vapor = p_gas / (1+mol_film/(B*mol_gas)); % Vapor pressure of film coolant at interface with combustion gas
    end

    delta = 1.3; % empirical entrainment augmentation factor. For injection parallel to core flow, 1.0 to 1.6; for swirl injection, maybe as low as 0.4
    X = (dens_gas)^0.5 * v_inj_gas * (T_gas/T_interface)^0.25 / surface_tension; % 1/N^0.5 - liquid entrainment correlation parameter CHECK IF THERE SHOULD BE A G
    A = -1.005245E-9*X^2 + 0.0003704722077*X + 1.52524723076; % 1/m - liquid entrainment parameter

    a = 10*A*(1+3*(surface_tension*X)^-0.8); % 1/m - heat-transfer augmentation factor for liquid surface roughness

    Reynolds = mdot_film./(pi*d_avg*viscosity); % Reynold's number of film PROBABLY INVALID AROUND THROAT AND NOZZLE

    if (Reynolds < 3000) % Really should be < 2300, but using laminar flow in transition region for conservatism
        lam_square = [1, 1.43, 2, 3, 4, 8, 100000000000; 3.61, 3.73, 4.12, 4.79, 5.33, 6.49, 8.23]; % Nusselt number for laminar flow in square pipes (Lavine, table 8.1)
        Nusselt = interp1(lam_square(1,:),lam_square(2,:), max(l_chamber./d_avg, d_avg/l_chamber)); % REVIEW, PROBABLY WRONG
    elseif (Prandtl >= 0.7 && Prandtl <= 2000 && Reynolds >= 3000 && Reynolds <  5E6)  
        f = (0.79 .* log(Reynolds) - 1.64).^(-2); % friction factor, smooth pipe approximation REVIEW IF THIS IS ACCURATE
        Nusselt = (f/8).*(Reynolds-1000).*Prandtl/(1+12.7*(f/8).^0.5.*(Prandtl.^(2/3)-1)); % Nusselt Number, Gnielenski correlation (turbulent)
    else
        error("No correlation for film coolant Reynold's and/or Prandtl number");
    end


    Stanton = Nusselt / (Reynolds*Prandtl);
    V = pi*d_avg*dens_gas*v_inj_gas*Stanton*B*a; % kg/m*s - surface vaporization rate
    length = log( 1 + A*mdot_film/V) / A;
end




% film_thickness = 0.376*x.*Re_film.^-0.2; % m - boundary layer thickness of the film coolant (fuel)

% film_eff_theory = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000; 1, 0.91, 0.77, 0.69, 0.64, 0.59, 0.56, 0.53, 0.52]; % Correlation between Reynold's and film efficiency (eta), 1969 paper
% film_eff = interp1(film_eff_theory(1,:), film_eff_theory(2,:), Re_film);

% film_length = film_eff*mdot_film*(cp_film_liq*(T_boil_film-T_film_initial) + heat_vaporization_film) ./ ( 2*pi*r1(1)*h_gas*(T_stagnation-T_boil_film));
