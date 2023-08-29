%% Finds the length of a film cooled region per NASA paper SP-8124

function length = film_cooled_length(mdot_film,  T_film, cp_film, dens_film, v_film, d_avg, h_exhaust, dens_exhaust, v_exhaust, T_exhaust, H_exhaust, p_exhaust, cp_exhaust, mol_exhaust)
    T_film = 500;
    surface_tension = py.CoolProp.CoolProp.PropsSI('I', 'T', T_film, 'Q', 0, 'Ethanol'); % N/m
    mol_film = py.CoolProp.CoolProp.PropsSI('M', 'T', T_film, 'Q', 0, 'Ethanol'); % N/m
    
    delta = 1.3; % empirical entrainment augmentation factor. For injection parallel to core flow, 1.0 to 1.6; for swirl injection, maybe as low as 0.4
    X = 0.47414 * (rho_gas/g)^0.5 * v_exhaust * (T_free/T_interface)^0.25 / surface_tension; % 1/N^0.5 - liquid entrainment correlation parameter
    A = -1005245E-9*X^2 + 0.0003704722077*X + 1.52524723076; % 1/m - liquid entrainment parameter

    H_exhaust_stag = H_exhaust + v_exhaust^2/2; % J - stagnation enthalpy of combustion gas
    H_film_saturation = py.CoolProp.CoolProp.PropsSI('I', 'P', p_vapor, 'Q', 0, 'Ethanol');
    H_film =  py.CoolProp.CoolProp.PropsSI('I', 'T', T_film, 'P', p_vapor, 'Ethanol');

    T0 = 0;
    T_if = T_film;
    while(abs(T_if-T0) > 1)
        H_exhaust_film = H_exhaust + cp_exhaust*(T_if-T_exhaust); % J - enthalpy of combustion gas at temperature of film coolant;
        B = (H_exhaust_stag - H_exhaust_film) / (H_film_saturation - H_film);
        p_vapor = p_exhaust / (1+mol_film/(B*mol_exhaust)); % Vapor pressure of film coolant at interface with combustion gas
        T0 = T_if;
        T_if = py.CoolProp.CoolProp.PropsSI('T', 'P', p_vapor, 'Q', 0, 'Ethanol'); % K - temperature at combustion gas / film coolant interface
    end

    a = 10*A*(1+3*(surface_tension*X)^-0.8); % 1/m - heat-transfer augmentation factor for liquid surface roughness

    Stanton = h_exhaust / (dens_film*v_film*cp_film); % Stanton number of film coolant
    V = pi*d_avg*dens_exhaust*v_exhaust*Stanton*B*a; % kg/m*s - surface vaporization rate
    length = log( 1 + A*mdot_film/V) / A;
end
