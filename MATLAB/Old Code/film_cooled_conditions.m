function T_ab_cooled = film_cooled_conditions()
    
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

    liquid_flow_ratio = (1+B)/(0.6*B) - 1;
    flow_ratio = (mdot_total-mdot_film)/mdot_film * (2*psi_L * xbar / r_injection * sqrt(1-liquid_flow_ratio) - (psi_L*xbar/r_injection)^2) + liquid_flow_ratio
    
    if  flow_ratio < (liquid_flow_ratio+0.6)
        shape_factor = 0.6 + 0.263*(flow_ratio-liuid_flow_ratio);
    else
        shape_factor = 0.758;
    end

    effectiveness = 1 / (shape_factor*(1+flow_ratio)); % effectiveness of film coolant downstream of injection point

    % Reactive Model
    H_aw = H_gas_stag - effectiveness*(H_gas_stag - H_film) - (1-Prandtl^(1/3))*(H_gas_stag-H_gas)
    vapor_fraction = exp(-epsilon*tau)

end