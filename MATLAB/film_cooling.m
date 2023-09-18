function [T_film, liquid] = film_cooling(h_gas, T_recovery, T_saturation, T_film_prev, x, r1, mdot_film, film_injection_x, cp_film, T_injection, dA, boiling_energy, T_wall)

    if cumtrapz(h_gas.*(T_recovery-T_saturation).*dA) <= boiling_energy
        % wall region is covered by liquid film coolant
		liquid = 1;
        T_film = T_saturation; % K
	else
		liquid = 0;
        dT = dx * 1.628 * (2*pi*r1*h_gas/cp_film) * ((v_gas/v_film) * (1/(2*pi*r1)) * (cp_film/k_film))^0.125 * (mdot_film)^-0.875 * (T_recovery-T_film_prev);
        T_film = T_film_prev + dT; % K
    end

end

