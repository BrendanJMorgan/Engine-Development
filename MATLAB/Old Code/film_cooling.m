%% Combustion Thermal Environment
T_free = Tc ./ (1 + (gamma_gas-1)/2 .* M_gas.^2); % Free-Stream Temperature

%% Wall Conduction

fin_thickness = gap_pipe; % m
fin_height = h_pipe; % m
fin_length = x_exit; % m

%% Film Flow

mdot_film = film_fraction*mdot_cc;
mdot_gas = (1-film_fraction)*mdot_cc;

T_cool = T_amb*ones(1,length(x));

film_cooled = zeros(1,length(x)); % boolean - true values signify there is liquid film cooling at that location
film_injection_x = [0, l_chamber]; % m - film cooling orifices around perimeter of injector and along bottom edge of chamber wall

% Compute Euclidean distances
delta_x = diff(x);
delta_r1 = diff(r1);
segment_lengths = sqrt(delta_x.^2 + delta_r1.^2);
total_path_lengths = [0, cumsum(segment_lengths)];
    
for j = 1:1:length(film_injection_x)
    i = floor(film_injection_x(j)/dx)+1;
    
    T_film = T_cool(i);
    Pr_film_liq = py.CoolProp.CoolProp.PropsSI('Prandtl', 'P', pc, 'Q', 0, 'Ethanol');
    visc_film_liq = py.CoolProp.CoolProp.PropsSI('V', 'P', pc, 'Q', 0, 'Ethanol'); % Pa-s
    d_avg = 2*r1(i); % m - CORRECT TO A BETTER AVERAGE
    
    v_gas(i)
    film_length(j) = film_cooled_length(0.5*mdot_film, T_film, Pr_film_liq, visc_film_liq, d_avg, l_chamber, dens_gas(i), v_gas(i), T_free(i), H_gas(i), p_gas(i), cp_gas(i), mol_gas(i));

    % Find where film injection starts
    injection_idx = find(x >= film_injection_x(j), 1);
    if isempty(injection_idx)
        error('Film injection point not in range.');
    else
        start_length = total_path_lengths(injection_idx);
        end_length = start_length + film_length(j);
    end
    
    % Create the cooled array
    cooled(j,:) = (total_path_lengths >= start_length) & (total_path_lengths <= end_length);

end

cooled = any(cooled);










