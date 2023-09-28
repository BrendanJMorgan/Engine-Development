%% Objective: finds the temperatue distributions of the combustion chamber and nozzle walls
%% Combustion Thermal Environment
T_free = Tc ./ (1 + (gamma_gas-1)/2 .* M_gas.^2); % Free-Stream Temperature

%% Wall Conduction

fin_thickness = gap_pipe; % m
fin_height = h_pipe; % m
fin_length = x_exit; % m

A_gas = 2*pi*r1*dx./n_pipe; % m2 - Area of dx per coolant channel
t_wall = r2 - r1; % Wall thickness (projected vertically)

%% Film Flow

mdot_film = [0.5, 0.5]*film_fraction*mdot_cc;
mdot_gas = (1-film_fraction)*mdot_cc;

T_injection = T_amb; % K - PLACEHOLDER UNTIL REPLACED WITH REGEN TEMPERATURE AT THAT LOCATION
heat_vap_film = py.CoolProp.CoolProp.PropsSI('H','P',pc,'Q',1,'Ethanol') - py.CoolProp.CoolProp.PropsSI('H','P',pc,'Q',0,'Ethanol'); % J/kg - heat of vaporization of film coolant
cp_film = py.CoolProp.CoolProp.PropsSI('C','P',pc,'Q',0,'Ethanol');
T_sat_film = py.CoolProp.CoolProp.PropsSI('T','P',pc,'Q',0,'Ethanol'); % THIS SHOULD BE DEPENDENT ON PRESSURE IN THROAT/NOZZLE
T_film = T_amb*ones(1,length(x));
film_injection_x = [0, l_chamber]; % m - film cooling orifices around perimeter of injector and along bottom edge of chamber wall

%% Thermal Balance


products = convertCharsToStrings(cea.output.eql.products);
fractions = convertCharsToStrings(cea.output.eql.fractions);
fractions = interp1( [0; x_combustor; x2_throat; x_exit], fractions', x )' ; % Interpolate

% Preallocate several arrays
T_wall_cold = zeros(1,length(x)); % K - wall temperature on coolant side
% T_wall_hot = T_amb*ones(1,length(x)); % K - wall temperature on chamber interior side
Tab = zeros(1,length(x));
T_cool = T_amb*ones(1,length(x));
T_ab = T_amb*ones(1,length(x));
cp_gas = zeros(1,length(x)); % J/kg - specific heat capacity (at constant pressure) of combustion gases
visc_gas = zeros(1,length(x)); % Pa-s - dynamic viscosity of combustion gases
cond_gas = zeros(1,length(x));  % W/m-K - thermal conductivity of combustion gases
h_gas = zeros(1,length(x)); % W / m2-K - convective heat transfer coefficient - MUST be initialized as all zeros for film_cooling() to work properly
q_gas = zeros(1,length(x));
Pr_gas = zeros(1,length(x));
Nu_cool = zeros(1,length(x));
h_cool = zeros(1,length(x));

warned = false;
start = true;

T_ref = T_free .* (1 + 0.032*M_gas.^2 + 0.58*(T_wall_hot./T_free-1)); % K - reference temperature for combustion gas thermophysical properties

for i = 1:1:length(x)
	[cp_gas(i), visc_gas(i), cond_gas(i)] = mixture(products, fractions(:,i), T_ref(i), p_gas(i)); % find properties of the mix of combustion gases at this location
end
Pr_gas = cp_gas.*visc_gas./cond_gas; % Prandtl Number

% Bartz relation
sigma = 1 ./ ( (0.5*(T_wall_hot./Tc).*(1+(gamma_gas-1)/2 .* M_gas.^2)+0.5).^0.68 .* (1+(gamma_gas-1)/2 .* M_gas.^2).^0.12 ); % Bartz Equation coefficient accounting for density & viscosity variation across boundary layer
h_gas = (0.026./d_throat.^0.2).*(visc_gas.^0.2 .* cp_gas ./ Pr_gas.^0.6 ) .* (p_gas./c_star).^0.8 .* (d_throat./rc_throat)^0.1 .* (0.5*d_throat./r1).^1.8 .* sigma; % W/m2-K - Convective Heat Transfer Coefficient - nonfilm case

T_ab = Tc.*(1+Pr_gas.^0.33.*(gamma_gas-1)/2.*M_gas.^2)./(1+(gamma_gas-1)/2.*M_gas.^2); % K - Adiabatic Wall Temperature. 0.33 ==> 0.5 IF GAS FLOW IS LAMINAR % TAB ONLY APPLIES TO CHAMBER, NOT NOZZLE??
T_recovery = c_star_eff^2 .* T_free + Pr_gas.^(1/3) .* (Tc - T_free); % K - recovery temperature of combustion gas (fluid state if it slowed to velocity of body it is contacting, considering friction and heat transfer).

film_phase = zeros(1,length(x)); % 0 means film coolant is a gas, and 1 means it is a liquid at that location

for i = 1:1:length(film_injection_x)
	% Liquid Film Coolant Distribution
	boiling_energy = injection_efficiency*mdot_film(i)*(cp_film*(T_sat_film - T_injection) + heat_vap_film); % J - energy required to boil all liquid film coolant into gas
	liquid_start = floor(film_injection_x(i)/dx)+1;
	accumulated_energy = cumtrapz( h_gas(1,liquid_start:end) .* ( T_recovery(1,liquid_start:end)-T_sat_film*ones(1,length(x)-liquid_start+1) ) ).*dA(liquid_start:end); % J - total heat energy outputted y combustion gas up to each point
	liquid_end = find(accumulated_energy >= boiling_energy, 1, 'first') + liquid_start - 1;
	if isempty(liquid_end)
		liquid_end = length(accumulated_energy) + liquid_start - 1;
	end
	film_phase(liquid_start:liquid_end) = ones(1, liquid_end-liquid_start+1);
	T_film(liquid_start:liquid_end) = T_sat_film;

	% Gaseous Film Coolant Distribution
	for j = liquid_end+1:1:length(x)
		cp_film = py.CoolProp.CoolProp.PropsSI('C','P',p_gas(j),'T',T_film(j-1),'Ethanol'); % J/kg-K
		k_film = py.CoolProp.CoolProp.PropsSI('L','P',p_gas(j),'T',T_film(j-1),'Ethanol'); % J/kg-K
		dT = dx * 1.628 * (2*pi*r1(j)*h_gas(j)/cp_film) * ((v_gas(j)/v_injection) * (1/(2*pi*r1(j))) * (cp_film/k_film))^0.125 * (mdot_film(i))^-0.875 * (T_recovery(j)-T_film(j-1)); % K
		T_film(j) = T_film(j-1) + dT;
	end
end

%%

if flow_direction == 1
    flow = 1:1:length(x);
elseif flow_direction == -1
    flow = length(x):-1:1;
end

T_input = T_amb;
for i = flow % Find heat and temperatures along engine contour
	%% Regenerative Coolant Properties
    cp_cool = py.CoolProp.CoolProp.PropsSI("C", "T", T_cool(i), "P", p_cool(i), "Ethanol"); % J/kg-K
    visc_cool = py.CoolProp.CoolProp.PropsSI("V", "T", T_cool(i), "P", p_cool(i), "Ethanol"); % Pa-s
    k_cool = py.CoolProp.CoolProp.PropsSI("L", "T", T_cool(i), "P", p_cool(i), "Ethanol"); % W/m-K
    Pr_cool = cp_cool*visc_cool/k_cool; % Prandtl Number

    if (Re_cool(i)<3000) % Really should be < 2300, but using laminar flow in transition region for conservatism
        lam_square = [1, 1.43, 2, 3, 4, 8, 100000000000; 3.61, 3.73, 4.12, 4.79, 5.33, 6.49, 8.23]; % Nusselt number for laminar flow in square pipes (Lavine, table 8.1)
        Nu_cool(i) = interp1(lam_square(1,:),lam_square(2,:), max(h_pipe./w_pipe(i), w_pipe(i)./h_pipe));
        if ~warned
            fprintf("Warning: regen coolant goes laminar at %g m from injector (%g m from exit plane)\n", x(i), x_exit-x(i));
            warned = true;
        end
    elseif (Pr_cool >= 0.7 && Pr_cool <= 2000 && Re_cool(i) >= 3000 && Re_cool(i) < 5E6)
        Nu_cool(i) = (f_cool(i)/8).*(Re_cool(i)-1000).*Pr_cool/(1+12.7*(f_cool(i)/8).^0.5.*(Pr_cool.^(2/3)-1)); % Nusselt Number, Gnielenski correlation (turbulent)
    else
        error("No correlation for coolant Reynold's and/or Prandtl number");
    end
    
    h_cool(i) = Nu_cool(i)*k_cool./d_hydraulic(i); % W/m2-K - Convective heat transfer coefficient
	m_coeff = sqrt(2*h_cool/(k_wall*fin_thickness)); % coefficient for fin efficency equation
    fin_eff = tanh(m_coeff*(h_pipe+fin_thickness/2))/(m_coeff*(h_pipe+fin_thickness/2)); % fin efficiency for each length step
    A_cool = (2*fin_height*fin_eff+w_pipe)*dx; % m2 - adjusted contact area of coolant on channel walls

    %% Heat Balancing

	if film_phase(i) == 1 % Liquid State
		q_wall(i) = ( T_film(i) - T_cool(i) ) / ( 1/(h_cool(i)*A_cool(i)) + 1/(k_wall*A_gas(i)/t_wall(i)) ); % W
		T_wall_cold(i) = q_wall(i)/(h_cool(i)*A_cool(i)) + T_cool(i); % K  
		T_wall_hot(i) = q_wall(i)/(k_wall*A_gas(i)/t_wall(i)) + T_wall_cold(i); % K
		if round(T_wall_hot(i)) ~= round(T_film(i))
			error("Wall temperature and film temperature do not match in liquid film cooled region at %g m from injector (%g m from exit plane)\n", x(i), x_exit-x(i));
		end
	elseif film_phase(i) == 0
		q_wall(i) = (T_film(i)-T_cool(i)) / ( 1/(h_gas(i)*A_gas(i)) + 1/(k_wall*A_gas(i)/t_wall(i)) + 1/(h_cool(i)*A_cool(i)) ); % W
		T_wall_cold(i) = q_wall(i)/(h_cool(i)*A_cool(i)) + T_cool(i); % K
		T_wall_hot(i)= q_wall(i)/(k_wall*A_gas(i)/t_wall(i)) + T_wall_cold(i); % K
	end

	dT = q_wall(i)*n_pipe(i) / (cp_cool*mdot_fuel_cc); % Coolant Temperature change at one channel MAYBE THIS SHOULD BE CV
    if i ~= flow(end)
        T_cool(i+flow_direction) = T_cool(i) + dT;
	end

    % Check if the coolant boils. Gaseous coolant will be extremely ineffective, and harder to analyze
	boiling_cool = py.CoolProp.CoolProp.PropsSI( 'T', 'P', p_cool(i), 'Q', 0, 'Ethanol'); % K - Ethanol Saturation Temperature
    if T_cool(i) > boiling_cool
        error("Coolant starts boiling at %g m from injector (%g m from exit plane)", x(i), x_exit-x(i));
	end

end