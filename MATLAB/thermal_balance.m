%% Combustion Thermal Environment
Tf = Tc ./ (1 + (gamma-1)/2 * M.^2); % Free-Stream Temperature

%% Wall Conduction

fin_thickness = gap_pipe; % m
fin_height = h_pipe; % m
fin_length  = x_exit; % m

%% Film Flow

mdot_film = film_fraction*mdot_cc;
mdot_gas = (1-film_fraction)*mdot_cc;

cp_film_liq = py.CoolProp.CoolProp.PropsSI('C', 'P', pc, 'Q', 0, 'Ethanol'); % J/kg-s % QUESTIONABLE REFERENCE TEMPERATURE
visc_film_liq = py.CoolProp.CoolProp.PropsSI('V', 'P', pc, 'Q', 0, 'Ethanol'); % Pa-s % QUESTIONABLE REFERENCE TEMPERATURE
dens_film_liq = py.CoolProp.CoolProp.PropsSI('D', 'P', pc, 'Q', 0, 'Ethanol'); % kg/m3 %  QUESTIONABLE REFERENCE TEMPERATURE
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

% Re_film = mdot_film./(pi*2*r1(1)*visc_film_liq); % Reynold's number of film PROBABLY INVALID AROUND THROAT AND NOZZLE
% film_thickness = 0.376*x.*Re_film.^-0.2; % m - boundary layer thickness of the film coolant (fuel)

% film_eff_theory = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000; 1, 0.91, 0.77, 0.69, 0.64, 0.59, 0.56, 0.53, 0.52]; % Correlation between Reynold's and film efficiency (eta), 1969 paper
% film_eff = interp1(film_eff_theory(1,:), film_eff_theory(2,:), Re_film);

% film_length = film_eff*mdot_film*(cp_film_liq*(T_boil_film-T_film_initial) + heat_vaporization_film) ./ ( 2*pi*r1(1)*h_gas*(T_stagnation-T_boil_film));

%% Thermal Balance

T_cool = T_amb*ones(1,length(x));
T_wall_cold = zeros(1,length(x)); % K - wall temperature on coolant side
T_wall_hot = T_amb*ones(1,length(x)); % K - wall temperature on chamber interior side
Tab = zeros(1,length(x));
Tref = zeros(1,length(x)); % K - reference tempeature for combustion gas thermophysical properties
h_film_gas = zeros(1,length(x));
h_film_liq = zeros(1,length(x));
T_film = zeros(1,length(x));

p_gas =1E5*cea.output.eql.pressure; % Pa - pressure of combustion gases
p_gas = interp1([0,x2_throat,x_exit], [p_gas(1),p_gas(2),p_gas(3)], x); % Interpolate

products = convertCharsToStrings(cea.output.eql.products);
fractions = convertCharsToStrings(cea.output.eql.fractions);
fractions = interp1( [0; x2_throat; x_exit], fractions', x )' ; % Interpolate

% Preallocate several arrays
cp_gas = zeros(1,length(x)); % J/kg - specific heat capacity (at constant pressure) of combustion gases
visc_gas = zeros(1,length(x)); % Pa-s - dynamic viscosity of combustion gases
cond_gas = zeros(1,length(x));  % W/m-K - thermal conductivity of combustion gases
h_gas = zeros(1,length(x));
q_gas = zeros(1,length(x));
Pr_gas = zeros(1,length(x));
Nu_cool = zeros(1,length(x));
h_cool = zeros(1,length(x));

if flow_direction == 1
    flow = 1:1:length(x);
elseif flow_direction == -1
    flow = length(x):-1:1;
end

warned = false;
start = true;

T_input = T_amb;
for i = flow % Find heat and temperatures along engine contour
    % Coolant Properties
    cp_cool = py.CoolProp.CoolProp.PropsSI("C", "T", T_cool(i), "P", p_cool(i), "Ethanol"); % J/kg-K
    visc_cool = py.CoolProp.CoolProp.PropsSI("V", "T", T_cool(i), "P", p_cool(i), "Ethanol"); % Pa-s
    k_cool = py.CoolProp.CoolProp.PropsSI("L", "T", T_cool(i), "P", p_cool(i), "Ethanol"); % W/m-K
    Pr_cool = cp_cool*visc_cool/k_cool; % Prandtl Number

    if (Re_cool(i)<3000) % Really should be < 2300, but using laminar flow in transition region for conservatism
        lam_square = [1, 1.43, 2, 3, 4, 8, 100000000000; 3.61, 3.73, 4.12, 4.79, 5.33, 6.49, 8.23]; % Nusselt number for laminar flow in square pipes (Lavine, table 8.1)
        Nu_cool(i) = interp1(lam_square(1,:),lam_square(2,:), max(h_pipe./w_pipe(i), w_pipe(i)./h_pipe));
        if ~warned
            fprintf("Warning: regen coolant goes laminar at %g m from injector (%g m before exit)\n", x(i), x_exit-x(i));
            warned = true;
        end
    elseif (Pr_cool >= 0.7 && Pr_cool <= 2000 && Re_cool(i) >= 3000 && Re_cool(i) < 5E6)        
        Nu_cool(i) = (f_cool(i)/8).*(Re_cool(i)-1000).*Pr_cool/(1+12.7*(f_cool(i)/8).^0.5.*(Pr_cool.^(2/3)-1)); % Nusselt Number, Gnielenski correlation (turbulent)
    else       
        x(i)
        Re_cool(i)
        Pr_cool
        error("No correlation for coolant Reynold's and/or Prandtl number");
    end
    
    h_cool(i) = Nu_cool(i)*k_cool./d_hydraulic(i); % W/m2-K - Convective heat transfer coefficient

    if ~start
        [cp_gas(i), visc_gas(i), cond_gas(i)] = mixture(products, fractions(:,i-flow_direction), Tref(i-flow_direction), p_gas(i-flow_direction)); % find properties of the mix of combustion gases at this location
    end

    % Heat Balancing
    residual = 2;

    while abs(residual) > 1 % Newton-Raphson Iteration

        rise = 0;
        run = 0;

        for j = -1:2:1 % Find two neighboring points so a local derivative can be taken

            Tref(i) = Tf(i) .* (1 + 0.032*M(i).^2 + 0.58*(T_input./Tf(i)-1)); % K - reference temperature for gas properties
            
            if start
                [cp_gas(i), visc_gas(i), cond_gas(i)] = mixture(products, fractions(:,i), Tref(i), p_gas(i)); % find properties of the mix of combustion gases at this location
            end

            Pr_gas(i) = cp_gas(i)*visc_gas(i)/cond_gas(i); % Prandtl Number
        
            % Bartz relation
            sigma = 1 ./ ( (0.5*(T_input/Tc)*(1+(gamma-1)/2 * M(i).^2)+0.5).^0.68 .* (1+(gamma-1)/2 * M(i).^2).^0.12 ); % Bartz Equation coefficent accounting for density & viscosity variation across boundary layer
            %h_gas(i) = (0.026./d_throat.^0.2).*(visc_gas(i).^0.2 .* cp_gas(i) ./ Pr_gas(i).^0.6 ) .* (pc./c_star).^0.8 .* (d_throat./rc_throat)^0.1 .* (0.5*d_throat./r1(i)).^1.8 .* sigma(i); % W/m2-K - Convective Heat Transfer Coefficient - nonfilm case
            h_film_gas(i) = (0.026./d_throat.^0.2).*(visc_film_avg.^0.2 .* cp_film_avg ./ Pr_film_avg.^(2/3) ) .* (pc/c_star).^0.8 .* (d_throat./rc_throat)^0.1 .* (0.5*d_throat./r1(i)).^1.8 .* sigma; % W/m2-K - heat transfer between exhaust and film AND gas film to wall  - film case

            m_coeff = sqrt(2*h_cool(i)/(k_wall*fin_thickness)); % coefficient for fin efficency equation
            fin_eff = tanh(m_coeff*(h_pipe+fin_thickness/2))/(m_coeff*(h_pipe+fin_thickness/2)); % fin efficiency for each length step
            A_cool = (2*fin_height*fin_eff+w_pipe(i))*dx; % m2 - adjusted contact area of coolant on channel walls 

            % TAB ONLY APPLIES TO CHAMBER, NOT NOZZLE
            Tab(i) = Tc.*(1+Pr_gas(i).^0.33.*(gamma-1)/2.*M(i).^2)./(1+(gamma-1)/2.*M(i).^2); % K - Adiabatic Wall Temperature. 0.33 --> 0.5 IF GAS FLOW IS LAMINAR
            boiling_fuel = py.CoolProp.CoolProp.PropsSI( 'T', 'P', p_cool(i), 'Q', 0, 'Ethanol' ); % K - Ethanol Saturation Temperature
            A_gas = 2*pi*r1(i)*dx/n_pipe(i); % m2 - Area of dx per coolant channel

            %% FILM COOLING
            % temporary placeholder variables
            v_boundary = 1; % m/s - Axial stream velocity of combustion gases at edge of boundary layer
            v_center = 100; % m/s - Average axial stream velocity of combustion gases
            v_avg = 1; % m/s - Axial stream velocity of combustion gases at the center line of the thrust chamber

            film_friction = 0.015; % Applicable friction coefficient for the two phase flow between combustion gases and liquid film coolant 
            film_eff = 0.5; % Film-cooling efficiency - ranges from 0.3 to 0. per Huzel and Huang
            
            singularity = Tc - heat_vaporization_film/cp_film_liq;
            delta = 1e-10; % small value to avoid singularity, adjust as needed
            
            if Tab(i) - (singularity+delta) > 0 % checks which side of the singularity the zero is on
                Tab(i);
                T_film(i) = fzero(@(T) film_function(T, film_eff, mdot_film, mdot_gas, r1(i), A_gas, n_pipe(i), v_boundary, film_friction, v_avg, v_center, cp_film_gas, cp_gas(i), Tab(i), cp_film_liq, Tc, heat_vaporization_film), [singularity + delta,10000]);
                T_film(i);
            else
                T_film(i) = fzero(@(T) film_function(T, film_eff, mdot_film, mdot_gas, r1(i), A_gas, n_pipe(i), v_boundary, film_friction, v_avg, v_center, cp_film_gas, cp_gas(i), Tab(i), cp_film_liq, Tc, heat_vaporization_film), [0,singularity - delta]);
            end

            q_gas(i) = h_film_gas(i)*A_gas*(T_film(i)-T_input); % W - Heat Flux per coolant channel, negative because it flows out of the gas

            %%

            % Check if the coolant boils. Gaseous coolant will be extremely ineffective, and harder to analyze
            if T_cool(i) > boiling_fuel
                error("Coolant starts boiling at %g m from injector (%g m from exit plane)", x(i), x_exit-x(i));
            end

            T_wall_cold(i) = q_gas(i) ./ (10*h_cool(i)*A_cool) + T_cool(i); % K - cold wall temperature, q_gas remains negative because it flows out of wall
            T_output = q_gas(i)*(r2(i)-r1(i))/(k_wall*A_gas) + T_wall_cold(i) - T_input; % K - temporary hot wall temperature DIFFERENCE for purposes of iteration, q_gas turned positive because it flows into wall
            
            rise = rise + T_output*j;
            run = run + T_input*j;

            T_input = T_input + 0.1; % 0.1 is some small lateral move so an approximate derivative can be taken
        end
        
        slope = rise/run;
        residual =  -T_output / slope;
        T_input = T_input + residual;

    end

    T_wall_hot(i) = T_input;
    start = false;

    if i ~= flow(end)
        T_wall_hot(i+flow_direction) = T_wall_hot(i); % Approximation, gets overwritten at every subsequent loop
    end
    
    dT = q_gas(i)*n_pipe(i) / (cp_cool*mdot_fuel_cc); % Coolant Temperature change at one channel MAYBE THIS SHOULD BE CV
    if i ~= flow(end)
        T_cool(i+flow_direction) = T_cool(i) + dT;
    end
end

% Film cooling function
function y = film_function(T, film_eff, mdot_film, mdot_gas, r1, A_gas, n_pipe, v_boundary, film_friction, v_avg, v_center, cp_film_gas, cp_gas, Taw, cp_film_liq, Tc, heat_vaporization_film)
    a = (2*v_boundary/(film_friction*v_avg)); % idk what this constant is
    b = v_center/v_boundary-1; % or this one
    H = cp_film_gas * (Taw - T) ./ (cp_film_liq*(T - Tc) + heat_vaporization_film); % film coolant enthalpy ratio
    y = film_eff * ( mdot_film/mdot_gas ) * ( pi*r1^2/(A_gas*n_pipe) ) * a * ( 1+ b^(cp_film_gas/cp_gas) )  -  H; % Huzel and Huang page 99
end





