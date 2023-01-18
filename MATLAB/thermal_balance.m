%% Coolant Thermal Properties
% Ethanol is the coolant
cp_cool = cp_fuel; % J/kg-K - specific heat
k_cool = k_fuel; % W/m-K - thermal conductivity

%% Combustion Thermal Environment

% WALL TEMPERATURE SHOULD NOT BE INDEPENDENT
%Twall = 500; % K - wall temperature, independent variable
Tf = Tc ./ (1 + (gamma-1)/2 * M.^2); % Free-Stream Temperature


%% Coolant Loop Flow

Pr_cool = cp_cool*density_cool*kin_visc_cool/k_cool; % Prandtl Number
 
Nu_cool = zeros(1,length(x));
lam_square = [1, 1.43, 2, 3, 4, 8, 100000000000; 3.61, 3.73, 4.12, 4.79, 5.33, 6.49, 8.23]; % Nusselt number for laminar flow in square pipes (Lavine, table 8.1)
for i = 1:1:length(x)    
    if (Re_cool(i)<3000) % Really should be < 2300, but using laminar flow in transition region for conservatism        
        Nu_cool(i) = interp1(lam_square(1,:),lam_square(2,:), max(h_pipe./w_pipe(i), w_pipe(i)./h_pipe));
    elseif (Pr_cool >= 0.7 && Pr_cool <= 2000 && Re_cool(i) >= 3000 && Re_cool(i) < 5E6)        
        Nu_cool(i) = (f_cool(i)/8).*(Re_cool(i)-1000).*Pr_cool/(1+12.7*(f_cool(i)/8).^0.5.*(Pr_cool.^(2/3)-1)); % Nusselt Number, Gnielenski correlation (turbulent)
    else       
        x(i)
        Re_cool(i)
        Pr_cool
        error("No correlation for coolant Reynold's and/or Prandtl number");
    end
end

h_cool = Nu_cool*k_cool./d_hydraulic; % W/m2-K - Convective heat transfer coefficient

%% Wall Conduction

fin_thickness = gap_pipe; % m
fin_height = h_pipe; % m
fin_length  = x_exit; % m

%% Thermal Balance

T_cool = zeros(1,length(x));
T_cool(1) = Tamb; % K - coolant inlet temperature
T_wall_cold = zeros(1,length(x)); % K - wall temperature on coolant side
T_wall_hot = zeros(1,length(x)); % K - wall temperature on chamber interior side
    T_wall_hot(1) = Tamb; % K - initial guess
Tab = zeros(1,length(x));
Tref = zeros(1,length(x));

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


for i = 1:1:length(x) % Find heat and temperatures along engine contour

    residual = 20;
    while abs(residual) > 1 % Newton-Raphson Iteration

        rise = 0;
        run = 0;

        T_input = T_wall_hot(i); % K
        for j = -1:2:1 % Find two neighboring points so a local derivative can be taken
            
            Tref(i) = Tf(i) .* (1 + 0.032*M(i).^2 + 0.58*(T_input./Tf(i)-1)); % K - reference temperature for gas properties
            
            [cp_gas(i), visc_gas(i), cond_gas(i)] = mixture(products, fractions(:,i), Tref(i), p_gas(i)); % find properties of the mix of combustion gases at this location
        
            Pr_gas(i) = cp_gas(i)*visc_gas(i)/cond_gas(i); % Prandtl Number
        
            % Bartz relation  
            sigma = 1 ./ ( (0.5*(T_input/Tc)*(1+(gamma-1)/2 * M(i).^2)+0.5).^0.68 .* (1+(gamma-1)/2 * M(i).^2).^0.12 ); % Meaningless (?) coefficient
            h_gas(i) = (0.026./d_throat.^0.2).*(visc_gas(i).^0.2 .* cp_gas(i) ./ Pr_gas(i).^0.6 ) .* (pc./c_star).^0.8 .* (d_throat./r_throat)^0.1 .* (0.5*d_throat./r1(i)).^1.8 .* sigma; % W/m2-K - Convective Heat Transfer Coefficient    
        
            m_coeff = sqrt(2*h_cool(i)/(k_al6061*fin_thickness)); % coefficient for fin efficency equation
            fin_eff = tanh(m_coeff*(h_pipe+fin_thickness/2))/(m_coeff*(h_pipe+fin_thickness/2)); % fin efficiency for each length step
            A_cool = (2*fin_height*fin_eff+w_pipe(i))*dx; % m2 - adjusted contact area of coolant on channel walls
        
            Tab(i) = Tc.*(1+Pr_gas(i).^0.33.*(gamma-1)/2.*M(i).^2)./(1+(gamma-1)/2.*M(i).^2); % Adiabatic Wall Temperature. 0.33 --> 0.5 IF GAS FLOW IS LAMINAR
            A_gas = 2*pi*r1(i)*dx/n_pipe(i); % m2 - Area of dx per coolant channel
            q_gas(i) = h_gas(i)*A_gas*(Tab(i)-T_input); % W - Heat Flux per coolant channel, negative because it flows out of the gas
        
            boiling_fuel = py.CoolProp.CoolProp.PropsSI( 'T', 'P', p_cool(i), 'Q', 0, 'Ethanol' ); % K - Ethanol Saturation Temperature
        
            % Check if the coolant boils. Gaseous coolant will be extremely ineffective, and harder to analyze
            if T_cool(i) > boiling_fuel
                error("Coolant starts boiling at %g m from injector (%g m before exit)", x(i), x_exit-x(i));
            end

            T_wall_cold(i) = q_gas(i) ./ (10*h_cool(i)*A_cool) + T_cool(i); % K = cold wall temperature, q_gas remains negative because it flows out of wall
            T_output = q_gas(i)*(r2(i)-r1(i))/(k_al6061*A_gas) + T_wall_cold(i) - T_input; % K - temporary hot wall temperature for purposes of iteration, q_gas turned positive because it flows into wall
            
            rise = rise + T_output*j;
            run = run + T_input*j;

            T_input = T_input + 0.1;

        end
        
        slope = rise/run;
        residual =  - T_output / slope;
        T_wall_hot(i) = T_input + residual;

    end

    if i<length(x)
        T_wall_hot(i+1) = T_wall_hot(i); % Approximation, gets overwritten at every subsequent loop
    end
    
    dT = q_gas(i)*n_pipe(i) / (cp_cool*mdot_fuel_cc); % Coolant Temperature change at one channel
    if i < length(x)
        T_cool(i+1) = T_cool(i) + dT;
    end
    
end







