%% Thermal Environment

% WALL TEMPERATURE SHOULD NOT BE INDEPENDENT
Twall = 500; % K - wall temperature, independent variable
Tf = Tc ./ (1 + (gamma-1)/2 * M.^2); % Free-Stream Temperature
Tref = Tf .* (1 + 0.032*M.^2 + 0.58*(Twall./Tf-1));

% UPDATE PROPERTIES TO BE FOR REFERENCE TEMPERATURE
cp_gas = cea.output.eql.cp; % Specific Heat, J/kg-K
cp_gas = interp1([0,x2_throat,x_exit], [cp_gas(1),cp_gas(2),cp_gas(3)], x); % Interpolate
visc = cea.output.eql.viscosity; % Dynamic Viscosity, Pa-s
visc = interp1([0,x2_throat,x_exit], [visc(1),visc(2),visc(3)], x); % Interpolate
Pr = struct2array(cea.output.eql.prandtl(:,1)); % Prandtl Number
Pr = interp1([0,x2_throat,x_exit], [Pr(1),Pr(2),Pr(3)], x); % Interpolate

% Bartz relation  
sigma = 1 ./ ( (0.5*(Twall/Tc)*(1+(gamma-1)/2 * M.^2)+0.5).^0.68 .* (1+(gamma-1)/2 * M.^2).^0.12 ); % Meaningless (?) coefficient
h_gas = (0.026./d_throat.^0.2).*(visc.^0.2 .* cp_gas ./ Pr.^0.6 ) .* (pc./c_star).^0.8 .* (d_throat./r_throat)^0.1 .* (0.5*d_throat./r1).^1.8 .* sigma; % W/m2-K - Convective Heat Transfer Coefficient

Tab = Tc.*(1+Pr.^0.33.*(gamma-1)/2.*M.^2)./(1+(gamma-1)/2.*M.^2); % Adiabatic Wall Temperature
A_gas = 2*pi*r1*dx/n_pipe; % Area of dx per coolant channel
q_gas = h_gas.*A_gas.*(Tab-Twall); % W - Heat Flux

%% Coolant Properties
% Ethanol is the coolant
density_cool =	density_fuel; % kg/m3 - ethanol at STP
cp_cool = cp_fuel; % J/kg-K - specific heat
kin_visc_cool = kin_visc_fuel; % m2/s - kinematic viscosity
k_cool = k_fuel; % W/m-K - thermal conductivity

%% Pressure Drop

p_cool = zeros(1,length(x));
v_cool = zeros(1,length(x)); % m/s - fluid bulk velocity
    v_cool(1) = mdot_fuel ./ (density_cool*n_pipe*w_pipe(1)*h_pipe); % m/s - fluid velocity
d_hydraulic = 4*w_pipe.*h_pipe./(2*w_pipe+2*h_pipe); % m - hydraulic diameter
Re_cool = zeros(1,length(x));
    Re_cool(1) = v_cool(1).*d_hydraulic(1)./kin_visc_cool; % Reynold's number

for i = 1:1:length(x)-1
    
    A1 = w_pipe(i)*h_pipe;
    A2 = w_pipe(i+1)*h_pipe;


    f_cool = (0.79 .* log(Re_cool) - 1.64).^(-2); % friction factor, smooth pipe approximation REVIEW IF THIS IS ACCURATE

    dp1(i) = 0.5*f_cool(i)*density_cool*v_cool(i)^2*dx./d_hydraulic(i); % Major Losses (channel wall roughness)

    d_ratio = d_hydraulic(i+1)/d_hydraulic(i);
    if A2 < A1
        K = 0.5 - 0.167*d_ratio - 0.125*d_ratio^2 - 0.208*d_ratio^3; % Contraction
    else
        K = (d_ratio^-2 - 1)^2; % Expansion
    end

    dp2(i) = 0.5*K*density_cool*v_cool(i)^2; % Minor Losses (diameter changes)
    dp3(i) = (2 / (n_pipe*(A1+A2)) ) * (1 / (density_cool*A1*n_pipe) - 1 / (density_cool*A2*n_pipe)); % Acceleration Losses

    dp(i) = dp1(i)+dp2(i)+dp3(i);

    p_cool(i+1) = p_cool(i) - dp(i);

    v_cool(i+1) = sqrt(v_cool(i)^2 - 2*dp(i)/density_cool); % THIS LINE BREAKS AT I = 367
    Re_cool(i+1) = v_cool(i+1)*d_hydraulic(i+1)/kin_visc_cool; % Reynold's 

    if imag(v_cool(i+1)) ~= 0
        error("Coolant Flow stalls out at %g m from injector (%g m before exit)", x(i), x_exit-x(i));
    end
end

p_cool = p_cool + pc - p_cool(end); % Pa - ending pressure EXPAND TO ACCOUNT FOR PLUMBING AND INJECTOR STIFFNESS 




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
m_coeff = sqrt(2*h_gas/(k_al6061*fin_thickness)); % coefficient for fin efficency equation
fin_eff = tanh(m_coeff*dx)/m_coeff*dx; % fin efficiency for each length step
area_cool = (2*fin_height*fin_eff+fin_thickness)*dx; % m2 - adjusted contact area of coolant on channel walls

%% Thermal Balance

T_cool = zeros(1,length(x));
T_cool(1) = Tamb; % K - coolant inlet temperature
T_wall_cold = zeros(1,length(x));
T_wall_hot = zeros(1,length(x));

for i = 1:1:length(x)

    if T_cool(i) > boiling_fuel
        error("Coolant starts boiling at %g m from injector (%g m before exit)", x(i+1), x_exit-x(i+1));
    end

    T_wall_cold(i) = q_gas(i) ./ (h_cool(i)*area_cool*n_pipe) + T_cool(i);
    T_wall_hot(i) = q_gas(i)*(r2(i)-r1(i))/(k_al6061*pi*2*r1(i)*dx) + T_wall_cold(i);
    
    dT = q_gas(i)*n_pipe / (cp_cool*mdot_fuel); % Coolant Temperature change
    if i < length(x)
        T_cool(i+1) = T_cool(i) + dT;
    end
    
end







