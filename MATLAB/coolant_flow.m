%% Coolant Bulk Properties
% Ethanol is the coolant
density_cool =	density_fuel; % kg/m3 - ethanol at STP % ADJUST WITH COOLPROP
kin_visc_cool = kin_visc_fuel; % m2/s - kinematic viscosity % ADJUST WITH COOLPROP


%% Pressure Drop

p_cool = zeros(1,length(x));
v_cool = zeros(1,length(x)); % m/s - fluid bulk velocity
    v_cool(1) = mdot_fuel_cc ./ (density_cool*n_pipe(1)*w_pipe(1)*h_pipe); % m/s - fluid velocity
d_hydraulic = 4*w_pipe.*h_pipe./(2*w_pipe+2*h_pipe); % m - hydraulic diameter
Re_cool = zeros(1,length(x));
    Re_cool(1) = v_cool(1).*d_hydraulic(1)./kin_visc_cool; % Reynold's number
f_cool = zeros(1,length(x));
f_cool(1) = (0.79 .* log(Re_cool(1)) - 1.64).^(-2); % friction factor, smooth pipe approximation REVIEW IF THIS IS ACCURATE

for i = 1:1:length(x)-1

    dp1 = 0.5*f_cool(i)*density_cool*v_cool(i)^2*dx./d_hydraulic(i); % Major Losses (channel wall roughness)
    
    if n_pipe(i+1) < n_pipe(i) % For when channels merge
        d_hydraulic_mod = 4*0.5*w_pipe(i+1).*h_pipe./(2*0.5*w_pipe(i+1)+2*h_pipe); % m - half-channel hydraulic diameter
        d_ratio = d_hydraulic_mod/d_hydraulic(i);
        A1 = w_pipe(i)*h_pipe;
        A2 = 0.5*w_pipe(i+1)*h_pipe;
    elseif n_pipe(i+1) > n_pipe(i) % For when channels divide
        d_hydraulic_mod = 4*0.5*w_pipe(i).*h_pipe./(2*w_pipe(i)+h_pipe); % m - half-channel hydraulic diameter
        d_ratio = d_hydraulic(i+1)/d_hydraulic_mod;
        A1 = 0.5*w_pipe(i)*h_pipe;
        A2 = w_pipe(i+1)*h_pipe;
    else % Everything else (vast majority of the time)
        d_ratio = d_hydraulic(i+1)/d_hydraulic(i);
        A1 = w_pipe(i)*h_pipe;
        A2 = w_pipe(i+1)*h_pipe;
    end

    if d_ratio < 1
        K = 0.5 - 0.167*d_ratio - 0.125*d_ratio^2 - 0.208*d_ratio^3; % Contraction
    else
        K = (d_ratio^-2 - 1)^2; % Expansion
    end

    dp2 = 0.5*K*density_cool*v_cool(i)^2; % Pa - Minor Losses (diameter changes)
    dp3 = (2 / (n_pipe(i)*(A1+A2)) ) * (1 / (density_cool*A1*n_pipe(i)) - 1 / (density_cool*A2*n_pipe(i))); % Pa - Acceleration Losses

    dp = dp1+dp2+dp3; % Pa - Total Losses

    p_cool(i+1) = p_cool(i) - dp;

    v_cool(i+1) = sqrt(v_cool(i)^2 - 2*dp/density_cool); % m/s - coolant velocity
    Re_cool(i+1) = v_cool(i+1)*d_hydraulic(i+1)/kin_visc_cool; % Coolant Reynold's Number
    f_cool(i+1) = (0.79 .* log(Re_cool(i+1)) - 1.64).^(-2); % friction factor, smooth pipe approximation REVIEW IF THIS IS ACCURATE

    if imag(v_cool(i+1)) ~= 0 % If v_cool drops to zero
        error("Coolant Flow stalls out at %g m from injector (%g m before exit)", x(i), x_exit-x(i));
    end
end

p_cool = p_cool + pc - p_cool(end); % Pa - ending pressure EXPAND TO ACCOUNT FOR PLUMBING AND INJECTOR STIFFNESS


