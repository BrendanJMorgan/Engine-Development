%% Coolant Bulk Properties
% Ethanol is the coolant
density_cool = PropsSI('D','T',T_amb,'P',p_amb,['Ethanol[',num2str(proof),']&Water[',num2str(1-proof),']']); % kg/m3 - ethanol at STP
kin_visc_cool = PropsSI('V','T',T_amb,'P',p_amb,['Ethanol[',num2str(proof),']&Water[',num2str(1-proof),']'])/density_cool; % m2/s - kinematic viscosity
mdot_cool = mdot_fuel_cc; % kg/s

%% Pressure Drop

p_cool = ones(1,length(x)); % Pa - coolant pressure

if flow_direction == 1
    flow = 1:1:length(x)-1;
    v_cool = mdot_fuel_cc ./ (density_cool*n_pipe.*w_pipe*h_pipe); % m/s - fluid bulk speed
    flow_exit_index = length(x);
elseif flow_direction == -1
    flow = length(x):-1:2;
    v_cool = mdot_fuel_cc ./ (density_cool*n_pipe.*w_pipe*h_pipe); % m/s - fluid bulk speed
        flow_exit_index = 1;
end

d_hydraulic = 4*w_pipe.*h_pipe./(2*w_pipe+2*h_pipe); % m - hydraulic diameter
Re_cool = ones(1,length(x))*v_cool(1).*d_hydraulic(1)./kin_visc_cool; % Reynold's number
f_cool = ones(1,length(x))*(0.79 .* log(Re_cool(1)) - 1.64).^(-2); % friction factor, smooth pipe approximation REVIEW IF THIS IS ACCURATE

for i = flow
    
    dp1 = 0.5*f_cool(i)*density_cool*v_cool(i)^2*dx./d_hydraulic(i); % Pa - Major Losses (channel wall roughness)
    
    if n_pipe(i+flow_direction) < n_pipe(i) % For when channels merge
        d_hydraulic_mod = 4*0.5*w_pipe(i+flow_direction).*h_pipe./(2*0.5*w_pipe(i+flow_direction)+2*h_pipe); % m - half-channel hydraulic diameter
        d_ratio = d_hydraulic_mod/d_hydraulic(i);
        A1 = w_pipe(i)*h_pipe;
        A2 = 0.5*w_pipe(i+flow_direction)*h_pipe;
    elseif n_pipe(i+flow_direction) > n_pipe(i) % For when channels divide
        d_hydraulic_mod = 4*0.5*w_pipe(i).*h_pipe./(2*w_pipe(i)+h_pipe); % m - half-channel hydraulic diameter
        d_ratio = d_hydraulic(i+flow_direction)/d_hydraulic_mod;
        A1 = 0.5*w_pipe(i)*h_pipe;
        A2 = w_pipe(i+flow_direction)*h_pipe;
    else % Everything else (vast majority of the time)
        d_ratio = d_hydraulic(i+flow_direction)/d_hydraulic(i);
        A1 = w_pipe(i)*h_pipe;
        A2 = w_pipe(i+flow_direction)*h_pipe;
    end

    if d_ratio < 1
        K = 0.5 - 0.167*d_ratio - 0.125*d_ratio^2 - 0.208*d_ratio^3; % Contraction
    else
        K = (d_ratio^-2 - 1)^2; % Expansion
    end

    dp2 = 0.5*K*density_cool*v_cool(i)^2; % Pa - Minor Losses (diameter changes)
    dp3 = (2 / (n_pipe(i)*(A1+A2)) ) * (1 / (density_cool*A1*n_pipe(i)) - 1 / (density_cool*A2*n_pipe(i))) * mdot_cool^2; % Pa - Acceleration Losses

    dp = dp1+dp2+dp3; % Pa - Total Losses

    p_cool(i+flow_direction) = p_cool(i) - dp;
    Re_cool(i+flow_direction) = v_cool(i+flow_direction)*d_hydraulic(i+flow_direction)/kin_visc_cool; % Coolant Reynold's Number
    f_cool(i+flow_direction) = (0.79 .* log(Re_cool(i+flow_direction)) - 1.64).^(-2); % friction factor, smooth pipe approximation REVIEW IF THIS IS ACCURATE

end


p_cool = p_cool + p_cc*(1+cc_stiffness) - p_cool(end); % Pa - ending pressure - this is a msitake but I'm rolling with it bc manaufacturing already started
% p_cool = p_cool + p_cc*(1+cc_stiffness) - p_cool(flow_exit_index); % Pa - ending pressure

