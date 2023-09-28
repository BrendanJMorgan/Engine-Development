%% Blades

% Leading and Trailing Edges
r_inlet = r_eye; % m - this assumption can be changed

% Surface of Revolution
rms_curve =  [shroud_curve(:,1) - 0.5*crosswise_gap.*cos(normal_angles), shroud_curve(:,2) - 0.5*crosswise_gap.*sin(normal_angles)]; % [m,m] - halfway between shroud and impeller

% Inlet Blade Angles
inlet_gap = crosswise_gap(find(rms_curve > r_inlet, 1));
u_inlet = shaft_speed*r_inlet; % m/s - rotational speed at the inlet
v_inlet = vdot_fuel/(2*pi*r_inlet*inlet_gap); % m/s - fluid velocity at the inlet, assumed to be entirely radial. Can be found with CFD instead for more accuracy
angle_inlet = atan(u_inlet/v_inlet); % rad - inlet angle of blades relative to straight radial direction

% Outlet Blade Angles
% slip_factor

%%
