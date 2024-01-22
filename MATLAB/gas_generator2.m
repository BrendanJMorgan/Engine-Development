clear cea_gg

cea_gg = CEA('problem','rocket', ...
     'equilibrium', ...  % equilibrium flow, infinite chamber, 
     'o/f',OF_gg, ...          % Define OF ratio
     'case','cea_ggM-rocket1', ...
     'p(psi)',p_gg/6894.76,'pi/p',p_gg/p_amb, ... % Define nozzle with pressure ratio
     'reactants', ... 
     'fuel','C2H5OH(L)','wt%',proof,'t(k)',T_amb,'H2O(L)','wt%',1-proof,'t(k)',T_amb, ... % 190-Proof Ethanol for Fuel, at ambient
     'oxid','O2(L)','wt%',100, ... % Liquid Oxygen for Oxidizer, at boiling point
     'output','short','massf','transport','mks', ...
     'end');

Ae_At_gg = cea_gg.output.eql.aeat(end); % Expansion / Area Ratio
Tt_gg = cea_gg.output.eql.temperature(1); % K - stagnation, gg chamber temperature
gamma_gg = cea_gg.output.eql.gamma(end); % unitless - ratio of specific heats of final exhaust

c_star_ideal_gg = cea_gg.output.eql.cstar(1); % m/s - Characteristic Velocity
c_star_gg = c_star_ideal_gg*c_star_eff_gg; % m/s

p_exit_gg = 1E5*cea_gg.output.eql.pressure(end); % Pa - exit pressure
c_tau_ideal_gg = sqrt( (2*gamma_gg^2/(gamma_gg-1) * (2/(gamma_gg+1))^((gamma_gg+1)/(gamma_gg-1)) * (1-(p_amb/p_gg)^((gamma_gg-1)/gamma_gg) ) ) ) + Ae_At*(p_exit-p_amb)/p_gg; % Ideal Thrust Coefficient
lambda_cone_gg = 0.5*(1+cos(diverge_angle_gg)); % Divergence correction factor for conical nozzle
c_tau_gg = c_tau_ideal_gg*lambda_cone_gg*c_tau_eff_gg; % Thrust Coefficient

v_exhaust_ideal_gg = cea_gg.output.eql.sonvel(end)*cea_gg.output.eql.mach(end); % m/s - ideal exhaust velocity
v_exhaust_gg = c_star_gg*c_tau_gg; % m/s - actual exhaust velocity