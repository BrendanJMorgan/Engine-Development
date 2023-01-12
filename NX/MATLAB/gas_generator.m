%%

OF_gg = 0.05;
p_gg = 100*6894.76; % Pa - chamber pressure inside gas generator
gg_fraction = 0.1;

mdot_gg = mdot_cc*gg_fraction/(1-gg_fraction);     % kg/s - Propellant Mass Flow Rate
mdot_fuel_gg = mdot_gg*(1/(1+OF_gg)); % kg/s - Fuel Mass Flow Rate
mdot_ox_gg = mdot_gg*(OF_gg/(1+OF_gg)); % kg/s - Oxidizer Mass Flow Rate

mdot

 cea_gg = CEA('problem','rocket', ...
     'equilibrium', ...  % equilibrium flow, infinite chamber, 
     'o/f',OF_gg, ...          % Define OF ratio
     'case','CEAM-rocket1', ...
     'p(psi)',p_gg/6894.76,'pi/p',p_gg/p_amb, ... % Define nozzle with pressure ratio
     'reactants', ... 
     'fuel','C2H5OH(L)','wt%',proof,'t(k)',Tamb,'H2O(L)','wt%',1-proof,'t(k)',Tamb, ... % 190-Proof Ethanol for Fuel, at ambient
     'oxid','O2(L)','wt%',100, ... % Liquid Oxygen for Oxidizer, at boiling point
     'output','short','massf','transport','mks', ...
     'end');

T_gg = cea_gg.output.eql.temperature(1); % Stagnation / Chamber Temperature in gas generator
gamma_gg = mean(cea_gg.output.eql.gamma); % Ratio of Specific Heats % MAKE X DEPENDENT

c_star_ideal_gg = cea.output.eql.cstar(3); % Characteristic Velocity
c_star_eff_gg = 0.75; % Char Vel Efficiency, experimental
c_star_gg = c_star_ideal_gg*c_star_eff_gg;

M_inj_gg  = cea_gg.output.eql.mach(1); % Mach at injector
M_comb_gg = 0.002; % Mach at start of convergence (combustor?)
M_throat_gg = cea_gg.output.eql.mach(2); % Mach at throat, should be 1 unless something is very wrong
if M_throat_gg < 0.99 || M_throat_gg > 1.01
    error("Sonic point is not at throat");
end
M_exit_gg = cea_gg.output.eql.mach(3); % Mach at exit, should be supersonic

c_exit_gg = cea_gg.output.eql.sonvel(3); % Sonic velocity at exit