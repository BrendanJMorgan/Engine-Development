%%

cea_gg = CEA('problem','rocket', ...
     'equilibrium', ...  % equilibrium flow, infinite chamber, 
     'o/f',OF_gg, ...          % Define OF ratio
     'case','CEAM-rocket1', ...
     'p(psi)',p_gg/6894.76,'pi/p',p_gg/p_amb, ... % Define nozzle with pressure ratio
     'reactants', ... 
     'fuel','C2H5OH(L)','wt%',proof,'t(k)',T_amb,'H2O(L)','wt%',1-proof,'t(k)',T_amb, ... % 190-Proof Ethanol for Fuel, at ambient
     'oxid','O2(L)','wt%',100, ... % Liquid Oxygen for Oxidizer, at boiling point
     'output','short','massf','transport','mks', ...
     'end');

Tt_gg = cea_gg.output.eql.temperature(1); % K - Stagnation / Chamber Temperature in gas generator
cp_gg = 1000*mean(cea_gg.output.eql.cp); % J/kg-K - average Specific Heat MODIFY TURBINE TO NOT BE CALORICALLY PERFECT?
gamma_gg = mean(cea_gg.output.eql.gamma); % Ratio of Specific Heats
R_gg = 8.3145 / (0.001*cea.output.eql.mw(3)); % J/kg-K - Specific Gas Constant, fully combusted GG PRODUCTS SHOULD FULLY COMBUST BEFORE NOZZLE PLATE

c_star_ideal_gg = cea.output.eql.cstar(3); % m/s - ideal characteristic velocity
c_star_gg = c_star_ideal_gg*c_star_eff_gg; % m/s - actual characteristic velocity
