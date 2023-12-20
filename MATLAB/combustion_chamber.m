%% CEA Run Rocket

 cea = CEA('problem','rocket', ...
     'equilibrium', ...  % equilibrium flow
     'fac','acat',pi*r1_chamber^2 / A_throat, ... % finite area combustor, contraction ratio
     'o/f',OF, ...          % Define OF ratio
     'case','CEAM-rocket1', ...
     'p(psi)',p_cc/6894.76,'pi/p',p_cc/p_amb, ... % Define nozzle with pressure ratio
     'reactants', ...    
     'fuel','C2H5OH(L)','wt%',proof,'t(k)',T_amb,'H2O(L)','wt%',1-proof,'t(k)',T_amb, ... % 190-Proof Ethanol for Fuel, at ambient
     'oxid','O2(L)','wt%',100, ... % Liquid Oxygen for Oxidizer, at boiling point
     'output','short','massf','transport','mks', ...
     'end');

 gamma_avg = mean([cea.output.eql.gamma(3), cea.output.eql.gamma(4)]); % Ratio of Specific Heats, estimate
 Tt_cc = cea.output.eql.temperature(1); % K