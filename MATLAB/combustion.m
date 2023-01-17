%% CEA Run Rocket

 cea = CEA('problem','rocket', ...
     'equilibrium', ...  % equilibrium flow, infinite chamber, 
     'o/f',OF, ...          % Define OF ratio
     'case','CEAM-rocket1', ...
     'p(psi)',pc/6894.76,'pi/p',pc/p_amb, ... % Define nozzle with pressure ratio
     'reactants', ... 
     'fuel','C2H5OH(L)','wt%',proof,'t(k)',Tamb,'H2O(L)','wt%',1-proof,'t(k)',Tamb, ... % 190-Proof Ethanol for Fuel, at ambient
     'oxid','O2(L)','wt%',100, ... % Liquid Oxygen for Oxidizer, at boiling point
     'output','short','massf','transport','mks', ...
     'end');

Tc = cea.output.eql.temperature(1); % Stagnation / Chamber Temperature
gamma = mean(cea.output.eql.gamma); % Ratio of Specific Heats % MAKE X DEPENDENT

M_inj = cea.output.eql.mach(1); % Mach at injector
M_comb = 0.002; % Mach at start of convergence (combustor?) UPDATE TO FINITE AREA COMBUSTOR
M_throat = cea.output.eql.mach(2); % Mach at throat, should be 1 unless something is very wrong
if M_throat < 0.99 || M_throat > 1.01
    error("Sonic point is not at throat");
end
M_exit = cea.output.eql.mach(3); % Mach at exit, should be supersonic