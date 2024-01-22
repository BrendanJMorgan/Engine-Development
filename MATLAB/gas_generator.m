
initial_p_gg = 100*6894.76; % Pa
% options = optimset('Display','iter','PlotFcns',@optimplotfval);
optimal_p_gg = fminsearch(@(p_gg) objectiveFunction(p_gg, norm(v_spouting), OF_gg, p_amb, T_amb, proof), initial_p_gg);

%%

cea_gg = CEA('problem','rocket', ...
     'equilibrium', ...  % equilibrium flow, infinite chamber, 
     'o/f',OF_gg, ...          % Define OF ratio
     'case','CEAM-rocket1', ...
     'p(psi)',optimal_p_gg/6894.76,'pi/p',p_gg/p_amb, ... % Define nozzle with pressure ratio
     'reactants', ... 
     'fuel','C2H5OH(L)','wt%',proof,'t(k)',T_amb,'H2O(L)','wt%',1-proof,'t(k)',T_amb, ... % 190-Proof Ethanol for Fuel, at ambient
     'oxid','O2(L)','wt%',100, ... % Liquid Oxygen for Oxidizer, at boiling point
     'output','short','massf','transport','mks', ...
     'end');

Tt_gg = cea_gg.output.eql.temperature(1);       % K - Stagnation / Chamber Temperature in gas generator
cp_gg = 1000*mean(cea_gg.output.eql.cp);        % J/kg-K - average Specific Heat MODIFY TURBINE TO NOT BE CALORICALLY PERFECT?
gamma_gg = mean(cea_gg.output.eql.gamma);       % Ratio of Specific Heats
R_gg = 8.3145 / (0.001*cea.output.eql.mw(3));   % J/kg-K - Specific Gas Constant, fully combusted GG PRODUCTS SHOULD FULLY COMBUST BEFORE NOZZLE PLATE
dens_gg = cea_gg.output.eql.density(1);         % kg/m3 - Stagnation/Chamber density in gas generator
sonic_gg = cea_gg.output.eql.sonvel(3);         % m/s - speed of sound at ambient pressure
c_star_ideal_gg = cea.output.eql.cstar(3);      % m/s - ideal characteristic velocity
c_star_gg = c_star_ideal_gg*c_star_eff_gg;      % m/s - actual characteristic velocity

optimal_p_gg/6894.76
%%

function error = objectiveFunction(p_gg, v_nozzle_abs, OF_gg, p_amb, T_amb, proof)
    v_nozzle_abs
    % p_gg/6894.76
    v_exit = calculateExitVelocity(p_gg, OF_gg, p_amb, T_amb, proof);
    error = -(v_nozzle_abs - v_exit);
end

function v_exit = calculateExitVelocity(p_gg, OF_gg, p_amb, T_amb, proof)
    p_gg/6894.76
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
    
    v_exit = cea_gg.output.eql.mach(3)*cea_gg.output.eql.sonvel(3) % m/s
end