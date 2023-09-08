%% Objective: find all relevant properties of the combustion gas flow through the combustion chamber and nozzle

Tc = cea.output.eql.temperature(1); % Stagnation / Chamber Temperature

n1 = floor(x_combustor/dx);
n2 = floor(x2_throat/dx) - n1;  
n3 = length(x) - n2 - n1;

%% Mach

M1 = cea.output.eql.mach(1); % Mach at injector
M2 = cea.output.eql.mach(2); % Mach at start of convergence (combustor?) UPDATE TO FINITE AREA COMBUSTOR
M3 = cea.output.eql.mach(3); % Mach at throat, should be 1 unless something is very wrong
M4 = cea.output.eql.mach(4); % Mach at exit, should be supersonic

% THE STRETCHING IS SUS
M_gas = zeros(1,length(x)); % Mach Number

% Chamber
M_gas(1:round(l_chamber/dx)) = linspace(M1, M2, round(l_chamber/dx));
dm = 0.0001;

% Converging
m = 0:dm:1;
for i = round(l_chamber/dx):1:round(x2_throat/dx)
    % Isentropic Area Relation
    [~,index] = min(abs( (1./m)*(0.5*(gamma_gas+1)).^(-0.5*(gamma_gas+1)./(gamma_gas-1)).*(1+0.5*(gamma_gas-1) .* m.^2).^((gamma_gas+1)./(2*(gamma_gas-1))) - (r1(i)*2./d_throat).^2 ) )   ;
    M_gas(i) = (-dm + dm * index);
end
factor1 = (1-M2)/(1-M_gas(round(l_chamber/dx)));
M_gas(round(l_chamber/dx):round(x2_throat/dx)) = 1 - ((1-M_gas(round(l_chamber/dx):round(x2_throat/dx))) * factor1); % Stretch to fit CEA data

% Diverging
m = 1:dm:3;
for i = round(x2_throat/dx):1:length(x)
    % Isentropic Area Relation
    [~,index] = min(abs( (1./m)*(0.5*(gamma_gas+1)).^(-0.5*(gamma_gas+1)./(gamma_gas-1)).*(1+0.5*(gamma_gas-1) .* m.^2).^((gamma_gas+1)./(2*(gamma_gas-1))) - (r1(i)*2./d_throat).^2 ) )   ;
    M_gas(i) = (1 -dm + dm * index) * (M4);
end
factor2 = M4/M_gas(length(x));
M_gas(round(x2_throat/dx):length(x)) = M_gas(round(x2_throat/dx):length(x)) * factor2; % Stretch to fit CEA data

%% Enthalpy
H1 = cea.output.eql.enthalpy(1);
H2 = cea.output.eql.enthalpy(2);
H3 = cea.output.eql.enthalpy(3);
H4 = cea.output.eql.enthalpy(4);
segment1 = linspace(H1, H2, n1);
segment2 = linspace(H2, H3, n2);
segment3 = linspace(H3, H4, n3);
H_gas = 1000 * [segment1, segment2, segment3]; % J/kg

%% Density
dens1 = cea.output.eql.density(1);
dens2 = cea.output.eql.density(2);
dens3 = cea.output.eql.density(3);
dens4 = cea.output.eql.density(4);
segment1 = linspace(dens1, dens2, n1);
segment2 = linspace(dens2, dens3, n2);
segment3 = linspace(dens3, dens4, n3);
dens_gas = [segment1, segment2, segment3]; % kg/m^3

%% Pressure
p1 = cea.output.eql.pressure(1);
p2 = cea.output.eql.pressure(2);
p3 = cea.output.eql.pressure(3);
p4 = cea.output.eql.pressure(4);
segment1 = linspace(p1, p2, n1);
segment2 = linspace(p2, p3, n2);
segment3 = linspace(p3, p4, n3);
p_gas = 1E5 * [segment1, segment2, segment3]; % Pa

%% Specific Heat, Constant Pressure
cp1 = cea.output.eql.cp(1);
cp2 = cea.output.eql.cp(2);
cp3 = cea.output.eql.cp(3);
cp4 = cea.output.eql.cp(4);
segment1 = linspace(cp1, cp2, n1);
segment2 = linspace(cp2, cp3, n2);
segment3 = linspace(cp3, cp4, n3);
cp_gas = [segment1, segment2, segment3]; % J/(kg*K)

%% Ratio of Specific Heats
gamma1 = cea.output.eql.gamma(1);
gamma2 = cea.output.eql.gamma(2);
gamma3 = cea.output.eql.gamma(3);
gamma4 = cea.output.eql.gamma(4);
segment1 = linspace(gamma1, gamma2, n1);
segment2 = linspace(gamma2, gamma3, n2);
segment3 = linspace(gamma3, gamma4, n3);
gamma_gas = [segment1, segment2, segment3]; % unitless

%% Molar Mass
mol1 = cea.output.eql.mw(1);
mol2 = cea.output.eql.mw(2);
mol3 = cea.output.eql.mw(3);
mol4 = cea.output.eql.mw(4);
segment1 = linspace(mol1, mol2, n1);
segment2 = linspace(mol2, mol3, n2);
segment3 = linspace(mol3, mol4, n3);
mol_gas = [segment1, segment2, segment3]; % g/mol

%% Velocity
cs_gas = sqrt(gamma_gas.*p_gas./dens_gas); % m/s - speed of sound
v_gas = M_gas.*cs_gas; % m/s



