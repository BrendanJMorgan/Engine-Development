%% Objective: find all relevant properties of the combustion gas flow through the combustion chamber and nozzle

Tc = cea.output.eql.temperature(1); % Stagnation / Chamber Temperature

n1 = floor(x_combustor/dx);
n2 = floor(x2_throat/dx) - n1;  
n3 = length(x) - n2 - n1;

%% Enthalpy
H_gas = cea.output.eql.enthalpy;
segment1 = linspace(H_gas(1), H_gas(2), n1);
segment2 = linspace(H_gas(2), H_gas(3), n2);
segment3 = linspace(H_gas(3), H_gas(4), n3);
H_gas = 1000 * [segment1, segment2, segment3]; % J/kg

%% Density
dens_gas = cea.output.eql.density;
segment1 = linspace(dens_gas(1), dens_gas(2), n1);
segment2 = linspace(dens_gas(2), dens_gas(3), n2);
segment3 = linspace(dens_gas(3), dens_gas(4), n3);
dens_gas = [segment1, segment2, segment3]; % kg/m^3

%% Pressure
p_gas = cea.output.eql.pressure;
segment1 = linspace(p_gas(1), p_gas(2), n1);
segment2 = linspace(p_gas(2), p_gas(3), n2);
segment3 = linspace(p_gas(3), p_gas(4), n3);
p_gas = 1E5 * [segment1, segment2, segment3]; % Pa

%% Specific Heat, Constant Pressure
cp_gas = cea.output.eql.cp;
segment1 = linspace(cp_gas(1), cp_gas(2), n1);
segment2 = linspace(cp_gas(2), cp_gas(3), n2);
segment3 = linspace(cp_gas(3), cp_gas(4), n3);
cp_gas = [segment1, segment2, segment3]; % J/(kg*K)

%% Ratio of Specific Heats
gamma_gas = cea.output.eql.gamma;
segment1 = linspace(gamma_gas(1), gamma_gas(2), n1);
segment2 = linspace(gamma_gas(2), gamma_gas(3), n2);
segment3 = linspace(gamma_gas(3), gamma_gas(4), n3);
gamma_gas = [segment1, segment2, segment3]; % unitless

%% Molar Mass
mol_gas = cea.output.eql.mw;
segment1 = linspace(mol_gas(1), mol_gas(2), n1);
segment2 = linspace(mol_gas(2), mol_gas(3), n2);
segment3 = linspace(mol_gas(3), mol_gas(4), n3);
mol_gas = [segment1, segment2, segment3]; % g/mol

%% Mach
cs_gas = sqrt(gamma_gas.*p_gas./dens_gas); % m/s - speed of sound
M1 = cea.output.eql.mach(1) + v_injection/cs_gas(1); % Mach at injector
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
    [~,index] = min(abs( (1./m)*(0.5*(gamma_avg+1)).^(-0.5*(gamma_avg+1)./(gamma_avg-1)).*(1+0.5*(gamma_avg-1) .* m.^2).^((gamma_avg+1)./(2*(gamma_avg-1))) - (r1(i)*2./d_throat).^2 ) )   ;
    M_gas(i) = (-dm + dm * index);
end
factor1 = (1-M2)/(1-M_gas(round(l_chamber/dx)));
M_gas(round(l_chamber/dx):round(x2_throat/dx)) = 1 - ((1-M_gas(round(l_chamber/dx):round(x2_throat/dx))) * factor1); % Stretch to fit CEA data

% Diverging
m = 1:dm:3;
for i = round(x2_throat/dx):1:length(x)
    % Isentropic Area Relation
    [~,index] = min(abs( (1./m)*(0.5*(gamma_avg+1)).^(-0.5*(gamma_avg+1)./(gamma_avg-1)).*(1+0.5*(gamma_avg-1) .* m.^2).^((gamma_avg+1)./(2*(gamma_avg-1))) - (r1(i)*2./d_throat).^2 ) )   ;
    M_gas(i) = (1 -dm + dm * index) * (M4);
end
factor2 = M4/M_gas(length(x));
M_gas(round(x2_throat/dx):length(x)) = M_gas(round(x2_throat/dx):length(x)) * factor2; % Stretch to fit CEA data

%% Velocity

v_gas = M_gas.*cs_gas; % m/s



