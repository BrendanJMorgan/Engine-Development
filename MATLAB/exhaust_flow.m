%% Exhaust Flow Properties

M = zeros(1,length(x)); % Mach Number

% Chamber
M(1:round(l_chamber/dx)) = linspace(M_inj, M_comb, round(l_chamber/dx));
dm = 0.0001;

% Converging
m = 0:dm:1;
for i = round(l_chamber/dx):1:round(x2_throat/dx)
    % Isentropic Area Relation
    [~,index] = min(abs( (1./m)*(0.5*(gamma+1)).^(-0.5*(gamma+1)./(gamma-1)).*(1+0.5*(gamma-1) .* m.^2).^((gamma+1)./(2*(gamma-1))) - (r1(i)*2./d_throat).^2 ) )   ;
    M(i) = (-dm + dm * index);
end
factor1 = (1-M_comb)/(1-M(round(l_chamber/dx)));
M(round(l_chamber/dx):round(x2_throat/dx)) = 1 - ((1-M(round(l_chamber/dx):round(x2_throat/dx))) * factor1); % Stretch to fit CEA data

% Diverging
m = 1:dm:3;
for i = round(x2_throat/dx):1:length(x)
    % Isentropic Area Relation
    [~,index] = min(abs( (1./m)*(0.5*(gamma+1)).^(-0.5*(gamma+1)./(gamma-1)).*(1+0.5*(gamma-1) .* m.^2).^((gamma+1)./(2*(gamma-1))) - (r1(i)*2./d_throat).^2 ) )   ;
    M(i) = (1 -dm + dm * index) * (M_exit);
end
factor2 = M_exit/M(length(x));
M(round(x2_throat/dx):length(x)) = M(round(x2_throat/dx):length(x)) * factor2; % Stretch to fit CEA data
