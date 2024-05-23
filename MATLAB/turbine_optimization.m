clear stress
clear base
clear pitchline
clear tip


d_throat_nozzle = 1/4*0.0254;   % m - diameter of each nozzle leading off the manifold - mdot_gg is a direct function of this and nozzle_number
nozzle_number = 4;              % unitless
blade_width_rotor = 1/2*0.0254; % m - real turbines seems to be about 0.5 inches or so


gears = [1, 65/60, 70/55, 75/50, 80/45, 85/40, 90/35, 95/30, 100/25];
d = 0.1:0.001:0.4;

figure(1); clf;
figure(2); clf;

for i = 1:1:length(gears)
    for j = 1:1:length(d)
        gear_ratio = gears(i);           % unitless - needs to add up to 125
        d_throat_nozzle = d(j)*0.0254;   % m - diameter of each nozzle leading off the manifold - mdot_gg is a direct function of this and nozzle_number
        turbine
    
        stress(i,j) = stress_radial;
        base(i,j) = r_rotor_base;
        pitchline(i,j) = r_rotor_pitchline;
        tip(i,j) = r_rotor_tip;
    end

end

figure(1); clf; hold on;
for i = 1:1:length(gears)
    plot(d, stress(i,:)/6.895e+6);
end
yscale('log')
legend('1', '65/60', '70/55', '75/50', '80/45', '85/40', '90/35', '95/30', '100/25')
hold off

figure(2); clf; hold on;
for i = 1:1:length(gears)
    plot(d,tip(i,:)/0.0254)
end
legend('1', '65/60', '70/55', '75/50', '80/45', '85/40', '90/35', '95/30', '100/25')
hold off

% figure(3)
% surf(d, gears, stress)
