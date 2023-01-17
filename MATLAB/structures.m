%% Check for hoop stress failure in combustion chamber
hoop = p_gas*0.5.*(r2+r1) ./ (r2-r1);

[max_hoop, max_hoop_loc] = max(hoop); 

% Aluminum 6061 Yield Strength Temperature Dependence
al6061_temp = ones(1,10) + [20,100,150,200,250,300,350,400,450,500]; % K - Approximated from https://firesciencereviews.springeropen.com/articles/10.1186/s40038-015-0007-5#:~:text=6061%2DT651%20yield%20strength%20exhibits,%C2%B0C%20(~240%20MPa).
al6061_yield = [320,310,300,280,180,100,60,35,15,10] * 1E6; % MPa
yield_cc = interp1(al6061_temp, al6061_yield, T_wall_hot);

if (yield_cc/max_hoop < FS_design)
    error("Peak hoop stress FSy of %g at %gm from injector", max_hoop/yield_al6061, x(max_hoop_loc));
end