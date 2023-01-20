%% Check for hoop stress failure in combustion chamber
hoop = p_gas*0.5.*(r2+r1) ./ (r2-r1);

[max_hoop, max_hoop_loc] = max(hoop); 

switch wall
    case "aluminum" % 6061
        % Aluminum 6061 Yield Strength Temperature Dependence
        wall_temp = 273*ones(1,10) + [20,100,150,200,250,300,350,400,450,500]; % K - Approximated from https://firesciencereviews.springeropen.com/articles/10.1186/s40038-015-0007-5#:~:text=6061%2DT651%20yield%20strength%20exhibits,%C2%B0C%20(~240%20MPa).
        wall_yield = [320,310,300,280,180,100,60,35,15,10] * 1E6; % Pa
        yield_cc = interp1(wall_temp, wall_yield, T_wall_hot);
    case "steel"
        wall_temp = [100,2000]; % K - fake
        wall_yield = [415,415] * 1E6; % Pa
        yield_cc = interp1(wall_temp, wall_yield, T_wall_hot);
    case "copper"

end

FSy = yield_cc/max_hoop;
if (FSy < FS_design)
    fprintf("Peak hoop stress FSy of %g at %gm from injector\n", max_hoop/yield_al6061, x(max_hoop_loc));
end