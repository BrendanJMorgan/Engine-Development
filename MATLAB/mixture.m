%% Find properties of gas mixtures
%
% mixture(compounds, fractions, Tref, pref, property_type)
%
% Tref - reference temeprature; pref - reference pressure
% property_type:
% C - specific heat
% V - viscosity
% L - thermal conductivity

function [cp, visc, cond, density] = mixture(compounds, fractions, Tref, pref)
    
    if length(compounds) ~= length(fractions)
        error("compounds do not match with compound fractions")
    end

    % proplist = zeros(1,length(compounds));
    % fractionlist = zeros(1,length(fractions));
    j = 1;

    for i = 1:1:length(compounds)
        % Convert CEA product names to CoolProp-valid fluid names
        switch compounds(i)
            case "*CO"
                name = "CarbonMonoxide";
                polar(j) = 1; % CO is a polar molecule
            case "*CO2"
                name = "CarbonDioxide";
            case "*H2"
                name = "Hydrogen";
            case "H2O"
                name = "Water";
                polar(j) = 1; % H2O is a polar molecule
            case "O2"
                name = "Oxygen";
            case "C2H4"
                name = "Ethylene";
            case "C2H6"
                name = "Ethane";
            otherwise
                continue; % Did not find a CoolProp-valid compound; skip to next
        end
        
        fractionlist(j) = fractions(i);

        denslist(j) = py.CoolProp.CoolProp.PropsSI("D", "T", Tref, "P", pref, name); % kg/m3
        cplist(j) = py.CoolProp.CoolProp.PropsSI("C", "T", Tref, "P", pref, name); % J/kg-K

        if name == "CarbonMonoxide" % CoolProp lacks viscosity and thermal conductivity data for carbon monoxide
            visclist(j) = 0.027098*Tref^(0.734156)*1E-5; % power law fitted to data from here: https://www.engineeringtoolbox.com/gases-absolute-dynamic-viscosity-d_1888.html
            condlist(j) = 0.227019*Tref^(0.828249)*0.001; % power law fitted to data from here: https://srd.nist.gov/jpcrdreprint/1.555827.pdf
        elseif name == "Ethylene"
            visclist(j) = (8E-10)*Tref^0.4845; % Terrible power law, fitted to 725 psi
        else
            visclist(j) = py.CoolProp.CoolProp.PropsSI("V", "T", Tref, "P", pref, name); % Pa-s
            condlist(j) = py.CoolProp.CoolProp.PropsSI("L", "T", Tref, "P", pref, name); % W/m-K
        end
        
        weights(j) = py.CoolProp.CoolProp.PropsSI("molemass", "T", Tref, "P", pref, name); % kg/kmol - Molar masses

        if name == "CarbonDioxide"
            boiling(j) = 194.65; % CO2 sublimates at 1 atmosphere; i.e. it has no saturation temperature
        else
            boiling(j) = py.CoolProp.CoolProp.PropsSI( 'T', 'P', 101325, 'Q', 0, name); % K - boiling points at one atmosphere
        end

         j = j+1;
    end

    if 1 - sum(fractionlist) > 0.05
        fprintf("Unaccounted for combustion products make up %g%% of exhaust\n", 100-100*sum(fractionlist));
    end  

    % Density
    density = sum(denslist.*fractionlist)/sum(fractionlist);

    % Specific Heat
    cp = sum(cplist.*fractionlist)/sum(fractionlist);
    
    % Viscosity - Wilke's Method
    phi = (sqrt(2)/4) * ( ones(length(visclist)) + ( visclist'*(1./visclist) ).^0.5 * ( (1./weights)'*weights ).^0.25 ).^2 ./ ( ones(length(visclist)) + weights'*(1./weights) ).^0.5; % Wilke's Coefficients
    phi = phi - diag(diag(phi)); % delete all j = i terms
    visc = dot( ( fractionlist.*visclist )', 1 ./ ( fractionlist' + phi*fractionlist' ) );
    
    % Thermal Conductivity - Lindsay and Bromley's Method
    sutherland = sqrt(2.25 * (boiling' * boiling)); % Sutherland Constants
    sutherland = ((0.733-1)*xor(polar',polar)+ones).*sutherland; % Adjust constants where one polar molecule is involved in the interaction
    A = 0.25 * ( 1 + (  visclist'*(1./visclist) .* ((1./weights)'*weights).^0.75  .* (ones(length(sutherland),1)  +  diag(sutherland)/Tref) * 1./(ones(1,length(sutherland))  +  diag(sutherland)'/Tref)  ).^0.5 ).^2 .* ( 1./(ones(length(sutherland))  +  diag(sutherland)/Tref ) * ( ones(length(sutherland)) + sutherland/Tref ) );
    cond = dot( condlist', 1 ./ ( 1./fractionlist' .* A*fractionlist' ) );

end



