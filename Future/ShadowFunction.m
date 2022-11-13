function [Psi] = ShadowFunction(r_I, s_I, d, height_atm)
%SHADOWFUNCTION Calculate Psi in range [0,1] to determine the amount of
%sunlight on the satellite, by using a perspective projection
 %   r_I = [m] satellite position in ECI frame
 %   s_I = [-] unit vector in ECI frame pointing towards sun
 %   d   = [m] distance earth centre to sun
 %   heigth_atm = [km] atmospheric height, if 0 than no atmospheric effects
if nargin < 4
    height_atm = 0;
end

%% [m] to [km];
d = d/(1E3); 
r_I = r_I/(1E3);

r_eq = 6378.1370; % [km], equatorial radius earth
r_po = 6356.7523; % [km], polar radius earth

a_atm = r_eq + height_atm;
b_atm = r_po*(a_atm/r_eq);

%% Calculate psi
if height_atm == 0 % No atmospheric effects
    [state2, Area_earth, ~, ~, r_sol] = ProjectionFunction(r_eq, r_po, r_I, s_I, d);
    Area_sol = pi*r_sol^2;
    
    if state2 == 0
        Psi = 1;
    elseif state2 == -1
        Psi = 0;
    else
        Psi = Area_earth/Area_sol;
    end
else
    %% Consider atmosphere
    [state1, Area_atm, dis1, dis0, r_sol] = ProjectionFunction(a_atm, b_atm, r_I, s_I, d);
    [state2, Area_earth, dis2, ~,~] = ProjectionFunction(r_eq, r_po, r_I, s_I, d);
    %% (Linear) absorption
    mu(1) = 0;
    mu(2) = 1;
    
%     a_linear = mu(1);
%     b_linear = mu(2);
%     
%     a_exp = mu(1);
%     b_exp = log(m(2)/mu(1));
%     
%     a_log = mu(1);
%     b_log = (mu(2)-mu(1))/log(2);
    
    Area_sol = pi*r_sol^2;

    thickness = dis1 - dis2;
    
    if state1 == 0 % Full phase
        Psi = 1;
    end
    %% Partly in atmosphere
    if (state1 == 1 || state1 == 2) && state2 == 0
        x(1) = 1;
        x(2) = thickness - (dis1-dis0)/thickness;
        u(1) = (mu(2)-mu(1))*x(1) + mu(1);
        u(2) = (mu(2)-mu(1))*x(2) + mu(1);
        
        coeff = (u(1)+u(2))/2;
        Psi = (Area_atm + (Area_sol-Area_atm)*coeff)/Area_sol;
    end
    %% Totally in atmosphere
    if  state1 == -1 && state2 == 0
        if not(dis2 == 0)
            x(1) = (dis0-dis2)/thickness;
            x(2) = (dis0-dis2 + 2*r_sol)/thickness;
            u(1) = (mu(2)-mu(1))*x(1) + mu(1);
            u(2) = (mu(2)-mu(1))*x(2) + mu(1);
        elseif abs(dis2)<1E-9 % No projection for solid earth, but sat is in atmos projection
            u(1) = 0.5;
            u(2) = 0.5;
        end
        Psi = (u(1) + u(2))/2;
    end
    %% Partly in the earth and atmosphere
    if state1 == -1 && (state2 == 1 || state2 == 2)
        x(1) = (2*r_sol - (dis2-dis0))/thickness;
        x(2) = 0;
        
        u(1) = (mu(2)-mu(1))*x(1) + mu(1);
        u(2) = (mu(2)-mu(1))*x(2) + mu(1);
        
        coeff = (u(1)+u(2))/2;
        Psi = (Area_earth*coeff)/Area_sol;
    end
    %% Totally in umbra
    if state2 == -1
        Psi = 0;
    end
    %% Sun bigger than atmospheric thickness
    if (state1 == 1 || state1 == 2) && (state2 == 1 || state2 == 2)
        area1 = Area_atm;
        area2 = Area_earth - Area_atm;
        x(1) = 1;
        x(2) = 0;
        
        u(1) = (mu(2)-mu(1))*x(1) + mu(1);
        u(2) = (mu(2)-mu(1))*x(2) + mu(1);
        
        Psi = (area1 + area2*(u(1)+u(2))/2)/Area_sol;
    end
    
    % Solves problem near transition from full phase to penumbra
    if Psi > 1
        Psi = 1;
    end
end

