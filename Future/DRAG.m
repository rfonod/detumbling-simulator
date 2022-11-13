function [T] = DRAG(n,v_B,acc)
%DRAG Summary of this function goes here
%   n = surface normal unit vector
%% Useful constants and values
global dim CoM parameters
v_B    = [v_B(1); v_B(2); v_B(3)]; % Make sure column vector
v      = norm(v_B);
v_unit = v_B/v;
rho    = parameters.rho;
Cdmin  = 2.2;
Cdmax  = 4;

%% Extract dimensions
d_x   = dim(1);
d_y   = dim(2);
d_z   = dim(3);
d_x_p = dim(4);
d_y_p = dim(5);
d_z_p = dim(6);

alpha = 0.2; % Fraction with increasing drag coefficient

%% Calculate Torques
T = [0; 0; 0]; % Preallocated torque

%% X-surface
if abs(n(1)) == 1 
    % Surface element dimensions and vectors
    dx = n(1)*1/2*d_x;
    dy = linspace( (1/2-alpha)*d_y - CoM(2), 1/2*d_y - CoM(2), acc);
    dz = d_z;
    % Surface elements have equal size
    dS = abs( dz*(dy(2)-dy(1)) );
 
    % Centroids of surface elements
    dy = 1/2*(dy(1:end-1) + dy(2:end));

    % Vector for increasing drag coefficient
    Cdvec = linspace(Cdmin, Cdmax, acc); % Linearly increasing Cd for 1/4th of the surface
    
    for i = 1:length(dy)
        dF     = Cdvec(i)*dS*v_unit;
        dT     = [dy(i)* dF(3); -dx*dF(3); -dy(i)*dF(1) + dx*dF(2)];
        T      = T + dT; 
    end
    
    % Uniform surface
    A = (1-alpha)*d_y*dz;
    F = Cdmin*A*v_unit;
    r = [dx; CoM(2)-alpha*d_y; 0];
    T = T + cross(r,F);
    
    % Plate
    A = d_y_p*d_z_p;
    F = Cdmin*A*v_unit;
    r = [n(1)*1/2*d_x_p; -CoM(2) + 1/2*(d_y + d_y_p); 0];
    T = T + cross(r,F);
    
    T = -1/2*rho*v^2*dot(n,v_unit)*T; % Take -1/2rho v^2 ndotv out of summation to reduce computation time
    
    
    
%% Y-surface
elseif n(2) > 0 % Drag on plate -> full symmetry: no torque
    
    F = -1/2*Cdmin*d_x_p*d_z_p*v^2*rho*dot(n,v_unit)*v_unit;
    r = [0; 1/2*d_y - CoM(2) + d_y_p; 0];
    
    T = cross(r,F);
    return
elseif n(2) < 0        % Drag on body and plate
    %% Body
    A = d_x*d_z;
    
    F = Cdmin*A*v_unit;
    r = [0; 1/2*d_y + CoM(2); 0];
    
    T1 = cross(r,F);
    
    %% plate
    xdir = sign(v_B(1));
    zdir = sign(v_B(2));

    Cd = linspace(Cdmax, Cdmin, acc-1);

    % Long side
    dx = xdir*linspace(1/2*d_x, 1/2*d_x_p, acc);
    dx = 1/2*(dx(1:end-1) + dx(2:end));     % Centroids
    dz = d_z_p;
    
    dS = abs( (dx(2)-dx(1))*dz ); 
    
    T2 = [0; 0; 0];
    for i = 1:length(dx)
        dF  = Cd(i)*dS*v_unit;
        dT2 = [(d_y - CoM(2))*dF(3); -dx(i)*dF(3); -(d_y - CoM(2))*dF(1) + dx(i)*dF(2)];
        T2  = T2 + dT2;
    end
    
    % Short side
    dz = zdir*linspace(1/2*d_z, 1/2*d_z_p, acc);
    dz = 1/2*(dz(1:end-1) + dz(2:end));     % Centroids
    dx = d_x;
    
    dS = abs( (dz(2)-dz(1))*dx ); 
    
    T3 = [0; 0; 0];
    
    for i = 1:length(dx)
        dF  = Cd(i)*dS*v_unit;
        dT3 = [-dz(i)*dF(2) + (d_y - CoM(2))*dF(3) 
                dz(i)*dF(1) 
               -(d_y - CoM(2))*dF(1)];
        T3  = T3 + dT3;
    end
        
T = T1 + T2 + T3;

T = -1/2*rho*v^2*dot(n,v_unit)*T; % Take -1/2rho v^2 ndotv out of summation to reduce computation time
    
%% Z-surface
elseif abs(n(3)) == 1 
    dx = d_x;
    dy = linspace( (1/2-alpha)*d_y - CoM(2), 1/2*d_y - CoM(2), acc);
    % Surface elements have equal size
    dS = abs( dx*(dy(2)-dy(1)) );
    
    dz = n(3)*1/2*d_z;
    
    % Centroids of surface elements
    dy = 1/2*(dy(1:end-1) + dy(2:end));

    % Matrix for increasing drag coefficient
    Cdvec = linspace(Cdmin, Cdmax, acc); % Linearly increasing Cd for 1/4th of the surface
    
    T = 0;
    for i = 1:length(dy)
            dF     = Cdvec(i)*dS*v_unit;
            dT     = [-dz*dF(2) + dy(i)*dF(3); dz*dF(1); -dy(i)*dF(1)];
            T      = T + dT;
    end
    
    
    % Uniform surface
    A = (1-alpha)*d_y*dx;
    F = Cdmin*A*v_unit;
    r = [0; CoM(2)-alpha*d_y; dz];
    T = T + cross(r,F);
    
    % Plate
    A = d_x_p*d_y_p;
    F = Cdmin*A*v_unit;
    r = [0; -CoM(2) + 1/2*(d_y + d_y_p); n(3)*1/2*d_z_p];
    T = T + cross(r,F);
    
    T = -1/2*rho*v^2*dot(n,v_unit)*T; % Take -1/2rho v^2 ndotv out of summation to reduce computation time
    
else
    error('No valid normal vector given. Provide in in either x, y or z direction (BODY frame)')

end
end
