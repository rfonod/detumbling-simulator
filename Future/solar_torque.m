function [T] = solar_torque(C_I2B, r_I, s_I, d)
%SOLAR_TORQUE function to compute the solar torque acting on the satellite.
%  C_I2B = [-], rotation matrix from ECI to BODY frame
%  r_I   = [m], position vector of satellite in ECI frame
%  s_I   = [-], unit solar pointing vector in ECI frame
%  d     = [m], distance Earth to sun

global k prevPsi Psi % Needs to be global in order to be able to use them in next iterations
%% Shadow function - determine if the satellite receives any sunlight
% Calculate Psi every 20 iterations if it was previously in full
% illumination or shadow
prevPsi = 0;

if (prevPsi == 1 || prevPsi == 0) && (mod(k,20) == 0 || k == 1 || k == 2)
    prevPsi = Psi;
    Psi = ShadowFunction(r_I, s_I, d, 50);
elseif mod(k,4) == 0 % Must be in penumbra because prevPsi =/= 0 or 1
    prevPsi = Psi;
    Psi = ShadowFunction(r_I, s_I, d, 50);
end

if Psi == 0 % In Earth's shadow, no use to do further computation
    T = [0; 0; 0];
    return
end

%% Convert solar unit vector to body frame
s_B   = C_I2B*s_I; % Points from spacecraft to sun

%% Material properties  -> Needs to be refined
SolarPanelArea = 0.8;       % Approximately 80% of the PQ is covered in solar panels
AlBeta  = 0.825;            % Specular Reflection coefficient for Aluminum
AlGamma = 0.64;             % Diffuse Reflection coefficient for Aluminum
SPBeta  = 0.05;
SPGamma = 0.05;

beta    = SolarPanelArea*SPBeta  + (1-SolarPanelArea)*AlBeta;
gamma   = SolarPanelArea*SPGamma + (1-SolarPanelArea)*AlGamma;

%% Check which side is illuminated
if s_B(1) > 0 % Then positive x side is illuminated
    nx = [1 0 0]';
else
    nx = [-1 0 0]';
end

if s_B(2) > 0 % Then positive y side (plate) is illuminated
    ny = [0 1 0]';
else
    ny = [0 -1 0]';
end

if s_B(3) > 0 % Then positive z side is illuminated
    nz = [0 0 1]';
else
    nz = [0 0 -1]';
end

%% Calculate Torques with SRT Function
Tx = SRT(nz, s_B, beta, gamma);
Ty = SRT(ny, s_B, beta, gamma);
Tz = SRT(nx, s_B, beta, gamma);

% Add them and scale with shadow function
T = Psi*(Tx+Ty+Tz);
end

function [T] = SRT(n, s_B, beta, gamma )
%SRT Calculates the Solar Radiation Torque on one surface
%   n = surface normal

%% Make sure vectors column vectors for cross product
n   = n(:);
s_B = s_B(:);
%% Variables
global dim CoM parameters
Theta = acos(s_B'*n/(norm(s_B)*norm(n))); % Angle between normal and sunrays
phiC  = parameters.phiC;

%% Extract dimensions
d_x   = dim(1);
d_y   = dim(2);
d_z   = dim(3);
d_x_p = dim(4);
d_y_p = dim(5);
d_z_p = dim(6);

%% Calculate torques T = r x F
if abs(n(1)) == 1 % x-surface
    r_B    = [n(1)*1/2*d_x; - CoM(2); 0];                     % Body surface centroid vector
    r_p    = [n(1)*1/2*d_x_p; - CoM(2)+1/2*d_y+1/2*d_y_p; 0]; % Plate surface centroid vector
    area_B = d_y*d_z;
    area_p = d_y_p*d_z_p;
    F_B    = -phiC*area_B*cos(Theta)*(2*(gamma/3 + beta*cos(Theta))*n + (1-beta)*s_B);
    F_p    = -phiC*area_p*cos(Theta)*(2*(gamma/3 + beta*cos(Theta))*n + (1-beta)*s_B);
    T = cross(r_B, F_B) + cross(r_p, F_p);
    return
    
elseif n(2) == 1 % y-surface body
    dirx    = sign(s_B(1));
    dirz    = sign(s_B(3));             % Determine which side of the plate is illuminated
    
    r_B     = [0; 1/2*d_y + CoM(2); 0];                    % Body surface centroid vector
    r_p1    = [dirx*1/2*(d_x+d_x_p); - CoM(2)+1/2*d_y; 0]; % Plate surface centroid vector
    r_p2    = [0; - CoM(2)+1/2*d_y; dirz*1/2*(d_z+d_z_p)]; % Plate surface centroid vector
    
    area_B  = d_y*d_z;
    area_p1 = 1/2*(d_x_p-d_x)*d_z_p;
    area_p2 = 1/2*(d_z_p-d_z)*d_x;
    
    F_B     = -phiC*area_B *cos(Theta)*(2*(gamma/3 + beta*cos(Theta))*n + (1-beta)*s_B);
    F_p1    = -phiC*area_p1*cos(Theta)*(2*(gamma/3 + beta*cos(Theta))*n + (1-beta)*s_B);
    F_p2    = -phiC*area_p2*cos(Theta)*(2*(gamma/3 + beta*cos(Theta))*n + (1-beta)*s_B);
    
    T = cross(r_B, F_B) + cross(r_p1, F_p1) + cross(r_p2, F_p2);
    return
    
elseif n(2) == -1 % y-surface plate
    r_p    = [0; - CoM(2) + 1/2*d_y + d_y_p; 0]; % Plate surface centroid vector
    area_p = d_x_p*d_z_p;
    F_p    = -phiC*area_p*cos(Theta)*(2*(gamma/3 + beta*cos(Theta))*n + (1-beta)*s_B);
    
    T = cross(r_p, F_p);
    return
    
elseif abs(n(3)) == 1 % z-surface
    r_B    = [0; - CoM(2); n(3)*1/2*d_z];                     % Body surface centroid vector
    r_p    = [0; - CoM(2)+1/2*d_y+1/2*d_y_p; n(3)*1/2*d_z_p]; % Plate surface centroid vector
    area_B = d_x*d_y;
    area_p = d_x_p*d_y_p;
    F_B    = -phiC*area_B*cos(Theta)*(2*(gamma/3 + beta*cos(Theta))*n + (1-beta)*s_B);
    F_p    = -phiC*area_p*cos(Theta)*(2*(gamma/3 + beta*cos(Theta))*n + (1-beta)*s_B);
    
    T = cross(r_B, F_B) + cross(r_p, F_p);
    return
else
    error('No valid normal vector given. Provide it in either x, y or z direction')
end
end




