function [T_aer] = dist_aer(v_I,q)
global rho areas
% v_I - satellite velocity in ECI [m/s]
% q   - attitude quaternion [-]

C_d = 2.1; % satellite's drag coeficient (rectangular box) [-]
CoP = [0.0054 0.0020 -0.0082]'; % worst-case CoP estimate from CoM [m]

% Unit surface area directions and satellite velocity in Body
I3 = eye(3);
A  = [I3 -I3];
C_I2B = q2dcm(q);
normA = C_I2B*A;
V = norm(v_I);
v_I_norm = v_I/V;
v_B_norm = C_I2B*v_I_norm;


% Compute aerodynamic forces on the satellite areas
F = 0;
F_const = 0.5*C_d*rho*V^2;
Theta = v_I_norm'*normA;
for i = 1:6
    if  Theta(i) > 0
        F = F - F_const*Theta(i)*areas(i);
    end
end

% Calculate torque
T_aer = F*cross(v_B_norm,CoP);
end