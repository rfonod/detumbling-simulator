function [T] = gg_torque(r_I)
% Gravity Gradient Torque
% Description TBC
% This function should by a more accuate model

% Edited: 26-2-19 (T)
% Added J2 gravity gradient model:
% https://ac.els-cdn.com/S1270963800010816/1-s2.0-S1270963800010816-main.pdf?_tid=2d77e6d7-42c1-4f5d-8fee-95831285f9dd&acdnat=1551205020_afbd68964bb275d41bb57877bc84f0fe

global I mu
% 
% r_I_mag  = norm(r_I);
% r_I_norm = r_I/r_I_mag;
% T  = 3*mu*cross(r_I_norm,I*r_I_norm)/(r_I_mag^3);

z = [0 0 1]';
r_I_mag  = norm(r_I);
J2 = 1.081874E-3; % Needs reference

T = 3*mu/r_I_mag^5*cross(r_I, I*r_I) + ...
    mu*J2*(-3*cross(z,I*z)/r_I_mag^5 + 15*r_I(3)/r_I_mag^7*(cross(z,I*r_I) + cross(r_I,I*z)) + ...
    15/2*(1/r_I_mag^7 - 7*r_I(3)^2/r_I_mag^9)*cross(r_I,I*r_I));

end
