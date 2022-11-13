function [r_out,v_out] = kepler2cart(a,e,i,Omega,omega,nu,mu)
% Kepler elements to Cartesian components (in ECI)
% Inputs are in [m] and [rad]
% Outputs are in [m] and [m/s]

% WGS84
if nargin <= 6
    mu = 3.9860044181e14; % Earth's standard gravitational parameter [m^3/s^2]
end

r      = a*(1-e^2)/(1+e*cos(nu));
nu_dot = sqrt(mu/a^3)*(1+e*cos(nu))^2/(sqrt((1-e^2)^3));
r_dot  = a*(1-e^2)*e*sin(nu)*nu_dot/((1+e*cos(nu))^2);

C1 = [cos(Omega)*cos(nu+omega)-sin(Omega)*cos(i)*sin(nu+omega)
      sin(Omega)*cos(nu+omega)+cos(Omega)*cos(i)*sin(nu+omega)
      sin(i)*sin(nu+omega)];

C2 = [-cos(Omega)*sin(nu+omega)-sin(Omega)*cos(i)*cos(nu+omega)
      -sin(Omega)*sin(nu+omega)+cos(Omega)*cos(i)*cos(nu+omega)
      sin(i)*cos(nu+omega)];

r_out = r*C1;
v_out = r_dot*C1 + r*nu_dot*C2;

end

