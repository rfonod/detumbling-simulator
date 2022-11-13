function [parameters] = SolarActivity(r_I, Sdate, activity, activityA, magnInd, t)
%SOLARACTIVITY Calculates different parameters dependant on the solar
%activity (low, med, high)
 %   Densities: [Alt. low med high]
 %   PhiC:      Constant in solar radiation torque formula, phi/C = K/R^2

%% 
if nargin < 2
    Sdate = [2018, 1, 4];
end
if nargin < 3
    activity = 'high';    Sdate = [2018, 1, 4];
end
parameters = struct('rho',0,'phiC',0');

%% Constants
c = 299792458;  % Speed of light

%% Solar Constant
dayNum   = day(datetime(Sdate), 'dayofyear');         % Day number
phi      = 1366.1.*(1+0.033*cosd(360*dayNum/365.25)); % Solar constant
SCchange = 0.1/100;                                   % Change in solar constant with high / low activity
rescaled = 2*activity - 1;      % [-1, 1] 
parameters.phiC = phi*(1+SCchange*rescaled)/c;

%% Density
parameters.rho = nrlmsise00(Sdate, r_I, activity, activityA, magnInd, t);

end