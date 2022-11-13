function [s,d] = SolUnitVec(Sdate,t)
%SOLARUNITVECTOR calculates ECI frame unit vector s pointing towards the
%sun, which is at distance d in [m].
% Ref: http://www.dept.aoe.vt.edu/~cdhall/courses/aoe4140/attde.pdf
%     Ldate = date in decimal year fraction
%     t     = time (seconds)
%%
Ldate = decyear(Sdate(1), Sdate(2), Sdate(3));
JD = DecDate2JulianDate(Ldate,t); % Julian date

TUT1       = (JD-2451545)/36525;
LambdaMSun = 280460618 + 36000.7705361*TUT1; % Mean longitude of sun
TTDB       = TUT1;
MSun       = 357.5277233 + 35999.05034*TTDB; % Mean anomaly of sun

% Ecliptic longitude sun:
LambdaEcliptic = LambdaMSun + 1.91466471*sind(MSun) + 0.918994643*sind(2*MSun);
Epsilon        = 23.439291 - 0.0130042*TTDB;

s = [cosd(LambdaEcliptic); cosd(Epsilon)*sind(LambdaEcliptic); sind(Epsilon)*sind(LambdaEcliptic)];
s = s./norm(s);
d = 1.000140612-0.016708617*cosd(MSun)-0.00139589*cosd(2*MSun); % Distance [AU]
d = 149597870700*d;                                             % [AU] to [m]
end

function [JD] = DecDate2JulianDate(Ldate,t)
% Check for leap year or not
if mod(floor(Ldate),4) == 0
    ndays = 366;
else
    ndays = 365;
end

time = t / (ndays*24*3600);
date = Ldate + time;
year = floor(date);
partialYear = mod(date,1);
date0 = datenum(num2str(year), 'yyyy');
date1 = datenum(num2str(year+1),'yyyy');
daysInYear = date1-date0;
num = date0 + partialYear*daysInYear;
date = datevec(num);
% Convert date to julian date according to source
year = date(1); month = date(2); day = date(3); hour = date(4); min = date(5); sec = date(6);
JD = juliandate(date);
end