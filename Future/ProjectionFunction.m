function [state, area_bright, dis_boundary, dis_circle, R0] = ProjectionFunction(a, b, r_I, s_I,dist_sun_earth)
%SHADOWFUNCTION Calculate Psi in range [0,1] to determine the amount of
%sunlight on the satellite. ECI frame with positive z-axis pointing toward
%the sun (s_I = z);
%   Detailed explanation goes here
%% Define parameters
r_eq = a;
r_po = b;
r = r_I;
s = s_I*dist_sun_earth;
ds_s = s - r;  % Vector pointing from sat to sun
dis_sat_sun   = norm(ds_s);
dis_sat_earth = norm(r);
rs = 696392;   % [km], radius sun

n = (r-s)/norm(r-s);    % Unit vector along ISF z-axis
A = diag([1/r_eq^2, 1/r_eq^2, 1/r_po^2]);
%% Determine full phase or not
s1 =  sqrt(1/(n'/A*n));
s2 = -s1;

t1 = (dot(r,n) - 1/s1)/(dot(n,n));
t2 = (dot(r,n) - 1/s2)/(dot(n,n));

xs1 = r - t1*n;
xs2 = r - t2*n;

ds1 = norm(xs1 - s);
ds2 = norm(xs2 - s);

if ds1<=ds2; temp = ds1; else; temp = ds2;
end
if ds1>=ds2; ds2  = ds1; else; ds2 = ds2;
end
ds1 = temp;

Delta = (s'*A*n)^2 - n'*A*n*(s'*A*s - 1);

if Delta > 0 % sun-sat line intersects the Earth ellipsoid
    t1 = (-2*s'*A*n + sqrt(Delta))/2/(n'*A*n);
    t2 = (-2*s'*A*n - sqrt(Delta))/2/(n'*A*n);
    
    if t1 < t2; ds1 = t2; else; ds1 = t2;
    end
else
    p  = cross(r, n);
    ns = cross(n, p);
    ns = ns/norm(ns);
    t  = ns'/A*ns;
    lam1 = sqrt(1/t);
    lam2 = -lam1;
    
    s1 = dot(ns, s) - lam1*t;
    s2 = dot(ns, s) - lam2*t;
    
    if abs(s1) < abs(s2)
        lam = lam1;
    else
        lam = lam2;
    end
    dis = lam*(n'/A*ns) - dot(n, s);
    
    ds1 = dis;
end

%% Check Positioin of Sat
if dis_sat_sun <= ds1 % Satellite between Earth and Sun
    state = 0;
    [area_bright, dis_boundary, dis_circle, R0] = deal(zeros(1,1));
    return
elseif dis_sat_sun > ds2 % Ellipse
    state = 2;
elseif dis_sat_sun > ds1 && dis_sat_sun <= ds2 % Hyperbola
    state = 1;
else
    error('Could not determine state!')
end

%% Build ISF Coordinate system
f = 1000;             % Some random number
o = r-f*n;            % Origin
u = n - 1/dot(r,n)*r; % Unit vector along ISF x-axis, does not exist when in full phase or umbra
u = u/norm(u);
v = cross(n,u);       % Unit vector along ISF y-axis;

% Projection of Earth's center PEC in ECEF
t_dis = sqrt( (f*dis_sat_earth/dot(r,n))^2 - f^2);
PEC = t_dis*u;
PEC_isf = [dot(u, PEC); dot(v,PEC); dot(n,PEC)]; % Projection Earth Center in ISF frame

%% Compute M and k0-k5
M = A*(r*r')*A' - (r'*A*r-1)*A;

k0 = u'*M*u*1E6;
k1 = u'*M*v*1E6;
k2 = v'*M*v*1E6;
k3 = -2*f*u'*M*n*1E6;
k4 = -2*f*v'*M*n*1E6;
k5 = f^2*n'*M*n*1E6;

k = [k0 k1 k2 k3 k4 k5];
Omega = (k(3)*k(4)^2-2*k(2)*k(4)*k(5)+k(1)*k(5)^2)/(4*(k(1)*k(3)-k(2)^2));
B_mat = [k0 k1; k1 k2];

temp = sqrt( (k0+k2)^2 - 4*(k0*k2 - k1^2) );
lambda(1) = (k0 + k2 + temp)/2;
lambda(2) = (k0 + k2 - temp)/2;

if abs(k1) < 1E-12
    r1 = [1 0]';
    r2 = [0 1]';
else
    r1 = [lambda(1) - k2; k1];
    r2 = [lambda(2) - k2; k1];
end
r1 = r1/norm(r1);
r2 = r2/norm(r2);
tx = (r1(1)*k(4) + r1(2)*k(5))/lambda(1);
ty = (r2(1)*k(4) + r2(2)*k(5))/lambda(2);

% Semi-major and semi-minor axis
AA = (Omega - k5)/lambda(1);
BB = (Omega - k5)/lambda(2);
R0 = f*rs/dis_sat_sun;

%% Solve the quartic
A = lambda(1)*(R0-1/2*tx)^2 + 1/4*lambda(2)*ty^2 - Omega + k(6);
B = 2*lambda(2)*R0*ty;
C = lambda(1)*(1/2*tx^2-2*R0^2) + lambda(2)*(1/2*ty^2+4*R0^2) - 2*Omega + 2*k(6);
D = 2*lambda(2)*R0*ty;
E = lambda(1)*(R0+1/2*tx)^2 + 1/4*lambda(2)*ty^2 - Omega + k(6);

ce = [B/A; C/A; D/A; E/A];

[numreal, x] = EtaFinder(ce(1), ce(2), ce(3), ce(4));
x = sort(x, 'descend');

% Find intersection sun and earth
if numreal == 3 || numreal == 4
    error('WARNING: ProjectionFunction: %d intersections!', numreal);
end

%% Calculate the sun projection centre Os in the transformed frame
Os = [0.5*tx; 0.5*ty];

PEC_new(1) = PEC_isf(1)*r1(1) + PEC_isf(2)*r1(2) + 1/2*tx;
PEC_new(2) = PEC_isf(1)*r2(1) + PEC_isf(2)*r2(2) + 1/2*ty;
PEC_new = PEC_new';


Q1 = [ (1-x(1)^2)/(1+x(1)^2)*R0 + 1/2*tx; 2*x(1)/(1+x(1)^2)*R0 + 1/2*ty ];
Q2 = [ (1-x(2)^2)/(1+x(2)^2)*R0 + 1/2*tx; 2*x(2)/(1+x(2)^2)*R0 + 1/2*ty ];

% Figure out the intersection between the line from PEC_new to OS and the
% conical curve
PEC_Os_len = norm(Os - PEC_new);
d  = (Os - PEC_new)/PEC_Os_len;

aa = lambda(1)*(d(1))^2 + lambda(2)*(d(2))^2;
bb = 2*(lambda(1)*PEC_new(1)*d(1) + lambda(2)*PEC_new(2)*d(2));
cc = lambda(1)*(PEC_new(1))^2     + lambda(2)*(PEC_new(2))^2 + k(6) - Omega;

t1 = (-bb + sqrt(bb^2-4*aa*cc))/(2*aa);
t2 = (-bb - sqrt(bb^2-4*aa*cc))/(2*aa);

if t1*t2 < 0
    if t1>0;  temp=t1; else; temp = t2;
    end
else
    if t1<t2; temp=t1; else; temp = t2;
    end
end
dis_boundary = temp;

t = [tx; ty];
% Figure out the intersection between the line form PEC_new to Os and the
% solar circle
aa = (norm(d))^2;
bb = 2*dot(PEC_new, d) - dot(d,t);
cc = (norm(PEC_new))^2 - dot(PEC_new, t) + 1/4*norm(t)^2 - R0^2;

t1 = (-bb + sqrt(bb^2-4*aa*cc))/aa/2;
t2 = (-bb - sqrt(bb^2-4*aa*cc))/aa/2;

if t1*t2 < 0
    if t1>0;  temp=t1; else; temp = t2;
    end
else
    if t1<t2; temp=t1; else; temp = t2;
    end
end
dis_circle = temp;


%% Check if ellipse
if det(B_mat) > 0  % Ellipse
    if Omega/(Omega-k(6)) < 1; in_out = 1; else; in_out = 0;
    end
    
    if in_out == 1 && numreal < 2 % Umbra
        state = -1;
        area_bright = 0;
        return
    end
    if in_out == 0 && numreal <2 % Full phase
        state = 0;
        area_bright = pi*R0^2;
        return
    end
    % Penumbra
    area_bright = area_ellipse(Q1, Q2, Os, R0, sqrt(AA), sqrt(BB), in_out);
    
elseif det(B_mat) < 0 % Hyperbola
    if Omega/(Omega-k(6))<=1; in_out = 0; else; in_out = 1;
    end
    
    if in_out == 1 && numreal < 2
        state = -1; % Umbra
        area_bright = 0;
        return
    end
    if in_out == 0 && numreal < 2
        state = 0;  % Full phase
        area_bright = pi*rs^2;
        return
    end
    % Penumbra
    if AA > 0
        a = sqrt(AA);
        x_axis = 1;
    else
        a = sqrt(-AA);
        x_axis = 0;
    end
    
    if BB > 0
        b = sqrt(BB);
        x_axis = false;
    else
        b = sqrt(-BB);
    end
    area_bright = area_hyperbola(Q1, Q2, Os, PEC_new, R0, a, b, in_out, x_axis);
    
end

end

function [num_real, x] = EtaFinder(a,b,c,d)
%ETAFINDER Finds the roots of the 'eta-function' x^4 + a*x^3 + b*x^2 + c*x + d = 0
%          and returns only the real roots.
%   a,b,c,d  = coefficients of 'eta-function'
%   num_real = number of real solutions (= number of intersections between
%              earth and sun projection)
%   x        = root of 'eta-function'

eps = 1E-5;   % Threshold for complex part
num_real = 0; % Number of real solutions

retval = roots([1,a,b,c,d]); % Find all (real and complex) roots
x = [0 0];                   % Initialize real roots, if none than returns zeros

%% Find real roots
for i = 1:4
    if abs(imag(retval(i))) < eps
        num_real = num_real + 1;
        x(num_real) = real(retval(i));
    end
end
end

function [area] = area_ellipse(Q1, Q2, Os, R0, a, b, in_out)
%AREA_ELLIPSE Used to calculated the area of the elliptical eclipse
%   Detailed explanation goes here
%% Define parameters
as = Q1 - Os;
bs = Q2 - Os;

ae = Q1;
be = Q2;

TQ1Q2Os = 0.5*abs(as(1)*bs(2) - bs(1)*as(2));
TQ1Q2Oe = 0.5*abs(ae(1)*be(2) - be(1)*ae(2));

cts = dot(as, bs)/(norm(as)*norm(bs));

SQ1Q2Os = 0.5*acos(cts)*R0^2;
SQ1Q2Oe = 0;

if ae(1) > 0 && ae(2) > 0 % First quadrant
    q1 = atan(ae(2)/ae(1));
elseif ae(1) < 0 && ae(2) > 0 % Second quadrant
    q1 = pi - atan(-ae(2)/ae(1));
elseif ae(1) < 0 && ae(2) < 0
    q1 = pi + atan(ae(2)/ae(1));
elseif ae(1) > 0 && ae(2) < 0
    q1 = 2*pi - atan(-ae(2)/ae(1));
end

if be(1) > 0 && be(2) > 0 % First quadrant
    q2 = atan(be(2)/be(1));
elseif be(1) < 0 && be(2) > 0 % Second quadrant
    q2 = pi - atan(-be(2)/be(1));
elseif be(1) < 0 && be(2) < 0
    q2 = pi + atan(be(2)/be(1));
elseif be(1) > 0 && be(2) < 0
    q2 = 2*pi - atan(-be(2)/be(1));
end

s_a = a*b/2*(q1 - atan2( (b-a)*sin(2*q1), (a+b) + (b-a)*cos(2*q1) ) );
s_b = a*b/2*(q2 - atan2( (b-a)*sin(2*q2), (a+b) + (b-a)*cos(2*q2) ) );

SQ1Q2Oe = abs(s_b - s_a);
S1 = SQ1Q2Oe - TQ1Q2Oe;

%% Calculate areas
if in_out == 1
    area = SQ1Q2Os - TQ1Q2Os - S1;
else
    area_shadow = SQ1Q2Os - TQ1Q2Os + S1;
    area = pi*R0^2 - area_shadow;
end



end

function [area] = area_hyperbola(Q1, Q2, Os, Oe, Rs, a, b, in_out, x_axis)
%AREA_HYPERBOLA Used to calculate the solar area of the hyperbolic eclipse
%   Detailed explanation goes here
%% Define parameters
as = Q1 - Os;
bs = Q2 - Os;

ae = Q1 - Oe;
be = Q2 - Oe;

TQ1Q2Os = 0.5*abs(as(1)*bs(2) - bs(1)*as(2));
TQ1Q2Oe = 0.5*abs(ae(1)*be(2) - be(1)*ae(2));

cts = dot(as, bs)/(norm(as)*norm(bs));
SQ1Q2Os = 0.5*acos(cts)*Rs^2;

%% Calculate area
if x_axis == 0
    tmp = Q1(1); Q1(1) = Q1(2); Q1(2) = -tmp;
    tmp = Q2(1); Q2(1) = Q2(2); Q2(2) = -tmp;
    tmp = a;
    a = b;
    b = tmp;
end

xQ1 = abs(Q1(1));
xQ2 = abs(Q2(1));

s1 = b*( xQ1*sqrt(  (xQ1*xQ1/a/a) -1.0 ) - a*acosh(xQ1/a) );
s2 = b*( xQ2*sqrt(  (xQ2*xQ2/a/a) -1.0 ) - a*acosh(xQ2/a) );

if isnan(s1) || isnan(s2)
    testc = 0;
    warning('s1 and s2 are NaN!')
end

if s1 <= s2; ss = s1; else; ss = s2;
end
if s1 >= s2; s1 = s1; else; s1 = s2;
end

if Q1(2)*Q2(2) < 0
    if abs(Q2(1) - Q1(1)) > 1E-10
        k = (Q2(2) - Q1(2))/(Q2(1) - Q1(1));
        m = Q1(2) - k*Q1(1);
        x_m = -m/k;
    else
        x_m = Q1(1);
    end
    
    if Q1(1) > Q2(1)
        S2 = ss + (s1-ss)/2 - 1/2*abs( (x_m-Q1(1))*Q1(2) ) + 1/2*abs( (Q2(1)-x_m)*Q2(2) );
    else
        S2 = ss + (s1-ss)/2 - 1/2*abs( (Q2(1)-x_m)*Q2(2) ) + 1/2*abs( (x_m-Q1(1))*Q1(2) );
    end
else
    s_trapizium = 0.5*abs( Q1(1)-Q2(1) )*( abs(Q1(2)) + abs(Q2(2)) );
    S2 = (s1 - ss -2*s_trapizium)/2;
end

if in_out == 1 % Centre of sun inside ellipse
    area = SQ1Q2Os - TQ1Q2Os - S2;
else % Centre of sun outside ellipse
    area_shadow = SQ1Q2Os - TQ1Q2Os + S2;
    area = pi*Rs^2 - area_shadow;
end

end
