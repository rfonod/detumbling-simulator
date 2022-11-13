function [data] = main(seed)
%--------------------------------------------------------------------------
%                             MAIN FUNCTION
%           Simplified ADCS Simulator for Delfi-PQ (Detumbling mode)
% Author: Robert Fonod (c)
% Version: 0.8
% Last Modified: 12 October 2018
% Note: muT == micro Tesla
%--------------------------------------------------------------------------
if nargin == 0
    tic, clc
    MC_on = 0;        % Monte Carlo {1 - yes, 0 - no}
    seed  = 1;        % noise seed assignment
else
    MC_on = 1;
end

%% Some Useful Constants and Global Parameters
global I Iinv rho areas
rng(seed,'v5normal')  % seed for randn
rng(seed,'v5uniform') % seed for rand
rng(seed,'twister');  % seed for randi
d2r = pi/180;         % deg to rad conversion
r2d = 1/d2r;          % rad to deg conversion

%% Simulation Related Parameters
N_orb = 10;           % number of orbits (1h day ~ 16 orb) per simulation [-]
q0 = [1 0 0 0]';      % initial attitude (BODY wrt ECI) quaternion [q_v,q_s] [-]
w0 = 30*[1;1;1]*d2r;  % initial angular velocity (BODY wrt ECI) [rad/s]
no_C  = 0;            % duration of no-control after release [s] /1800 = 30 min/
AirDens  = 'low';     % air density {low, medium, high}
att_solv = 'RK4';     % attitude propagation solver {'RK4','ODE45'}
%att_solv = 'ODE45';  % attitude propagation solver {'RK4','ODE45'}

savePlots = 0;        % save plots {1 - yes, 0 - no}
saveExcel = 0;        % save data to excel file for implementation testing

%% Delfi-PQ Related Parameters

% A - Magnetorques (MTQs)
m_rise  = .01;  % magnetorquers rise/fall time [s]
m_x_max = .002; % max magnetic dipole moment along x-axis [A.m^2]
m_y_max = .002; % max magnetic dipole moment along y-axis [A.m^2]
m_z_max = .002; % max magnetic dipole moment along z-axis [A.m^2]
if MC_on        % adding 15% uncertainty truncated at 3sigma
    m_all_unc = truncatedGaussian(m_x_max/15,3*m_x_max/15,3);
    m_x_max = m_x_max + m_all_unc(1);
    m_y_max = m_y_max + m_all_unc(2);
    m_z_max = m_z_max + m_all_unc(3);
end
m_max = [m_x_max m_y_max m_z_max]';

m_res_mag = .0001;   % res. mag. dipole moment magnitude [A.m^2]
if MC_on             % res. mag. dipole moment direction (unit vector) [-]
    az = 2*pi*rand; el = 2*pi*rand;
    m_res_dir = [sin(az)*cos(el); sin(az)*sin(el); cos(az)];
    m_res_dir = m_res_dir/norm(m_res_dir);
    m_res_mag = m_res_mag + truncatedGaussian(m_res_mag/10,3*m_res_mag/10);
else
    m_res_dir = [-0.6045; -0.4251; -0.6737]; % arbitrary direction - sample run
end
m_res = m_res_mag*m_res_dir; % res. mag. dipole moment vector [A.m^2]

m_x_pol = 1;  % magnetorquer polarity in x-axis [-]
m_y_pol = 1;  % magnetorquer polarity in y-axis [-]
m_z_pol = 1;  % magnetorquer polarity in z-axis [-]
m_pol   = [m_x_pol m_y_pol m_z_pol]';

% B - Magnetometers (MTMs)
mag1_rms = .5;       % MTM 1 noise (rms) [muT]
mag2_rms = .5;       % MTM 2 noise (rms) [muT]
mag1_res = .3;       % MTM 1 resolution [muT/LSb] (BMX055=.3, MAG3110=.1)
mag2_res = .3;       % MTM 2 resolution [muT/LSb] (BMX055=.3, MAG3110=.1)
mag1_S2B = [1 0 0;0 1 0;0 0 -1]; % MTM 1 frame to B frame rotation
mag2_S2B = [1 0 0;0 1 0;0 0 -1]; % MTM 2 frame to B frame rotation

mag1_bias_mag = .4;  % MTM 1 bias magnitude [muT]
mag2_bias_mag = .4;  % MTM 2 bias magnitude [muT]
if MC_on             % MTM 1/2 bias directions (unit vectors) [-]
    az1 = 2*pi*rand; el1 = 2*pi*rand;
    az2 = 2*pi*rand; el2 = 2*pi*rand;
    mag1_bias_dir = [sin(az1)*cos(el1); sin(az1)*sin(el1); cos(az1)];
    mag2_bias_dir = [sin(az2)*cos(el2); sin(az2)*sin(el2); cos(az2)];
    mag1_bias_dir = mag1_bias_dir/norm(mag1_bias_dir);
    mag2_bias_dir = mag2_bias_dir/norm(mag2_bias_dir);
else
    mag1_bias_dir = [-0.9416; -0.2857; -0.1781]; % arbitrary dir. for sample run
    mag2_bias_dir = [-0.2574;  0.0577;  0.9646]; % arbitrary dir. for sample run
end
mag1_bias = mag1_bias_mag*mag1_bias_dir; % MTM 1 bias vector [muT]
mag2_bias = mag2_bias_mag*mag2_bias_dir; % MTM 2 bias vector [muT]

w1 = .5; % weight factor of MTM 1
w2 = .5; % weight factor of MTM 2

% C - Mass & Dimension
mass = .6;       % nominal mass  [kg]
mass_uncer = .1; % +/- uncertainty on mass (for MC) [kg]

if MC_on
    mass = mass + truncatedGaussian(mass_uncer,3*mass_uncer);
end

d_x   = .05;   % main body x-axis length [m]
d_y   = .05;   % main body y-axis length [m]
d_z   = .178;  % main body z-axis length [m]

d_x_p = .0640; % plate x-axis length [m]
d_y_p = .0016; % plate y-axis length [m]
d_z_p = .1920; % plate z-axis length [m]

area_x = (d_y+d_y_p)*d_z+(d_z_p-d_z)*d_y_p; % cross-sectional area x face [m^2]
area_y = d_x_p*d_z_p;                       % cross-sectional area y face [m^2]
area_z = d_x*d_y+(d_x_p-d_x)*d_y_p;         % cross-sectional area z face [m^2]
areas = [area_x area_y area_z area_x area_y area_z]; %(x,y,z,-x,-y,-z)
plate_pos = [0;d_y+d_y_p;0]/2; % plate position from main body centre [m]

mbody_vol = d_x*d_y*d_z;       % main body volume [m^3]
plate_vol = d_x_p*d_y_p*d_z_p; % plate volume [m^3]

density = mass/(mbody_vol + plate_vol); % satellite density [kg/m^3]

mbody_mass = density*mbody_vol; % main body mass [kg]
plate_mass = density*plate_vol; % plate body mass [kg]

% D - CoM & Inertia
CoM = plate_mass*plate_pos/mass;

% Mass moment of inertia (with the plate) along x,y,z-axis [kg.m^2]
I_x = mbody_mass*(d_z^2 + d_y^2)/12 + plate_mass*(d_z_p^2 + d_y_p^2)/12 + ...
    mbody_mass*CoM(2)^2 + plate_mass*(plate_pos(2)-CoM(2))^2;
I_y = mbody_mass*(d_z^2 + d_x^2)/12 + plate_mass*(d_z_p^2 + d_x_p^2)/12;
I_z = mbody_mass*(d_y^2 + d_x^2)/12 + plate_mass*(d_y_p^2 + d_x_p^2)/12 + ...
    mbody_mass*CoM(2)^2 + plate_mass*(plate_pos(2)-CoM(2))^2;

% Mass moment of inertia (without the plate) along x,y,z-axis [kg.m^2]
% I_x = mass*(d_z^2 + d_y^2)/12; % mass moment of inertia along x-axis [kg.m^2]
% I_y = mass*(d_z^2 + d_x^2)/12; % mass moment of inertia along y-axis [kg.m^2]
% I_z = mass*(d_y^2 + d_x^2)/12; % mass moment of inertia along z-axis [kg.m^2]

I_nom = [I_x,I_y,I_z];

if MC_on % adding 5% additional uncertainty (16.67% due to mass uncertainty)
    I_x = I_x + truncatedGaussian(I_x/20,3*I_x/20); 
    I_y = I_y + truncatedGaussian(I_y/20,3*I_y/20); 
    I_z = I_z + truncatedGaussian(I_z/20,3*I_z/20); 
end

I    = diag([I_x,I_y,I_z]);    % inertia matrix [kg.m^2]
Iinv = inv(I);                 % inertia matrix inverse [kg^-1.m^-2]

%% Orbit Related Parameters
alt    = 350e3;           % altitude (LEO) [m]
re     = 6378.137e3;      % Earth (equatorial) radius [m]
a      = re + alt;        % semi-major axis (assuming circular orbit) [m]
ecc    = 0;               % eccentricity [-]
inc    = 96.84895*d2r;    % orbital inclination [rad]
omega  = 0;               % argument of perigee (periapsis) [rad]
OMEGA  = 210*d2r;         % RAAN [rad]
M0     = 60*d2r;          % mean anomaly at epoch (t0) [rad]
M      = 5.97218e24;      % mass of the Earth [kg]
G      = 6.67408e-11;     % gravitational constant [m^3/(kg.s^2)]
mu     = 3.9860044181e14; % Earth's standard gravitational parameter [m^3/s^2]
n      = sqrt(mu/a^3);    % mean motion [rad/s]
T_orb  = 2*pi/n;          % orbital period [s]
I2E0   = 0;               % initial rotation of the ECEF wrt ECI frame [rad]
we     = 7.292115e-5;     % Earth's rotational rate (GWS84) [rad/s]
Ldate  = decyear(2018,3,31);    % date for IGRF 12
T_srd  = 24*60^2*365.24/366.24; % sidereal day [s]
C_I2Ev = [0 1 0;-1 0 0;0 0 0];  % for conversion of v_I to v_E

% Air density [kg/m^3], see Ref. [http://www.braeunig.us/space/atmos.htm]
switch AirDens
    case 'low'
        rho = 0.5*(2.50 + 1.510)*1e-12; % lows solar activity [kg/m^3]
    case 'medium'
        rho = 0.5*(1.16 + 0.799)*1e-11; % medium solar activity [kg/m^3]
    case 'high'
        rho = 0.5*(1.03 + 0.805)*1e-10; % extremely high solar activity [kg/m^3]
end

%% Detuble Algorithm Related Parameters
f_c   = 4;            % control/sensing loop frequency [Hz]
T_c   = 1/f_c;        % control/sensing loop sampling rate [s]
delta = .6;           % MTQs duty cycle [fraction of T_c]
w_des = 5*d2r;        % desired detumbling (stopping) rate for all axes [rad/s]
alpha = 1/200;        % filter costant (detumbling parameter) [-]
p_tmb = 2*[1;1;1];    % initial tumbling parameter [muT/s]
p_bar_l = .075;       % lower tumble parameter threshold [muT/s]
p_bar_u = .085;       % upper tumble parameter threshold [muT/s]
C_det = 0;            % detumbled counter [-]
C_tum = 0;            % tumbling counter [-]
t_bar_det = 1*60*60;  % confirmation time for detumbled state [s]
t_bar_tum = 1*2*60;   % confirmation time for tumbling state [s]
g_m   = 11.44*d2r;    % Earth geom. field tild angle wrt the polar axis [rad]
xi_m  = acos(cos(inc)*cos(g_m)+sin(inc)*sin(g_m)); % Ref. [Avanzi & Giulietti]
k_w_T = 2*n*(1+sin(xi_m))*min(I_nom);              % B-dot gain [A.m^2.s/T]
k_w   = k_w_T*1e6;                                 % B-dot gain [A.m^2.s/muT]
%% Load Orbital Data (R_I V_I R_E V_E B_E) if Existent
if ~exist('Existing_Orbits','dir')
    mkdir('Existing_Orbits')
end
addpath('Existing_Orbits')
for i=0:31
    orbit_str = sprintf('Orbit_%d_%d_%d.mat',alt/1e3,f_c,N_orb+i);
    if exist(orbit_str,'file')
        load(orbit_str,'B_E','R_E','R_I','V_I')
        orbit_exists = 1; break;
    end
end
if ~exist('orbit_exists','var')
    orbit_exists = 0;
    orbit_str = sprintf('Orbit_%d_%d_%d.mat',alt/1e3,f_c,N_orb);
    disp('Orbit not yet simulated. Simulation On!')
end

%% Simulation Initialization
ADCS_on  = 0;                   % initial ADCS mode {1 - actuating, 0 - sensing}
q        = q0;                  % initial quaternion [-]
w        = w0;                  % initial angular rate [rad/s]
t_on_cnt = [0;0;0];             % sum of the opening times (per axis) [s]
t0       = 0;                   % initial time [s]
t_tot    = ceil(N_orb*T_orb);   % total simulation time [s]
t        = t0:T_c:t0+t_tot;     % time vector [s]
t_len    = length(t);           % total number of sample/time points [-]
k_orb    = 0;                   % number of orbits passed [-]
[t_det_w,t_det_p] = deal(nan);  % stopping times (real and estimated) [s]
[t_on_sum_w,t_on_sum_p] = deal(nan(3,1)); % sum of t_on times - whole period [s]
t_on_orb = nan(3,N_orb);        % sum of t_on times - per orbit basis [s]

%% Vector Pre-allocation
if ~MC_on
    if ~orbit_exists
        [R_I,V_I,R_E,V_E,B_E] = deal(zeros(3,t_len));
    end
    if saveExcel
        [B_B_meas1,B_B_meas2] = deal(zeros(3,t_len));
        [C_tum_vec,C_det_vec] = deal(nan(t_len,1));
    end
    [W,B_dot,B_B,B_B_meas,M_del,T_on,S_on,P_tmb] = deal(zeros(3,t_len));
    Q = zeros(4,t_len);
end


%% Start Simulation
for k = 1:t_len
    
    %---------------------------------------------------------------------------
    % Environment Simulation
    %---------------------------------------------------------------------------
    
    % Attitude Propagation
    if k>1
        tspan     = [t(k-1),t(k)];
        inp.m_del = m_del;
        inp.m_res = m_res;
        inp.b_I   = b_I;
        inp.v_I   = v_I;
        r_I_norm = r_I/norm(r_I);
        inp.T_gg  = 3*mu*cross(r_I_norm,I*r_I_norm)/(norm(r_I)^3);
        [q, w]    = propag_att(tspan,[q; w],inp,att_solv);
    end
    
    % Rotation matrix from ECI to ECEF
    phi   = mod(I2E0+2*pi*t(k)/T_srd,2*pi);
    C_I2E = rot_mat(phi,3);
    
    % Rotation matrix from ECI to BODY
    C_I2B = q2dcm(q);
    
    if orbit_exists
        b_E = B_E(:,k);
        r_I = R_I(:,k);
        v_I = V_I(:,k);
    else
        % Satellite position and velocity in ECI
        trueAnom = M0 + t(k)*n;
        [r_I,v_I] = kepler2cart(a,ecc,inc,OMEGA,omega,trueAnom,mu);
        
        % Satellite position and velocity in ECEF
        r_E = C_I2E*r_I;
        v_E = C_I2E*(v_I+we*C_I2Ev*r_I);
        
        % Satellite position in Geodetic coordinates
        LLA = ecef2lla(r_E');
        latDeg = LLA(1); lonDeg = LLA(2); altGeo = LLA(3);
        
        % Earth's geomegnetic field in NED (two implementations are equivalent)
        b_NED = igrfmagm(altGeo,latDeg,lonDeg,Ldate,12)'*1e-3;     % [muT]
        %b_NED = wrldmagm(altGeo,latDeg,lonDeg,Ldate,'2015')*1e-3; % [muT]
        
        % Earth's geomagnetic field NED to ECEF conversion
        [b_Ex,b_Ey,b_Ez] = ned2ecefv(b_NED(1),b_NED(2),b_NED(3),latDeg,lonDeg);
        b_E = [b_Ex;b_Ey;b_Ez];
    end
    
    % Magnetic field in ECI and BODY
    b_I = C_I2E'*b_E;
    b_B = C_I2B*b_I;
    
    %---------------------------------------------------------------------------
    % ADCS
    %---------------------------------------------------------------------------
    
    % Sensing (frame transf. + noise + bias + sensor resolution error) [muT]
    b_S_raw1 = round((mag1_S2B'*(b_B + mag1_rms*randn(3,1) + mag1_bias))/mag1_res)*mag1_res;
    b_S_raw2 = round((mag2_S2B'*(b_B + mag2_rms*randn(3,1) + mag2_bias))/mag2_res)*mag2_res;
    
    b_B_meas1 = mag1_S2B*b_S_raw1;
    b_B_meas2 = mag2_S2B*b_S_raw2;
    
    b_B_meas = w1*b_B_meas1 + w2*b_B_meas2;  % weighted average (un-normalized)
    b_B_meas_norm = b_B_meas/norm(b_B_meas); % weighted average (normalized)
    
    if exist('b_B_meas_prev_norm','var')
        b_dot      = (b_B_meas - b_B_meas_prev)/T_c;           % not used
        b_dot_norm = (b_B_meas_norm - b_B_meas_prev_norm)/T_c; % for B-dot
        
        % Tumble parameter update 
        p_tmb = alpha*abs(b_dot_norm) + (1-alpha)*p_tmb;
    else
        b_dot      = zeros(3,1);
        b_dot_norm = zeros(3,1);
    end
    b_B_meas_prev = b_B_meas;
    b_B_meas_prev_norm = b_B_meas_norm;
    
    % Counter update for the "detumbled" state
    if all(p_tmb<=p_bar_l)
        if C_det < 65535 % for a 16 bit (unsigned) inteeger only
            C_det = C_det + 1;
        end
    else
        C_det = 0;
    end
    
    % Counter update for the "tumbling" state
    if any(p_tmb>=p_bar_u)
        if C_tum < 65535 % for a 16 bit (unsigned) inteeger only
            C_tum = C_tum + 1;
        end
    else
        C_tum = 0;
    end
    
    %  Check detumbling stoping criteria - true
    if all(abs(w)<=w_des) && isnan(t_det_w)
        t_det_w = t(k);
        t_on_sum_w = t_on_cnt;
    end
    
    %  Check "detumbled" state stoping criteria - estimate
    if C_det*T_c >= t_bar_det
        if isnan(t_det_p)
            t_det_p = t(k);
            t_on_sum_p = t_on_cnt;
        end
        ADCS_on = 0;
    end
    
    %  Check "tumbling" state start criteria - estimate
    if C_tum*T_c >= t_bar_tum
        ADCS_on = 1;
    end
    
    % Commands computation
    if (t(k) > no_C) && ADCS_on
        m_des = -k_w*b_dot_norm/norm(b_B_meas);
        t_on  = delta*T_c*min(1,abs(m_des)./m_max);
        s_on  = m_pol.*sign(m_des);
        m_del = s_on.*(t_on-m_rise).*m_max/T_c;
    else
        [t_on, s_on, m_del] = deal(zeros(3,1));
    end
    t_on_cnt = t_on_cnt + t_on;
    
    if MC_on && ~isnan(t_det_w) && ~isnan(t_det_p)
        break
    end
    
    %---------------------------------------------------------------------------
    % House Keeping
    %---------------------------------------------------------------------------
    if ~MC_on
        if ~orbit_exists
            R_I(:,k) = r_I;
            V_I(:,k) = v_I;
            R_E(:,k) = r_E;
            V_E(:,k) = v_E;
            B_E(:,k) = b_E;
        end
        if saveExcel
            B_S_raw1(:,k) = b_S_raw1;
            B_S_raw2(:,k) = b_S_raw2;
            C_det_vec(k)  = C_det;
            C_tum_vec(k)  = C_tum;
        end
        Q(:,k) = q;
        W(:,k) = w;
        B_B(:,k) = b_B;
        B_B_meas(:,k)  = b_B_meas;
        B_dot(:,k) = b_dot;
        M_del(:,k) = m_del;
        T_on(:,k)  = t_on;
        S_on(:,k)  = s_on;
        P_tmb(:,k) = p_tmb;
    end
    
    %---------------------------------------------------------------------------
    % Orbit Count Display
    %---------------------------------------------------------------------------
    if  t(k)>=(k_orb+1)*T_orb
        k_orb = k_orb + 1;
        if k_orb == 1
            t_on_orb(:,k_orb) =  t_on_cnt(:);
        else
            t_on_orb(:,k_orb) =  t_on_cnt(:) - sum(t_on_orb(:,1:k_orb-1),2);
        end
        if ~MC_on, fprintf('Orbit %d/%d passed\n',k_orb,N_orb); end
    end 
end

%% Save Orbital Data (R_I V_I R_E V_E B_E) if Non-existent
if ~orbit_exists && ~MC_on
    cd('Existing_Orbits')
    save(orbit_str,'R_I','V_I','R_E','V_E','B_E')
    cd('..')
end

%% Plot Results
if saveExcel
    T = table(B_S_raw1',B_S_raw2',T_on',S_on',P_tmb',C_tum_vec,C_det_vec);
    filename = 'data4test_03_in.csv';
    writetable(T,filename,'Delimiter',',')
end

if MC_on
    data(1,1)     = t_det_w;
    data(2,1)     = t_det_p;
    data(3,1)     = mass;
    data(4:6,1)   = diag(I);
    data(7:9,1)   = m_max;
    data(10:12,1) = m_res;
    data(13:15,1) = p_tmb;
    data(16:18,1) = w0;
    data(19:21,1) = t_on_sum_w;
    data(22:24,1) = t_on_sum_p;
    data(25:25+N_orb-1,1) = t_on_orb;
else
    toc
    all_var = who;
    for i=1:size(all_var,1)
        putvar(all_var{i})
    end
    plots
end

end