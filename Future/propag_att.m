function [q_out,w_out] = propag_att(tspan,xk,inp,att_solv)

switch att_solv
    case 'ODE45'
        %opts = odeset('RelTol',1e-15,'AbsTol',1e-15); %,'InitialStep',1e-4,'Maxstep',1e-4);
        opts = [];
        [~,x_out] = ode45(@(t,x) att_dyn_kin(x,inp),tspan,xk,opts);
        q_out = x_out(end,1:4)'./norm(x_out(end,1:4)); % normalize quaternion
        w_out = x_out(end,5:7)';
    case 'RK4'
        dt = diff(tspan);
        rk_step = ceil(dt/0.125); 
        % cFreq = 1 [Hz] -> rk_step = 20 & dt = 0.125
        % cFreq = 2 [Hz] -> rk_step = 10 & dt = 0.125
        % cFreq = 4 [Hz] -> rk_step = 5  & dt = 0.125
        
        dt = dt/rk_step;
        
        for i = 1:rk_step
            xk = RK4(@att_dyn_kin,dt,xk,inp);
            xk(1:4) = xk(1:4)/norm(xk(1:4));
        end
        
        q_out = xk(1:4);
        w_out = xk(5:7);
        
    otherwise
        error('Wrong Attitude Propagation Solver Selected')
end
%disp(norm(q_out)-1) % check the normalization

end

function dx = att_dyn_kin(x,inp)
global I Iinv
q  = x(1:4);
w  = x(5:7);

%% Attitude Dynamics

% Magnetic field in BODY at the given instant of time [T]
C_I2B = q2dcm(q);
b_B   = C_I2B*inp.b_I*1e-6;    % ECI to B conversion anc conversion to [T]

acc = 100;
T_mtq = cross(inp.m_del,b_B);  % torques due to the MTQs [Nm]
T_res = cross(inp.m_res,b_B);  % dist. torques due to residual mag. dipole [Nm]
T_gg  = gg_torque(inp.r_I);    % dist. torques due to gravity gradient [Nm]
T_aer = drag_torque(C_I2B,inp.v_I,acc);% dist. torque due to atmospheric drag [Nm]
T_sr  = solar_torque(C_I2B,inp.r_I,inp.s_I,inp.d); % dist. torque due to solar radiation [Nm]

T_all = T_mtq + T_res + T_gg + T_sr + T_aer;

dw = Iinv*(T_all - cross(w,I*w));

%% For plotting (Comment / delete this during Monte Carlo)
global k Tres Tgg Taer Tsr t_len
Tres(1:3,k) = T_res;
Tgg(1:3,k)  = T_gg;
Taer(1:3,k) = T_aer;
Tsr(1:3,k)  = T_sr;
% Save them for quicker access
if k == t_len
    fprintf('Saving Torques')
    save('Torques.mat', 'Tres', 'Tgg', 'Taer', 'Tsr')
end
%% Attitude Kinematics (the two equations are equivalent, but the 1st one executes slightly faster)
dq = 0.5*W_matrix(w)*q;
%dq = 0.5*Q_matrix(q)*w;

dx = [dq; dw];
end

function W = W_matrix(w)
W = [0     w(3) -w(2)  w(1)
    -w(3)  0     w(1)  w(2)
     w(2) -w(1)  0     w(3)
    -w(1) -w(2) -w(3)  0   ];
end

function Q = Q_matrix(q)
Q = [-q(2) -q(3) -q(4)
    q(1) -q(4)  q(3)
    q(4)  q(1) -q(2)
    -q(3)  q(2)  q(1) ];
end

function S = S_matrix(w)
% Skew symetric Matrix
S = [0    -w(3)  w(2)
    w(3)  0    -w(1)
    -w(2)  w(1)  0   ];
end

function [x_out] = RK4(odefun,dt,x0,inp)
% Runge Kutta 4
k1 = dt*feval(odefun,x0,inp);
k2 = dt*feval(odefun,x0+0.5*k1,inp);
k3 = dt*feval(odefun,x0+0.5*k2,inp);
k4 = dt*feval(odefun,x0+k3,inp);

% Output
x_out = x0 + (k1+2*(k2+k3)+k4)/6;
end