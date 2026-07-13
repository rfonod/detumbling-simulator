function tests = test_dynamics
% Unit tests for the attitude propagator and the aerodynamic torque model.
% Both read the simulator's global parameters, which setupOnce initializes.
% Run with: runtests('tests')
tests = functiontests(localfunctions);
end

function setupOnce(~)
global I Iinv rho areas %#ok<GVMIS>
I     = diag([1.2e-3, 2.0e-3, 0.6e-3]); % representative picosat inertia [kg.m^2]
Iinv  = inv(I);
rho   = 0;                              % vacuum unless a test sets otherwise
areas = [8.9e-3 1.2e-2 2.6e-3 8.9e-3 1.2e-2 2.6e-3]; % (x,y,z,-x,-y,-z) [m^2]
end

function test_quaternion_stays_normalized(testCase)
global rho %#ok<GVMIS>
rho = 0;
inp = struct('m_del',zeros(3,1),'m_res',zeros(3,1),'b_I',zeros(3,1), ...
    'v_I',[7.7e3;0;0],'T_gg',zeros(3,1));
q = [0.2; -0.4; 0.1; 0.887]; q = q/norm(q);
w = deg2rad([180; -120; 60]);
for k = 1:40
    [q,w] = propag_att([0 0.25],[q; w],inp,'RK4');
end
verifyEqual(testCase, norm(q), 1, 'AbsTol', 1e-10);
end

function test_torque_free_conserves_momentum_and_energy(testCase)
% With no external torque, Euler's equations conserve the angular momentum
% magnitude and the rotational kinetic energy. Checked at the 30 deg/s rate
% of the default run, where the RK4 truncation error is negligible.
global I rho %#ok<GVMIS>
rho = 0;
inp = struct('m_del',zeros(3,1),'m_res',zeros(3,1),'b_I',zeros(3,1), ...
    'v_I',[7.7e3;0;0],'T_gg',zeros(3,1));
q = [0; 0; 0; 1];
w = deg2rad([30; 30; 30]);
H0 = norm(I*w);
E0 = 0.5*w'*I*w;
for k = 1:100
    [q,w] = propag_att([0 0.25],[q; w],inp,'RK4');
end
verifyEqual(testCase, norm(I*w), H0, 'RelTol', 1e-4);
verifyEqual(testCase, 0.5*w'*I*w, E0, 'RelTol', 1e-4);
end

function test_fast_tumble_energy_drift_is_integrator_truncation(testCase)
% At the paper's 180 deg/s tumble, the RK4 sub-step of 0.125 s that propag_att
% picks for a 4 Hz control loop loses a few percent of the rotational energy
% over 25 s of torque-free motion. Refining the sub-step by a factor of ten
% recovers it, which shows the drift is integrator truncation and not a
% modelling error. Fidelity studies at fast tumble rates should use ODE45 or a
% shorter propagation span.
global I rho %#ok<GVMIS>
rho = 0;
inp = struct('m_del',zeros(3,1),'m_res',zeros(3,1),'b_I',zeros(3,1), ...
    'v_I',[7.7e3;0;0],'T_gg',zeros(3,1));
w0 = deg2rad([180; 180; 180]);
E0 = 0.5*w0'*I*w0;

% Coarse: the 0.125 s sub-step used by the simulator at f_c = 4 Hz
q = [0; 0; 0; 1]; w = w0;
for k = 1:100
    [q,w] = propag_att([0 0.25],[q; w],inp,'RK4');
end
drift_coarse = abs(0.5*w'*I*w - E0)/E0;

% Fine: the same 25 s, propagated in 0.0125 s steps
q = [0; 0; 0; 1]; w = w0;
for k = 1:2000
    [q,w] = propag_att([0 0.0125],[q; w],inp,'RK4');
end
drift_fine = abs(0.5*w'*I*w - E0)/E0;

verifyLessThan(testCase, drift_coarse, 0.05);   % a few percent, not a blow-up
verifyLessThan(testCase, drift_fine, 1e-4);     % vanishes as the step shrinks
verifyLessThan(testCase, drift_fine, drift_coarse);
end

function test_rk4_and_ode45_agree(testCase)
% The two solvers must propagate the same trajectory to within their tolerance
global rho %#ok<GVMIS>
rho = 0;
inp = struct('m_del',zeros(3,1),'m_res',[1e-4;-1e-4;5e-5],'b_I',[20;-15;30], ...
    'v_I',[7.7e3;0;0],'T_gg',zeros(3,1));
x0 = [0; 0; 0; 1; deg2rad([30; -20; 10])];
[q_rk,w_rk] = propag_att([0 0.25],x0,inp,'RK4');
[q_od,w_od] = propag_att([0 0.25],x0,inp,'ODE45');
% Tolerance is set by ODE45's default RelTol (1e-3), not by RK4
verifyEqual(testCase, q_rk, q_od, 'AbsTol', 1e-3);
verifyEqual(testCase, w_rk, w_od, 'AbsTol', 1e-3);
end

function test_solver_name_is_validated(testCase)
inp = struct('m_del',zeros(3,1),'m_res',zeros(3,1),'b_I',zeros(3,1), ...
    'v_I',[7.7e3;0;0],'T_gg',zeros(3,1));
verifyError(testCase, ...
    @() propag_att([0 0.25],[0;0;0;1;0;0;0],inp,'EULER'), ?MException);
end

function test_magnetorquer_torque_spins_up_the_satellite(testCase)
% A dipole perpendicular to the field must produce a non-zero torque
global rho %#ok<GVMIS>
rho = 0;
inp = struct('m_del',[0.002;0;0],'m_res',zeros(3,1),'b_I',[0;0;30], ...
    'v_I',[7.7e3;0;0],'T_gg',zeros(3,1));
[~,w] = propag_att([0 0.25],[0;0;0;1;0;0;0],inp,'RK4');
verifyGreaterThan(testCase, norm(w), 0);
end

function test_no_aero_torque_in_vacuum(testCase)
global rho %#ok<GVMIS>
rho = 0;
verifyEqual(testCase, dist_aer([7.7e3;0;0],[0;0;0;1]), zeros(3,1), 'AbsTol', 1e-15);
end

function test_aero_torque_scales_with_density(testCase)
% The aerodynamic torque is linear in the air density
global rho %#ok<GVMIS>
v_I = [7.7e3; 1.0e2; -50];
q   = [0.1; 0.2; 0.3; 0.927]; q = q/norm(q);
rho = 1e-12; T1 = dist_aer(v_I,q);
rho = 2e-12; T2 = dist_aer(v_I,q);
verifyEqual(testCase, T2, 2*T1, 'RelTol', 1e-10);
verifyGreaterThan(testCase, norm(T1), 0);
end
