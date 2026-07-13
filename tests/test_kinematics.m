function tests = test_kinematics
% Unit tests for the attitude and orbit kinematics helpers.
% Run with: runtests('tests')
tests = functiontests(localfunctions);
end

function test_q2dcm_identity(testCase)
% The identity quaternion [q_v; q_s] = [0;0;0;1] must map to the identity DCM
verifyEqual(testCase, q2dcm([0;0;0;1]), eye(3), 'AbsTol', 1e-12);
end

function test_q2dcm_orthonormal(testCase)
% Any unit quaternion must map to a proper rotation matrix
rng(0,'twister');
for k = 1:20
    q = randn(4,1); q = q/norm(q);
    C = q2dcm(q);
    verifyEqual(testCase, C*C', eye(3), 'AbsTol', 1e-12);
    verifyEqual(testCase, det(C), 1, 'AbsTol', 1e-12);
end
end

function test_q2dcm_x_rotation(testCase)
% A 90 deg rotation about the x-axis, as a quaternion and via rot_mat
q = [sin(pi/4); 0; 0; cos(pi/4)];
verifyEqual(testCase, q2dcm(q), rot_mat(pi/2,1), 'AbsTol', 1e-12);
end

function test_rot_mat_zero_angle(testCase)
for axis = 1:3
    verifyEqual(testCase, rot_mat(0,axis), eye(3), 'AbsTol', 1e-12);
end
end

function test_rot_mat_orthonormal(testCase)
for axis = 1:3
    C = rot_mat(0.7,axis);
    verifyEqual(testCase, C*C', eye(3), 'AbsTol', 1e-12);
    verifyEqual(testCase, det(C), 1, 'AbsTol', 1e-12);
end
end

function test_rot_mat_bad_axis_falls_back_to_identity(testCase)
% An unknown axis warns (unidentified warning) and returns the identity
w = warning('off','all');
restore = onCleanup(@() warning(w)); %#ok<NASGU>
verifyEqual(testCase, rot_mat(0.1,4), eye(3), 'AbsTol', 1e-12);
end

function test_kepler2cart_circular_orbit(testCase)
% On a circular orbit the radius is constant, the speed is sqrt(mu/a),
% and the position and velocity vectors are perpendicular
mu = 3.9860044181e14;
a  = 6378.137e3 + 350e3;
for nu = linspace(0,2*pi,9)
    [r,v] = kepler2cart(a,0,deg2rad(96.8),deg2rad(210),0,nu,mu);
    verifyEqual(testCase, norm(r), a, 'RelTol', 1e-10);
    verifyEqual(testCase, norm(v), sqrt(mu/a), 'RelTol', 1e-10);
    verifyEqual(testCase, dot(r,v)/(norm(r)*norm(v)), 0, 'AbsTol', 1e-10);
end
end

function test_kepler2cart_inclination(testCase)
% The orbit normal must be inclined by exactly the requested inclination
mu  = 3.9860044181e14;
a   = 6378.137e3 + 350e3;
inc = deg2rad(96.84895);
[r,v] = kepler2cart(a,0,inc,deg2rad(210),0,deg2rad(60),mu);
h = cross(r,v);
verifyEqual(testCase, acos(h(3)/norm(h)), inc, 'AbsTol', 1e-9);
end
