function [T] = drag_torque(C_I2B,v_I,acc)
%DRAG_TORQUE function to compute the aerodynamic drag torque acting on the satellite.
%   Detailed explanation goes here
if nargin < 5 % No accuracy specified
    acc = 100;
end

%% Determine velocity in BODY frame
v_B    = C_I2B*v_I; %C_I2B*(v_I + wx*C_I2B*r_I);

%% Determine which side in unshield
if v_B(1) > 0       % Moves in positive BODY x-direction -> drag on positive x-surface
    nx = [1 0 0]';
else                % Moves in negative BODY x-direction -> drag on negative x-surface
    nx = [-1 0 0];
end

if v_B(2) > 0       % Moves in positive BODY y-direction -> drag on positive y-surface (plate)
    ny = [0 1 0]';
else                % Moves in negative BODY y-direction -> drag on positive y-surface (plate and body)
    ny = [0 -1 0];
end

if v_B(3) > 0       % Moves in positive BODY z-direction -> drag on positive z-surface
    nz = [0 0 1]';
else                % Moves in negative BODY z-direction -> drag on negative z-surface
    nz = [0 0 -1]';
end

%% Calculate torque vectors
Tx = DRAG(nx, v_B, acc);  % Torque due to drag on x-surface
Ty = DRAG(ny, v_B, acc);  % Torque due to drag on y-surface
Tz = DRAG(nz, v_B, acc);  % Torque due to drag on z-surface

T = Tx + Ty + Tz;
end

