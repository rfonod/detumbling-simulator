function [C] = rot_mat(theta,axis)
% Creates rotation matrix for a rotation of theta around axis

ct = cos(theta);
st = sin(theta);

switch axis
    case 1
        C = [1  0  0 ;
             0  ct st;
             0 -st ct;];
    case 2
        C = [ct 0 -st;
             0  1  0 ;
             st 0  ct;];
    case 3
        C = [ ct st 0;
             -st ct 0;
              0  0  1;];  
    otherwise
        warning('Only 1, 2 or 3 possible for axis');
        C = eye(3);
end


end