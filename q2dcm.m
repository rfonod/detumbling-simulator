function DCM = q2dcm(q)
% Quaternion 2 DCM

q12=q(1)^2;
q22=q(2)^2;
q32=q(3)^2;

q1q2 = q(1)*q(2);
q1q3 = q(1)*q(3);
q2q4 = q(2)*q(4);
q1q4 = q(1)*q(4);
q2q3 = q(2)*q(3);
q3q4 = q(3)*q(4);

DCM = [1-2*(q22+q32),2*(q1q2+q3q4),2*(q1q3-q2q4)
    2*(q1q2-q3q4),1-2*(q12+q32),2*(q2q3+q1q4)
    2*(q1q3+q2q4),2*(q2q3-q1q4),1-2*(q12+q22)];

end