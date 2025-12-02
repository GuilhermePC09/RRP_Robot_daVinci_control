function [ IK_joints ] = IKin()
%IKIN Robot Inverse Kinematics
%   EndEff_p = pe (3x1)

syms Pe_x Pe_y Pe_z a2 real

% Solving for q1
q1 = atan2(Pe_y, Pe_x);

% Solving for q2 and q3
q3_1 = sqrt(Pe_x^2 + Pe_y^2 + Pe_z^2 - a2^2);
q3_2 = - sqrt(Pe_x^2 + Pe_y^2 + Pe_z^2 - a2^2);

alpha = atan2(Pe_z, sqrt(Pe_x^2 + Pe_y^2));
beta = atan2(q3_1, a2);

q2_1 = alpha + beta;
q2_2 = alpha - beta;

IK_joints = [q1, q2_1, q3_1; q1, q2_2, q3_2];

end