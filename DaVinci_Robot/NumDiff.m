function J = NumDiff( Robot , Robot_T )
%NUMDIFF Numerial Differentiation of the Direct Kinematics 
% Numerical Diff to calculate the Jacobian 

delta = 1e-6;  % Small perturbation
num_joints = size(Robot,1);      % number of robot joints
robot_q = symvar(Robot);         % names of robot coordinates

% Forward kinematics
o0 = Robot_T(1:3, 4);
R0 = Robot_T(1:3, 1:3);

% Initialize Jacobian
Jv = sym(zeros(3, num_joints));
Jomega = sym(zeros(3, num_joints));

for i = 1:num_joints
    q_perturb = robot_q;
    q_perturb(i) = robot_q(i) + delta;
    
    % Compute perturbed forward kinematics
    T_perturb = subs(Robot_T, robot_q, q_perturb);
    o_perturb = T_perturb(1:3, 4);
    R_perturb = T_perturb(1:3, 1:3);
    
    % Compute linear velocity (position difference / delta)
    Jv(:, i) = (o_perturb - o0) / delta;
    
    % Compute angular velocity (rotation matrix difference / delta)
    R_diff = R_perturb * R0.' - eye(3);
    Jomega(:, i) = [R_diff(3, 2); R_diff(1, 3); R_diff(2, 1)] / delta;

end

J = [Jv; Jomega];

end

