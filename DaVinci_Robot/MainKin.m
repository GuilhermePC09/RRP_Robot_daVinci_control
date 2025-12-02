%% Generate Matlab and Simulink Code for Robot
%  Create Simulink Library file Robot_Lib.mdl
% new_system('Robot_daVinci_Lib','Library');
open_system('Robot_daVinci_Lib');

%% Cleaning
clear;
clc;

%% Direct Kinematics

[ Robot ] = Robot_dV;
[RobotdV_T] = DKin(Robot);
RobotdV_R = RobotdV_T(1:3,1:3);
RobotdV_p = RobotdV_T(1:3,4);

disp('Rotation Matrix (R):')
pretty(RobotdV_R)
disp('Position Matrix (p):')
pretty(RobotdV_p)

% Value of a2
a2_val = 0.1;

%% Inverse Kinematics

IK_joints = IKin();
disp('Joint Configurations:')
pretty(IK_joints)

%% Generate optimized embeded Matlab function blocks for Simulink

syms a2 real
RobotdV_R_subs = subs(RobotdV_R, a2, a2_val);
RobotdV_p_subs = subs(RobotdV_p, a2, a2_val);
IK_joints_subs = subs(IK_joints, a2, a2_val);

matlabFunctionBlock('Robot_daVinci_Lib/Robot_daVinci_Direct_Kinematics',RobotdV_R_subs, RobotdV_p_subs);
matlabFunctionBlock('Robot_daVinci_Lib/Robot_daVinci_Inverse_Kinematics',IK_joints_subs);


%% Geometric Jacobian
J = GeoJac(Robot_dV);

J_subs = subs(J, a2, a2_val);

% Generate Matlab function block 
matlabFunctionBlock('Robot_daVinci_Lib/Robot_daVinci_Geometric_Jacobian',J_subs);

%% Numerical Differentiation and Jacobian Validation

J_num = NumDiff(Robot_dV, RobotdV_T);

J_num_subs = subs(J_num, a2, a2_val);

% Generate Matlab function block 
matlabFunctionBlock('Robot_daVinci_Lib/Robot_daVinci_NumericalDiff_Jacobian',J_num_subs);


%% Singularities analysis

% Position Jacobian
J_pos = J(1:3, :);

detJ = det(J_pos);
detJ_simple = simplify(detJ);

disp('Jacobian determinant (det(J)):');
pretty(detJ_simple);

%% Plot Robot

PlotRobot_dV([0; 0; 0.1], a2_val);

% q_start = [-pi/2; -pi/4; 0.05];
% q_end = [pi/2; pi/4; 0.2];
% 
% steps = 100;
% for i = 1:steps
%     % Interpolate joints
%     q_current = q_start + (q_end - q_start) * (i/steps);
% 
%     % Call the all-in-one plotting function
%     PlotRobot_dV(q_current, a2_val);
% 
%     drawnow; % Force MATLAB to update the plot
%     pause(0.005); % Pause to make the animation visible
% end


%% Save library in current Directory
save_system('Robot_daVinci_Lib');
% close_system('Robot_daVinci_Lib');

