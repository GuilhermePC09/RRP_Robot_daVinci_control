function [A_num, B_num, C_num, D_num, sym_list, num_list] = NumSS(A_lin, B_lin, C_lin, D_lin)
%NUMSS Substitutes numeric values on the state space model.
%
% Inputs:
%   A_lin, B_lin, C_lin, D_lin - The four symbolic linearized state-space matrices
%
% Outputs:
%   A_num, B_num, C_num, D_num - The four numeric linearized state-space matrices
%

syms a2 q3_bar g m1 m2 m3 lc2 lc3 bv1 bv2 bv3 K_env3
syms Ix1 Iy1 Iz1 Ix2 Iy2 Iz2 Ix3 Iy3 Iz3

% System Parameters
g_val = 9.81;
a2_val = 0.4; % 0.4
q3_bar_val = 0.2; % 0.2

% Link 1 (Cylinder)
m1_val = 4.3175;
r1_val = 0.05;
h1_val = 0.2;
% (Assuming CoM(1,:) = [0 0 0], so lc1 = 0)
Ix1_val = (1/12)*m1_val*(3*r1_val^2 + h1_val^2); 
Iy1_val = (1/2)*m1_val*r1_val^2;
Iz1_val = Ix1_val;

% Link 2 (Cylinder)
m2_val = 8.635;
r2_val = 0.05;
h2_val = 0.4;
lc2_val = h2_val / 2; % CoM at center (0.2m)
Ix2_val = (1/2)*m2_val*r2_val^2;
Iy2_val = (1/12)*m2_val*(3*r2_val^2 + h2_val^2);
Iz2_val = Iy2_val;

% Link 3 (Block)
m3_val = 1.32;
L3_val = 0.2; % x-axis
h3_val = 0.04; % y-axis
b3_val = 0.06; % z-axis
lc3_val = L3_val / 2; % CoM at center (0.1m)
Ix3_val = (1/12)*m3_val*(h3_val^2 + b3_val^2); 
Iy3_val = (1/12)*m3_val*(L3_val^2 + b3_val^2);
Iz3_val = (1/12)*m3_val*(L3_val^2 + h3_val^2);

% Damping (Viscosity) Parameters (Assumed small values)
bv1_val = 2; 
bv2_val = 2;
bv3_val = 5;

% Stiffness Parameter 
K_env3_val = 500;

% Substitution Lists
sym_list = [a2, q3_bar, g, m1, m2, m3, lc2, lc3, ...
            Ix1, Iy1, Iz1, Ix2, Iy2, Iz2, Ix3, Iy3, Iz3, ...
            bv1, bv2, bv3, K_env3];
        
num_list = [a2_val, q3_bar_val, g_val, m1_val, m2_val, m3_val, lc2_val, lc3_val, ...
            Ix1_val, Iy1_val, Iz1_val, Ix2_val, Iy2_val, Iz2_val, ...
            Ix3_val, Iy3_val, Iz3_val, ...
            bv1_val, bv2_val, bv3_val, K_env3_val];


% disp('Substituting symbolic matrices to numeric...');
A_num = double(subs(A_lin, sym_list, num_list));
B_num = double(subs(B_lin, sym_list, num_list));
C_num = double(subs(C_lin, sym_list, num_list));
D_num = double(subs(D_lin, sym_list, num_list));

end

