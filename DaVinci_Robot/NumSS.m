function [A_num, B_num, C_num, D_num, E_num, sym_list, num_list] = NumSS(A_lin, B_lin, C_lin, D_lin, E_lin)
%NUMSS Substitutes numeric values on the state space model.
%
% Inputs:
%   A_lin, B_lin, C_lin, D_lin, E_lin - The symbolic linearized state-space matrices
%
% Outputs:
%   A_num, B_num, C_num, D_num, E_num - The numeric linearized state-space matrices
%

syms a2 q3_bar g m1 m2 m3 lc2 lc3 bv1 bv2 bv3 K_env3
syms Ix [3, 1] real
syms Ixy [3, 1] real
syms Ixz [3, 1] real
syms Iy [3, 1] real
syms Iyz [3, 1] real
syms Iz [3, 1] real

% Functions for coordinates rotation
A_rad = @(deg) deg * pi / 180;
R_x = @(deg) [1, 0, 0; 0, cos(A_rad(deg)), -sin(A_rad(deg)); 0, sin(A_rad(deg)), cos(A_rad(deg))];
R_y = @(deg) [cos(A_rad(deg)), 0, sin(A_rad(deg)); 0, 1, 0; -sin(A_rad(deg)), 0, cos(A_rad(deg))];
R_z = @(deg) [cos(A_rad(deg)), -sin(A_rad(deg)), 0; sin(A_rad(deg)), cos(A_rad(deg)), 0; 0, 0, 1];


% System Parameters
g_val = 9.81;
a2_val = 0.45;
q3_bar_val = 0.36; %0.36


% LINK 1
m1_val = 3.679; % kg

I1_CAD = [0.022, 0.00,  0.00;
          0.00,  0.022, 0.00;
          0.00,  0.00,  0.013]; % kg.m^2, CM


% 90 deg rotation on X
R1 = R_x(90);
I1_DH = R1' * I1_CAD * R1;
J_ref1 = 0.024; 
I1_DH = I1_DH + diag([0, J_ref1, 0]);

Ix1_val  = I1_DH(1,1);
Ixy1_val = I1_DH(1,2);
Ixz1_val = I1_DH(1,3);
Iy1_val  = I1_DH(2,2);
Iyz1_val = I1_DH(2,3);
Iz1_val  = I1_DH(3,3);


% LINK 2
m2_val = 8.056; % kg
lc2_val = 0.35; % m


I2_CAD = [ 0.11,      -3.889e-7, -3.165e-6;
          -3.889e-7,   0.095,    -0.022;
          -3.165e-6,  -0.022,     0.04 ]; % kg.m^2, CM

% 90 deg rotation on Y and 180 on new X 
R2 = R_y(90) * R_x(180); 
I2_DH = R2' * I2_CAD * R2;
J_ref2 = 0.169; 
I2_DH = I2_DH + diag([0, J_ref2, 0]);

Ix2_val  = I2_DH(1,1);
Ixy2_val = I2_DH(1,2);
Ixz2_val = I2_DH(1,3);
Iy2_val  = I2_DH(2,2);
Iyz2_val = I2_DH(2,3);
Iz2_val  = I2_DH(3,3);


% LINK 3

m3_val = 1.949; % kg
lc3_val = 0.3; % m

I3_CAD = [ 0.003,      1.742e-8,  0.00;
           1.742e-8,   0.024,    -1.681e-9;
           0.00,      -1.681e-9,  0.022 ]; % kg.m^2, CM

% Same as link 2
R3 = R2;
I3_DH = R3' * I3_CAD * R3;

Ix3_val  = I3_DH(1,1);
Ixy3_val = I3_DH(1,2);
Ixz3_val = I3_DH(1,3);
Iy3_val  = I3_DH(2,2);
Iyz3_val = I3_DH(2,3);
Iz3_val  = I3_DH(3,3);




% Damping (Viscosity) Parameters (Assumed small values)
bv1_val = 2; 
bv2_val = 2;
bv3_val = 5;

% Stiffness Parameter 
K_env3_val = 500;

% Substitution Lists
sym_list = [a2, q3_bar, g, m1, m2, m3, lc2, lc3, ...
            Ix1, Ixy1, Ixz1, Iy1, Iyz1, Iz1, ...
            Ix2, Ixy2, Ixz2, Iy2, Iyz2, Iz2, ...
            Ix3, Ixy3, Ixz3, Iy3, Iyz3, Iz3, ...
            bv1, bv2, bv3, K_env3];


        
num_list = [a2_val, q3_bar_val, g_val, m1_val, m2_val, m3_val, lc2_val, lc3_val, ...
            Ix1_val, Ixy1_val, Ixz1_val, Iy1_val, Iyz1_val, Iz1_val, ...
            Ix2_val, Ixy2_val, Ixz2_val, Iy2_val, Iyz2_val, Iz2_val, ...
            Ix3_val, Ixy3_val, Ixz3_val, Iy3_val, Iyz3_val, Iz3_val, ...
            bv1_val, bv2_val, bv3_val, K_env3_val];


% disp('Substituting symbolic matrices to numeric...');
A_num = double(subs(A_lin, sym_list, num_list));
B_num = double(subs(B_lin, sym_list, num_list));
C_num = double(subs(C_lin, sym_list, num_list));
D_num = double(subs(D_lin, sym_list, num_list));
E_num = double(subs(E_lin, sym_list, num_list));

end

