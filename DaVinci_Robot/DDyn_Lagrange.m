function [ tau, B, phi, G ] = DDyn_Lagrange(Robot, M, CoM, I, g_vec)
% DDyn_Lagrange Calculates the robot dynamics (tau) using the
%               Lagrangian formulation (L = K - U).
%
%   Inputs:
%       Robot - Robot D-H table (symbolic)
%       M     - [n x 1] vector of masses
%       CoM   - [n x 3] matrix of CoM vectors 
%       I     - {n x 1} cell array of 3x3 inertia tensors
%       g_vec - [3 x 1] gravity vector 
%
%
%   Output:
%       tau   - [n x 1] symbolic vector of joint torques/forces
%       B     - [n x n] mass matrix
%       phi   - [n x 1] coriolis/centrifugal vector
%       G     - [n x 1] gravity vector

n = size(Robot, 1);

syms dq1 dq2 dq3 real
syms ddq1 ddq2 ddq3 real

dq = [dq1; dq2; dq3];
ddq = [ddq1; ddq2; ddq3];

% Create Symbolic q vector (positions)
q = sym(zeros(n, 1));
for i = 1:n
    p = Robot(i,:);
    if ~isempty(symvar(p(2))) % Revolute joint
        q(i) = p(2);
    else % Prismatic joint
        q(i) = p(1);
    end
end

% --- Calculate total kinetic (K) and potential (U) energy ---
% Iterate through each link, find its Jacobians (Jv, Jw) for its CoM, and sum its energies

K_total = sym(0);
U_total = sym(0);

% Recursive kinematics
T_prev = eye(4);
p_origins = sym(zeros(3, n+1)); % p_origins(:,j) = p_{j-1}^0
z_axes = sym(zeros(3, n+1));    % z_axes(:,j) = z_{j-1}^0
z_axes(:, 1) = [0; 0; 1];

% Cell arrays to store all T_i^0 and R_i^0
T_globals = cell(n, 1);
R_globals = cell(n, 1);

% Calculate all T_i^0 and R_i^0
for i = 1:n
    T_i_prev = DHTransf(Robot(i,:));
    T_globals{i} = T_prev * T_i_prev;
    R_globals{i} = T_globals{i}(1:3, 1:3);
    
    % Store origins and z-axes for Jacobian loop
    p_origins(:, i+1) = T_globals{i}(1:3, 4);
    z_axes(:, i+1) = T_globals{i}(1:3, 3);
    
    T_prev = T_globals{i};
end

% Calculate K_i and U_i for each link
for i = 1:n
    m_i = M(i);
    p_ci_local = CoM(i, :)'; % r_{i,Ci}^i
    I_ci_local = I{i};       % I_Ci^i
    
    R_i_0 = R_globals{i};
    p_i_0 = p_origins(:, i+1);
    
    % Global CoM position
    p_ci_0 = p_i_0 + R_i_0 * p_ci_local;
    
    % Jacobians for the CoM
    Jv_ci = sym(zeros(3, n)); % Jv for CoM_i
    Jw_i  = sym(zeros(3, n)); % Jw for CoM_i (same as for origin_i)
    
    for j = 1:i
        p_j_minus_1 = p_origins(:, j);
        z_j_minus_1 = z_axes(:, j);
        
        p_j_link_type = Robot(j,:);
        
        if ~isempty(symvar(p_j_link_type(2))) % Revolute joint j
            Jw_i(:, j) = z_j_minus_1;
            Jv_ci(:, j) = cross(z_j_minus_1, p_ci_0 - p_j_minus_1);
        else % Prismatic joint j
            Jw_i(:, j) = sym(0);
            Jv_ci(:, j) = z_j_minus_1;
        end
    end
    
    % Calculate velocities
    v_ci = Jv_ci * dq;
    w_i  = Jw_i * dq;
    
    % Calculate K_i and U_i
    I_i_0 = R_i_0 * I_ci_local * R_i_0';
    
    K_i = 1/2 * m_i * (v_ci' * v_ci) + 1/2 * w_i' * I_i_0 * w_i;
    U_i = -m_i * g_vec' * p_ci_0;
    
    K_total = K_total + K_i;
    U_total = U_total + U_i;
end

K_total = simplify(K_total);
U_total = simplify(U_total);

% --- Calculate G (gravity vector) ---
G = simplify(jacobian(U_total, q)');

% Calculate B (mass matrix) ---
% B(i,j) = d^2(K) / (d(dq_i) * d(dq_j))
B = simplify(hessian(K_total, dq));

% --- Calculate phi (coriolis/centrifugal vector) ---
% phi = C(q,dq)*dq, derived from Christoffel symbols
% C_ijk = 1/2 * (dB_ij/dq_k + dB_ik/dq_j - dB_jk/dq_i)
phi = sym(zeros(n, 1));
for i = 1:n
    sum_C_ijk_dq_j_dq_k = sym(0);
    for j = 1:n
        for k = 1:n
            % Calculate Christoffel symbol C_ijk
            dB_ij_dq_k = diff(B(i,j), q(k));
            dB_ik_dq_j = diff(B(i,k), q(j));
            dB_jk_dq_i = diff(B(j,k), q(i));
            
            C_ijk = 1/2 * (dB_ij_dq_k + dB_ik_dq_j - dB_jk_dq_i);
            
            sum_C_ijk_dq_j_dq_k = sum_C_ijk_dq_j_dq_k + C_ijk * dq(j) * dq(k);
        end
    end
    phi(i) = sum_C_ijk_dq_j_dq_k;
end
phi = simplify(phi);


tau = B*ddq + phi + G;

tau = simplify(tau);


end