% CLEANUP AND INITIALIZATION

clear;
clc;
close all;

% Define the name of the full system matrices workspace file
file_name = 'DV_model.mat';

% Check for file existence before loading to prevent unhandled runtime errors
if exist(file_name, 'file')
    load(file_name);
    fprintf('=> daVinci full model successfully loaded from "%s"!\n', file_name);
else
    error('Error: File "%s" not found. Please run the "MainDyn" script first.', file_name);
end

% Assign state-space matrices from the loaded numeric models
A = A_num;
B = B_num;
C = C_num;
D = D_num;
E = E_num;

n_states  = size(A, 1);
n_inputs  = size(B, 2);
n_outputs = size(C, 1);

%% LQR AND REDUCED-ORDER OBSERVER SYNTHESIS

% --------- LQR Control Configuration ---------
err_q1_max  = 0.001;
err_q2_max  = 0.001;
err_q3_max  = 0.0001;
err_dq1_max = 0.5;
err_dq2_max = 0.5;
err_dq3_max = 0.5;

tau1_max = 5.3;  
tau2_max = 41; 
f3_max   = 420;

% State penalty matrix Q (Expanded 6x6 layout for full model tracking)
Q_lqr = diag([1/err_q1_max^2,  1/err_q2_max^2,  1/err_q3_max^2, ...
          1/err_dq1_max^2, 1/err_dq2_max^2, 1/err_dq3_max^2]);

% Actuator penalty matrix R (3x3 mapping for the 3 active motors)
R_lqr = diag([tau1_max*10, tau2_max*5, f3_max/100]);

[K_lqr, P_ss, poles_cl] = lqr(A, B, Q_lqr, R_lqr);

disp('========================================');
disp('   CONTROL PARAMETERS - LQR REGULATOR   ');
disp('========================================');
disp('Penalty matrix Q:'); disp(Q_lqr);
disp('Penalty matrix R:'); disp(R_lqr);
disp('Gain matrix K:');    disp(K_lqr);
disp('Closed-loop poles:'); disp(poles_cl);

% --------- Reduced-Order Observer Configuration ---------
% V_mat transformation isolating the 3 joint velocities at the bottom of the vector
V_mat = [zeros(n_states - n_outputs, n_outputs), eye(n_states - n_outputs)];
T_transf = [C; V_mat];
Ti = inv(T_transf);

M = Ti(:,1:n_outputs); % Dimension: 6 x 3
N = Ti(:,n_outputs+1:end); % Dimension: 6 x 3

% Re-assemble Friedland subspace partitions using full 6x6 system maps
A11 = C * A * M;   
A12 = C * A * N;   
A21 = V_mat * A * M;   
A22 = V_mat * A * N;   
B1  = C * B;       
B2  = V_mat * B;    

disp('===================================================');
disp('   ESTIMATION PARAMETERS - REDUCED-ORDER OBSERVER  ');
disp('===================================================');

Qe_red = 1000 * eye(n_states - n_outputs);
Re_red = eye(n_outputs);
J = lqr(A22', A12', Qe_red, Re_red).';

disp('Gain matrix J:');    disp(J);

F_red = A22 - J * A12;
G_red = A21 - J * A11 + F_red * J;
H_red = B2  - J * B1;
S_red = M + N * J;

fprintf('Reduced observer poles:\n'); disp(eig(F_red));


%% SIMULATION SCENARIOS LOOP

% Fixed geometric and central parameters
q3_bar_val = 0.2;
center_y   = 0; 
center_z   = -q3_bar_val; 
a2_val     = 0.4;

% Tracking performance cost function weights
Q_seg = Q;
R_seg = R;
Q1    = 1 * Q_seg;
t_f   = 5;

for scenario_id = ["A", "B", "C"]
    
    % 1. Dynamic Time Horizon and Resolution Allocation per Scenario
    if strcmp(scenario_id, "C")
        T_task = 10.0;           % Extended horizon for the 3D helix trajectory
        num_points = 2000;       % Higher temporal resolution
        freq_mov = 2 * pi * 0.5; % 0.5 Hz frequency mapping (1 loop every 2 seconds)
    else
        T_task = 4.0;            % Standard horizon for planar circles
        num_points = 1000;
        freq_mov = 2 * pi * (1 / T_task); % 1 loop every 4 seconds
    end
    
    % Generate baseline and extended time arrays to eliminate LQT edge-transient shocks
    t_sim = linspace(0, T_task, num_points); 
    dt = t_sim(2) - t_sim(1);
    
    T_extended = T_task + 0.5;         % Added 0.5 seconds padding to the backward horizon
    num_points_ext = num_points + 250; % Proportional extension resolution
    t_sim_ext = linspace(0, T_extended, num_points_ext);
    
    % Operating equilibrium point vector for the full 6-state configuration
    x_eq = [0; 0; 0.36; 0; 0; 0]; 
    u_eq = Getnvector(x_eq, 0);
    
    % 2. Cartesian Space Task Trajectory References Generation (vec_ref)
    switch scenario_id
        case "A"
            disp('==== Running Scenario A: YZ-Planar Circle (4s) ====');
            radius = 0.04; % 4 cm
            x_rf = a2_val * ones(size(t_sim)); 
            y_rf = center_y + radius * sin(freq_mov * t_sim);
            z_rf = center_z - radius + radius * cos(freq_mov * t_sim);
            
            x_rf_ext = a2_val * ones(size(t_sim_ext)); 
            y_rf_ext = center_y + radius * sin(freq_mov * t_sim_ext);
            z_rf_ext = center_z - radius + radius * cos(freq_mov * t_sim_ext);
            
        case "B"
            disp('==== Running Scenario B: XY-Planar Circle (4s) ====');
            radius = 0.04; % 4 cm
            x_rf = a2_val + radius * cos(freq_mov * t_sim); 
            y_rf = center_y + radius * sin(freq_mov * t_sim);
            z_rf = center_z * ones(size(t_sim)); 
            
            x_rf_ext = a2_val + radius * cos(freq_mov * t_sim_ext); 
            y_rf_ext = center_y + radius * sin(freq_mov * t_sim_ext);
            z_rf_ext = center_z * ones(size(t_sim_ext));
            
        case "C"
            disp('==== Running Scenario C: 3D Helical Trajectory (10s) ====');
            radius = 0.03; % 3 cm
            x_start = a2_val;      
            x_end = a2_val + 0.10; % 10 cm linear advance along the X axis
            
            x_rf = linspace(x_start, x_end, num_points);
            y_rf = center_y + radius * sin(freq_mov * t_sim);
            z_rf = center_z - radius + radius * cos(freq_mov * t_sim); 
            
            x_rf_ext = linspace(a2_val, a2_val + 0.10 * (T_extended/T_task), num_points_ext);
            y_rf_ext = center_y + radius * sin(freq_mov * t_sim_ext);
            z_rf_ext = center_z - radius + radius * cos(freq_mov * t_sim_ext);
    end
    vec_ref = [x_rf; y_rf; z_rf];
    vec_ref_ext = [x_rf_ext; y_rf_ext; z_rf_ext];
    
    % 3. Execute IDyn Mapping
    [q_ref, dq_ref, ddq_ref, tau_feedforward] = IDyn(num_points, vec_ref, dt);
    [q_ref_ext, dq_ref_ext, ~, ~] = IDyn(num_points_ext, vec_ref_ext, t_sim_ext(2)-t_sim_ext(1));
    
    % Construct comprehensive state trajectory references and time derivatives
    x_ref_total     = [q_ref; dq_ref];
    x_ref_ext_total = [q_ref_ext; dq_ref_ext];
    dx_ref_total    = [dq_ref; ddq_ref];
    
    % Pre-compute optimal linear Feedforward input profiles via model inversion
    delta_u_traj = zeros(3, num_points);
    for i = 1:num_points
        delta_u_traj(:, i) = pinv(B) * (dx_ref_total(:, i) - A * x_ref_total(:, i));
    end
    
    % 4. Pre-calculate 6x6 Finite-Horizon LQT Matrix Backward Integrations
    [t_P, P_flat_sol, t_eta, eta_sol] = compute_backward_lqt(A, B, Q_seg, R_seg, t_sim_ext, x_ref_ext_total, P_ss);
    
    % --- DISTURBANCE CONDITIONS ITERATION LOOP (CLEAN TRACKING VS RESPIRATORY NOISE) ---
    for pert_id = ["No Disturbance", "With Disturbance"]
        
        if strcmp(pert_id, "With Disturbance")
            w_amp = 0.005; % Inject active 5 mm sinusoidal breathing disturbance profile
        else
            w_amp = 0.0;   % Ideal disturbance-free closed-loop environment
        end
        
        % Dynamic initialization of the unmeasured velocity subspace (z0 = w_hat0 - J * C * x0)
        build_z0 = @(w_hat0, x0) w_hat0 - J * C * x0;
        
        x0_plant = x_ref_total(:, 1); 
        z0_init  = build_z0(zeros(3, 1), x0_plant);
        X0_global = [x0_plant; z0_init]; % Unified 9-state solver vector (6 real + 3 observer)
        
        options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
        
        % Solve Tracking Strategy 1: Finite-Horizon LQT Simulation
        [~, X_out1] = ode45(@(t, X) rhs_finite_lqt(t, X, A, B, C, E, F_red, G_red, H_red, S_red, N, R_seg, t_P, P_flat_sol, t_eta, eta_sol, t_sim, x_ref_total, w_amp), t_sim, X0_global, options);
        
        % Solve Tracking Strategy 2: Assumed Model / Stationary Feedforward Simulation
        [~, X_out2] = ode45(@(t, X) rhs_feedforward(t, X, A, B, C, E, F_red, G_red, H_red, S_red, N, K_lqr, t_sim, x_ref_total, delta_u_traj, w_amp), t_sim, X0_global, options);
        
        % Post-simulation algebraic reconstruction of real-time control efforts (u)
        u_seg1 = zeros(num_points, 3); 
        u_seg2 = zeros(num_points, 3); 
        
        Ri = inv(R_seg);
        
        for i = 1:num_points
            t = t_sim(i);
            
            % --- Reconstruction Profile 1: Finite-Horizon LQT Atuations ---
            x_lqt = X_out1(i, 1:6)';
            z_lqt = X_out1(i, 7:9)';
            xhat_lqt = S_red * (C * x_lqt) + N * z_lqt;
            
            P_flat = interp1(t_P, P_flat_sol, t, 'linear', 'extrap')';
            P = reshape(P_flat, 6, 6);
            eta = interp1(t_eta, eta_sol, t, 'linear', 'extrap')';
            
            u_seg1(i, :) = (-Ri * B' * P * xhat_lqt + Ri * B' * eta)';
            
            % --- Reconstruction Profile 2: Assumed Model Atuations ---
            x_ff = X_out2(i, 1:6)';
            z_ff = X_out2(i, 7:9)';
            xhat_ff = S_red * (C * x_ff) + N * z_ff;
            
            xr = x_ref_total(:, i);
            du = delta_u_traj(:, i);
            
            u_seg2(i, :) = (-K_lqr * (xhat_ff - xr) + du)';
        end
        
        % Plots and Animations
        PlotSimControl(t_sim, vec_ref, q_ref, X_out1, X_out2, u_seg1, u_seg2, scenario_id, pert_id);
        AnimSimControl(t_sim, X_out2, vec_ref, scenario_id, pert_id);
    end
end




% =====================================================================
%  LOCAL DYNAMIC HELPER FUNCTIONS
% =====================================================================

function dX = rhs_finite_lqt(t, X, A, B, C, E, F_red, G_red, H_red, S_red, N_obs, R_seg, t_P, P_flat_sol, t_eta, eta_sol, t_sim, x_ref_total, w_amp)
    % rhs_finite_lqt - Evaluates state derivative steps for the time-varying LQT tracker
    x = X(1:6);   % 6 real plant states
    z = X(7:9);   % 3 reduced-order estimator states
    y = C * x;
    
    xhat = S_red * y + N_obs * z;
    
    % Interpolate time-varying 6x6 Riccati kernels and co-state maps
    P_flat = interp1(t_P, P_flat_sol, t, 'linear', 'extrap')';
    P = reshape(P_flat, 6, 6);
    eta = interp1(t_eta, eta_sol, t, 'linear', 'extrap')';
    
    Ri = inv(R_seg);
    u = -Ri * B' * P * xhat + Ri * B' * eta;
    
    w = w_amp * sin(2 * pi * 0.5 * t);
    
    dx = A * x + B * u + E * w;
    dz = F_red * z + G_red * y + H_red * u;
    dX = [dx; dz];
end

function dX = rhs_feedforward(t, X, A, B, C, E, F_red, G_red, H_red, S_red, N_obs, K_lqr, t_sim, x_ref_total, delta_u_traj, w_amp)
    % rhs_feedforward - Evaluates state derivative steps for the Assumed Model controller
    x = X(1:6);   
    z = X(7:9);   
    y = C * x;
    
    xhat = S_red * y + N_obs * z;
    
    xr = interp1(t_sim, x_ref_total', t, 'linear', 'extrap')';
    du = interp1(t_sim, delta_u_traj', t, 'linear', 'extrap')';
    
    u = -K_lqr * (xhat - xr) + du;
    
    w = w_amp * sin(2 * pi * 0.5 * t);
    
    dx = A * x + B * u + E * w;
    dz = F_red * z + G_red * y + H_red * u;
    dX = [dx; dz];
end

function [t_P, P_flat_sol, t_eta, eta_sol] = compute_backward_lqt(A, B, Q, R, t_sim, x_ref_total, P_ss)
    % compute_backward_lqt - Dispatches backward temporal EDO solvers for 6x6 Riccati structures
    t_f = t_sim(end);
    Ri = inv(R); 
    BRiBt = B * Ri * B';
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    
    % Anchor final boundary conditions to the infinite-horizon steady-state profile
    P1_flat = P_ss(:);
    [t_P, P_flat_sol] = ode45(@(t, P_flat) riccati_rhs(t, P_flat, A, BRiBt, Q), [t_f, 0], P1_flat, options);
    t_P = flipud(t_P); 
    P_flat_sol = flipud(P_flat_sol);
    
    eta1 = P_ss * x_ref_total(:, end);
    [t_eta, eta_sol] = ode45(@(t, eta) rhs_interpolated_eta(t, eta, A, BRiBt, Q, t_P, P_flat_sol, t_sim, x_ref_total), [t_f, 0], eta1, options);
    t_eta = flipud(t_eta); 
    eta_sol = flipud(eta_sol);
end

function dP_flat = riccati_rhs(~, P_flat, A, BRiBt, Q)
    % riccati_rhs - Formulates differential step matrix operations for 6x6 Riccati arrays
    P = reshape(P_flat, 6, 6);
    dP = -P * A - A' * P + P * BRiBt * P - Q;
    dP_flat = dP(:);
end

function deta = rhs_interpolated_eta(t, eta, A, BRiBt, Q, t_P, P_flat_sol, t_sim, x_ref_total)
    % rhs_interpolated_eta - Formulates differential step mapping for tracking co-states
    P_flat = interp1(t_P, P_flat_sol, t, 'linear', 'extrap')';
    P = reshape(P_flat, 6, 6);
    xr = interp1(t_sim, x_ref_total', t, 'linear', 'extrap')';
    deta = -(A' - P * BRiBt) * eta - Q * xr;
end