%  CLEANUP AND INITIALIZATION

clear;
clc;
close all;

file_name = 'DV_reduced_model.mat';

% Check for file existence before loading
if exist(file_name, 'file')
    load(file_name);
    fprintf('=> Model successfully loaded from "%s"!\n', file_name);
else
    error('Error: File "%s" not found. Please run the "MainDyn" script first.', file_name);
end

%% 0) SYSTEM ANALYSIS

% Reduced system
A = A_reduced;
B = B_reduced;
C = C_reduced;
D = D_reduced;
E = E_reduced;

n_states = size(A, 1);
n_inputs = size(B, 2);
n_outputs = size(C, 1);

% ControlAnalysis(A, B, C, D);

%% 1) LQR CONTROL

% Maximum tolerable limits (Bryson's approach)
err_q2_max = 0.001;  % 0.001 rad of tolerable error
err_q3_max = 0.0001;  % 0.0001 m (0.1 mm) of tolerable error

err_dq1_max = 0.5; % 0.5 rad/s tolerable
err_dq2_max = 0.5;
err_dq3_max = 0.5;

% Actuator max effort 
tau1_max = 5.3;  
tau2_max = 41; 
f3_max   = 420;

% Q matrix (5x5 for the reduced system)
Q_lqr = diag([1/err_q2_max^2, 1/err_q3_max^2, ...
              1/err_dq1_max^2, 1/err_dq2_max^2, 1/err_dq3_max^2]);

% R matrix (3x3 for the 3 actuators)
R_lqr = diag([tau1_max*10, tau2_max*5, f3_max/100]);

disp('========================================');
disp('   CONTROL PARAMETERS - LQR REGULATOR   ');
disp('========================================');
disp('Penalty matrix Q:');
disp(Q_lqr);
disp('Penalty matrix R:');
disp(R_lqr);

% Synthesize the optimal LQR gain matrix
[K_lqr, ~, poles_cl] = lqr(A, B, Q_lqr, R_lqr);

disp('Gain matrix K:');
disp(K_lqr);
disp('Closed-loop poles:');
disp(poles_cl);

% disp('Processing closed-loop analyses for LQR');
% AnalyzeCL(A, B, C, D, E, K_lqr, 'LQR');

%% 2) POLE PLACEMENT CONTROL

% Strategy based on LQR poles
% ctrl_poles = poles_cl';
ctrl_poles = [-5, -7 + 4i, -7 - 4i, -40 + 20i, -40 - 20i];

disp('=======================================');
disp('  CONTROL PARAMETERS - PP REGULATOR    ');
disp('=======================================');
disp('Desired closed-loop poles:');
disp(ctrl_poles');

% Calculation of the gain matrix using the Place algorithm
F_pp = place(A, B, ctrl_poles);
% F_pp = F_pp/1.2;

disp('Calculated gain matrix F:');
disp(F_pp);

% disp('Processing closed-loop analyses for Pole Placement');
% AnalyzeCL(A, B, C, D, E, F_pp, 'PP');

%% 3) REGULATORS COMPARATIVE ANALYSIS

disp('Processing comparative analysis between LQR and Pole Placement...');
CompReg(A, B, C, D, K_lqr, F_pp);

%% 4) FULL-ORDER OBSERVER - IDENTITY
% LQR dual approach (Kalman-Bucy filter tuning structure)
Wc = 1000 * eye(n_states); % "Process noise" covariance (increase = faster observer dynamics)
Vc = eye(n_outputs);       % "Measurement noise" covariance (increase = slower observer dynamics)
L = lqr(A', C', Wc, Vc).';

% Pole placement approach alternative
% poles_obs_id = [-30+8i, -30-8i, -20, -25+2i, -25-2i];
% L = place(A', C', poles_obs_id)';

disp('================================================');
disp('   ESTIMATION PARAMETERS - IDENTITY OBSERVER    ');
disp('================================================');
fprintf('Observer poles (eigenvalues of (A - L*C)):\n'); 
disp(eig(A - L * C));
fprintf('Observer gain matrix (L):\n'); 
disp(L);

% Isolated observer system formulation (inputs = [u; y], output = x_hat)
A_obs = A - L * C;
B_obs = [B, L];
C_obs = eye(n_states);
D_obs = zeros(n_states, n_inputs + n_outputs);
sys_obs = ss(A_obs, B_obs, C_obs, D_obs);

%% 5) REDUCED-ORDER OBSERVER

% Friedland state transformation and subspace partitioning
V_mat = [zeros(n_states - n_outputs, n_outputs), eye(n_states - n_outputs)];
T     = [C; V_mat];

if abs(det(T)) < 1e-10
    error('Error: Transformation matrix T is (nearly) singular. Select an alternative V_mat partition.');
end

Ti = inv(T);
M  = Ti(:, 1:n_outputs); % Dimension: n x m
N  = Ti(:, n_outputs+1:end); % Dimension: n x (n-m)

% Subspace block partitioning
A11 = C * A * M; % m x m
A12 = C * A * N;   % m x (n-m)
A21 = V_mat * A * M; % (n-m) x m
A22 = V_mat * A * N; % (n-m) x (n-m)
B1  = C * B; % m x n_inputs
B2  = V_mat * B; % (n-m) x n_inputs

disp('===================================================');
disp('   ESTIMATION PARAMETERS - REDUCED-ORDER OBSERVER  ');
disp('===================================================');

% Observability and numerical conditioning validation of the reduced pair
rank_obs_red = rank(obsv(A22, A12));
fprintf('Reduced observability rank = %d (Required n-m = %d)\n', rank_obs_red, n_states - n_outputs);

if rank_obs_red < (n_states - n_outputs)
    error('Error: Reduced-order pair (A22, A12) is NOT observable. Synthesis cannot proceed.');
end
fprintf('Condition number cond(obsv(A22,A12)) = %.3e\n\n', cond(obsv(A22, A12)));

% Observer gain (J) synthesis via LQE dual approach
Qe_red = 1000 * eye(n_states - n_outputs); % Process noise covariance (increase = faster estimation)
Re_red = eye(n_outputs); % Measurement noise covariance (increase = slower estimation)
J = lqr(A22', A12', Qe_red, Re_red).';

% Alternative Pole Placement approach
% J = place(A22', A12', 3*ctrl_poles).';

% Compute reduced-order observer reconstruction matrices
F_red = A22 - J * A12;
G_red = A21 - J * A11 + F_red * J;
H_red = B2  - J * B1;
S_red = M   + N * J;

fprintf('Reduced observer poles (eigenvalues of F):\n'); disp(eig(F_red));
fprintf('Gain matrix norm ||J|| = %.3e\n\n', norm(J));

% Open-loop augmented state-space assembly [x; z]
A_aug_red = [A,         zeros(n_states, n_states - n_outputs);
             G_red * C, F_red];
B_aug_red = [B; H_red];
E_aug_red = [E; zeros(n_states - n_outputs, size(E, 2))];
C_aug_red = [C, zeros(n_outputs, n_states - n_outputs)];

% Mapping matrix to extract x_hat from the augmented state: x_hat = S_red*C*x + N*z
recover_xhat = [S_red * C, N];

% --- Dynamic Observer Initialization Function ---
% Maps initial velocity estimates (w_hat0) and initial physical states (x0) 
% to the unmeasured state estimator space (z0 = w_hat0 - J * C * x0)
build_z0 = @(w_hat0, x0) w_hat0 - J * C * x0;

%% 6) OBSERVERS COMPARATIVE ANALYSIS

% Define initial conditions for convergence validation
t_sim = 0 : 1e-3 : 2;
x0_real = [0.05; 0.02; 1; 3; 2];
x0_hat  = zeros(n_states, 1);

obs_red.F_red = F_red;
obs_red.A_aug_red = A_aug_red;
obs_red.recover_xhat = recover_xhat;
obs_red.J = J;

state_labels = {'Joint q_2', 'Joint q_3', ...
                'Velocity dq_1', 'Velocity dq_2', 'Velocity dq_3'};

disp('Processing comparative analysis between Identity and Reduced-Order observers...');
CompObsv(A, C, L, obs_red, poles_cl, x0_real, t_sim, state_labels, x0_hat);


%% 7) CLOSED-LOOP SIMULATION (LQR + REDUCED OBSERVER)

% Real initial state vector of the plant
x0_real = [0.05; 0.02; 1; 3; 2]; 

% Run the complete closed-loop simulation and generate all performance plots
disp('Running full closed-loop simulation (LQR + Reduced-Order Observer)...');
SimCLTotal(A, B, C, E, K_lqr, F_red, G_red, H_red, S_red, N, x0_real, t_sim);


%% 8) FINITE-HORIZON LQ TRACKER (LQT) WITH REDUCED OBSERVER
% Cost function weights and time horizon configurations
Q_seg = Q_lqr;
R_seg = R_lqr;
Q1 = 1 * Q_seg;
t_f = 5;

t_span_prog = 0 : 0.001 : t_f; 
t_span_back = [t_f, 0];    % Time interval for backward integration

% Pre-compute input coupling matrices
Ri = inv(R_seg);
BRiBt = B * Ri * B';       

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% 1. Backward Integration: Differential Riccati Equation Matrix P(t)
disp('=> Solving backward differential Riccati equation...');
P1_flat = Q1(:); % Flatten terminal 5x5 matrix into a 25x1 vector
[t_P, P_flat_sol] = ode45(@(t, P_flat) riccati_rhs(t, P_flat, A, BRiBt, Q_seg), t_span_back, P1_flat, options);
t_P = flipud(t_P); P_flat_sol = flipud(P_flat_sol);

% 2. Backward Integration: Co-State Vector eta(t)
disp('=> Solving backward co-state differential equation...');
eta1 = Q1 * daVinci_ref(t_f); % Terminal co-state condition
[t_eta, eta_sol] = ode45(@(t, eta) eta_rhs(t, eta, A, BRiBt, Q_seg, t_P, P_flat_sol, n_states), t_span_back, eta1, options);
t_eta = flipud(t_eta); eta_sol = flipud(eta_sol);

% 3. Forward Closed-Loop Simulation (Physical Plant + Reduced Observer)
disp('=> Launching forward tracking simulation with estimated state feedback...');
x0_real = zeros(n_states, 1);
z0_inicial = build_z0(zeros(n_states - n_outputs, 1), x0_real);
X0_global = [x0_real; z0_inicial]; 

[t_sim, X_total] = ode45(@(t, X) forward_sim_rhs(t, X, A, B, C, E, F_red, G_red, H_red, S_red, N, Ri, t_P, P_flat_sol, t_eta, eta_sol, n_states), t_span_prog, X0_global, options);


% Plot and Statistics
P_flat_ss = P_flat_sol(1, :); 
P_ss = reshape(P_flat_ss, n_states, n_states);
poles_cl_lqt = eig(A - B * Ri * B' * P_ss);
poles_obs_red = eig(F_red);

PlotLQTResults(t_sim, X_total, A, B, C, S_red, N, Ri, t_P, P_flat_sol, t_eta, eta_sol, poles_cl_lqt, poles_obs_red);


% =====================================================================
%  LOCAL DYNAMIC HELPER FUNCTIONS
% =====================================================================

function xr = daVinci_ref(t)
    % daVinci_ref - Generates smooth exponential tracking references
    alpha_regime = 10;
    q2_f = 5 * pi / 180; % Target of 5 degrees converted to radians
    q3_f = 0.05;         % Target of 5 centimeters converted to meters
    
    % Trajectory position profiles
    q2_r = q2_f * (1 - exp(-alpha_regime * t));
    q3_r = q3_f * (1 - exp(-alpha_regime * t));
    
    % Trajectory velocity profiles (First derivatives)
    dq2_r = q2_f * alpha_regime * exp(-alpha_regime * t);
    dq3_r = q3_f * alpha_regime * exp(-alpha_regime * t);
    
    xr = [q2_r; q3_r; 0; dq2_r; dq3_r]; 
end

function dP_flat = riccati_rhs(~, P_flat, A, BRiBt, Q)
    % riccati_rhs - Evaluates the differential continuous Riccati matrix step
    n = sqrt(length(P_flat));
    P = reshape(P_flat, n, n);
    dP = -P * A - A' * P + P * BRiBt * P - Q;
    dP_flat = dP(:);
end

function deta = eta_rhs(t, eta, A, BRiBt, Q, t_P, P_flat_sol, n_states)
    % eta_rhs - Evaluates the time-varying backward tracking co-state vector step
    P_flat = interp1(t_P, P_flat_sol, t, 'linear', 'extrap')';
    P = reshape(P_flat, n_states, n_states);
    deta = -(A' - P * BRiBt) * eta - Q * daVinci_ref(t);
end

function dX = forward_sim_rhs(t, X, A, B, C, E, F_red, G_red, H_red, S_red, N, Ri, t_P, P_flat_sol, t_eta, eta_sol, n_states)
    % forward_sim_rhs - Couples the physical plant and the observer dynamics forward in time
    x = X(1:n_states);
    z = X(n_states+1:end);
    y = C * x;
    
    % Reconstruct current state estimation vector
    x_hat = S_red * y + N * z;
    
    % Recover time-varying matrices along the integrated optimal trajectory
    P_flat = interp1(t_P, P_flat_sol, t, 'linear', 'extrap')';
    P = reshape(P_flat, n_states, n_states);
    eta = interp1(t_eta, eta_sol, t, 'linear', 'extrap')';
    
    % Compute time-varying controller parameters
    K_t = Ri * B' * P;
    u_tilde = Ri * B' * eta;
    u = -K_t * x_hat + u_tilde;
    
    % Inject respiratory simulation disturbance profile through matrix E
    w = 0.005 * sin(2 * pi * 0.5 * t);
    
    % Compute differential vector steps
    dx = A * x + B * u + E * w;
    dz = F_red * z + G_red * y + H_red * u;
    
    dX = [dx; dz];
end


%% 9) ASSUMED MODEL TRACKER
% Stabilization decay parameter (5 / alpha = 0.5s settling time target)
alpha_regime = 10; 

% Construct the autonomous reference generator state matrix A_r (5x5)
A_r = [ 0,  0,  0,   1,   0;   % d(q2_r)/dt  = dq2_r
        0,  0,  0,   0,   1;   % d(q3_r)/dt  = dq3_r
        0,  0,  0,   0,   0;   % d(dq1_r)/dt = 0
        0,  0,  0, -alpha_regime,  0;   % d(dq2_r)/dt
        0,  0,  0,   0, -alpha_regime]; % d(dq3_r)/dt

% Feedforward gain (G_r) synthesis based on steady-state model inversion
M_sel = C; 
A_cl_inv = inv(A - B * K_lqr);
N_ff = pinv(M_sel * A_cl_inv * B) * M_sel * A_cl_inv;
G_r = N_ff * (A - A_r);

% Assemble the global 13-state augmented system: [x (5); z (3); x_r (5)]
A_ex = [ (A - B * K_lqr * S_red * C),            (-B * K_lqr * N),         (B * (K_lqr - G_r));
         (G_red * C - H_red * K_lqr * S_red * C), (F_red - H_red * K_lqr * N), (H_red * (K_lqr - G_r));
         zeros(5, 5),                             zeros(5, 3),              A_r ];

% Map respiratory disturbance input solely to the physical plant channel
E_ex = [ E; 
         zeros(n_states - n_outputs, size(E, 2));
         zeros(5, size(E, 2)) ];

sys_ex = ss(A_ex, E_ex, eye(13), zeros(13, 1));

% Simulation time configurations
t_sim = 0 : 0.001 : 5; 

% Initial conditions vector configuration
x0_plant = zeros(n_states, 1); 
z0_init  = build_z0(zeros(3, 1), x0_plant);

% Calculate reference generator shooting velocities for a 5-deg and 5-cm target
v2_0 = alpha_regime * (5 * pi / 180); 
v3_0 = alpha_regime * 0.05;          
x0_reference = [0; 0; 0; v2_0; v3_0]; 

% Unified 13x1 initial state vector
X0_global = [x0_plant; z0_init; x0_reference];

% Execute the linear time simulation under continuous respiratory disturbance
w_breathing = 0.005 * sin(2 * pi * 0.5 * t_sim);
[X_out, ~] = lsim(sys_ex, w_breathing, t_sim, X0_global);


% Plot and Statistics
poles_obs_red = eig(F_red);

PlotAMResults(t_sim, X_out, A, B, C, S_red, N, K_lqr, G_r, poles_cl, poles_obs_red);