%% ----- Cleaning -----
clear;
clc;
close all;

%% ----- Direct Dynamics (Lagrange) -----

[Robot, M, CoM, I] = Robot_dV();

syms g real;
g_vec = [0; 0; -g];

[ tau, B, phi, G ] = DDyn_Lagrange(Robot, M, CoM, I, g_vec);
n = phi + G;

% disp('Mass Matrix (B):')
% pretty(B)
% disp('Coriolis Vector (phi):')
% pretty(phi)
% disp('Gravity Vector (G):')
% pretty(G)


%% ----- Non-linear State Space Model -----

% Adding viscosity on each joint
syms dq1 dq2 dq3 real
dq_sym = [dq1; dq2; dq3];

syms bv1 bv2 bv3 real
Bv = diag([bv1, bv2, bv3]);

tau_friction = Bv * dq_sym;

% Adding an environment "stiffness" (spring)
syms q1 q2 q3 real
q_sym = [q1; q2; q3];

syms K_env3 real
K_env_matrix = diag([0, 0, K_env3]);

syms q3_bar real
q_bar_env = [0; 0; q3_bar];

% Adding pysiological perturbation
syms w real
w_env = [0; 0; w];

tau_env = K_env_matrix * (q_sym - (q_bar_env + w_env));


% New n = phi + G + tau_b + tau_spring
n = n + tau_friction + tau_env;

% disp('Total n vector (phi + G + Bv*dq + tau_ext):');
% pretty(n)

[f_ss, x_ss, u_ss, w_ss] = StateSpaceFunc(B, n);

% disp('State Vector (x):');
% pretty(x_ss)
% disp('Input Vector (u):');
% pretty(u_ss)
% disp('Perturbation Variable (w):');
% pretty(w_ss)
% disp('State Function f(x,u) = x_dot:');
% pretty(f_ss)

% Torque needed to mantain equilibrium
% tau_eq = subs(n, [q_sym; dq_sym], [0; 0; q3_bar; 0; 0; 0]);
% disp(tau_eq);


%% ----- Linearizarion (Taylor expansion) -----

[A_lin, B_lin, C_lin, D_lin, E_lin] = LinModel(f_ss, x_ss, u_ss, w_ss, G);

% disp('Linearized State Matrix A:');
% pretty(A_lin)
% disp('Linearized Input Matrix B:');
% pretty(B_lin)
% disp('Linearized Output Matrix C:');
% pretty(C_lin)
% disp('Linearized Feedthrough Matrix D:');
% pretty(D_lin)
% disp('Linearized Perturbation Matrix E:');
% pretty(E_lin)


%% ----- Numeric parameters substitution -----

[A_num, B_num, C_num, D_num, E_num, sym_list, num_list] = NumSS(A_lin, B_lin, C_lin, D_lin, E_lin);

disp('--- Numeric Matrices (6x6) ---');
disp('Numeric State Matrix A:');
disp(A_num);
disp('Numeric Input Matrix B:');
disp(B_num);
disp('Numeric Output Matrix C:');
disp(C_num);
disp('Numeric Feedthrough Matrix D:');
disp(D_num);
disp('Numeric Perturbation Matrix E:');
disp(E_num);

% disp(eig(A_num));


% EXPORTING MATRICES FOR MAIN CONTROL
disp('Saving original matrices for control analysis...');
save('DV_model.mat', 'A_num', 'B_num', 'C_num', 'D_num', 'E_num');
disp('File "DV_model.mat" generated!');


%% ----- Reduced system -----

% Define the states we want to keep
states_of_interest_idx = [2, 3, 4, 5, 6];

% Reduced A and B matrices
A_reduced = A_num(states_of_interest_idx, states_of_interest_idx);
B_reduced = B_num(states_of_interest_idx, :);

% Reduced C and D matrices (for the output y = [q2; q3])
% (Original C was 3x6, new C is 2x5)
output_of_interest_idx = [2, 3];
C_reduced = C_num(output_of_interest_idx, states_of_interest_idx);
D_reduced = D_num(output_of_interest_idx, :);

% Reduced E matrix
E_reduced = E_num(states_of_interest_idx, :);

disp('--- Reduced Matrices (5x5) ---');
disp('Reduced State Matrix A (5x5):');
disp(A_reduced);
disp('Reduced Input Matrix B (5x3):');
disp(B_reduced);

% EXPORTING MATRICES FOR MAIN CONTROL
disp('Saving reduced matrices for control analysis...');
save('DV_reduced_model.mat', 'A_reduced', 'B_reduced', 'C_reduced', 'D_reduced', 'E_reduced');
disp('File "DV_reduced_model.mat" generated!');

%% ----- Model analysis (FT, poles, zeros, routh table) -----

sys_reduced = ss(A_reduced, B_reduced, C_reduced, D_reduced);
sys_reduced.InputName = {'Tau_1', 'Tau_2', 'Force_3'};
sys_reduced.OutputName = {'q_2', 'q_3'};

% TRANSFER FUNCTIONS

% Convert SS to TF
TFs = tf(sys_reduced);

% Extract individual Transfer Functions for the main diagonal
% (The direct effect of actuator i on joint i)

disp('--- Transfer Functions ---');
% G22(s) = q2(s) / Tau2(s)
disp('Transfer Function G22 (Arm - q2/tau2):');
G22 = TFs(1,2); % Output 1 (q2), Input 2 (Tau2) from reduced system
G22 = minreal(G22, 1e-4); % Cancels close poles/zeros  
display(G22);

disp('Transfer Function G33 (Prismatic - q3/F3):');
G33 = TFs(2,3); % Output 2 (q3), Input 3 (F3)
G33 = minreal(G33, 1e-4);
display(G33);

% POLES AND ZEROS

% Poles
poles = pole(sys_reduced);
disp('Poles of the System:');
disp(poles);

% Transmission zeros (MIMO zeros)
zeros_mimo = tzero(sys_reduced);
disp('Transmission Zeros (MIMO):');
disp(zeros_mimo);


% ROUTH-HURWITZ

% Get the characteristic polynomial: det(sI - A) = 0
poly_char = poly(A_reduced); 
disp('Characteristic Polynomial:');
s = tf('s');
poly_display = 0;
for i = 1:length(poly_char)
    poly_display = poly_display + poly_char(i) * s^(length(poly_char)-i);
end
display(poly_display);

% Routh table
disp('Routh Table:');
rh_table = RouthTab(poly_char);
disp(rh_table);


%% ----- Frequency response (bode plots) -----

% Recall sys_reduced inputs: {Tau_1, Tau_2, Force_3}
% Recall sys_reduced outputs: {q_2, q_3}

% G_arm: Input 2 (Tau_2) -> Output 1 (q_2)
G_arm = sys_reduced(1, 2); 

% G_prism: Input 3 (Force_3) -> Output 2 (q_3)
G_prism = sys_reduced(2, 3);

figure('Name', 'Frequency response: Bode plots', 'Color', 'w');
subplot(1, 2, 1);
bode(G_arm);
grid on;
subplot(1, 2, 2);
bode(G_prism);
grid on;

% Performance metrics
disp('--- Frequency Domain Metrics ---');

[Gm_2, Pm_2, Wcg_2, Wcp_2] = margin(G_arm);
bw_2 = bandwidth(G_arm);
fprintf('\nArm Link (q2):\n');
fprintf('  - Bandwidth:      %.2f rad/s (System speed)\n', bw_2);
fprintf('  - Phase Margin:   %.2f deg   (Damping/Stability)\n', Pm_2);
if ~isinf(Gm_2)
    fprintf('  - Gain Margin:    %.2f dB\n', 20*log10(Gm_2));
else
    fprintf('  - Gain Margin:    Infinite\n');
end

[Gm_3, Pm_3, Wcg_3, Wcp_3] = margin(G_prism);
bw_3 = bandwidth(G_prism);
fprintf('\nPrismatic Link (q3):\n');
fprintf('  - Bandwidth:      %.2f rad/s\n', bw_3);
fprintf('  - Phase Margin:   %.2f deg\n', Pm_3);


%% ----- Time domain respose (step input) -----

% Simulation Parameters
t_sim = 0:0.01:10; % 10 seconds
num_steps = length(t_sim);

% PLOT 1: Step on joint 2
u_step_q2 = zeros(num_steps, 3);
u_step_q2(:, 2) = 1; % 1 Nm

% Simulate
[y_q2, t_q2] = lsim(sys_reduced, u_step_q2, t_sim);

figure('Name', 'Time Response (Joint 2)', 'Color', 'w');
subplot(2,1,1);
plot(t_q2, u_step_q2(:,2), 'r--', 'LineWidth', 1.5);
title('Input: Torque on Joint 2');
ylabel('Torque (N.m)'); grid on;
subplot(2,1,2);
plot(t_q2, y_q2(:,1), 'b-', 'LineWidth', 2); % Output 1 is q2
title('Output: Angular Displacement of Joint 2');
ylabel('Angle (rad)'); xlabel('Time (s)'); grid on;


% PLOT 2: Step on joint 3
u_step_q3 = zeros(num_steps, 3);
u_step_q3(:, 3) = 1; % 1 N

% Simulate
[y_q3, t_q3] = lsim(sys_reduced, u_step_q3, t_sim);

figure('Name', 'Time Response (Joint 3)', 'Color', 'w');
subplot(2,1,1);
plot(t_q3, u_step_q3(:,3), 'm--', 'LineWidth', 1.5);
title('Input: Force on Joint 3');
ylabel('Force (N)'); grid on;
xlim([0, 5]);
subplot(2,1,2);
plot(t_q3, y_q3(:,2), 'k-', 'LineWidth', 2); % Output 2 is q3
title('Output: Linear Displacement of Joint 3');
ylabel('Displacement (m)'); xlabel('Time (s)'); grid on;
xlim([0, 5]);

% PLOT 3: Initial condition recovery (disturbance)
% The robot was bumped 5cm away from equilibrium -> see it return to zero (stability check)

% Initial State: q3 displaced by 0.05 m (5 cm)
% x_reduced state order: [q2, q3, dq1, dq2, dq3]
x0 = [0; 0.05; 0; 0; 0]; 

% 2. Simulate (Zero input)
[y_init, t_init] = initial(sys_reduced, x0, t_sim);

figure('Name', 'Initial Condition Response', 'Color', 'w');
plot(t_init, y_init(:,2), 'g-', 'LineWidth', 2);
title('Disturbance Rejection: Initial Displacement q_3 = 5cm');
ylabel('Displacement q_3 (m)'); xlabel('Time (s)'); 
grid on;
xlim([0, 5]);
yline(0, 'k--'); % Equilibrium line

% DC gains calculation (Steady state)
dc_gains = dcgain(sys_reduced);

fprintf('\n--- Steady State Analysis (Step Inputs) ---\n');
fprintf('Joint 2: For 1 N.m input, displaces %.4f rad\n', dc_gains(1,2));
fprintf('   -> Equivalent gravitational stiffness: %.2f N.m/rad\n', 1/dc_gains(1,2));

fprintf('Joint 3: For 1 N input, displaces %.4f m\n', dc_gains(2,3));
fprintf('   -> Environment stiffness (K_env): %.2f N/m\n', 1/dc_gains(2,3));



%% ----- Time domain respose (realistic input) -----
disp('--- Time Domain Simulations ---')

% Generate matlab functions to handle non-linear B and n functions
B_numeric_expression = subs(B, sym_list, num_list);
n_numeric_expression = subs(n, sym_list, num_list);

vars_state = [q1; q2; q3; dq1; dq2; dq3];

matlabFunction(B_numeric_expression, 'File', 'GetBmatrix', 'Vars', {vars_state});
matlabFunction(n_numeric_expression, 'File', 'Getnvector', 'Vars', {vars_state, w});

% Dynamics function
function dxdt = robot_dynamics(t, x, u_applied)
    dq = x(4:6);
    dxdt = [dq; GetBmatrix(x) \ (u_applied - Getnvector(x, 0))];
end

% Ode options
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

% Torque needed no mantain equilibirum pose 
x_eq = [0; 0; 0.36; 0; 0; 0]; % q3 = q3_bar = 0.36
u_eq = Getnvector(x_eq, 0);
fprintf('\nEquilibrium Torque: [%.2f; %.2f; %.2f]\n', u_eq(1), u_eq(2), u_eq(3));

% --------------------------------------------------

%% SCENARIO 1A: STEP ON JOINT 2 (0.1 N.m)
func_value = 0.1;
sim_time = 5;
fprintf('\nRunning Scenario 1A: %.2fNm step on joint 2\n', func_value);

% Linear Simulation
t_vec = linspace(0, sim_time, 500);
u_step_lin = zeros(length(t_vec), 3);
u_step_lin(:, 2) = func_value; % Step on Tau2
[y_lin, t_lin] = lsim(sys_reduced, u_step_lin, t_vec);
q3_lin = y_lin(:,1) + x_eq(2); 

% Nonlinear Simulation
delta_u = [0; func_value; 0];
total_u = u_eq + delta_u;
[t_non, x_non] = ode45(@(t,x) robot_dynamics(t, x, total_u), [0 sim_time], x_eq, opts);

% Plot Comparison
figure('Name', 'Scenario 1A: Comparison', 'Color', 'w');
plot(t_lin, q3_lin, 'b--', 'LineWidth', 2); hold on;
plot(t_non, x_non(:,2), 'r-', 'LineWidth', 1.5);
title('Scenario 1A: Step Response on joint 2 (0.1Nm)');
xlabel('Time (s)'); ylabel('Joint 2 Position (rad)');
legend('Linear Model', 'Nonlinear Model');
grid on;


%% SCENARIO 1B: STEP ON JOINT 3 (100 N)
func_value = -100;
sim_time = 5; 
fprintf('Running Scenario 1B: %.2fNm step on joint 3\n', func_value);

% Linear Simulation
t_vec = linspace(0, sim_time, 500);
u_step_lin = zeros(length(t_vec), 3);
u_step_lin(:, 3) = func_value; % Step on F3
[y_lin, t_lin] = lsim(sys_reduced, u_step_lin, t_vec);
q3_lin = y_lin(:,2) + x_eq(3); 

% Nonlinear Simulation
delta_u = [0; 0; func_value];
total_u = u_eq + delta_u;
[t_non, x_non] = ode45(@(t,x) robot_dynamics(t, x, total_u), [0 sim_time], x_eq, opts);

% Plot Comparison
figure('Name', 'Scenario 1B: Comparison', 'Color', 'w');
plot(t_lin, q3_lin, 'b--', 'LineWidth', 2); hold on;
plot(t_non, x_non(:,3), 'r-', 'LineWidth', 1.5);
title('Scenario 1B: Step Response on joint 3 (-100N)');
xlabel('Time (s)'); ylabel('Joint 3 Position (m)');
legend('Linear Model', 'Nonlinear Model');
grid on;


%% SCENARIO 2A: SIN ON JOINT 2
sim_time = 20; 
t_vec = linspace(0, sim_time, 500);
func_value = 1;
A = func_value;
omega = 5;

fprintf('Running Scenario 2A: %dsin(%dt) on joint 2\n', A, omega);

% Linear Simulation
u_sin_lin = zeros(length(t_vec), 3);
u_sin_lin(:, 2) = A * sin(omega * t_vec); % sine applied to torque Tau2
[y_lin, t_lin] = lsim(sys_reduced, u_sin_lin, t_vec);
q2_lin = y_lin(:,1) + x_eq(2); 

% Nonlinear Simulation
u_sin = @(t) u_eq + [0;
                     A*sin(omega*t);
                     0];

[t_non, x_non] = ode45(@(t,x) robot_dynamics(t, x, u_sin(t)), [0 sim_time], x_eq, opts);

% Plot Comparison
figure('Name', 'Scenario 2A: Comparison', 'Color', 'w');
plot(t_lin, q2_lin, 'b--', 'LineWidth', 2); hold on;
plot(t_non, x_non(:,2), 'r-', 'LineWidth', 1.5);
title('Scenario 2A: Sin Response on joint 2 (sin(5t))');
xlabel('Time (s)'); ylabel('Joint 2 Position (rad)');
legend('Linear Model', 'Nonlinear Model');
grid on;


%% SCENARIO 2B: SIN ON JOINT 3 
sim_time = 20; 
t_vec = linspace(0, sim_time, 500);
func_value = 50;
A = func_value;  
omega = 2;

fprintf('Running Scenario 2B: %dsin(%dt) on joint 3\n', A, omega);

% Linear Simulation
u_sin_lin = zeros(length(t_vec), 3);
u_sin_lin(:, 3) = A * sin(omega * t_vec); % sine applied to torque Tau3
[y_lin, t_lin] = lsim(sys_reduced, u_sin_lin, t_vec);
q3_lin = y_lin(:,2) + x_eq(3); 

% Nonlinear Simulation
u_sin = @(t) u_eq + [0;
                     0;
                     A*sin(omega*t)];

[t_non, x_non] = ode45(@(t,x) robot_dynamics(t, x, u_sin(t)), [0 sim_time], x_eq, opts);

% Plot Comparison
figure('Name', 'Scenario 2B: Comparison', 'Color', 'w');
plot(t_lin, q3_lin, 'b--', 'LineWidth', 2); hold on;
plot(t_non, x_non(:,3), 'r-', 'LineWidth', 1.5);
title('Scenario 2B: Sin Response on joint 3 (50sin(2t))');
xlabel('Time (s)'); ylabel('Joint 3 Position (m)');
legend('Linear Model', 'Nonlinear Model');
grid on;

% --------------------------------------------------
%% SCENARIO 3: "SUTURING" TASK SIMULATION

% Generate IK function
IK_sym = IKin(); 
syms a2 real 
IK_subs = subs(IK_sym, a2, 0.4);

% Inputs for IKin are Pe_x, Pe_y, Pe_z
syms Pe_x Pe_y Pe_z real
matlabFunction(IK_subs, 'File', 'GetIKnumeric', 'Vars', {Pe_x, Pe_y, Pe_z});

% Full 6x6 model
sys_full = ss(A_num, B_num, C_num, D_num);

% --------------------------------------------------
% SCENARIO 3A (YZ plane circle)
disp('Running Scenario 3A: YZ plane circle');

% Define Trajectoy
T_task = 4.0; 
num_points = 1000;
t_traj = linspace(0, T_task, num_points); 
dt = t_traj(2) - t_traj(1);

% Center of operation (from equilibrium q = [0; 0; 0.2])
% x ~ a2, y = 0, z = -q3
q3_bar_val = 0.2;
center_x = 0.4; center_y = 0; center_z = -q3_bar_val; 

% Stitch parameters
radius = 0.04; % 4 cm
freq = 2 * pi * 1/T_task;

% Desired Path (Circle in YZ plane)
x_ref = center_x * ones(size(t_traj)); 
y_ref = center_y + radius * sin(freq * t_traj);
z_ref = center_z - radius + radius * cos(freq * t_traj);
vec_ref = [x_ref; y_ref; z_ref];

% Inverse Dynamics
[q_ref, dq_ref, ~, tau_feedforward] = IDyn(num_points, vec_ref, dt);

% Simulation Inputs
% Delta u (for linear model)
delta_u = tau_feedforward - u_eq; 

% Initial condition
x0_task = [q_ref(:,1); dq_ref(:,1)];

% Linear Simulation
[~, ~, x_lin] = lsim(sys_full, delta_u,t_traj, x0_task-x_eq);

% Recover total state: x_linear = delta_x + x_equilibrium
x_lin(:,1) = x_lin(:,1) + x_eq(1);
x_lin(:,2) = x_lin(:,2) + x_eq(2);
x_lin(:,3) = x_lin(:,3) + x_eq(3);

% Nonlinear Simulation
% Interpolator for the torque input
u_interp = @(t) [interp1(t_traj, tau_feedforward(1,:), t, 'linear', 'extrap');
                 interp1(t_traj, tau_feedforward(2,:), t, 'linear', 'extrap');
                 interp1(t_traj, tau_feedforward(3,:), t, 'linear', 'extrap')];


[t_non, x_non] = ode45(@(t,x) robot_dynamics(t, x, u_interp(t)), t_traj, x0_task, opts);

% Plot 
PlotSim(t_traj, vec_ref, q_ref, tau_feedforward, x_lin, t_non, x_non);

% Animation 
AnimSim(t_non, x_non, vec_ref);


%% SCENARIO 3B (XY plane circle)
disp('Running Scenario 3B: XY plane circle');

% Desired Path (Circle in XY plane)
% x varies (reach), y varies (base rotation), z is constant
x_ref = center_x + radius * cos(freq * t_traj); 
y_ref = center_y + radius * sin(freq * t_traj);
z_ref = center_z * ones(size(t_traj)); 
vec_ref = [x_ref; y_ref; z_ref];

% Inverse Dynamics
[q_ref, dq_ref, ~, tau_feedforward] = IDyn(num_points, vec_ref, dt);

% Simulation Inputs
% Delta u (for linear model)
delta_u = tau_feedforward - u_eq; 

% Initial condition
x0_task = [q_ref(:,1); dq_ref(:,1)];

% Linear Simulation
[~, ~, x_lin] = lsim(sys_full, delta_u, t_traj, x0_task-x_eq);

% Recover total state: x_linear = delta_x + x_equilibrium
x_lin(:,1) = x_lin(:,1) + x_eq(1);
x_lin(:,2) = x_lin(:,2) + x_eq(2);
x_lin(:,3) = x_lin(:,3) + x_eq(3);

% Nonlinear Simulation
% Interpolator for the torque input
u_interp = @(t) [interp1(t_traj, tau_feedforward(1,:), t, 'linear', 'extrap');
                 interp1(t_traj, tau_feedforward(2,:), t, 'linear', 'extrap');
                 interp1(t_traj, tau_feedforward(3,:), t, 'linear', 'extrap')];

[t_non, x_non] = ode45(@(t,x) robot_dynamics(t, x, u_interp(t)), t_traj, x0_task, opts);

% Plot 
PlotSim(t_traj, vec_ref, q_ref, tau_feedforward, x_lin, t_non, x_non);

% Animation 
AnimSim(t_non, x_non, vec_ref);

%% Scenario 3C (Helical trajectory, circular in YZ advancing in X)
disp('Running Scenario 3C: Helical trajectory');

% Trajectory Setup
T_task = 10.0; % Longer duration for multiple turns
num_points = 2000;
t_traj = linspace(0, T_task, num_points); 
dt = t_traj(2) - t_traj(1);

% Parameters
q3_bar_val = 0.2;
center_y = 0; 
center_z = -q3_bar_val; 

% Spiral Parameters
a2_val = 0.4;
radius = 0.03; % 3 cm radius
freq = 2 * pi * 0.5; % 0.5 Hz (1 turn every 2 seconds)
x_start = a2_val;      % Start at "elbow" distance
x_end = a2_val + 0.10; % Advance 10 cm
x_velocity = (x_end - x_start) / T_task;

% 2. Desired Path (Helix along X)
% X advances linearly
x_ref = linspace(x_start, x_end, num_points);
% Y and Z circle around the center
y_ref = center_y + radius * sin(freq * t_traj);
z_ref = center_z - radius + radius * cos(freq * t_traj); 
vec_ref = [x_ref; y_ref; z_ref];

% Inverse Dynamics
[q_ref, dq_ref, ~, tau_feedforward] = IDyn(num_points, vec_ref, dt);

% Simulation Inputs
% Delta u (for linear model)
delta_u = tau_feedforward - u_eq; 

% Initial condition
x0_task = [q_ref(:,1); dq_ref(:,1)];

% Linear Simulation
[y_lin, t_lin, x_lin] = lsim(sys_full, delta_u, t_traj, x0_task-x_eq);

% Recover total state: x_linear = delta_x + x_equilibrium
x_lin(:,1) = x_lin(:,1) + x_eq(1);
x_lin(:,2) = x_lin(:,2) + x_eq(2);
x_lin(:,3) = x_lin(:,3) + x_eq(3);


% Nonlinear Simulation
% Interpolator for the torque input
u_interp = @(t) [interp1(t_traj, tau_feedforward(1,:), t, 'linear', 'extrap');
                 interp1(t_traj, tau_feedforward(2,:), t, 'linear', 'extrap');
                 interp1(t_traj, tau_feedforward(3,:), t, 'linear', 'extrap')];

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_non, x_non] = ode45(@(t,x) robot_dynamics(t, x, u_interp(t)), t_traj, x0_task, opts);

% Plot 
PlotSim(t_traj, vec_ref, q_ref, tau_feedforward, x_lin, t_non, x_non);

% Animation 
AnimSim(t_non, x_non, vec_ref);










