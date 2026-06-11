% CLEANUP AND INITIALIZATION

clear;
clc;
close all;

% Loading reduced model matrices
nome_arquivo = 'DV_reduced_model.mat';

if exist(nome_arquivo, 'file')
    load(nome_arquivo);
    fprintf('=> Modelo do daVinci carregado com sucesso de "%s"!\n', nome_arquivo);
else
    error('Erro: Arquivo "%s" não encontrado. Rode o script MainDyn primeiro.', nome_arquivo);
end


%% SYSTEM ANALYSIS

% Reduced system
A = A_reduced;
B = B_reduced;
C = C_reduced;
D = D_reduced;
E = E_reduced;

n_states = size(A, 1);
n_inputs = size(B, 2);
n_outputs = size(C, 1);

%ControlAnalysis(A, B, C, D);

%% LQR CONTROL

% Definition of maximum tolerable limits (Bryson's approach)
err_q2_max = 0.001;  % 0.001 rad of tolerable error
err_q3_max = 0.001;  % 0.001 m (1 mm) of tolerable error

err_dq1_max = 0.1;   % 0.1 rad/s tolerable
err_dq2_max = 0.1;   % 0.1 rad/s tolerable
err_dq3_max = 0.01;   % 0.1 m/s tolerable

% Definition of actuator effort (Motor Capacity)
tau1_max = 5.3;  
tau2_max = 41; 
f3_max   = 420;

% Q matrix (5x5 for the reduced system)
Q_lqr = diag([1/err_q2_max^2, 1/err_q3_max^2, ...
              1/err_dq1_max^2, 1/err_dq2_max^2, 1/err_dq3_max^2]);

% R matrix (3x3 for the 3 actuators)
R_lqr = diag([tau1_max*10, tau2_max*5, f3_max/10]);
% R_lqr = diag([tau1_max*10, tau2_max*10, f3_max/10]);
% R_lqr = diag([1e2, 1e2, 1e1]);

disp('==============================');
disp(' PARÂMETROS DE CONTROLE - LQR ');
disp('==============================');

disp('Matriz de penalidade Q:');
disp(Q_lqr);
disp('Matriz de penalidade R:');
disp(R_lqr);

[K_lqr, ~, poles_cl] = lqr(A, B, Q_lqr, R_lqr);

disp('Matriz de ganhos K:');
disp(K_lqr);

disp('Polos em malha fechada:');
disp(poles_cl);

disp('Processando análises em malha fechada para LQR...');
% AnalyzeCL(A, B, C, D, E, K_lqr, 'LQR');

%% POLE PLACEMENT CONTROL

% Strategy based on LQR poles
% ctrl_poles = poles_cl';

% Surgical damping strategy (Butterworth in tissue + Real at joints)
ctrl_poles = [-5, -7 + 3i, -7 - 3i, -18 + 12i, -18 - 12i];

disp('============================================');
disp(' PARÂMETROS DE CONTROLE - ALOCAÇÃO DE POLOS ');
disp('============================================');
disp('Polos desejados para a malha fechada:');
disp(ctrl_poles');

% Calculation of the gain matrix using the Place algorithm
F_pp = place(A, B, ctrl_poles);
% F_pp = F_pp/1.2;

disp('Matriz de ganhos F calculada:');
disp(F_pp);

disp('Processando análises em malha fechada para alocação de polos...');
% AnalyzeCL(A, B, C, D, E, F_pp, 'PP');


%% FULL-ORDER OBSERVER - IDENTITY

% LQR approach
Wc = 5000 * eye(n_states); % "process noise"  — increase = faster observer
Vc = eye(n_outputs);    % "measurement noise"  — increase = slower observer
L = lqr(A', C', Wc, Vc).';

%Pole Placement approach
% poles_obs_id = [-30+8i,-30+8i, -20, -25+2i, -25-2i];
% L = place(A',C',poles_obs_id)';

fprintf('Ganho do observador L:\n'); disp(L);
fprintf('Autovalores de (A - L*C):\n'); disp(eig(A - L*C).');

% Isolated observer system (inputs = [u; y], output = x_hat)
A_obs = A - L * C;
B_obs = [B, L];
C_obs = eye(n_states);
D_obs = zeros(n_states, n_inputs + n_outputs);
sys_obs = ss(A_obs, B_obs, C_obs, D_obs);

% pzmap(sys_obs)

%% REDUCED-ORDER OBSERVER

% Construction of transformation T and partitioning
V_mat = [zeros(n_states - n_outputs, n_outputs), eye(n_states - n_outputs)];
T  = [C; V_mat];

if abs(det(T)) < 1e-10
    error('T = [C; V_mat] e (quase) singular: escolha outra particao V_mat.');
end
Ti = inv(T);

M = Ti(:, 1:n_outputs);            % n x m
N = Ti(:, n_outputs+1:end);        % n x (n-m)

A11 = C     * A * M;   % m     x m
A12 = C     * A * N;   % m     x (n-m)
A21 = V_mat * A * M;   % (n-m) x m
A22 = V_mat * A * N;   % (n-m) x (n-m)
B1  = C     * B;       % m     x n_inputs
B2  = V_mat * B;       % (n-m) x n_inputs

% Observability and conditioning of the reduced pair (A22, A12)
rank_obs_red = rank(obsv(A22, A12));
fprintf('Posto da observabilidade reduzida = %d (necessario n-m = %d)\n', ...
        rank_obs_red, n_states - n_outputs);
if rank_obs_red < n_states - n_outputs
    error('Par (A22, A12) NAO observavel: observador reduzido nao realizavel.');
end
fprintf('cond(obsv(A22,A12)) = %.3e   (se >> 1e6, o ganho J tende a crescer)\n', ...
        cond(obsv(A22, A12)));

% Gain J via LQE in the fictitious system (A22', A12')
Qe_red = 8000 * eye(n_states - n_outputs); % "process noise"  — increase = faster
Re_red = eye(n_outputs);              % "measurement noise"  — increase = slower
J = lqr(A22', A12', Qe_red, Re_red).';

% Gain J via PP in the fictitious system (A22', A12')
% J = place(A22', A12', 3*ctrl_poles).';

% Observer matrices
F_red = A22 - J * A12;
G_red = A21 - J * A11 + F_red * J;
H_red = B2  - J * B1;
S_red = M   + N * J;

fprintf('Polos do observador reduzido (eig F):\n'); disp(eig(F_red).');
fprintf('norm(J) = %.3e\n', norm(J));

% Augmented system [x; z] in open loop

% [dx]   [ A       0 ] [x]   [ B ]      [ E ]
% [dz] = [ G*C     F ] [z] + [ H ] u  + [ 0 ] w
A_aug_red = [A,         zeros(n_states, n_states - n_outputs);
             G_red * C, F_red];
B_aug_red = [B; H_red];
E_aug_red = [E; zeros(n_states - n_outputs, size(E, 2))];
C_aug_red = [C, zeros(n_outputs, n_states - n_outputs)];

% Reconstruction of x_hat from the augmented state:

% x_hat = S*C*x + N*z  =  recover_xhat * [x; z]
recover_xhat = [S_red * C, N];

% --- Correct initialization of z (use in block 3) ---------------------
% The reduced observer reads y exactly; it only has freedom in the unmeasured part (w).
% Given an initial estimate of velocities w_hat0 and the real state x0:
%   w_hat = z + J*y   =>   z0 = w_hat0 - J*C*x0
% w_hat0 = 0  => observer "unaware" of velocities.
build_z0 = @(w_hat0, x0) w_hat0 - J * C * x0;

%% OBSERVER ERROR CONVERGENCE (OPEN LOOP)

% Simulation time
t_sim = 0 : 1e-3 : 2;

% REAL initial state of the plant
x0_real = [0.05; 0.02; 1; 3; 2];

% Initial estimated state
x0_hat = zeros(n_states, 1);

% Initial error
e0 = x0_real - x0_hat;

% Reduced-order observer packing
obs_red.F_red        = F_red;
obs_red.A_aug_red    = A_aug_red;
obs_red.recover_xhat = recover_xhat;
obs_red.J            = J;

% State names
state_names = {'q_2','q_3','dq_1','dq_2','dq_3'};

Comp_Obsv(A, C, L, obs_red, x0_real, t_sim, state_names, x0_hat)

%% CONVERGENCE ANALYSIS AND POLE SEPARATION

obs_red.F_red        = F_red;
obs_red.A_aug_red    = A_aug_red;
obs_red.recover_xhat = recover_xhat;
obs_red.J            = J;

state_labels  = {'Posicao q2 (rad)','Posicao q3 (m)', ...
           'Velocidade dq1 (rad/s)','Velocidade dq2 (rad/s)','Velocidade dq3 (m/s)'};
x0_real = [0.05; 0.02; 1; 3; 2];
x0_hat  = [0.05; 0.02; 0; 0; 0];

AnalyseCL_Obsv(A, C, L, obs_red, poles_cl, x0_real, t_sim, state_labels, x0_hat)

%% FULL CLOSED-LOOP SIMULATION (LQR + REDUCED OBSERVER)

disp('==================================================');
disp(' INICIANDO SIMULAÇÃO EM MALHA FECHADA COMPLETA... ');
disp('==================================================');

% 1. Construction of the Augmented System Matrices in Closed Loop

% Combining the real plant dynamics with the control law u = -K_lqr * x_hat
A_cl_total = [ (A - B * K_lqr * S_red * C),            (-B * K_lqr * N);
               (G_red * C - H_red * K_lqr * S_red * C), (F_red - H_red * K_lqr * N) ];

% Disturbance input: the signal w enters only the physical plant
B_cl_total = [ E; 
               zeros(n_states - n_outputs, size(E, 2)) ];

% Auxiliary C matrix to extract all state variables from the ss block
C_cl_total = eye(2 * n_states - n_outputs);
D_cl_total = zeros(2 * n_states - n_outputs, 1);

% Instantiation of the global Closed-Loop system
sys_cl_total = ss(A_cl_total, B_cl_total, C_cl_total, D_cl_total);

% 2. Time Configuration and Sinusoidal Respiratory Disturbance
t_sim2 = 0 : 0.001 : 5;
f_respiracao = 0.5;    % 30 breaths per minute
omega_w = 2 * pi * f_respiracao;
w_senoidal = 0.005 * sin(omega_w * t_sim2); % 5 mm breathing amplitude

% 3. Correct Initial Conditions (Misaligned plant and Observer in the dark)
x0_real = [0.05; 0.02; 1; 3; 2]; 
z0_correto = build_z0(zeros(n_states - n_outputs, 1), x0_real);
x0_cl_total = [x0_real; z0_correto];

% 4. Execution of the Combined Time Simulation (Jolt + Continuous Breathing)
[X_cl_total, ~] = lsim(sys_cl_total, w_senoidal, t_sim2, x0_cl_total);

% Separation of Real States and Observer States
x_real_sim2 = X_cl_total(:, 1:n_states);
z_sim2      = X_cl_total(:, n_states+1:end);

% Reconstruction of estimated variables (x_hat) and control efforts (u) over time
x_hat_sim2 = zeros(length(t_sim2), n_states);
u_sim2     = zeros(length(t_sim2), n_inputs);

for i = 1:length(t_sim2)
    % Reconstruct the x_hat vector using the output equation of the reduced-order observer
    x_hat_sim2(i, :) = (S_red * C * x_real_sim2(i, :)' + N * z_sim2(i, :)')';
    
    % Calculate the actual torque injected at the joints: u = -K * x_hat
    u_sim2(i, :) = (-K_lqr * x_hat_sim2(i, :)')';
end

%  PLOT GENERATION

% --- PLOT 1: JOINT POSITION BEHAVIOR ---
figure('Name', 'Malha Fechada Completa: Posicoes das Juntas', 'Color', 'w');

subplot(2,1,1);
plot(t_sim2, x_real_sim2(:,1), 'b', 'LineWidth', 2); grid on;
title('Resposta da Junta q2 (Elevacao) com Controle via Estado Estimado');
ylabel('Posicao Real (rad)');

subplot(2,1,2);
plot(t_sim2, x_real_sim2(:,2), 'b', 'LineWidth', 2); grid on;
title('Resposta da Junta q3 (Insercao) com Controle via Estado Estimado');
xlabel('Tempo (s)'); ylabel('Posicao Real (m)');

% --- PLOT 2: REALISTIC CONTROL EFFORT BEHAVIOR ---
figure('Name', 'Malha Fechada Completa: Esforcos de Controle', 'Color', 'w');
plot(t_sim2, u_sim2(:,1), 'b', 'LineWidth', 2); hold on;
plot(t_sim2, u_sim2(:,2), 'r', 'LineWidth', 2);
plot(t_sim2, u_sim2(:,3), 'g', 'LineWidth', 2); grid on;
title('Esforco de Controle dos Motores (Realimentacao por Observador Reduzido)');
xlabel('Tempo (s)'); ylabel('Torque (Nm) / Forca (N)');
legend('tau1 (Base)', 'tau2 (Elevacao)', 'F3 (Insercao)', 'Location', 'best');

disp('=> Simulação do Cenário 2 concluída e gráficos gerados!');