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
err_q3_max = 0.0001;  % 0.0001 m (0.1 mm) of tolerable error

err_dq1_max = 0.5;   % 0.5 rad/s tolerable
err_dq2_max = 0.5;   % 0.5 rad/s tolerable
err_dq3_max = 0.5;   % 0.5 m/s tolerable

% Definition of actuator effort (Motor Capacity)
tau1_max = 5.3;  
tau2_max = 41; 
f3_max   = 420;

% Q matrix (5x5 for the reduced system)
Q_lqr = diag([1/err_q2_max^2, 1/err_q3_max^2, ...
              1/err_dq1_max^2, 1/err_dq2_max^2, 1/err_dq3_max^2]);

% R matrix (3x3 for the 3 actuators)
R_lqr = diag([tau1_max*10, tau2_max*5, f3_max/100]);

disp('========================================');
disp(' PARÂMETROS DE CONTROLE - REGULADOR LQR ');
disp('========================================');

disp('Matriz de penalidade Q:');
disp(Q_lqr);
disp('Matriz de penalidade R:');
disp(R_lqr);

[K_lqr, ~, poles_cl] = lqr(A, B, Q_lqr, R_lqr);

disp('Matriz de ganhos K:');
disp(K_lqr);

disp('Polos em malha fechada:');
disp(poles_cl);

% disp('Processando análises em malha fechada para LQR...');
% AnalyzeCL(A, B, C, D, E, K_lqr, 'LQR');

%% POLE PLACEMENT CONTROL

% Strategy based on LQR poles
% ctrl_poles = poles_cl';

% Surgical damping strategy (Butterworth in tissue + Real at joints)
ctrl_poles = [-5, -7 + 4i, -7 - 4i, -40 + 20i, -40 - 20i];

disp('=======================================');
disp(' PARÂMETROS DE CONTROLE - REGULADOR PP');
disp('=======================================');
disp('Polos desejados para a malha fechada:');
disp(ctrl_poles');

% Calculation of the gain matrix using the Place algorithm
F_pp = place(A, B, ctrl_poles);
% F_pp = F_pp/1.2;

disp('Matriz de ganhos F calculada:');
disp(F_pp);

% disp('Processando análises em malha fechada para alocação de polos...');
% AnalyzeCL(A, B, C, D, E, F_pp, 'PP');

%% RESPOSTA EM FREQUÊNCIA (BODE / SIGMA)

% Sistemas em malha aberta e malha fechada para comparação
sys_ol     = ss(A,            B, C, D);
sys_cl_lqr = ss(A - B*K_lqr, B, C, D);
sys_cl_pp  = ss(A - B*F_pp,  B, C, D);

% --- Sigma plot: visão MIMO global ---
figure('Name', 'Resposta em Frequencia: Valores Singulares (MIMO)', 'Color', 'w');
sigma(sys_ol,     'k--', sys_cl_lqr, 'b-', sys_cl_pp, 'r-');
grid on;
title('Resposta em Frequencia — Valores Singulares (Comparativo)');
xlabel('Frequencia (rad/s)'); ylabel('Ganho (dB)');
legend('Malha Aberta', 'Malha Fechada LQR', 'Malha Fechada PP', 'Location', 'best');

% --- Bode por canal de saída ---
bode_out_names = {'q_2 (Elevacao)', 'q_3 (Insercao)'};
bode_in_names  = {'tau1 (Base)', 'tau2 (Elevacao)', 'F3 (Insercao)'};

for out_idx = 1:n_outputs
    figure('Name', ['Bode — Saida ' bode_out_names{out_idx}], 'Color', 'w');
    for in_idx = 1:n_inputs
        subplot(n_inputs, 1, in_idx);
        bode(sys_ol(out_idx, in_idx),     'k--', ...
             sys_cl_lqr(out_idx, in_idx), 'b-',  ...
             sys_cl_pp(out_idx, in_idx),  'r-');
        grid on;
        title(['Saida: ' bode_out_names{out_idx} '   |   Entrada: ' bode_in_names{in_idx}]);
        if in_idx == 1
            legend('Malha Aberta', 'LQR', 'PP', 'Location', 'best');
        end
    end
end

%% FULL-ORDER OBSERVER - IDENTITY

% LQR approach
Wc = 1000 * eye(n_states); % "process noise"  — increase = faster observer
Vc = eye(n_outputs);    % "measurement noise"  — increase = slower observer
L = lqr(A', C', Wc, Vc).';

%Pole Placement approach
% poles_obs_id = [-30+8i,-30+8i, -20, -25+2i, -25-2i];
% L = place(A',C',poles_obs_id)';

disp('================================================');
disp(' PARÂMETROS DE CONTROLE - OBSERVADOR IDENTIDADE ');
disp('================================================');

fprintf('Polos do observador (autovalores de (A - L*C)):\n'); disp(eig(A - L*C));
fprintf('Matriz de ganhos do observador (L):\n'); disp(L);


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

disp('===================================================');
disp(' PARÂMETROS DE CONTROLE - OBSERVADOR ORD. REDUZIDA ');
disp('===================================================');

% Observability and conditioning of the reduced pair (A22, A12)
rank_obs_red = rank(obsv(A22, A12));
fprintf('Posto da observabilidade reduzida = %d (necessario n-m = %d)\n', ...
        rank_obs_red, n_states - n_outputs);
if rank_obs_red < n_states - n_outputs
    error('Par (A22, A12) NAO observavel: observador reduzido nao realizavel.');
end
fprintf('cond(obsv(A22,A12)) = %.3e   (se >> 1e6, o ganho J tende a crescer)\n\n', ...
        cond(obsv(A22, A12)));

% Gain J via LQE in the fictitious system (A22', A12')
Qe_red = 1000 * eye(n_states - n_outputs); % "process noise"  — increase = faster
Re_red = eye(n_outputs);              % "measurement noise"  — increase = slower
J = lqr(A22', A12', Qe_red, Re_red).';

% Gain J via PP in the fictitious system (A22', A12')
% J = place(A22', A12', 3*ctrl_poles).';

% Observer matrices
F_red = A22 - J * A12;
G_red = A21 - J * A11 + F_red * J;
H_red = B2  - J * B1;
S_red = M   + N * J;

fprintf('Polos do observador reduzido (eig F):\n'); disp(eig(F_red));
fprintf('norm(J) = %.3e\n\n', norm(J));

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
x0_hat  = [0; 0; 0; 0; 0];

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
t_sim2 = 0 : 0.001 : 2;
f_respiracao = 0.5;    % 30 breaths per minute
omega_w = 2 * pi * f_respiracao;
w_senoidal = 0.005 * sin(omega_w * t_sim2); % 5 mm breathing amplitude

% 3. Correct Initial Conditions (Misaligned plant and Observer in the dark)
x0_real = [0.05; 0.02; 1; 3; 2]; 
z0 = build_z0(zeros(n_states - n_outputs, 1), x0_real);
x0_cl_total = [x0_real; z0];

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

% --- PLOT 3: ESTADO REAL vs ESTIMADO — POSIÇÕES ---
figure('Name', 'Regulador: Estados Reais vs Estimados (Posicoes)', 'Color', 'w');

subplot(2,1,1);
plot(t_sim2, x_real_sim2(:,1), 'b-',  'LineWidth', 2); hold on;
plot(t_sim2, x_hat_sim2(:,1),  'r--', 'LineWidth', 1.5); grid on;
title('q_2 (Elevacao) — Real vs Estimado pelo Obs. Reduzido');
ylabel('Posicao (rad)');
legend('Real', 'Estimado', 'Location', 'best');

subplot(2,1,2);
plot(t_sim2, x_real_sim2(:,2), 'b-',  'LineWidth', 2); hold on;
plot(t_sim2, x_hat_sim2(:,2),  'r--', 'LineWidth', 1.5); grid on;
title('q_3 (Insercao) — Real vs Estimado pelo Obs. Reduzido');
xlabel('Tempo (s)'); ylabel('Posicao (m)');

% --- PLOT 4: ESTADO REAL vs ESTIMADO — VELOCIDADES ---
figure('Name', 'Regulador: Estados Reais vs Estimados (Velocidades)', 'Color', 'w');
vel_names_cl   = {'dq_1 (rad/s)', 'dq_2 (rad/s)', 'dq_3 (m/s)'};
vel_ylabels_cl = {'Vel. (rad/s)', 'Vel. (rad/s)', 'Vel. (m/s)'};

for k = 3:5
    subplot(3,1,k-2);
    plot(t_sim2, x_real_sim2(:,k), 'b-',  'LineWidth', 2); hold on;
    plot(t_sim2, x_hat_sim2(:,k),  'r--', 'LineWidth', 1.5); grid on;
    title([vel_names_cl{k-2} ' — Real vs Estimado pelo Obs. Reduzido']);
    ylabel(vel_ylabels_cl{k-2});
    if k == 5, xlabel('Tempo (s)'); end
    if k == 3, legend('Real', 'Estimado', 'Location', 'best'); end
end

% --- PLOT 5: DERIVADA DO ESFORÇO DE CONTROLE (du/dt) ---
dt_sim2 = t_sim2(2) - t_sim2(1);
du_sim2 = diff(u_sim2) / dt_sim2;
t_du2   = t_sim2(1:end-1);

figure('Name', 'Regulador: Taxa de Variacao do Esforco de Controle', 'Color', 'w');

subplot(3,1,1);
plot(t_du2, du_sim2(:,1), 'b', 'LineWidth', 1.5); grid on;
title('du/dt — tau1 (Base)');
ylabel('d\tau_1/dt  (Nm/s)');

subplot(3,1,2);
plot(t_du2, du_sim2(:,2), 'r', 'LineWidth', 1.5); grid on;
title('du/dt — tau2 (Elevacao)');
ylabel('d\tau_2/dt  (Nm/s)');

subplot(3,1,3);
plot(t_du2, du_sim2(:,3), 'g', 'LineWidth', 1.5); grid on;
title('du/dt — F3 (Insercao)');
xlabel('Tempo (s)'); ylabel('dF_3/dt  (N/s)');

disp('=> Simulação do Cenário 2 concluída e gráficos gerados!');


%% SEGUIDOR VIA LQ (HORIZONTE FINITO)

% Configurações Iniciais de Custo e Tempo
Q_seg = Q_lqr;             % Peso de rastreamento (mesmo do LQR)
R_seg = R_lqr;             % Peso dos motores (mesmo do LQR)
Q1 = 1 * Q_seg;           % Penalidade do estado terminal (igualmente escalada)

t_f = 5.0;                 % Horizonte de tempo final da cirurgia (2 segundos)
t_span_prog = 0 : 0.001 : t_f; 
t_span_back = [t_f, 0];    % Intervalo para integração retroativa

% Pré-cálculo da matriz de acoplamento de entrada
Ri = inv(R_seg);
BRiBt = B * Ri * B';       % Dimensão (n_states x n_states)

% Integração Retroativa da EDO de Riccati (P)
disp('=> Resolvendo EDO de Riccati retroativa...');
P1_flat = Q1(:); % Achata a matriz terminal 5x5 em um vetor 25x1
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

[t_P, P_flat_sol] = ode45(@(t, P_flat) riccati_rhs(t, P_flat, A, BRiBt, Q_seg), t_span_back, P1_flat, options);

% Inverte os resultados para a ordem cronológica correta (0 -> t_f)
t_P = flipud(t_P);
P_flat_sol = flipud(P_flat_sol);

% Integração Retroativa da EDO do Co-Estado (eta)
disp('=> Resolvendo EDO do co-estado retroativa...');
x_ref_tf = daVinci_ref(t_f);
eta1 = Q1 * x_ref_tf; % Condição final do co-estado

[t_eta, eta_sol] = ode45(@(t, eta) eta_rhs(t, eta, A, BRiBt, Q_seg, t_P, P_flat_sol, n_states), t_span_back, eta1, options);

t_eta = flipud(t_eta);
eta_sol = flipud(eta_sol);

% Simulação Progressiva em Malha Fechada Completa (Planta + Obs Reduzido)
disp('=> Iniciando simulação progressiva com realimentação por estimativa...');

% Condições iniciais: planta desalinhada e observador no escuro
x0_real = zeros(1, 5)';

z0_inicial = build_z0(zeros(n_states - n_outputs, 1), x0_real);
X0_global = [x0_real; z0_inicial]; % Vetor de estados combinados (5 reais + 3 do observador)

% Executa a simulação temporal acoplada
[t_sim, X_total] = ode45(@(t, X) forward_sim_rhs(t, X, A, B, C, E, F_red, G_red, H_red, S_red, N, Ri, t_P, P_flat_sol, t_eta, eta_sol, n_states), t_span_prog, X0_global, options);

% Pós-Processamento e Reconstrução dos Sinais
x_real = X_total(:, 1:n_states);
z_obs  = X_total(:, n_states+1:end);

x_hat = zeros(length(t_sim), n_states);
u_seg = zeros(length(t_sim), n_inputs);
y_ref = zeros(length(t_sim), n_outputs);

for i = 1:length(t_sim)
    t = t_sim(i);
    y_atual = C * x_real(i, :)';
    
    % Reconstrução do estado estimado pelo observador reduzido
    x_hat(i, :) = (S_red * y_atual + N * z_obs(i, :)')';
    
    % Recuperação das matrizes variantes no tempo por interpolação
    P_flat = interp1(t_P, P_flat_sol, t, 'linear', 'extrap')';
    P = reshape(P_flat, n_states, n_states);
    eta = interp1(t_eta, eta_sol, t, 'linear', 'extrap')';
    
    % Ganhos variantes no tempo
    K_t = Ri * B' * P;
    u_til = Ri * B' * eta;
    
    % Lei de controle ótima: u(t) = -K(t)*x_hat(t) + u_til(t)
    u_seg(i, :) = (-K_t * x_hat(i, :)' + u_til)';
    
    % Registra as referências de posição para o gráfico
    ref_estado = daVinci_ref(t);
    y_ref(i, :) = (C * ref_estado)';
end

% PLOTS DOS GRÁFICOS DE DESEMPENHO

figure('Name', 'Seguidor LQ: Rastreamento de Trajetoria', 'Color', 'w');

subplot(2,1,1);
plot(t_sim, y_ref(:,1), 'r--', 'LineWidth', 1.5); hold on;
plot(t_sim, x_real(:,1), 'b', 'LineWidth', 2); grid on;
title('Rastreamento da Junta q2 (Elevacao) - Seguidor LQ Finito');
ylabel('Posicao (rad)');
legend('Referencia', 'Resposta Real', 'Location', 'best');

subplot(2,1,2);
plot(t_sim, y_ref(:,2), 'r--', 'LineWidth', 1.5); hold on;
plot(t_sim, x_real(:,2), 'b', 'LineWidth', 2); grid on;
title('Rastreamento da Junta q3 (Insercao) - Seguidor LQ Finito');
xlabel('Tempo (s)'); ylabel('Posicao (m)');

figure('Name', 'Seguidor LQ: Esforco de Controle dos Motores', 'Color', 'w');
yyaxis left
plot(t_sim, u_seg(:,1), 'b', 'LineWidth', 2); hold on;
plot(t_sim, u_seg(:,2), 'r', 'LineWidth', 2);
ylabel('Torque dos Motores (Nm)');
set(gca, 'YColor', 'k');
yyaxis right
plot(t_sim, u_seg(:,3), 'g', 'LineWidth', 2); grid on;
ylabel('Forca de Insercao (N)');
set(gca, 'YColor', 'k');
grid on;
title('Seguidor LQ: Esforco de Controle dos Motores');
xlabel('Tempo (s)');
legend('tau1 (Base)', 'tau2 (Elevacao)', 'F3 (Insercao)', 'Location', 'best');

% Derivada do esforço — Seguidor LQ
dt_lqt   = t_sim(2) - t_sim(1);
du_lqt   = diff(u_seg) / dt_lqt;
t_du_lqt = t_sim(1:end-1);

figure('Name', 'Seguidor LQ: Taxa de Variacao do Esforco de Controle', 'Color', 'w');

subplot(3,1,1);
plot(t_du_lqt, du_lqt(:,1), 'b', 'LineWidth', 1.5); grid on;
title('du/dt — tau1 (Base) — Seguidor LQ');
ylabel('d\tau_1/dt  (Nm/s)');

subplot(3,1,2);
plot(t_du_lqt, du_lqt(:,2), 'r', 'LineWidth', 1.5); grid on;
title('du/dt — tau2 (Elevacao) — Seguidor LQ');
ylabel('d\tau_2/dt  (Nm/s)');

subplot(3,1,3);
plot(t_du_lqt, du_lqt(:,3), 'g', 'LineWidth', 1.5); grid on;
title('du/dt — F3 (Insercao) — Seguidor LQ');
xlabel('Tempo (s)'); ylabel('dF_3/dt  (N/s)');

% Estado real vs estimado — Seguidor LQ (posições)
figure('Name', 'Seguidor LQ: Estados Reais vs Estimados (Posicoes)', 'Color', 'w');

subplot(2,1,1);
plot(t_sim, y_ref(:,1),    'r--', 'LineWidth', 1.5); hold on;
plot(t_sim, x_real(:,1),   'b-',  'LineWidth', 2);
plot(t_sim, x_hat(:,1),    'g-',  'LineWidth', 1.2); grid on;
title('q_2 (Elevacao) — Referencia, Real e Estimado — Seguidor LQ');
ylabel('Posicao (rad)');
legend('Referencia', 'Real', 'Estimado', 'Location', 'best');

subplot(2,1,2);
plot(t_sim, y_ref(:,2),    'r--', 'LineWidth', 1.5); hold on;
plot(t_sim, x_real(:,2),   'b-',  'LineWidth', 2);
plot(t_sim, x_hat(:,2),    'g-',  'LineWidth', 1.2); grid on;
title('q_3 (Insercao) — Referencia, Real e Estimado — Seguidor LQ');
xlabel('Tempo (s)'); ylabel('Posicao (m)');


% ------ FUNÇÕES AUXILIARES DA DINÂMICA ------


% Definição da Trajetória de Referência da Cirurgia
function xr = daVinci_ref(t)
    % Parâmetros da mesma curva exponencial suavizada
    alpha_regime = 10;
    q2_f = 5 * pi / 180; % Alvo de 5 graus (rad)
    q3_f = 0.05;         % Alvo de 5 cm (m)

    % Posições exponenciais: q(t) = q_f * (1 - e^(-alpha * t))
    q2_r = q2_f * (1 - exp(-alpha_regime * t));
    q3_r = q3_f * (1 - exp(-alpha_regime * t));

    % Velocidades exponenciais (Derivadas): dq(t) = q_f * alpha * e^(-alpha * t)
    dq2_r = q2_f * alpha_regime * exp(-alpha_regime * t);
    dq3_r = q3_f * alpha_regime * exp(-alpha_regime * t);

    % Montagem do vetor de estados de referência: [q2; q3; dq1; dq2; dq3]
    xr = [q2_r; q3_r; 0; dq2_r; dq3_r]; 
end

% Lado Direito da EDO de Riccati (dP/dt)
function dP_flat = riccati_rhs(~, P_flat, A, BRiBt, Q)
    n = sqrt(length(P_flat));
    P = reshape(P_flat, n, n);
    dP = -P * A - A' * P + P * BRiBt * P - Q;
    dP_flat = dP(:);
end

% Lado Direito da EDO do Co-Estado (deta/dt)
function deta = eta_rhs(t, eta, A, BRiBt, Q, t_P, P_flat_sol, n_states)
    P_flat = interp1(t_P, P_flat_sol, t, 'linear', 'extrap')';
    P = reshape(P_flat, n_states, n_states);
    deta = -(A' - P * BRiBt) * eta - Q * daVinci_ref(t);
end

% Dinâmica Acoplada Progressiva (Planta Real + Observador Reduzido)
function dX = forward_sim_rhs(t, X, A, B, C, E, F_red, G_red, H_red, S_red, N, Ri, t_P, P_flat_sol, t_eta, eta_sol, n_states)
    % Separa os estados
    x = X(1:n_states);
    z = X(n_states+1:end);
    y = C * x;
    
    % Reconstrói x_hat
    x_hat = S_red * y + N * z;
    
    % Interpola as matrizes de ganho da trajetória do LQT
    P_flat = interp1(t_P, P_flat_sol, t, 'linear', 'extrap')';
    P = reshape(P_flat, n_states, n_states);
    eta = interp1(t_eta, eta_sol, t, 'linear', 'extrap')';
    
    % Calcula a lei de controle variantes no tempo
    K_t = Ri * B' * P;
    u_til = Ri * B' * eta;
    u = -K_t * x_hat + u_til;
    
    % Perturbação respiratória ativa entra pela matriz E
    w = 0.005 * sin(2 * pi * 0.5 * t);
    
    % Derivadas
    dx = A * x + B * u + E * w;
    dz = F_red * z + G_red * y + H_red * u;
    
    dX = [dx; dz];
end

%% ANÁLISE ESTATÍSTICA E POLOS DE MALHA FECHADA

disp('==================================================');
disp(' COMPILANDO RELATÓRIO ESTATÍSTICO DO SEGUIDOR... ');
disp('==================================================');

% 1. Extração dos Polos Assintóticos do Controlador (t = 0 / Regime Estacionário)
P_flat_ss = P_flat_sol(1, :); 
P_ss = reshape(P_flat_ss, n_states, n_states);
K_ss = Ri * B'; 
poles_cl_ss = eig(A - B * K_ss * P_ss);

% Recuperação dos polos do observador reduzido
poles_obs_red = eig(F_red);

% 2. Cálculo das Métricas de Esforço de Controle
max_u = max(abs(u_seg));
mean_u = mean(abs(u_seg));
rms_u = rms(u_seg);

% 3. Cálculo dos Erros de Rastreamento em Regime (Apenas para t >= 0.5s)
idx_regime = (t_sim >= 0.5);
erro_q2 = y_ref(idx_regime, 1) - x_real(idx_regime, 1);
erro_q3 = y_ref(idx_regime, 2) - x_real(idx_regime, 2);

max_error_regime = [max(abs(erro_q2)), max(abs(erro_q3))];
rms_error_regime = [rms(erro_q2), rms(erro_q3)];

% --- IMPRESSÃO DOS RESULTADOS NO TERMINAL ---
fprintf('\n===================================================================\n');
fprintf('        TABELA DE POLOS ATUALIZADA (SINTONIA FINA) \n');
fprintf('===================================================================\n');
fprintf(' Modo |   Controlador (Assintotico)  |   Observador Reduzido   \n');
fprintf('-------------------------------------------------------------------\n');
for i = 1:5
    str_ctrl = sprintf('%+6.2f %+6.2fi', real(poles_cl_ss(i)), imag(poles_cl_ss(i)));
    if i <= 3
        str_red = sprintf('%+6.2f %+6.2fi', real(poles_obs_red(i)), imag(poles_obs_red(i)));
        fprintf('  %d   | %-28s | %-20s \n', i, str_ctrl, str_red);
    else
        fprintf('  %d   | %-28s | %-20s \n', i, str_ctrl, '         -         ');
    end
end
fprintf('===================================================================\n');

fprintf('\n===================================================================\n');
fprintf('            MÉTRICAS DOS ATUADORES (ESFORÇO DE CONTROLE) \n');
fprintf('===================================================================\n');
fprintf(' Atuador       |   Valor Maximo   |   Valor Medio   |   Valor RMS     \n');
fprintf('-------------------------------------------------------------------\n');
fprintf(' tau1 (Base)   | %12.4f Nm | %11.4f Nm | %11.4f Nm \n', max_u(1), mean_u(1), rms_u(1));
fprintf(' tau2 (Elev.)  | %12.4f Nm | %11.4f Nm | %11.4f Nm \n', max_u(2), mean_u(2), rms_u(2));
fprintf(' F3 (Insercao) | %12.4f N  | %11.4f N  | %11.4f N  \n', max_u(3), mean_u(3), rms_u(3));
fprintf('===================================================================\n');

fprintf('\n===================================================================\n');
fprintf('          ERROS DE RASTREAMENTO EM REGIME (t >= 0.5 s) \n');
fprintf('===================================================================\n');
fprintf(' Junta         |     Erro Maximo Absoluto    |       Erro RMS        \n');
fprintf('-------------------------------------------------------------------\n');
fprintf(' q2 (Elevacao) | %18.6f rad       | %16.6f rad  \n', max_error_regime(1), rms_error_regime(1));
fprintf(' q3 (Insercao) | %18.6f m         | %16.6f m    \n', max_error_regime(2), rms_error_regime(2));
fprintf('===================================================================\n\n');



%% SEGUIDOR POR MODELO ASSUMIDO

disp('==================================================');
disp(' INICIANDO SEGUIDOR POR MODELO ASSUMIDO...        ');
disp('==================================================');

% Parâmetro de decaimento para estabilizar em 0.5s (5 / alpha = 0.5)
alpha_regime = 10; 

% Construção da Matriz Ar Autônoma (5x5)
A_r = [ 0,  0,  0,   1,   0;   % d(q2_r)/dt = dq2_r
        0,  0,  0,   0,   1;   % d(q3_r)/dt = dq3_r
        0,  0,  0,   0,   0;   % d(dq1_r)/dt = 0
        0,  0,  0, -alpha_regime,  0;   % d(dq2_r)/dt
        0,  0,  0,   0, -alpha_regime]; % d(dq3_r)/dt

% Matriz de seleção (saídas que devem seguir a referência)
M = C; 

% Cálculo do Ganho de Pré-alimentação Gr do Professor
A_cl_inv = inv(A - B * K_lqr);
N_ff = pinv(M * A_cl_inv * B) * M * A_cl_inv;
G_r = N_ff * (A - A_r);

% Montagem do Sistema Global de 13 Estados [x (5); z (3); x_r (5)]
A_ex = [ (A - B * K_lqr * S_red * C),            (-B * K_lqr * N),         (B * (K_lqr - G_r));
         (G_red * C - H_red * K_lqr * S_red * C), (F_red - H_red * K_lqr * N), (H_red * (K_lqr - G_r));
         zeros(5, 5),                             zeros(5, 3),              A_r ];

% Perturbação respiratória senoidal entra apenas na planta física
E_ex = [ E; 
         zeros(n_states - n_outputs, size(E, 2));
         zeros(5, size(E, 2)) ];

sys_ex = ss(A_ex, E_ex, eye(13), zeros(13, 1));

% Configuração das Condições Iniciais do Vetor Global
t_sim = 0 : 0.001 : 5; 

% Planta real parte do repouso no zero
x0_planta = [0; 0; 0; 0; 0]; 
z0_inicial = build_z0(zeros(3, 1), x0_planta);

% O gerador Ar nasce posicionado em 0, mas com a "velocidade de disparo" calculada
v2_0 = alpha_regime * (5 * pi / 180); % Velocidade inicial para atingir 5 graus
v3_0 = alpha_regime * 0.05;          % Velocidade inicial para atingir 5 cm
x0_referencia = [0; 0; 0; v2_0; v3_0]; 

% Vetor de estado inicial unificado (13 x 1)
X0_global = [x0_planta; z0_inicial; x0_referencia];

% Execução da Simulação (LQR + Obs + Gerador de Trajetória + Respiração)
w_respiracao = 0.005 * sin(2 * pi * 0.5 * t_sim);
[X_out, ~] = lsim(sys_ex, w_respiracao, t_sim, X0_global);

% Extração dos Sinais para Análise
q2_real = X_out(:, 1); q3_real = X_out(:, 2);
q2_ref  = X_out(:, 9); q3_ref  = X_out(:, 10);

u_seg = zeros(length(t_sim), n_inputs);
for i = 1:length(t_sim)
    xhat = S_red * C * X_out(i, 1:5)' + N * X_out(i, 6:8)';
    xr   = X_out(i, 9:13)';
    u_seg(i, :) = (-K_lqr * xhat + (K_lqr - G_r) * xr)';
end

% ----- PLOTS DOS GRÁFICOS DE DESEMPENHO -----
figure('Name', 'Rastreamento por Gerador Exógeno Autônomo', 'Color', 'w');

subplot(2,1,1);
plot(t_sim, q2_ref, 'r--', 'LineWidth', 1.5); hold on;
plot(t_sim, q2_real, 'b', 'LineWidth', 2); grid on;
title('Rastreamento da Junta q2 (Elevacao) - Modelo Assumido');
ylabel('Posicao (rad)');

subplot(2,1,2);
plot(t_sim, q3_ref, 'r--', 'LineWidth', 1.5); hold on;
plot(t_sim, q3_real, 'b', 'LineWidth', 2); grid on;
title('Rastreamento da Junta q3 (Insercao) - Modelo Assumido');
xlabel('Tempo (s)'); ylabel('Posicao (m)');
legend('Referencia da EDO', 'Resposta Real', 'Location', 'best');

% --- GRÁFICO DE ESFORÇO ---
figure('Name', 'Modelo Assumido: Esforco de Controle dos Motores', 'Color', 'w');
yyaxis left
plot(t_sim, u_seg(:,1), 'b', 'LineWidth', 2); hold on;
plot(t_sim, u_seg(:,2), 'r', 'LineWidth', 2);
ylabel('Torque dos Motores (Nm)'); set(gca, 'YColor', 'k');

yyaxis right
plot(t_sim, u_seg(:,3), 'g', 'LineWidth', 2);
ylabel('Forca de Insercao (N)'); set(gca, 'YColor', 'k');
grid on; xlabel('Tempo (s)');
title('Modelo Assumido: Esforco de Controle dos Motores');
legend('tau1 (Base)', 'tau2 (Elevacao)', 'F3 (Insercao)', 'Location', 'best');

% Derivada do esforço — Modelo Assumido
dt_ma   = t_sim(2) - t_sim(1);
du_ma   = diff(u_seg) / dt_ma;
t_du_ma = t_sim(1:end-1);

figure('Name', 'Modelo Assumido: Taxa de Variacao do Esforco de Controle', 'Color', 'w');

subplot(3,1,1);
plot(t_du_ma, du_ma(:,1), 'b', 'LineWidth', 1.5); grid on;
title('du/dt — tau1 (Base) — Modelo Assumido');
ylabel('d\tau_1/dt  (Nm/s)');

subplot(3,1,2);
plot(t_du_ma, du_ma(:,2), 'r', 'LineWidth', 1.5); grid on;
title('du/dt — tau2 (Elevacao) — Modelo Assumido');
ylabel('d\tau_2/dt  (Nm/s)');

subplot(3,1,3);
plot(t_du_ma, du_ma(:,3), 'g', 'LineWidth', 1.5); grid on;
title('du/dt — F3 (Insercao) — Modelo Assumido');
xlabel('Tempo (s)'); ylabel('dF_3/dt  (N/s)');

% Estado real vs estimado — Modelo Assumido (posições)
x_hat_ma = zeros(length(t_sim), n_states);
for i = 1:length(t_sim)
    x_hat_ma(i, :) = (S_red * C * X_out(i, 1:5)' + N * X_out(i, 6:8)')';
end

figure('Name', 'Modelo Assumido: Estados Reais vs Estimados (Posicoes)', 'Color', 'w');

subplot(2,1,1);
plot(t_sim, q2_ref,          'r--', 'LineWidth', 1.5); hold on;
plot(t_sim, q2_real,         'b-',  'LineWidth', 2);
plot(t_sim, x_hat_ma(:,1),   'g-',  'LineWidth', 1.2); grid on;
title('q_2 (Elevacao) — Referencia, Real e Estimado — Modelo Assumido');
ylabel('Posicao (rad)');
legend('Referencia', 'Real', 'Estimado', 'Location', 'best');

subplot(2,1,2);
plot(t_sim, q3_ref,          'r--', 'LineWidth', 1.5); hold on;
plot(t_sim, q3_real,         'b-',  'LineWidth', 2);
plot(t_sim, x_hat_ma(:,2),   'g-',  'LineWidth', 1.2); grid on;
title('q_3 (Insercao) — Referencia, Real e Estimado — Modelo Assumido');
xlabel('Tempo (s)'); ylabel('Posicao (m)');
