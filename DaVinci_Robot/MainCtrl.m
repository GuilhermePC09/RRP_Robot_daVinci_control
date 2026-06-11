% LIMPEZA E INÍCIO

clear;
clc;
close all;

% Carregando matrizes do modelo reduzido
nome_arquivo = 'DV_reduced_model.mat';

if exist(nome_arquivo, 'file')
    load(nome_arquivo);
    fprintf('=> Modelo do daVinci carregado com sucesso de "%s"!\n', nome_arquivo);
else
    error('Erro: Arquivo "%s" não encontrado. Rode o script MainDyn primeiro.', nome_arquivo);
end


%% ANÁLSIE DO SISTEMA

% Sistema reduzido
A = A_reduced;
B = B_reduced;
C = C_reduced;
D = D_reduced;
E = E_reduced;

n_states = size(A, 1);
n_inputs = size(B, 2);
n_outputs = size(C, 1);

%ControlAnalysis(A, B, C, D);

%% CONTROLE POR LQR

% Definição dos limites máximos toleráveis (Abordagem de Bryson)
err_q2_max = 0.001;  % 0.001 rad de erro tolerável
err_q3_max = 0.001;  % 0.001 m (1 mm) de erro tolerável

err_dq1_max = 0.1;   % 0.1 rad/s tolerável
err_dq2_max = 0.1;   % 0.1 rad/s tolerável
err_dq3_max = 0.01;   % 0.1 m/s tolerável

% Definição do esforço dos atuadores (Capacidade dos Motores)
tau1_max = 5.3;  
tau2_max = 41; 
f3_max   = 420;

% Matriz Q (5x5 para o sistema reduzido)
Q_lqr = diag([1/err_q2_max^2, 1/err_q3_max^2, ...
              1/err_dq1_max^2, 1/err_dq2_max^2, 1/err_dq3_max^2]);

% Matriz R (3x3 para os 3 atuadores)
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

%% CONTROLE POR ALOCAÇÃO DE POLOS

% Estratégia baseada nos polos do LQR
% ctrl_poles = poles_cl';

% Estratégia de amortecimento cirúrgico (Butterworth no tecido + Real nas juntas)
ctrl_poles = [-5, -7 + 3i, -7 - 3i, -18 + 12i, -18 - 12i];

disp('============================================');
disp(' PARÂMETROS DE CONTROLE - ALOCAÇÃO DE POLOS ');
disp('============================================');
disp('Polos desejados para a malha fechada:');
disp(ctrl_poles');

% Cálculo da matriz de ganhos usando o algoritmo de Place
F_pp = place(A, B, ctrl_poles);
% F_pp = F_pp/1.2;

disp('Matriz de ganhos F calculada:');
disp(F_pp);

disp('Processando análises em malha fechada para alocação de polos...');
% AnalyzeCL(A, B, C, D, E, F_pp, 'PP');


%% OBSERVADOR DE ORDEM COMPLETA - IDENTIDADE
   
% Abordagem LQR
Wc = 5000 * eye(n_states); % "ruído de processo"  — subir = observador mais rápido
Vc = eye(n_outputs);    % "ruído de medição"   — subir = observador mais lento
L = lqr(A', C', Wc, Vc).';

%Abordagem Alocação de Polos
% poles_obs_id = [-30+8i,-30+8i, -20, -25+2i, -25-2i];
% L = place(A',C',poles_obs_id)';

fprintf('Ganho do observador L:\n'); disp(L);
fprintf('Autovalores de (A - L*C):\n'); disp(eig(A - L*C).');

% Sistema do observador isolado (entradas = [u; y], saída = x_hat)
A_obs = A - L * C;
B_obs = [B, L];
C_obs = eye(n_states);
D_obs = zeros(n_states, n_inputs + n_outputs);
sys_obs = ss(A_obs, B_obs, C_obs, D_obs);

% pzmap(sys_obs)

%% OBSERVADOR DE ORDEM REDUZIDA

% Construção da transformação T e particionamento
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

% Observabilidade e condicionamento do par reduzido (A22, A12)
rank_obs_red = rank(obsv(A22, A12));
fprintf('Posto da observabilidade reduzida = %d (necessario n-m = %d)\n', ...
        rank_obs_red, n_states - n_outputs);
if rank_obs_red < n_states - n_outputs
    error('Par (A22, A12) NAO observavel: observador reduzido nao realizavel.');
end
fprintf('cond(obsv(A22,A12)) = %.3e   (se >> 1e6, o ganho J tende a crescer)\n', ...
        cond(obsv(A22, A12)));

% Ganho J via LQE no sistema fictício (A22', A12')
Qe_red = 8000 * eye(n_states - n_outputs); % "ruído de processo"  — subir = mais rápido
Re_red = eye(n_outputs);              % "ruído de medição"   — subir = mais lento
J = lqr(A22', A12', Qe_red, Re_red).';

% Ganho J via PP no sistema fictício (A22', A12')
% J = place(A22', A12', 3*ctrl_poles).';

% Matrizes do observador 
F_red = A22 - J * A12;
G_red = A21 - J * A11 + F_red * J;
H_red = B2  - J * B1;
S_red = M   + N * J;

fprintf('Polos do observador reduzido (eig F):\n'); disp(eig(F_red).');
fprintf('norm(J) = %.3e\n', norm(J));

% Sistema aumentado [x; z] em malha aberta 

% [dx]   [ A       0 ] [x]   [ B ]      [ E ]
% [dz] = [ G*C     F ] [z] + [ H ] u  + [ 0 ] w
A_aug_red = [A,         zeros(n_states, n_states - n_outputs);
             G_red * C, F_red];
B_aug_red = [B; H_red];
E_aug_red = [E; zeros(n_states - n_outputs, size(E, 2))];
C_aug_red = [C, zeros(n_outputs, n_states - n_outputs)];

% Reconstrução de x_hat a partir do estado aumentado:

% x_hat = S*C*x + N*z  =  recover_xhat * [x; z]
recover_xhat = [S_red * C, N];

% --- Inicialização correta de z (usar no bloco 3) ---------------------
% O reduzido lê y exatamente; só tem liberdade na parte nao-medida (w).
% Dada uma estimativa inicial das velocidades w_hat0 e o estado real x0:
%   w_hat = z + J*y   =>   z0 = w_hat0 - J*C*x0
% w_hat0 = 0  => observador "ignorante" das velocidades.
build_z0 = @(w_hat0, x0) w_hat0 - J * C * x0;

%% CONVERGENCIA DO ERRO DOS OBSERVADORES (MALHA ABERTA)

% Tempo de simulação
t_sim = 0 : 1e-3 : 2;

% Estado inicial REAL da planta
x0_real = [0.05; 0.02; 1; 3; 2];

% Estado inicial estimado 
x0_hat = zeros(n_states, 1);

% Erro inicial
e0 = x0_real - x0_hat;

% Compactação observador reduzido
obs_red.F_red        = F_red;
obs_red.A_aug_red    = A_aug_red;
obs_red.recover_xhat = recover_xhat;
obs_red.J            = J;

% Nomes dos estados
state_names = {'q_2','q_3','dq_1','dq_2','dq_3'};

Comp_Obsv(A, C, L, obs_red, x0_real, t_sim, state_names, x0_hat)

%% ANÁLISE DE CONVERGÊNCIA E SEPARAÇÃO DE POLOS

obs_red.F_red        = F_red;
obs_red.A_aug_red    = A_aug_red;
obs_red.recover_xhat = recover_xhat;
obs_red.J            = J;

state_labels  = {'Posicao q2 (rad)','Posicao q3 (m)', ...
           'Velocidade dq1 (rad/s)','Velocidade dq2 (rad/s)','Velocidade dq3 (m/s)'};
x0_real = [0.05; 0.02; 1; 3; 2];
x0_hat  = [0.05; 0.02; 0; 0; 0];

AnalyseCL_Obsv(A, C, L, obs_red, poles_cl, x0_real, t_sim, state_labels, x0_hat)

%% SIMULAÇÃO EM MALHA FECHADA COMPLETA (LQR + OBS REDUZIDO)

disp('==================================================');
disp(' INICIANDO SIMULAÇÃO EM MALHA FECHADA COMPLETA... ');
disp('==================================================');

% 1. Construção das Matrizes do Sistema Aumentado em Malha Fechada

% Combinando a dinâmica da planta real com a lei de controle u = -K_lqr * x_hat
A_cl_total = [ (A - B * K_lqr * S_red * C),            (-B * K_lqr * N);
               (G_red * C - H_red * K_lqr * S_red * C), (F_red - H_red * K_lqr * N) ];

% Entrada de perturbação: o sinal w entra apenas na planta física
B_cl_total = [ E; 
               zeros(n_states - n_outputs, size(E, 2)) ];

% Matriz C auxiliar para extrair todas as variáveis de estado do bloco ss
C_cl_total = eye(2 * n_states - n_outputs);
D_cl_total = zeros(2 * n_states - n_outputs, 1);

% Instanciação do sistema global de Malha Fechada
sys_cl_total = ss(A_cl_total, B_cl_total, C_cl_total, D_cl_total);

% 2. Configuração do Tempo e do Distúrbio Respiratório Senoidal
t_sim2 = 0 : 0.001 : 5;
f_respiracao = 0.5;    % 30 respirações por minuto
omega_w = 2 * pi * f_respiracao;
w_senoidal = 0.005 * sin(omega_w * t_sim2); % 5 mm de amplitude da respiração

% 3. Condições Iniciais Corretas (Planta desalinhada e Observador no escuro)
x0_real = [0.05; 0.02; 1; 3; 2]; 
z0_correto = build_z0(zeros(n_states - n_outputs, 1), x0_real);
x0_cl_total = [x0_real; z0_correto];

% 4. Execução da Simulação Temporal Combinada (Tranco + Respiração Contínua)
[X_cl_total, ~] = lsim(sys_cl_total, w_senoidal, t_sim2, x0_cl_total);

% Separação dos Estados Reais e dos Estados do Observador
x_real_sim2 = X_cl_total(:, 1:n_states);
z_sim2      = X_cl_total(:, n_states+1:end);

% Reconstrução das variáveis estimadas (x_hat) e dos esforços (u) no tempo
x_hat_sim2 = zeros(length(t_sim2), n_states);
u_sim2     = zeros(length(t_sim2), n_inputs);

for i = 1:length(t_sim2)
    % Reconstrói o vetor x_hat usando a equação de saída do observador reduzido
    x_hat_sim2(i, :) = (S_red * C * x_real_sim2(i, :)' + N * z_sim2(i, :)')';
    
    % Calcula o torque real injetado nas juntas: u = -K * x_hat
    u_sim2(i, :) = (-K_lqr * x_hat_sim2(i, :)')';
end

%  GERAÇÃO DOS GRÁFICOS

% --- GRÁFICO 1: COMPORTAMENTO DAS POSIÇÕES DAS JUNTAS ---
figure('Name', 'Malha Fechada Completa: Posicoes das Juntas', 'Color', 'w');

subplot(2,1,1);
plot(t_sim2, x_real_sim2(:,1), 'b', 'LineWidth', 2); grid on;
title('Resposta da Junta q2 (Elevacao) com Controle via Estado Estimado');
ylabel('Posicao Real (rad)');

subplot(2,1,2);
plot(t_sim2, x_real_sim2(:,2), 'b', 'LineWidth', 2); grid on;
title('Resposta da Junta q3 (Insercao) com Controle via Estado Estimado');
xlabel('Tempo (s)'); ylabel('Posicao Real (m)');

% --- GRÁFICO 2: COMPORTAMENTO DOS ESFORÇOS DE CONTROLE REALISTAS ---
figure('Name', 'Malha Fechada Completa: Esforcos de Controle', 'Color', 'w');
plot(t_sim2, u_sim2(:,1), 'b', 'LineWidth', 2); hold on;
plot(t_sim2, u_sim2(:,2), 'r', 'LineWidth', 2);
plot(t_sim2, u_sim2(:,3), 'g', 'LineWidth', 2); grid on;
title('Esforco de Controle dos Motores (Realimentacao por Observador Reduzido)');
xlabel('Tempo (s)'); ylabel('Torque (Nm) / Forca (N)');
legend('tau1 (Base)', 'tau2 (Elevacao)', 'F3 (Insercao)', 'Location', 'best');

disp('=> Simulação do Cenário 2 concluída e gráficos gerados!');