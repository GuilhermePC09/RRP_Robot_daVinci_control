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

%% ======================
%  0) ANÁLISE DO SISTEMA
%  ======================

% Sistema reduzido
A = A_reduced;
B = B_reduced;
C = C_reduced;
D = D_reduced;
E = E_reduced;

n_states  = size(A, 1);
n_inputs  = size(B, 2);
n_outputs = size(C, 1);

% ControlAnalysis(A, B, C, D);

%% =====================================================================
%  1) OBSERVADOR DE ORDEM COMPLETA - IDENTIDADE
%  =====================================================================

% Wc = 35 * eye(n_states);     % "ruído de processo"  — subir = observador mais rápido
Wc = 5000 * eye(n_states); 
Vc = eye(n_outputs);    % "ruído de medição"   — subir = observador mais lento

L = lqr(A', C', Wc, Vc).';

fprintf('Ganho do observador L_pp:\n'); disp(L);
fprintf('Autovalores de (A - L_pp*C):\n'); disp(eig(A - L*C).');

% Sistema do observador isolado (entradas = [u; y], saída = x_hat)
A_obs = A - L * C;
B_obs = [B, L];
C_obs = eye(n_states);
D_obs = zeros(n_states, n_inputs + n_outputs);
sys_obs = ss(A_obs, B_obs, C_obs, D_obs);

% pzmap(sys_obs)
%% =====================================================================
%  2) OBSERVADOR DE ORDEM REDUZIDA (Método de Friedland)
%  =====================================================================
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
Qe_red = 35 * eye(n_states - n_outputs);   % "ruído de processo"  — subir = mais rápido
Qe_red = 8000 * eye(n_states - n_outputs);
Re_red = eye(n_outputs);              % "ruído de medição"   — subir = mais lento
J = lqr(A22', A12', Qe_red, Re_red).';

% --- Matrizes do observador --------------------------------------------
F_red = A22 - J * A12;
G_red = A21 - J * A11 + F_red * J;
H_red = B2  - J * B1;
S_red = M   + N * J;

fprintf('Polos do observador reduzido (eig F):\n'); disp(eig(F_red).');
fprintf('norm(J) = %.3e\n', norm(J));

% --- Sistema aumentado [x; z] em malha aberta -------------------------
% [dx]   [ A       0 ] [x]   [ B ]      [ E ]
% [dz] = [ G*C     F ] [z] + [ H ] u  + [ 0 ] w
A_aug_red = [A,         zeros(n_states, n_states - n_outputs);
             G_red * C, F_red];
B_aug_red = [B; H_red];
E_aug_red = [E; zeros(n_states - n_outputs, size(E, 2))];
C_aug_red = [C, zeros(n_outputs, n_states - n_outputs)];

% Reconstrução de x_hat a partir do estado aumentado:
%   x_hat = S*C*x + N*z  =  recover_xhat * [x; z]
recover_xhat = [S_red * C, N];

% --- Inicialização correta de z (usar no bloco 3) ---------------------
% O reduzido lê y exatamente; só tem liberdade na parte nao-medida (w).
% Dada uma estimativa inicial das velocidades w_hat0 e o estado real x0:
%   w_hat = z + J*y   =>   z0 = w_hat0 - J*C*x0
% w_hat0 = 0  => observador "ignorante" das velocidades.
build_z0 = @(w_hat0, x0) w_hat0 - J * C * x0;

%% =====================================================================
%  3) CONVERGENCIA DO ERRO DOS OBSERVADORES (MALHA ABERTA)
%  =====================================================================

% Tempo de simulação
t_sim = 0 : 1e-3 : 2;

% Estado inicial REAL da planta (paciente deslocado do equilíbrio)
x0_real = [0.05; 0.02; 1; 3; 2];

% Estado inicial estimado (sabemos só q2 e q3 pelas medidas - resto = 0)
x0_hat = zeros(n_states, 1);

% Erro inicial
e0 = x0_real - x0_hat;


% ---------------------------------------------------------------------
% TESTE 3.1: Convergência do erro com a planta em MALHA ABERTA
% (foca em comparar APENAS a dinâmica de erro dos observadores)
% ---------------------------------------------------------------------
% Dinâmica do erro: de/dt = (A - L*C) e   para observador de ordem completa
%                   dz_e/dt = F * z_e     para observador reduzido (em coord. z)

sys_err_pp = ss(A_obs, zeros(n_states,1), eye(n_states), zeros(n_states,1));

[e_pp, ~] = initial(sys_err_pp, e0, t_sim);

% Para o reduzido, o erro "verdadeiro" em x é:
%   e_x = x - x_hat = x - (S*C*x + N*z) = (I - S*C) x - N z
% Iniciando com z(0) = -V*x_hat(0) = 0 e x(0) = x0_real, simulamos a
% planta + observador reduzido em paralelo e comparamos no espaço dos
% estados originais.

% Sistema aumentado em MALHA ABERTA: u = 0, w = 0, partindo de
% [x0_real; z0]  com z0 = V*x0_hat - J*C*x0_hat = 0 (pois x0_hat = 0)
sys_aug_open = ss(A_aug_red, zeros(2*n_states - n_outputs, 1), ...
                  eye(2*n_states - n_outputs), zeros(2*n_states - n_outputs, 1));

z0 = build_z0(zeros(n_states - n_outputs, 1), x0_real);
[X_aug, ~] = initial(sys_aug_open, [x0_real; z0], t_sim);

% [X_aug, ~] = initial(sys_aug_open, [x0_real; zeros(n_states - n_outputs, 1)], t_sim);

x_real   = X_aug(:, 1:n_states);
z_red    = X_aug(:, n_states+1:end);
xhat_red = (recover_xhat * X_aug.').';
e_red    = x_real - xhat_red;

state_names = {'q_2', 'q_3', 'dq_1', 'dq_2', 'dq_3'};

for k = 1:n_states
    figure('Name', ['Erro de estimacao - ' state_names{k}], 'Color', 'w');
    plot(t_sim, e_pp(:, k),  'b-',  'LineWidth', 1.6); hold on;
    plot(t_sim, e_red(:, k), 'g-.', 'LineWidth', 1.6);
    grid on;
    xlabel('Tempo [s]');
    ylabel(['e_{' state_names{k} '}']);
    title(['Erro de estimacao: ' state_names{k} ' (malha aberta, sem entrada)']);
    legend('PP (ordem completa)', 'Reduzido', 'Location', 'best');
end


% ---------------------------------------------------------------------
% TESTE 3.2: Mapa de polos dos observadores
% ---------------------------------------------------------------------
figure('Name', 'Polos dos observadores', 'Color', 'w');
hold on; grid on;
plot(real(eig(A - L*C)), imag(eig(A - L*C)), 'bx', ...
     'MarkerSize', 12, 'LineWidth', 2);
plot(real(eig(F_red)),      imag(eig(F_red)),      'gs', ...
     'MarkerSize', 12, 'LineWidth', 2);
xline(0, 'k--'); yline(0, 'k--');
xlabel('Re'); ylabel('Im');
title('Polos: observador PP vs observador reduzido');
legend('Obs identidade (n=5)', 'Obs Reduzido (n=3)', 'Location', 'best');
axis equal;



%% =====================================================================
%  3) ANÁLISE DE CONVERGÊNCIA E SEPARAÇÃO DE POLOS - v2
%  =====================================================================
t_sim = 0 : 1e-3 : 2; % Janela de simulação de 2 segundos

% --- 3.1) TABELA COMPARATIVA DE POLOS (PRINCÍPIO DA SEPARAÇÃO) ---
poles_control = poles_cl; % Seus polos de malha fechada vindos do LQR/PP
poles_obs_full = eig(A - L*C);
poles_obs_red  = eig(F_red);

fprintf('\n===================================================================\n');
fprintf('        TABELA COMPARATIVA DE POLOS DO SISTEMA GLOBAL \n');
fprintf('===================================================================\n');
fprintf(' Modo |   Controlador (MF)    |   Obs. Completo   |   Obs. Reduzido   \n');
fprintf('-------------------------------------------------------------------\n');
for i = 1:5
    str_ctrl = sprintf('%+6.2f %+6.2fi', real(poles_control(i)), imag(poles_control(i)));
    str_full = sprintf('%+6.2f %+6.2fi', real(poles_obs_full(i)), imag(poles_obs_full(i)));
    
    if i <= 3
        str_red = sprintf('%+6.2f %+6.2fi', real(poles_obs_red(i)), imag(poles_obs_red(i)));
        fprintf('  %d   | %-21s | %-17s | %-17s \n', i, str_ctrl, str_full, str_red);
    else
        fprintf('  %d   | %-21s | %-17s | %-17s \n', i, str_ctrl, str_full, '       -       ');
    end
end
fprintf('===================================================================\n\n');


% --- 3.2) CONVERGÊNCIA DOS ESTADOS: PLANTA VS ESTIMADO ---
% Estado inicial REAL da planta (paciente fora do equilíbrio com velocidades)
x0_real = [0.05; 0.02; 1; 3; 2];

% Estado inicial estimado pelo computador (só conhece as posições pelos sensores)
x0_hat = [0.05; 0.02; 0; 0; 0]; 
e0 = x0_real - x0_hat;

% Simulação do erro do observador cheio para reconstruir seu xhat
sys_err_pp = ss(A_obs, zeros(n_states,1), eye(n_states), zeros(n_states,1));
[e_pp, ~] = initial(sys_err_pp, e0, t_sim);

% Correção da inicialização do estado z do observador reduzido
z0 = build_z0(zeros(n_states - n_outputs, 1), x0_real);

% Simulação do sistema aumentado em Malha Aberta (u=0, w=0)
sys_aug_open = ss(A_aug_red, zeros(2*n_states - n_outputs, 1), ...
                  eye(2*n_states - n_outputs), zeros(2*n_states - n_outputs, 1));
[X_aug, ~] = initial(sys_aug_open, [x0_real; z0], t_sim);

% Separação dos estados reais e estimados
x_real   = X_aug(:, 1:n_states);
xhat_red = (recover_xhat * X_aug.').';
xhat_pp  = x_real - e_pp; % Reconstrução matemática do xhat do observador cheio

% Nomes para formatação dos gráficos dinâmicos
state_labels = {'Posicao q2 (rad)', 'Posicao q3 (m)', 'Velocidade dq1 (rad/s)', 'Velocidade dq2 (rad/s)', 'Velocidade dq3 (m/s)'};
short_names  = {'q2', 'q3', 'dq1', 'dq2', 'dq3'};

for k = 1:n_states
    figure('Name', ['Convergencia de Estado - ' short_names{k}], 'Color', 'w');
    plot(t_sim, x_real(:, k),   'r-.',  'LineWidth', 2.5); hold on;
    plot(t_sim, xhat_pp(:, k),  'b--', 'LineWidth', 1.5);
    plot(t_sim, xhat_red(:, k), 'g-', 'LineWidth', 1.5);
    grid on;
    xlabel('Tempo (s)'); ylabel(state_labels{k});
    title(['Rastreamento e Convergencia do Estado: ' short_names{k}]);
    legend('Planta Real', 'Obs Ordem Completa', 'Obs Reduzido', 'Location', 'best');
end


% --- 3.3) MAPA GLOBAL DE POLOS EXPANDIDO ---
figure('Name', 'Mapa Global de Polos do Sistema', 'Color', 'w');
hold on; grid on;
plot(real(poles_control), imag(poles_control), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(poles_obs_full), imag(poles_obs_full), 'bs', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(real(poles_obs_red), imag(poles_obs_red), 'g+', 'MarkerSize', 8, 'LineWidth', 1.5);
xline(0, 'k--'); yline(0, 'k--');
xlabel('Eixo Real (1/s)'); ylabel('Eixo Imaginario (1/s)');
title('Mapa de Polos Global: Princípio da Separação');
legend('Controlador (Malha Fechada)', 'Obs Ordem Completa (n=5)', 'Obs Reduzido (n=3)', 'Location', 'best');




%% =====================================================================
%  CENÁRIO: SIMULAÇÃO EM MALHA FECHADA COMPLETA (LQR + OBS REDUZIDO)
%  =====================================================================
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

% =====================================================================
%  GERAÇÃO DOS GRÁFICOS
%  =====================================================================

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