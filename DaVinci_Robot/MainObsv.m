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
%  1) OBSERVADOR DE ORDEM COMPLETA
%  =====================================================================

Wc = 35 * eye(n_states);     % "ruído de processo"  — subir = observador mais rápido
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

pzmap(sys_obs)
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
[X_aug, ~] = initial(sys_aug_open, [x0_real; zeros(n_states - n_outputs, 1)], t_sim);

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


