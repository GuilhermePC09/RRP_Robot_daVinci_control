% LIMPEZA 
clc;


%% =====================================================================
%  0) ANÁLSIE DO SISTEMA
%  =====================================================================

% Sistema reduzido
A = A_reduced;
B = B_reduced;
C = C_reduced;
D = D_reduced;
E = E_reduced;

n = size(A, 1);     % numero de estados
m = size(B, 2);     % numero de entradas
p = size(C, 1);     % numero de saidas

ControlAnalysis(A, B, C, D);

%% =====================================================================
%  1) CONTROLE POR LQR
%  =====================================================================

% --- Definição dos limites máximos toleráveis (Abordagem de Bryson) ---
err_q2_max = 0.001;  % 0.001 rad de erro tolerável
err_q3_max = 0.001;  % 0.001 m (1 mm) de erro tolerável

err_dq1_max = 0.1;   % 0.1 rad/s tolerável
err_dq2_max = 0.1;   % 0.1 rad/s tolerável
err_dq3_max = 0.1;   % 0.1 m/s tolerável

% Matriz Q (5x5 para o sistema reduzido)
Q_lqr = diag([1/err_q2_max^2, 1/err_q3_max^2, ...
              1/err_dq1_max^2, 1/err_dq2_max^2, 1/err_dq3_max^2]);

% --- Definição do esforço dos atuadores (Capacidade dos Motores) ---
tau1_max = 15;  % Exemplo: Motor aguenta 10 Nm
tau2_max = 37;  % Exemplo: Motor aguenta 10 Nm
f3_max   = 424;  % Exemplo: Motor linear aguenta 50 N

% Matriz R (3x3 para os 3 atuadores)
% R_lqr = diag([1/tau1_max^2, 1/tau2_max^2, 1/f3_max^2]);
R_lqr = diag([100, 100, 10]);

[K_lqr, P_matrix, poles_cl] = lqr(A, B, Q_lqr, R_lqr);

disp('Matriz de Ganhos K calculada com sucesso:');
disp(K_lqr);

disp('Novos Polos em Malha Fechada:');
disp(poles_cl);

G_lqr = pinv(-C * inv(A - B * K_lqr) * B);
sys_lqr_bf = ss(A - B * K_lqr, B * G_lqr, C, zeros(p, p));


%  VISUALIZAÇÃO E COMPARAÇÃO DE DESEMPENHO (MIMO STEP)
t = 0:0.001:0.5; 
[y, t] = step(sys_lqr_bf, t); 

figure('Name', 'Rastreamento de Referência em Malha Fechada (LQR)', 'Color', 'w');

% Saída q2 devido à Ref q2
subplot(2,2,1);
plot(t, y(:,1,1), 'b', 'LineWidth', 2); grid on;
title('Ref q_2 \rightarrow Posição q_2'); 
ylabel('Ângulo (rad)');

% Saída q3 devido à Ref q2
subplot(2,2,3);
plot(t, y(:,2,1), 'b', 'LineWidth', 2); grid on;
title('Ref q_2 \rightarrow Posição q_3'); 
ylabel('Deslocamento (m)'); xlabel('Tempo (s)');

% Saída q2 devido à Ref q3
subplot(2,2,2);
plot(t, y(:,1,2), 'r', 'LineWidth', 2); grid on;
title('Ref q_3 \rightarrow Posição q_2');

% Saída q3 devido à Ref q3
subplot(2,2,4);
plot(t, y(:,2,2), 'r', 'LineWidth', 2); grid on;
title('Ref q_3 \rightarrow Posição q_3'); 
xlabel('Tempo (s)');

% Plot comparativo de rejeição de perturbação
x0 = [0.05; 0.02; 0; 0; 0]; % Robô deslocado do alvo
[y_ol, t_ol] = initial(ss(A, B, C, D), x0, 4);   % Resposta em Malha Aberta
[y_cl, t_cl] = initial(sys_lqr_bf, x0, 4);     % Resposta com LQR

figure('Name', 'Rejeição de Perturbação: Malha Aberta vs LQR', 'Color', 'w');
subplot(2,1,1);
plot(t_ol, y_ol(:,1), 'r--', 'LineWidth', 1.5); hold on;
plot(t_cl, y_cl(:,1), 'b', 'LineWidth', 2);
grid on; legend('Malha Aberta', 'LQR Controlado');
title('Recuperação da Junta q_2 após distúrbio'); ylabel('Ângulo (rad)');

subplot(2,1,2);
plot(t_ol, y_ol(:,2), 'r--', 'LineWidth', 1.5); hold on;
plot(t_cl, y_cl(:,2), 'b', 'LineWidth', 2);
grid on;
title('Recuperação da Junta q_3 após distúrbio'); ylabel('Deslocamento (m)');
xlabel('Tempo (s)');


%  CASO ESCALONAMENTO REALISTA

% % Q focado na precisão
% Q_lqr = diag([1e6, 1e6, 100, 100, 100]);
% 
% % 2. AUMENTAR R (Penaliza severamente o uso excessivo de torque/força)
% % Subir esses valores faz o LQR "economizar" os motores, reduzindo o pico de torque
% R_lqr = diag([100, 100, 10]);
% 
% % Re-calculando o ganho ótimo e o pré-filtro
% [K_lqr, ~, ~] = lqr(A, B, Q_lqr, R_lqr);
% G_lqr = pinv(-C * inv(A - B * K_lqr) * B);

% 3. DEFINIÇÃO DE AMPLITUDES CIRÚRGICAS REAIS (Substituindo o degrau de 1 metro)
amp_q2 = 5 * pi / 180;  % 5 graus de rotação (em radianos)
amp_q3 = 0.01;          % 1 centímetro de penetração (em metros)

% Matriz de escala para aplicar as amplitudes reais nas entradas do sistema MIMO
Matriz_Escala = diag([amp_q2, amp_q3]);
sys_lqr_bf = ss(A - B * K_lqr, B * G_lqr * Matriz_Escala, C, zeros(2,2));


%  EXTRAÇÃO E PLOT DO ESFORÇO DE CONTROLE

t = 0:0.001:1;
[~, ~, x_traj] = step(sys_lqr_bf, t); 

u_degrau_q2 = zeros(length(t), m);
u_degrau_q3 = zeros(length(t), m);

for idx = 1:length(t)
    % Cenário 1: Cirurgião moveu 5 graus em q2 (Ref = [amp_q2; 0])
    x_t_q2 = squeeze(x_traj(idx, :, 1))';
    u_degrau_q2(idx, :) = (-K_lqr * x_t_q2 + G_lqr * [amp_q2; 0])';
    
    % Cenário 2: Cirurgião afundou 1 cm em q3 (Ref = [0; amp_q3])
    x_t_q3 = squeeze(x_traj(idx, :, 2))';
    u_degrau_q3(idx, :) = (-K_lqr * x_t_q3 + G_lqr * [0; amp_q3])';
end

% --- GERANDO O PLOT ---
figure('Name', 'Esforco de Controle Realista (LQR Suave)', 'Color', 'w');

subplot(2,1,1);
plot(t, u_degrau_q2(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t, u_degrau_q2(:,2), 'r-', 'LineWidth', 2);
plot(t, u_degrau_q2(:,3), 'g-', 'LineWidth', 2); grid on;
title('Esforço de Controle: Comando de 5\circ em q_2');
ylabel('Torque (Nm) / Força (N)');
legend('\tau_1 (Base)', '\tau_2 (Elevação)', 'F_3 (Inserção)');

subplot(2,1,2);
plot(t, u_degrau_q3(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t, u_degrau_q3(:,2), 'r-', 'LineWidth', 2);
plot(t, u_degrau_q3(:,3), 'g-', 'LineWidth', 2); grid on;
title('Esforço de Controle: Comando de 1cm em q_3');
xlabel('Tempo (s)');
ylabel('Torque (Nm) / Força (N)');

%% =====================================================================
%  2) CONTROLE POR ALOCACAO DE POLOS
%  =====================================================================
% NOTA: place exige multiplicidade <= m. Com todos em -1 e m<n, vai
%       falhar -- ajuste os polos antes de rodar.
ctrl_poles = -1 * ones(1, n);

F_pp = place(A, B, ctrl_poles);
G_pp = inv(-C * inv(A - B*F_pp) * B);

sys_pp_bf = ss(A - B*F_pp, B*G_pp, C, zeros(p, m));

%% =====================================================================
%  3) OBSERVADOR POR ALOCACAO DE POLOS
%  =====================================================================
% NOTA: mesma restricao de multiplicidade (limite = p).
obs_poles = -1 * ones(1, n);

L_pp = place(A', C', obs_poles)';

A_obs_pp = A - L_pp * C;
B_obs_pp = [B, L_pp];
C_obs_pp = eye(n);
D_obs_pp = zeros(n, m + p);

sys_obs_pp = ss(A_obs_pp, B_obs_pp, C_obs_pp, D_obs_pp);

%% =====================================================================
%  4) OBSERVADOR LQR
%  =====================================================================
W_obs = eye(n);
V_obs = eye(p);

L_lqr = lqr(A', C', W_obs, V_obs)';

A_obs_lqr = A - L_lqr * C;
B_obs_lqr = [B, L_lqr];
C_obs_lqr = eye(n);
D_obs_lqr = zeros(n, m + p);

sys_obs_lqr = ss(A_obs_lqr, B_obs_lqr, C_obs_lqr, D_obs_lqr);

%% =====================================================================
%  5) CONTROLE COM ACAO INTEGRADORA -- LQR
%  =====================================================================
%  xa = [x; xi],  xi_dot = r - y
%  Ae = [ A   0 ]   Be = [ B ]   Ce = [ C  0 ]
%       [-C   0 ]        [ 0 ]
%  u = -F xhat + H xi,  [F  -H] = lqr(Ae, Be, Qe, Re)

Ae = [A,   zeros(n, p);
      -C,  zeros(p, p)];
Be = [B; zeros(p, m)];
Ce = [C, zeros(p, p)];

Qe_lqr = eye(n + p);
Re_lqr = eye(m);

Fe_lqr = lqr(Ae, Be, Qe_lqr, Re_lqr);
F_i_lqr = Fe_lqr(:, 1:n);
H_i_lqr = -Fe_lqr(:, n+1:end);

A_cl_int_lqr = [A - B*F_i_lqr,  B*H_i_lqr;
                -C,              zeros(p, p)];
B_cl_int_lqr = [zeros(n, p); eye(p)];
C_cl_int_lqr = [C, zeros(p, p)];
D_cl_int_lqr = zeros(p, p);

sys_int_lqr = ss(A_cl_int_lqr, B_cl_int_lqr, C_cl_int_lqr, D_cl_int_lqr);

%% =====================================================================
%  6) CONTROLE COM ACAO INTEGRADORA -- ALOCACAO DE POLOS
%  =====================================================================
% NOTA: place no sistema aumentado limita multiplicidade a m.
int_poles = -1 * ones(1, n + p);

Fe_pp = place(Ae, Be, int_poles);
F_i_pp = Fe_pp(:, 1:n);
H_i_pp = -Fe_pp(:, n+1:end);

A_cl_int_pp = [A - B*F_i_pp,  B*H_i_pp;
               -C,             zeros(p, p)];
B_cl_int_pp = [zeros(n, p); eye(p)];
C_cl_int_pp = [C, zeros(p, p)];
D_cl_int_pp = zeros(p, p);

sys_int_pp = ss(A_cl_int_pp, B_cl_int_pp, C_cl_int_pp, D_cl_int_pp);

%% =====================================================================
%  7) OBSERVACAO COM ACAO INTEGRADORA (planta + observador + integrador)
%  =====================================================================
%  xc = [x; xhat; xi]  -> 2n + p estados, entrada = r

F_use = F_i_lqr;     % ou F_i_pp
H_use = H_i_lqr;     % ou H_i_pp
L_use = L_lqr;       % ou L_pp

A_full = [ A,         -B*F_use,                B*H_use;
           L_use*C,    A - B*F_use - L_use*C,  B*H_use;
          -C,          zeros(p, n),            zeros(p, p)];
B_full = [zeros(n, p); zeros(n, p); eye(p)];
C_full = [C, zeros(p, n), zeros(p, p)];
D_full = zeros(p, p);

sys_full = ss(A_full, B_full, C_full, D_full);

%% =====================================================================
%  VERIFICACOES RAPIDAS
%  =====================================================================
fprintf('\n--- Polos em malha fechada ---\n');
fprintf('LQR              : '); disp(eig(A - B*F_lqr).');
fprintf('Pole placement   : '); disp(eig(A - B*F_pp).');
fprintf('Observador PP    : '); disp(eig(A - L_pp*C).');
fprintf('Observador LQR   : '); disp(eig(A - L_lqr*C).');
fprintf('Integ. LQR (aum) : '); disp(eig(Ae - Be*Fe_lqr).');
fprintf('Integ. PP  (aum) : '); disp(eig(Ae - Be*Fe_pp).');
fprintf('Sistema completo : '); disp(eig(A_full).');


%% =====================================================================
%  8) SEGUIDOR LQ COM PRÉ-ALIMENTAÇÃO (Módulo 8 - Seção 1)
%  Lei de controle: u(t) = K*e(t) + R^{-1}*B'*(eta - P*xr)
%  onde e(t) = xr(t) - x(t)  e  K = R^{-1}*B'*P  (mesmo K do LQR)
%  eta satisfaz: d(eta)/dt = -(A - B*K)'*eta - Q*xr,  eta(t1) = Q1*xr(t1)
%  =====================================================================
%  Reutiliza F_lqr e as matrizes Q_lqr, R_lqr do item 1.
%  Implementação: integração BACKWARD de eta, depois FORWARD de x.

Q_seg   = Q_lqr;          % mesma ponderação de estados
Q1_seg  = Q_lqr;          % penalização terminal (pode ajustar)
R_seg   = R_lqr;          % mesma ponderação de entradas

% --- Ganho do seguidor (MESMO K do regulador LQR) ---
[F_seg, P_seg, ~] = lqr(A, B, Q_seg, R_seg);
% F_seg == F_lqr  (confirmação)

% --- Horizonte de simulação ---
t1   = 10;                         % tempo final [s]
dt   = 1e-3;                       % passo de integração [s]
t_fw = 0 : dt : t1;                % vetor de tempo forward
N_t  = length(t_fw);

% --- Referência: degrau unitário em todos os n estados ---
%     Ajuste xr_func conforme o sinal de referência desejado
xr_func = @(t) ones(n, 1);        % referência constante (degrau)

% --- Integração BACKWARD de eta ---
%  d(eta)/dt = -(A - B*F_seg)'*eta - Q_seg*xr,  eta(t1) = Q1_seg*xr(t1)
A_cl_seg = A - B * F_seg;
eta      = zeros(n, N_t);
eta(:, end) = Q1_seg * xr_func(t1);

for k = N_t-1 : -1 : 1
    t_k       = t_fw(k+1);
    xr_k      = xr_func(t_k);
    deta      = -A_cl_seg' * eta(:, k+1) - Q_seg * xr_k;
    eta(:, k) = eta(:, k+1) - dt * deta;   % Euler backward
end

% --- Integração FORWARD de x ---
%  dx/dt = A*x + B*u,  u = K*e + R^{-1}*B'*(eta - P*xr)
x_seg = zeros(n, N_t);
u_seg = zeros(m, N_t);
x_seg(:, 1) = zeros(n, 1);        % condição inicial (ajuste se necessário)

for k = 1 : N_t-1
    xr_k      = xr_func(t_fw(k));
    e_k       = xr_k - x_seg(:, k);
    ff_k      = R_seg \ (B' * (eta(:, k) - P_seg * xr_k));  % pré-alimentação
    u_seg(:,k) = F_seg * e_k + ff_k;
    dx        = A * x_seg(:, k) + B * u_seg(:, k);
    x_seg(:, k+1) = x_seg(:, k) + dt * dx;
end
u_seg(:, end) = u_seg(:, end-1);  % repete último ponto

fprintf('\n--- Seguidor LQ (Módulo 8 - Seção 1) ---\n');
fprintf('Polos malha fechada seguidor : ');
disp(eig(A_cl_seg).');

%% =====================================================================
%  9) SEGUIDOR COM MODELO DE VARIÁVEIS EXÓGENAS (Módulo 8 - Seção 2)
%  Referência: dxr/dt = Ar*xr  |  perturbação: dw/dt = Aw*w
%  Lei de controle: u = K*e - N*F_exo*xo
%    com F_exo = [A - Ar, E],  xo = [xr; w]
%    e   N = inv(M * inv(A) * B) * M * inv(A),  M escolhida pelo usuário
%  =====================================================================
%  Parâmetros do modelo exógeno — ajuste conforme o seu sistema.
%  Aqui: referência constante (Ar = 0) e sem perturbação modelada.

Ar  = zeros(n, n);             % modelo da referência (0 → degrau)
E   = zeros(n, m);             % matriz de entrada da perturbação
Aw  = zeros(m, m);             % modelo da perturbação (0 → constante)

% Variável exógena aumentada: xo = [xr; w]
r_xo = n + m;                  % dimensão de xo

% Matriz F_exo: agrupa (A - Ar) e E
F_exo = [A - Ar, E];           % n × (n + m)

% Escolha de M: para seguir saídas medidas, use M = C.
%   M deve ser r×n com r = número de referências rastreáveis (r ≤ m).
%   Como m = 3 e p = 2, usamos as p primeiras linhas de C.
M = C;                         % p × n  (ajuste se quiser rastrear outros sinais)

% Ganho de pré-alimentação N = inv(M*inv(A)*B) * M*inv(A)
MiAB = M * (A \ B);            % p × m
if rank(MiAB) < size(MiAB, 1)
    warning(['M*inv(A)*B nao e invertivel. ' ...
             'Nao e possivel rastrear todos os sinais de referencia com os atuadores disponiveis.']);
    N_ff = pinv(MiAB) * M / A; % pseudo-inversa como fallback
else
    N_ff = MiAB \ (M / A);     % p × n
end

% Lei de controle: u = K*e - N_ff * F_exo * xo
% Reutiliza F_lqr como ganho de realimentação K
K_exo = F_lqr;                 % 3 × n

fprintf('\n--- Seguidor c/ variáveis exógenas (Módulo 8 - Seção 2) ---\n');
fprintf('Dimensão de N_ff        : %d × %d\n', size(N_ff));
fprintf('Dimensão de F_exo       : %d × %d\n', size(F_exo));
fprintf('Polos malha fechada     : ');
disp(eig(A - B * K_exo).');

% Verificação do erro em regime permanente
%   e_inf = inv(A) * (F_exo - B*N_ff*F_exo) * xo_inf
%   Para xo constante (Ar=0, Aw=0), xo_inf = [xr_inf; w_inf]
fprintf('Verificando cancelamento do erro em regime permanente...\n');
residuo = M / A * (F_exo - B * N_ff * F_exo);
fprintf('||M*inv(A)*(F_exo - B*N*F_exo)|| = %.2e  (esperado ~0)\n', norm(residuo));