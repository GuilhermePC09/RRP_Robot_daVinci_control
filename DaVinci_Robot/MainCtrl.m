clc;
%% Sistema reduzido
A = A_reduced;
B = B_reduced;
C = C_reduced;
D = D_reduced;

n = size(A, 1);     % numero de estados
m = size(B, 2);     % numero de entradas
p = size(C, 1);     % numero de saidas

ControlAnalysis(A, B, C, D, sys_reduced);

%% =====================================================================
%  1) CONTROLE POR LQR
%  =====================================================================
Q_lqr = eye(n);
R_lqr = eye(m);

F_lqr = lqr(A, B, Q_lqr, R_lqr);
G_lqr = inv(-C * inv(A - B*F_lqr) * B);

sys_lqr_bf = ss(A - B*F_lqr, B*G_lqr, C, zeros(p, m));

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