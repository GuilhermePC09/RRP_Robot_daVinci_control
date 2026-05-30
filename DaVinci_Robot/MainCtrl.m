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
%  8) VERIFICACOES RAPIDAS
%  =====================================================================
fprintf('\n--- Polos em malha fechada ---\n');
fprintf('LQR              : '); disp(eig(A - B*F_lqr).');
fprintf('Pole placement   : '); disp(eig(A - B*F_pp).');
fprintf('Observador PP    : '); disp(eig(A - L_pp*C).');
fprintf('Observador LQR   : '); disp(eig(A - L_lqr*C).');
fprintf('Integ. LQR (aum) : '); disp(eig(Ae - Be*Fe_lqr).');
fprintf('Integ. PP  (aum) : '); disp(eig(Ae - Be*Fe_pp).');
fprintf('Sistema completo : '); disp(eig(A_full).');
