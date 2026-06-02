function ControlAnalysis(A, B, C, D)
% ControlAnalysis Performs a comprehensive open-loop modern control analysis
%
% Inputs: A, B, C, D - State-space numeric matrices

n_states  = size(A, 1);
n_inputs  = size(B, 2);
n_outputs = size(C, 1);
tol = 1e-8;

% Instantiating the state-space system object internally to prevent mismatches
sys_control = ss(A, B, C, D);

fprintf('=========================================================\n');
fprintf(' MODERN CONTROL ANALYSIS  -  daVinci ROBOT (linearized)  \n');
fprintf('=========================================================\n');
fprintf(' States  (n) = %d  -> [q2, q3, dq1, dq2, dq3]\n', n_states);
fprintf(' Inputs  (m) = %d  -> [Tau_1, Tau_2, Force_3]\n', n_inputs);
fprintf(' Outputs (p) = %d  -> [q2, q3]\n\n',              n_outputs);

disp('A ='); disp(A);
disp('B ='); disp(B);
disp('C ='); disp(C);


%% 1) STABILITY ANALYSIS
fprintf('---------------------------------------------------------\n');
fprintf(' 1) STABILITY ANALYSIS\n');
fprintf('---------------------------------------------------------\n');

eig_A  = eig(A);
re_eig = real(eig_A);

poles_reduced = pole(sys_control);
zeros_reduced = tzero(sys_control);

disp('Poles of the system:'); disp(poles_reduced);
% disp('Zeros of the system:'); disp(zeros_reduced);

if all(re_eig < -tol)
    asym_stable = true;
    fprintf('=> ASYMPTOTICALLY STABLE.\n\n');
elseif any(re_eig > tol)
    asym_stable = false;
    fprintf('=> UNSTABLE.\n\n');
else
    asym_stable = false;
    fprintf('=> MARGINALLY STABLE.\n\n');
end

% Lyapunov check execution (only if strictly asymptotically stable)
if asym_stable
    Q = eye(n_states);
    P = lyap(A, Q);
    if all(eig(P) > 0)
        fprintf('Lyapunov check: P = P^T > 0 found (consistent with asymptotic stability).\n\n');
    else
        fprintf('Lyapunov check: P is not positive definite (numerical issue).\n\n');
    end
end

figure('Name', 'Open-loop poles and zeros', 'Color', 'w');
pzmap(sys_control); grid on;


%% 2) CONTROLLABILITY
fprintf('---------------------------------------------------------\n');
fprintf(' 2) CONTROLLABILITY\n');
fprintf('---------------------------------------------------------\n');

Mc = ctrb(A, B);
rank_Mc = rank(Mc);

fprintf('Rank of Controllability Matrix (Mc) = %d  (Required n = %d)\n', rank_Mc, n_states);

if rank_Mc == n_states
    fully_controllable = true;
    fprintf('=> COMPLETELY CONTROLLABLE.\n\n');
else
    fully_controllable = false;
    fprintf('=> NOT completely controllable. Uncontrollable modes: %d\n\n', n_states - rank_Mc);
end

fprintf('PBH test:\n');
fprintf('  %-28s %-10s %-15s\n', 'Eigenvalue', 'rank', 'status');
for k = 1:length(eig_A)
    lam = eig_A(k);
    r   = rank([lam*eye(n_states) - A, B]);
    status = 'controllable';
    if r ~= n_states, status = 'UNCONTROLLABLE'; end
    fprintf('  %-28s %-10d %-15s\n', sprintf('%+8.4f %+8.4fi', real(lam), imag(lam)), r, status);
end
fprintf('\n');


%% 3) CONTROLLABLE / UNCONTROLLABLE DECOMPOSITION
fprintf('---------------------------------------------------------\n');
fprintf(' 3) CONTROLLABILITY DECOMPOSITION\n');
fprintf('---------------------------------------------------------\n');

[Abar_c, Bbar_c, Cbar_c, ~, k_c] = ctrbf(A, B, C);
n_ctrb   = sum(k_c);
n_unctrb = n_states - n_ctrb;

fprintf('Controllable states  : %d\n',   n_ctrb);
fprintf('Uncontrollable states: %d\n\n', n_unctrb);

disp('Abar:'); disp(Abar_c);
disp('Bbar:'); disp(Bbar_c);
disp('Cbar:'); disp(Cbar_c);

if n_unctrb > 0
    Anc = Abar_c(1:n_unctrb,     1:n_unctrb);
    Ac  = Abar_c(n_unctrb+1:end, n_unctrb+1:end);
    Bc  = Bbar_c(n_unctrb+1:end, :);
    fprintf('Eigenvalues of Anc:\n'); disp(eig(Anc));
    fprintf('Eigenvalues of Ac:\n');  disp(eig(Ac));
else
    fprintf('All eigenvalues are controllable.\n\n');
end


%% 4) STABILIZABILITY
fprintf('---------------------------------------------------------\n');
fprintf(' 4) STABILIZABILITY\n');
fprintf('---------------------------------------------------------\n');

if n_unctrb == 0
    stabilizable = true;
    fprintf('=> STABILIZABLE (fully controllable).\n\n');
else
    eig_unctrb = eig(Anc);
    if all(real(eig_unctrb) < -tol)
        stabilizable = true;
        fprintf('=> STABILIZABLE. Uncontrollable modes:\n'); disp(eig_unctrb);
    else
        stabilizable = false;
        fprintf('=> NOT STABILIZABLE. Unstable uncontrollable modes:\n');
        disp(eig_unctrb(real(eig_unctrb) >= -tol));
    end
end


%% 5) OBSERVABILITY
fprintf('---------------------------------------------------------\n');
fprintf(' 5) OBSERVABILITY\n');
fprintf('---------------------------------------------------------\n');

Mo = obsv(A, C);
rank_Mo = rank(Mo);

fprintf('Rank of Observability Matrix (Mo) = %d  (Required n = %d)\n', rank_Mo, n_states);

if rank_Mo == n_states
    fully_observable = true;
    fprintf('=> COMPLETELY OBSERVABLE.\n\n');
else
    fully_observable = false;
    fprintf('=> NOT completely observable. Unobservable modes: %d\n\n', n_states - rank_Mo);
end

fprintf('PBH test:\n');
fprintf('  %-28s %-10s %-15s\n', 'Eigenvalue', 'rank', 'status');
for k = 1:length(eig_A)
    lam = eig_A(k);
    r   = rank([lam*eye(n_states) - A; C]);
    status = 'observable';
    if r ~= n_states, status = 'UNOBSERVABLE'; end
    fprintf('  %-28s %-10d %-15s\n', sprintf('%+8.4f %+8.4fi', real(lam), imag(lam)), r, status);
end
fprintf('\n');


%% 6) OBSERVABLE / UNOBSERVABLE DECOMPOSITION
fprintf('---------------------------------------------------------\n');
fprintf(' 6) OBSERVABILITY DECOMPOSITION\n');
fprintf('---------------------------------------------------------\n');

[Abar_o, Bbar_o, Cbar_o, ~, k_o] = obsvf(A, B, C);

n_obs   = sum(k_o);
n_unobs = n_states - n_obs;

fprintf('Observable states  : %d\n',   n_obs);
fprintf('Unobservable states: %d\n\n', n_unobs);

disp('Abar:'); disp(Abar_o);
disp('Bbar:'); disp(Bbar_o);
disp('Cbar:'); disp(Cbar_o);

if n_unobs > 0
    Ano = Abar_o(1:n_unobs,     1:n_unobs);
    Ao  = Abar_o(n_unobs+1:end, n_unobs+1:end);
    Co  = Cbar_o(:,             n_unobs+1:end);
    fprintf('Eigenvalues of Ano:\n'); disp(eig(Ano));
    fprintf('Eigenvalues of Ao:\n');  disp(eig(Ao));
else
    fprintf('All eigenvalues are observable.\n\n');
end


%% 7) DETECTABILITY
fprintf('---------------------------------------------------------\n');
fprintf(' 7) DETECTABILITY\n');
fprintf('---------------------------------------------------------\n');

if n_unobs == 0
    detectable = true;
    fprintf('=> DETECTABLE (fully observable).\n\n');
else
    eig_unobs = eig(Ano);
    if all(real(eig_unobs) < -tol)
        detectable = true;
        fprintf('=> DETECTABLE. Unobservable modes:\n'); disp(eig_unobs);
    else
        detectable = false;
        fprintf('=> NOT DETECTABLE. Unstable unobservable modes:\n');
        disp(eig_unobs(real(eig_unobs) >= -tol));
    end
end


%% SUMMARY
fprintf('=========================================================\n');
fprintf(' SUMMARY\n');
fprintf('=========================================================\n');
fprintf('  Asymptotically stable : %s\n',                         yesno(asym_stable));
fprintf('  Controllable          : %s   (rank Mc = %d / %d)\n',  yesno(fully_controllable), rank_Mc, n_states);
fprintf('  Stabilizable          : %s\n',                         yesno(stabilizable));
fprintf('  Observable            : %s   (rank Mo = %d / %d)\n',  yesno(fully_observable),   rank_Mo, n_states);
fprintf('  Detectable            : %s\n',                         yesno(detectable));
fprintf('=========================================================\n');

end


%% Helper
function s = yesno(cond)
    if cond, s = 'YES'; else, s = 'NO '; end
end
