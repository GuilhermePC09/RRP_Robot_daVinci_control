function AnalyzeCL(A, B, C, D, E, K, controller_name)
% AnalyzeCL Performs a refined closed-loop analysis of a regulator
%
% Inputs: A, B, C, D, E   - State-space matrices of the system
%         K               - Feedback gain matrix (u = -Kx)
%         controller_name - String with the controller name ('LQR' or 'PP')

n_states  = size(A, 1);
n_inputs  = size(B, 2);
n_outputs = size(C, 1);



% 1) POLE MAP: NEW POLES COMPARED TO OLD POLES

poles_ma = eig(A);
poles_cl = eig(A - B * K);

figure('Name', ['Comparativo de Polos - ' controller_name], 'Color', 'w');
plot(real(poles_ma), imag(poles_ma), 'b+', 'MarkerSize', 8, 'LineWidth', 2); hold on;
plot(real(poles_cl), imag(poles_cl), 'rx', 'MarkerSize', 8, 'LineWidth', 2);
grid on; xline(0, 'k--'); yline(0, 'k--');
title(['Mapeamento de polos: MA vs MF (' controller_name ')']);
xlabel('Eixo Real (1/s)');
ylabel('Eixo Imaginario (1/s)');
legend('MA', ['MF (' controller_name ')'], 'Location', 'best');



% 2) POSITION REGULATOR SCENARIO (INITIAL CONDITION RESPONSE)

% Loop configuration for initial condition
sys_ma_init = ss(A, zeros(size(B)), C, zeros(size(D)));
sys_cl_init = ss(A - B * K, zeros(size(B)), C, zeros(size(D)));

t = 0:0.001:2;  % Time vector
x0 = [0.05; 0.02; 0; 0; 0]; % Initial deviation from equilibrium

[y_ma_init, ~, ~]        = initial(sys_ma_init, x0, t);
[y_cl_init, ~, x_cl_init] = initial(sys_cl_init, x0, t);

% 2.1) Position plot (comparing CL with OL)
figure('Name', ['Regulacao de Posicao - Reg. Inicial - ' controller_name], 'Color', 'w');

subplot(2,1,1);
plot(t, y_ma_init(:,1), 'r--', 'LineWidth', 1.5); hold on;
plot(t, y_cl_init(:,1), 'b', 'LineWidth', 2); grid on;
title(['Regulação da junta q_2 | Resposta à condição inicial (' controller_name ')']);
ylabel('Desvio angular (rad)');
legend('MA', 'MF', 'Location', 'best');

subplot(2,1,2);
plot(t, y_ma_init(:,2), 'r--', 'LineWidth', 1.5); hold on;
plot(t, y_cl_init(:,2), 'b', 'LineWidth', 2); grid on;
title(['Regulação da junta q_3 | Resposta à condição inicial (' controller_name ')']);
xlabel('Tempo (s)'); ylabel('Desvio linear (m)');

% 2.2) Control effort plot (u = -Kx)
u_init = zeros(length(t), n_inputs);
for idx = 1:length(t)
    u_init(idx, :) = (-K * x_cl_init(idx, :)')';
end

figure('Name', ['Esforco de Controle - Reg. Inicial - ' controller_name], 'Color', 'w');
yyaxis left
plot(t, u_init(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t, u_init(:,2), 'r-', 'LineWidth', 2);
ylabel('Torque dos motores (Nm)');
set(gca, 'YColor', 'k');
yyaxis right
plot(t, u_init(:,3), 'g-', 'LineWidth', 2); grid on;
ylabel('Força de inserção (N)');
set(gca, 'YColor', 'k');
title(['Esforço de controle do regulador (' controller_name ')']);
xlabel('Tempo (s)');
legend('\tau_1', '\tau_2', 'F_3', 'Location', 'best');


%  3) SINUSOIDAL DISTURBANCE SCENARIO
t = 0:0.001:5;  % Time vector

% Definition of the breathing signal: 0.5 Hz, amplitude of 5 mm (0.005 m)
f_respiracao = 0.5;
omega_w = 2 * pi * f_respiracao;
w_senoidal = 0.005 * sin(omega_w * t);

% Systems for disturbance analysis (Input: w, Output: y)
sys_dist_ma = ss(A, E, C, zeros(n_outputs, 1));
sys_dist_cl = ss(A - B * K, E, C, zeros(n_outputs, 1));

% Simulação temporal das posições
y_dist_ma = lsim(sys_dist_ma, w_senoidal, t);
[y_dist_cl, ~, x_dist_cl] = lsim(sys_dist_cl, w_senoidal, t);

% 3.1) Position plot (comparing CL with OL)
figure('Name', ['Regulacao de Posicao - Dist. Respiratorio - ' controller_name], 'Color', 'w');

subplot(2,1,1);
plot(t, y_dist_ma(:,1), 'r--', 'LineWidth', 1.5); hold on;
plot(t, y_dist_cl(:,1), 'b', 'LineWidth', 2); grid on;
title(['Desvio da junta q_2 sob perturbação (' controller_name ')']);
ylabel('Erro angular (rad)');
legend('MA', 'MF');

subplot(2,1,2);
plot(t, y_dist_ma(:,2), 'r--', 'LineWidth', 1.5); hold on;
plot(t, y_dist_cl(:,2), 'b', 'LineWidth', 2); grid on;
title(['Desvio da junta q_3 sob perturbação (' controller_name ')']);
xlabel('Tempo (s)'); ylabel('Erro linear (m)');

% 3.2) Control effort plot
u_dist = zeros(length(t), n_inputs);
for idx = 1:length(t)
    u_dist(idx, :) = (-K * x_dist_cl(idx, :)')';
end

figure('Name', ['Esforco de Controle - Dist. Respiratorio - ' controller_name], 'Color', 'w');
yyaxis left
plot(t, u_dist(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t, u_dist(:,2), 'r-', 'LineWidth', 2);
ylabel('Torque dos motores (Nm)');
set(gca, 'YColor', 'k');
yyaxis right
plot(t, u_dist(:,3), 'g-', 'LineWidth', 2); grid on;
ylabel('Força de inserção (N)');
set(gca, 'YColor', 'k');
title(['Esforço de controle sob perturbação (' controller_name ')']);
xlabel('Tempo (s)');
legend('\tau_1', '\tau_2', 'F_3', 'Location', 'best');


%  4) CLOSED-LOOP BODE DIAGRAM

% Evaluates the frequency response of the disturbance channel (w) to the joints (y)
figure('Name', ['Diagrama de Bode de Malha Fechada - ' controller_name], 'Color', 'w');
bode(sys_dist_cl); grid on;
title(['Diagrama de Bode de MF: Distúrbio w -> juntas medidas (' controller_name ')']);

end