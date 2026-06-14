function AnalyzeCL(A, B, C, D, E, K, controller_name)
% AnalyzeCL - Performs a refined closed-loop analysis of a regulator
%
% Inputs: A, B, C, D, E   - State-space matrices of the system
%         K               - Feedback gain matrix (u = -Kx)
%         controller_name - String with the controller name ('LQR' or 'PP')

n_states  = size(A, 1);
n_inputs  = size(B, 2);
n_outputs = size(C, 1);

% Conversion factors for improved visual analysis
rad_to_deg = 180 / pi;
m_to_mm    = 1000;

% Define and create the destination folder structure for saving figures
folder_path = fullfile('Figures', 'Regulators');
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

% 1) POLE MAP: NEW POLES COMPARED TO OLD POLES
poles_ma = eig(A);
poles_cl = eig(A - B * K);

fig1 = figure('Name', ['Pole Comparison - ' controller_name], 'Color', 'w');
plot(real(poles_ma), imag(poles_ma), 'b+', 'MarkerSize', 8, 'LineWidth', 2); hold on;
plot(real(poles_cl), imag(poles_cl), 'rx', 'MarkerSize', 8, 'LineWidth', 2);
grid on; xline(0, 'k--'); yline(0, 'k--');

title(['Pole Mapping: OL vs CL (' controller_name ')']);
xlabel('Real Axis (1/s)');
ylabel('Imaginary Axis (1/s)');
legend('OL (Open-Loop)', ['CL (' controller_name ')'], 'Location', 'best');
hold off;

% Save Figure 1
saveas(fig1, fullfile(folder_path, ['Pole_Mapping_' controller_name '.png']));

% 2) POSITION REGULATOR SCENARIO (INITIAL CONDITION RESPONSE)
sys_ma_init = ss(A, zeros(size(B)), C, zeros(size(D)));
sys_cl_init = ss(A - B * K, zeros(size(B)), C, zeros(size(D)));

t = 0:0.001:2;              % Time vector
x0 = [0.05; 0.02; 0; 0; 0]; % Initial state deviation from equilibrium

[y_ma_init, ~, ~]         = initial(sys_ma_init, x0, t);
[y_cl_init, ~, x_cl_init] = initial(sys_cl_init, x0, t);

% 2.1) Position plot with unit conversions (rad -> deg, m -> mm)
fig2 = figure('Name', ['Position Regulation - Initial Condition - ' controller_name], 'Color', 'w');
subplot(2,1,1);
plot(t, y_ma_init(:,1) * rad_to_deg, 'r--', 'LineWidth', 1.5); hold on;
plot(t, y_cl_init(:,1) * rad_to_deg, 'b', 'LineWidth', 2); grid on;
title(['Joint q_2 Regulation | Initial Condition Response (' controller_name ')']);
ylabel('Angular Deviation (deg)');
legend('OL', 'CL', 'Location', 'best');

subplot(2,1,2);
plot(t, y_ma_init(:,2) * m_to_mm, 'r--', 'LineWidth', 1.5); hold on;
plot(t, y_cl_init(:,2) * m_to_mm, 'b', 'LineWidth', 2); grid on;
title(['Joint q_3 Regulation | Initial Condition Response (' controller_name ')']);
xlabel('Time (s)'); ylabel('Linear Deviation (mm)');
hold off;

% Save Figure 2
saveas(fig2, fullfile(folder_path, ['Position_Regulation_Initial_' controller_name '.png']));

% 2.2) Control effort plot (u = -Kx)
u_init = zeros(length(t), n_inputs);
for idx = 1:length(t)
    u_init(idx, :) = (-K * x_cl_init(idx, :)')';
end

fig3 = figure('Name', ['Control Effort - Initial Condition - ' controller_name], 'Color', 'w');
yyaxis left
plot(t, u_init(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t, u_init(:,2), 'r-', 'LineWidth', 2);
ylabel('Motor Torques (Nm)');
set(gca, 'YColor', 'k');

yyaxis right
plot(t, u_init(:,3), 'g-', 'LineWidth', 2); grid on;
ylabel('Insertion Force (N)');
set(gca, 'YColor', 'k');

title(['Regulator Control Effort (' controller_name ')']);
xlabel('Time (s)');
legend('\tau_1', '\tau_2', 'F_3', 'Location', 'best');
hold off;

% Save Figure 3
saveas(fig3, fullfile(folder_path, ['Control_Effort_Initial_' controller_name '.png']));

% 3) SINUSOIDAL DISTURBANCE SCENARIO
t = 0:0.001:5;  % Time vector

f_breathing = 0.5;
omega_w = 2 * pi * f_breathing;
w_sinusoidal = 0.005 * sin(omega_w * t);

sys_dist_ma = ss(A, E, C, zeros(n_outputs, 1));
sys_dist_cl = ss(A - B * K, E, C, zeros(n_outputs, 1));

y_dist_ma = lsim(sys_dist_ma, w_sinusoidal, t);
[y_dist_cl, ~, x_dist_cl] = lsim(sys_dist_cl, w_sinusoidal, t);

% 3.1) Position plot with unit conversions (rad -> deg, m -> mm)
fig4 = figure('Name', ['Position Regulation - Respiratory Disturbance - ' controller_name], 'Color', 'w');
subplot(2,1,1);
plot(t, y_dist_ma(:,1) * rad_to_deg, 'r--', 'LineWidth', 1.5); hold on;
plot(t, y_dist_cl(:,1) * rad_to_deg, 'b', 'LineWidth', 2); grid on;
title(['Joint q_2 Deviation under Disturbance (' controller_name ')']);
ylabel('Angular Error (deg)');
legend('OL', 'CL');

subplot(2,1,2);
plot(t, y_dist_ma(:,2) * m_to_mm, 'r--', 'LineWidth', 1.5); hold on;
plot(t, y_dist_cl(:,2) * m_to_mm, 'b', 'LineWidth', 2); grid on;
title(['Joint q_3 Deviation under Disturbance (' controller_name ')']);
xlabel('Time (s)'); ylabel('Linear Error (mm)');
hold off;

% Save Figure 4
saveas(fig4, fullfile(folder_path, ['Position_Regulation_Disturbance_' controller_name '.png']));

% 3.2) Control effort plot under disturbance
u_dist = zeros(length(t), n_inputs);
for idx = 1:length(t)
    u_dist(idx, :) = (-K * x_dist_cl(idx, :)')';
end

fig5 = figure('Name', ['Control Effort - Respiratory Disturbance - ' controller_name], 'Color', 'w');
yyaxis left
plot(t, u_dist(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t, u_dist(:,2), 'r-', 'LineWidth', 2);
ylabel('Motor Torques (Nm)');
set(gca, 'YColor', 'k');

yyaxis right
plot(t, u_dist(:,3), 'g-', 'LineWidth', 2); grid on;
ylabel('Insertion Force (N)');
set(gca, 'YColor', 'k');

title(['Control Effort under Disturbance (' controller_name ')']);
xlabel('Time (s)');
legend('\tau_1', '\tau_2', 'F_3', 'Location', 'best');
hold off;

% Save Figure 5
saveas(fig5, fullfile(folder_path, ['Control_Effort_Disturbance_' controller_name '.png']));

% 4) CLOSED-LOOP BODE DIAGRAM
fig6 = figure('Name', ['Closed-Loop Bode Diagram - ' controller_name], 'Color', 'w');
bode(sys_dist_cl); grid on;
title(['Closed-Loop Bode: Disturbance w -> Measured Joints (' controller_name ')']);

% Save Figure 6
saveas(fig6, fullfile(folder_path, ['Closed_Loop_Bode_' controller_name '.png']));

end