function SimCLTotal(A, B, C, E, K_lqr, F_red, G_red, H_red, S_red, N_obs, x0_real, t_sim)
% SIMCLTOTAL - Performs a full closed-loop simulation of the combined system
% (LQR regulator + Friedland Reduced-Order Observer) under respiratory disturbance.
%
% Inputs: A, B, C, E     - State-space matrices of the reduced system (5x5)
%         K_lqr          - LQR feedback gain matrix
%         F_red, G_red, H_red, S_red, N_obs - Reduced-order observer matrices
%         x0_real        - Real initial state vector of the plant
%         t_sim2         - Simulation time vector

n_states  = size(A, 1);
n_inputs  = size(B, 2);
n_outputs = size(C, 1);

% Conversion factors
rad_to_deg = 180 / pi;
m_to_mm = 1000;

% Define and verify destination folder structure for closed-loop figures
folder_path = fullfile('Figures', 'ClosedLoop');
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end


%  1) AUGMENTED CLOSED-LOOP SYSTEM ASSEMBLY [x; z]

% Combining plant dynamics and observer estimation law: u = -K_lqr * x_hat
A_cl_total = [ (A - B * K_lqr * S_red * C),            (-B * K_lqr * N_obs);
               (G_red * C - H_red * K_lqr * S_red * C), (F_red - H_red * K_lqr * N_obs) ];

% Disturbance mapping: Breathing signal w affects only the physical plant channel
B_cl_total = [ E; 
               zeros(n_states - n_outputs, size(E, 2)) ];

C_cl_total = eye(2 * n_states - n_outputs);
D_cl_total = zeros(2 * n_states - n_outputs, 1);

sys_cl_total = ss(A_cl_total, B_cl_total, C_cl_total, D_cl_total);


%  2) DISTURBANCE INPUT & INITIAL CONDITIONS CONFIGURATION

f_breathing = 0.5; % 0.5 Hz frequency mapping
omega_w = 2 * pi * f_breathing;
w_sinusoidal = 0.005 * sin(omega_w * t_sim); % 5 mm breathing amplitude

% Dynamic observer initialization (z0 = w_hat0 - J * C * x0_real)
% Computing directly from the state reconstruction compatibility error at t=0:
z0 = pinv(N_obs) * (zeros(n_states, 1) - S_red * C * x0_real);
x0_cl_total = [x0_real; z0];


%  3) TEMPORAL SIMULATION EXECUTION

[X_cl_total, ~] = lsim(sys_cl_total, w_sinusoidal, t_sim, x0_cl_total);

% Separate physical plant states from the unmeasured observer states
x_real_sim = X_cl_total(:, 1:n_states);
z_sim = X_cl_total(:, n_states+1:end);

% Reconstruct estimated variables (x_hat) and control actions (u) point by point
x_hat_sim = zeros(length(t_sim), n_states);
u_sim = zeros(length(t_sim), n_inputs);

for i = 1:length(t_sim)
    x_hat_sim(i, :) = (S_red * C * x_real_sim(i, :)' + N_obs * z_sim(i, :)')';
    u_sim(i, :)     = (-K_lqr * x_hat_sim(i, :)')';
end


%  4) FIGURE GENERATION & EXPORT

% --- PLOT 1: REAL VS ESTIMATED STATES - POSITIONS ---
fig3 = figure('Name', 'Observer Validation: Joint Positions', 'Color', 'w');
subplot(2,1,1);
plot(t_sim, x_real_sim(:,1) * rad_to_deg, 'b-',  'LineWidth', 2); hold on;
plot(t_sim, x_hat_sim(:,1) * rad_to_deg,  'r--', 'LineWidth', 1.5); grid on;
title('Joint q_2: Real State vs Reduced-Order Estimate');
ylabel('Position (deg)'); legend('Real', 'Estimated', 'Location', 'best');

subplot(2,1,2);
plot(t_sim, x_real_sim(:,2) * m_to_mm, 'b-',  'LineWidth', 2); hold on;
plot(t_sim, x_hat_sim(:,2) * m_to_mm,  'r--', 'LineWidth', 1.5); grid on;
title('Joint q_3: Real State vs Reduced-Order Estimate');
xlabel('Time (s)'); ylabel('Position (mm)');
saveas(fig3, fullfile(folder_path, 'Tracking_Real_vs_Estimated_Positions.png'));

% --- PLOT 2: REAL VS ESTIMATED STATES - VELOCITIES ---
fig4 = figure('Name', 'Observer Validation: Joint Velocities', 'Color', 'w');
vel_names_cl   = {'dq_1', 'dq_2', 'dq_3 '};
vel_ylabels_cl = {'Velocity (rad/s)', 'Velocity (rad/s)', 'Velocity (m/s)'};
for k = 3:5
    subplot(3,1,k-2);
    plot(t_sim, x_real_sim(:,k), 'b-',  'LineWidth', 2); hold on;
    plot(t_sim, x_hat_sim(:,k),  'r--', 'LineWidth', 1.5); grid on;
    title([vel_names_cl{k-2} ': Real State vs Reduced-Order Estimate']);
    ylabel(vel_ylabels_cl{k-2});
    if k == 5, xlabel('Time (s)'); end
    if k == 3, legend('Real', 'Estimated', 'Location', 'best'); end
end
saveas(fig4, fullfile(folder_path, 'Tracking_Real_vs_Estimated_Velocities.png'));

% --- PLOT 3: CONTROL EFFORT ---
fig2 = figure('Name', 'Closed Loop: Actuator Control Effort', 'Color', 'w');
yyaxis left
plot(t_sim, u_sim(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t_sim, u_sim(:,2), 'r-', 'LineWidth', 2);
ylabel('Motor Torques (Nm)'); set(gca, 'YColor', 'k');

yyaxis right
plot(t_sim, u_sim(:,3), 'g-', 'LineWidth', 2);
ylabel('Insertion Force F_3 (N)'); set(gca, 'YColor', 'k');
grid on;

title('Actuator Control Effort');
xlabel('Time (s)');
legend('\tau_1', '\tau_2', 'F_3', 'Location', 'best');
saveas(fig2, fullfile(folder_path, 'Closed_Loop_Control_Effort.png'));


% --- PLOT 5: CONTROL EFFORT DERIVATIVE (du/dt) ---
dt_sim2 = t_sim(2) - t_sim(1);
du_sim2 = diff(u_sim) / dt_sim2;
t_du2   = t_sim(1:end-1);

fig5 = figure('Name', 'Control Effort Rate of Change', 'Color', 'w');
subplot(3,1,1);
plot(t_du2, du_sim2(:,1), 'b', 'LineWidth', 1.5); grid on;
title('Control Rate: d\tau_1/dt'); ylabel('Rate (Nm/s)');

subplot(3,1,2);
plot(t_du2, du_sim2(:,2), 'r', 'LineWidth', 1.5); grid on;
title('Control Rate: d\tau_2/dt'); ylabel('Rate (Nm/s)');

subplot(3,1,3);
plot(t_du2, du_sim2(:,3), 'g', 'LineWidth', 1.5); grid on;
title('Control Rate: dF_3/dt'); xlabel('Time (s)'); ylabel('Rate (N/s)');
saveas(fig5, fullfile(folder_path, 'Closed_Loop_Control_Rates.png'));

fprintf('=> Closed-loop simulation figures successfully saved to "%s"\n\n', folder_path);
end