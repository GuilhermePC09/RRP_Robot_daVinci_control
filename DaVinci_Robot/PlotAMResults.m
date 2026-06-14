function PlotAssumedModelResults(t_sim, X_out, A, B, C, S_red, N_obs, K_lqr, G_r, poles_cl, poles_obs_red)
% PLOTAMRESULTS - Processes raw 13-state simulation data from the 
% Assumed Model Tracker

n_states  = size(A, 1);
n_inputs  = size(B, 2);
n_outputs = size(C, 1);

% Conversion factors
rad_to_deg = 180 / pi;
m_to_mm    = 1000;

% Define and verify destination folder structure for tracking plots
folder_path = fullfile('Figures', 'Tracking');
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

% Extract specific state channels from the raw 13-column output matrix
q2_real = X_out(:, 1); q3_real = X_out(:, 2);
q2_ref  = X_out(:, 9); q3_ref  = X_out(:, 10);

num_points = length(t_sim);
u_seg = zeros(num_points, n_inputs);
x_hat_ma = zeros(num_points, n_states);

% Point-by-point algebraic reconstruction of control actions and estimates
for i = 1:num_points
    xhat_t = S_red * C * X_out(i, 1:5)' + N_obs * X_out(i, 6:8)';
    xr_t   = X_out(i, 9:13)';
    
    % Control law implementation: u(t) = -K*(x_hat - x_r) + u_ff
    u_seg(i, :) = (-K_lqr * xhat_t + (K_lqr - G_r) * xr_t)';
    x_hat_ma(i, :) = xhat_t';
end


%  PLOTS

% --- PLOT 1: TRAJECTORY TRACKING PERFORMANCE ---
fig1 = figure('Name', 'Model Follower: Trajectory Tracking', 'Color', 'w');
subplot(2,1,1);
plot(t_sim, q2_ref * rad_to_deg, 'r--', 'LineWidth', 1.5); hold on;
plot(t_sim, q2_real * rad_to_deg, 'b', 'LineWidth', 2); grid on;
title('Joint q_2 Trajectory Tracking | Assumed Model');
ylabel('Position (deg)'); legend('Reference', 'Assumed Model Tacker', 'Location', 'best');

subplot(2,1,2);
plot(t_sim, q3_ref * m_to_mm, 'r--', 'LineWidth', 1.5); hold on;
plot(t_sim, q3_real * m_to_mm, 'b', 'LineWidth', 2); grid on;
title('Joint q_3 Trajectory Tracking | Assumed Model');
xlabel('Time (s)'); ylabel('Position (mm)');
saveas(fig1, fullfile(folder_path, 'Assumed_Model_Trajectory_Tracking.png'));

% --- PLOT 2: ACTUATOR EFFORT (DUAL INDEPENDENT AXES) ---
fig2 = figure('Name', 'Model Follower: Actuator Control Effort', 'Color', 'w');
yyaxis left
plot(t_sim, u_seg(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t_sim, u_seg(:,2), 'r-', 'LineWidth', 2);
ylabel('Motor Torques (Nm)'); set(gca, 'YColor', 'k');

yyaxis right
plot(t_sim, u_seg(:,3), 'g-', 'LineWidth', 2);
ylabel('Insertion Force F_3 (N)'); set(gca, 'YColor', 'k');
grid on;

title('Actuator Control Effort | Assumed Model');
xlabel('Time (s)');
legend('\tau_1', '\tau_2', 'F_3', 'Location', 'best');
saveas(fig2, fullfile(folder_path, 'Assumed_Model_Control_Effort.png'));

% --- PLOT 3: CONTROL EFFORT VARIATION RATE (du/dt) ---
dt_ma   = t_sim(2) - t_sim(1);
du_ma   = diff(u_seg) / dt_ma;
t_du_ma = t_sim(1:end-1);

fig3 = figure('Name', 'Model Tracker: Control Effort Rate of Change', 'Color', 'w');
subplot(3,1,1);
plot(t_du_ma, du_ma(:,1), 'b', 'LineWidth', 1.5); grid on;
title('Control Rate: d\tau_1/dt | Assumed Model'); ylabel('Rate (Nm/s)');

subplot(3,1,2);
plot(t_du_ma, du_ma(:,2), 'r', 'LineWidth', 1.5); grid on;
title('Control Rate: d\tau_2/dt | Assumed Model'); ylabel('Rate (Nm/s)');

subplot(3,1,3);
plot(t_du_ma, du_ma(:,3), 'g', 'LineWidth', 1.5); grid on;
title('Control Rate: dF_3/dt | Assumed Model'); xlabel('Time (s)'); ylabel('Rate (N/s)');
saveas(fig3, fullfile(folder_path, 'Assumed_Model_Control_Rates.png'));

fprintf('=> Assumed Model simulation figures successfully saved to "%s".\n\n', folder_path);

% Statistics
TrackingStatReport(poles_cl, poles_obs_red, u_seg, t_sim, [q2_ref, q3_ref], [q2_real, q3_real]);

end