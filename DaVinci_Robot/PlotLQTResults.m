function PlotLQTResults(t_sim, X_total, A, B, C, S_red, N_obs, Ri, t_P, P_flat_sol, t_eta, eta_sol, poles_cl_lqt, poles_obs_red)
% PLOTLQTRESULTS - Processes simulation data from the finite-horizon LQT

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

% Pre-allocate arrays for reconstruction loop
num_points = length(t_sim);
x_real = X_total(:, 1:n_states);
z_obs = X_total(:, n_states+1:end);

x_hat = zeros(num_points, n_states);
u_seg = zeros(num_points, n_inputs);
y_ref = zeros(num_points, n_outputs);

% Point-by-point algebraic reconstruction
for i = 1:num_points
    t = t_sim(i);
    y_current = C * x_real(i, :)';
    
    % Reconstruct estimated state using the reduced-order observer mapping
    x_hat(i, :) = (S_red * y_current + N_obs * z_obs(i, :)')';
    
    % Interpolate time-varying differential matrices
    P_flat = interp1(t_P, P_flat_sol, t, 'linear', 'extrap')';
    P = reshape(P_flat, n_states, n_states);
    eta = interp1(t_eta, eta_sol, t, 'linear', 'extrap')';
    
    % Compute time-varying optimal tracking gains
    K_t = Ri * B' * P;
    u_tilde = Ri * B' * eta;
    
    % Optimal control law application: u(t) = -K(t)*x_hat(t) + u_tilde(t)
    u_seg(i, :) = (-K_t * x_hat(i, :)' + u_tilde)';
    
    % Extract trajectory reference mapping at the current instant
    ref_state = daVinci_ref(t);
    y_ref(i, :) = (C * ref_state)';
end

%  PLOTS

% --- PLOT 1: TRAJECTORY TRACKING PERFORMANCE ---
fig1 = figure('Name', 'LQT: Trajectory Tracking Performance', 'Color', 'w');
subplot(2,1,1);
plot(t_sim, y_ref(:,1) * rad_to_deg, 'r--', 'LineWidth', 1.5); hold on;
plot(t_sim, x_real(:,1) * rad_to_deg, 'b', 'LineWidth', 2); grid on;
title('Joint q_2 Trajectory Tracking | LQT');
ylabel('Position (deg)'); legend('Reference', 'LQT', 'Location', 'best');

subplot(2,1,2);
plot(t_sim, y_ref(:,2) * m_to_mm, 'r--', 'LineWidth', 1.5); hold on;
plot(t_sim, x_real(:,2) * m_to_mm, 'b', 'LineWidth', 2); grid on;
title('Joint q_3 Trajectory Tracking | LQT');
xlabel('Time (s)'); ylabel('Position (mm)');
saveas(fig1, fullfile(folder_path, 'LQT_Trajectory_Tracking.png'));

% --- PLOT 2: ACTUATOR EFFORT ---
fig2 = figure('Name', 'LQT: Actuator Control Effort', 'Color', 'w');
yyaxis left
plot(t_sim, u_seg(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t_sim, u_seg(:,2), 'r-', 'LineWidth', 2);
ylabel('Motor Torques (Nm)'); set(gca, 'YColor', 'k');

yyaxis right
plot(t_sim, u_seg(:,3), 'g-', 'LineWidth', 2);
ylabel('Insertion Force F_3 (N)'); set(gca, 'YColor', 'k');
grid on;

title('Actuator Control Effort | LQT');
xlabel('Time (s)');
legend('\tau_1', '\tau_2', 'F_3', 'Location', 'best');
saveas(fig2, fullfile(folder_path, 'LQT_Control_Effort.png'));

% --- PLOT 3: CONTROL EFFORT VARIATION RATE (du/dt) ---
dt_lqt   = t_sim(2) - t_sim(1);
du_lqt   = diff(u_seg) / dt_lqt;
t_du_lqt = t_sim(1:end-1);

fig3 = figure('Name', 'LQT: Control Effort Rate of Change', 'Color', 'w');
subplot(3,1,1);
plot(t_du_lqt, du_lqt(:,1), 'b', 'LineWidth', 1.5); grid on;
title('Control Rate: d\tau_1/dt | LQT'); ylabel('Rate (Nm/s)');

subplot(3,1,2);
plot(t_du_lqt, du_lqt(:,2), 'r', 'LineWidth', 1.5); grid on;
title('Control Rate: d\tau_2/dt | LQT'); ylabel('Rate (Nm/s)');

subplot(3,1,3);
plot(t_du_lqt, du_lqt(:,3), 'g', 'LineWidth', 1.5); grid on;
title('Control Rate: dF_3/dt | LQT'); xlabel('Time (s)'); ylabel('Rate (N/s)');
saveas(fig3, fullfile(folder_path, 'LQT_Control_Rates.png'));

fprintf('=> LQT tracking simulation figures successfully saved to "%s"\n\n', folder_path);

% Statistics
TrackingStatReport(poles_cl_lqt, poles_obs_red, u_seg, t_sim, y_ref, x_real(:, 1:2));

end



% =====================================================================
%  LOCAL DYNAMIC HELPER FUNCTIONS
% =====================================================================

function xr = daVinci_ref(t)
    % daVinci_ref - Generates smooth exponential tracking references
    alpha_regime = 10;
    q2_f = 5 * pi / 180; % Target of 5 degrees converted to radians
    q3_f = 0.05;         % Target of 5 centimeters converted to meters
    
    % Trajectory position profiles
    q2_r = q2_f * (1 - exp(-alpha_regime * t));
    q3_r = q3_f * (1 - exp(-alpha_regime * t));
    
    % Trajectory velocity profiles (First derivatives)
    dq2_r = q2_f * alpha_regime * exp(-alpha_regime * t);
    dq3_r = q3_f * alpha_regime * exp(-alpha_regime * t);
    
    xr = [q2_r; q3_r; 0; dq2_r; dq3_r]; 
end