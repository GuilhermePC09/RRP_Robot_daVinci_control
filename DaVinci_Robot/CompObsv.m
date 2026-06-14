function CompObsv(A, C, L, obs_red, poles_control, x0_real, t_sim, state_labels, x0_hat)
% COMPOBSV - Performs a unified comparative analysis between Full-Order (Identity)
% and Reduced-Order (Friedland) state observers in both frequency and time domains.
%
% Inputs: A, C          - Plant state and output matrices (5x5 configuration)
%         L              - Full-order observer feedback gain matrix
%         obs_red        - Struct containing reduced-order observer matrices
%         poles_control  - Closed-loop poles of the controller (for separation principle table)
%         x0_real        - Real initial state vector of the plant
%         t_sim          - Simulation time vector
%         state_labels   - Cell array with string labels for each state
%         x0_hat         - Initial state estimate vector

n_states  = size(A, 1);
n_outputs = size(C, 1);
n_red     = n_states - n_outputs;
n_aug     = 2*n_states - n_outputs;

% Unit conversion multipliers for scaling plots
unit_scales = [180/pi, 1000, 1, 1, 1]; 
unit_strings = {'(deg)', '(mm)', '(rad/s)', '(rad/s)', '(m/s)'};

% Define and verify destination folder structure for saving figures
folder_path = fullfile('Figures', 'Observers');
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

F_red        = obs_red.F_red;
A_aug_red    = obs_red.A_aug_red;
recover_xhat = obs_red.recover_xhat;
J            = obs_red.J;

A_obs          = A - L*C;
poles_obs_full = eig(A_obs);
poles_obs_red  = eig(F_red);
poles_control  = poles_control(:);


%  1) COMPARATIVE POLE TABLE

fmt = @(p) sprintf('%+6.2f %+6.2fi', real(p), imag(p));
dash = '       -       ';

fprintf('\n===================================================================\n');
fprintf('        GLOBAL SYSTEM POLE COMPARISON (SEPARATION PRINCIPLE)     \n');
fprintf('===================================================================\n');
fprintf(' Mode |    Controller (CL)    |   Full-Order Obs.   |   Reduced-Order Obs. \n');
fprintf('-------------------------------------------------------------------\n');
for i = 1:n_states
    if i <= numel(poles_control),  s_ctrl = fmt(poles_control(i));  else, s_ctrl = dash; end
    if i <= numel(poles_obs_full), s_full = fmt(poles_obs_full(i)); else, s_full = dash; end
    if i <= numel(poles_obs_red),  s_red  = fmt(poles_obs_red(i));  else, s_red  = dash; end
    fprintf('  %d   | %-21s | %-17s | %-17s \n', i, s_ctrl, s_full, s_red);
end
fprintf('===================================================================\n\n');


%  2) MATHEMATICAL SIMULATION & ESTIMATION RECONSTRUCTION

e0 = x0_real - x0_hat;

% Reconstruct identity states from error dynamics: xhat = x_real - e
sys_err_full = ss(A_obs, zeros(n_states,1), eye(n_states), zeros(n_states,1));
[e_full, ~]  = initial(sys_err_full, e0, t_sim);

% Reconstruct reduced observer states via augmented parallel simulation
z0 = -J * C * x0_real;
sys_aug_open = ss(A_aug_red, zeros(n_aug,1), eye(n_aug), zeros(n_aug,1));
[X_aug, ~]   = initial(sys_aug_open, [x0_real; z0], t_sim);

x_real_sim   = X_aug(:, 1:n_states);
xhat_red_sim = (recover_xhat * X_aug.').';
xhat_full_sim = x_real_sim - e_full;

e_red_sim = x_real_sim - xhat_red_sim;


%  3) ANALYSIS PLOTS

% --- PLOT A: ESTIMATION ERROR CONVERGENCE ---
fig_error = figure('Name', 'Observer Dynamics: Estimation Error', 'Color', 'w', 'Position', [100 100 1100 850]);
for k = 1:n_states
    subplot(3, 2, k);
    plot(t_sim, e_full(:, k) * unit_scales(k), 'b-', 'LineWidth', 1.8); hold on;
    plot(t_sim, e_red_sim(:, k) * unit_scales(k), 'g-.', 'LineWidth', 1.8);
    grid on;
    
    ylabel([ 'Error ' unit_strings{k} ]);
    title([ 'Estimation Error: ' state_labels{k} ], 'FontSize', 9);
    if k >= 4, xlabel('Time (s)'); end
end
% Smart layout allocation for shared legend in the 6th slot
subplot(3,2,6); axis off;
legend(subplot(3,2,1), 'Full-Order Observer', 'Reduced-Order Observer', 'Location', 'best', 'FontSize', 11);
saveas(fig_error, fullfile(folder_path, 'Observer_Error_Convergence.png'));

% --- PLOT B: STATE TRACKING AND CONVERGENCE ---
fig_track = figure('Name', 'Observer Dynamics: State Tracking', 'Color', 'w', 'Position', [150 120 1100 850]);
for k = 1:n_states
    subplot(3, 2, k);
    plot(t_sim, x_real_sim(:, k) * unit_scales(k), 'r-', 'LineWidth', 2.2); hold on;
    plot(t_sim, xhat_full_sim(:, k) * unit_scales(k), 'b--', 'LineWidth', 1.5);
    plot(t_sim, xhat_red_sim(:, k) * unit_scales(k), 'g-.', 'LineWidth', 1.5);
    grid on;
    
    ylabel(unit_strings{k});
    title([ 'Tracking: ' state_labels{k} ], 'FontSize', 9);
    if k >= 4, xlabel('Time (s)'); end
end
subplot(3,2,6); axis off;
legend(subplot(3,2,1), 'Real Plant State', 'Full-Order Estimate', 'Reduced-Order Estimate', 'Location', 'best', 'FontSize', 10);
saveas(fig_track, fullfile(folder_path, 'Observer_State_Tracking.png'));

% --- PLOT C: EXPANDED GLOBAL POLE MAP ---
fig_map = figure('Name', 'Global Pole Map: Separation Principle', 'Color', 'w', 'Position', [200 200 700 550]);
plot(real(poles_control),  imag(poles_control),  'rx', 'MarkerSize', 10, 'LineWidth', 2); hold on;
plot(real(poles_obs_full), imag(poles_obs_full), 'bs', 'MarkerSize', 8,  'LineWidth', 1.5);
plot(real(poles_obs_red),  imag(poles_obs_red),  'g+', 'MarkerSize', 9,  'LineWidth', 1.5);
grid on; xline(0, 'k--'); yline(0, 'k--');

xlabel('Real Axis (1/s)'); ylabel('Imaginary Axis (1/s)');
title('Global Pole Map: Separation Principle Verification');
legend('Controller (Closed-Loop)', ...
       sprintf('Identity Observer (n=%d)', n_states), ...
       sprintf('Reduced-Order Observer (n=%d)', n_red), 'Location', 'best');
hold off;
saveas(fig_map, fullfile(folder_path, 'Global_Pole_Mapping.png'));

fprintf('=> Observer comparison analysis successfully generated and saved to "%s".\n\n', folder_path);
end