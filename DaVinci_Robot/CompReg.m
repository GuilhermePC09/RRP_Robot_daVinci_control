function CompReg(A, B, C, D, K_lqr, F_pp)
% COMPREG - Performs a comprehensive frequency and time-domain comparison
% between Open-Loop (OL), LQR, and Pole Placement (PP) regulators.
%
% Inputs: A, B, C, D - State-space matrices of the system (5x5 configuration)
%         K_lqr      - Optimal LQR feedback gain matrix
%         F_pp       - Pole Placement feedback gain matrix

n_states  = size(A, 1);
n_inputs  = size(B, 2);
n_outputs = size(C, 1);

% Conversion factors
rad_to_deg = 180 / pi;
m_to_mm    = 1000;

% Define and verify destination folder structure for saving comparison plots
folder_path = fullfile('Figures', 'Regulators');
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

% Instantiate the state-space system objects
sys_ol     = ss(A,            B, C, D);
sys_cl_lqr = ss(A - B*K_lqr, B, C, D);
sys_cl_pp  = ss(A - B*F_pp,  B, C, D);


%  1) FREQUENCY-DOMAIN COMPARISON: MIMO SINGULAR VALUES (SIGMA)

fig_sigma = figure('Name', 'Frequency Response: MIMO Singular Values', 'Color', 'w');
sigma(sys_ol, 'k--', sys_cl_lqr, 'b-', sys_cl_pp, 'r-');
grid on;

title('MIMO Frequency Response Comparison (Singular Values)');
xlabel('Frequency (rad/s)'); 
ylabel('Singular Values / Gain (dB)');
legend('Open-Loop (OL)', 'LQR Regulator', 'PP Regulator', 'Location', 'best');

% Save Sigma Plot
saveas(fig_sigma, fullfile(folder_path, 'Comparison_MIMO_Sigma.png'));


%  2) TIME-DOMAIN COMPARISON: INITIAL CONDITION RESPONSE

t_time = 0:0.001:2;
x0 = [0.05; 0.02; 0; 0; 0];

% Systems configured for initial condition excitation
sys_ol_init  = ss(A, zeros(size(B)), C, zeros(size(D)));
sys_lqr_init = ss(A - B*K_lqr, zeros(size(B)), C, zeros(size(D)));
sys_pp_init  = ss(A - B*F_pp, zeros(size(B)), C, zeros(size(D)));

[y_ol, ~]  = initial(sys_ol_init,  x0, t_time);
[y_lqr, ~] = initial(sys_lqr_init, x0, t_time);
[y_pp, ~]  = initial(sys_pp_init,  x0, t_time);

fig_time = figure('Name', 'Time Response: Initial Condition Regulation', 'Color', 'w', 'Position', [150 150 900 600]);

% Subplot 1: Joint q2 (Elevation) Regulation
subplot(2,1,1);
plot(t_time, y_ol(:,1) * rad_to_deg, 'k--', 'LineWidth', 1.5); hold on;
plot(t_time, y_lqr(:,1) * rad_to_deg, 'b-', 'LineWidth', 2);
plot(t_time, y_pp(:,1) * rad_to_deg, 'r-.', 'LineWidth', 1.8);
grid on;
title('Joint q_2 Regulation Transient Response');
ylabel('Deviation (deg)');
legend('Open-Loop (OL)', 'LQR Regulator', 'PP Regulator', 'Location', 'best');

% Subplot 2: Joint q3 (Insertion) Regulation
subplot(2,1,2);
plot(t_time, y_ol(:,2) * m_to_mm, 'k--', 'LineWidth', 1.5); hold on;
plot(t_time, y_lqr(:,2) * m_to_mm, 'b-', 'LineWidth', 2);
plot(t_time, y_pp(:,2) * m_to_mm, 'r-.', 'LineWidth', 1.8);
grid on;
title('Joint q_3 Regulation Transient Response');
xlabel('Time (s)'); ylabel('Deviation (mm)');
hold off;

% Save Time Response Plot
saveas(fig_time, fullfile(folder_path, 'Comparison_Time_Regulation.png'));


%  3) REFINED BODE MAGNITUDE COMPARISON (HIGH-QUALITY FORMATTING)

bode_out_names = {'Joint q_2', 'Joint q_3'};
bode_in_names  = {'Input \tau_1', 'Input \tau_2', 'Input F_3'};

w_vec = logspace(-1, 3, 400); 

for out_idx = 1:n_outputs
    fig_bode = figure('Name', ['Bode Magnitude Channel - Output ' num2str(out_idx)], 'Color', 'w', 'Position', [200 100 750 850]);
    
    for in_idx = 1:n_inputs
        subplot(n_inputs, 1, in_idx);
        
        % Numerically extract magnitude responses for the current SISO channel
        [mag_ol, ~, ~]  = bode(sys_ol(out_idx, in_idx), w_vec);
        [mag_lqr, ~, ~] = bode(sys_cl_lqr(out_idx, in_idx), w_vec);
        [mag_pp, ~, ~]  = bode(sys_cl_pp(out_idx, in_idx), w_vec);
        
        mag_ol_dB  = 20 * log10(squeeze(mag_ol));
        mag_lqr_dB = 20 * log10(squeeze(mag_lqr));
        mag_pp_dB  = 20 * log10(squeeze(mag_pp));
        
        semilogx(w_vec, mag_ol_dB, 'k--', 'LineWidth', 1.5); hold on;
        semilogx(w_vec, mag_lqr_dB, 'b-', 'LineWidth', 2);
        semilogx(w_vec, mag_pp_dB, 'r-.', 'LineWidth', 1.5);
        grid on;
        
        title(sprintf('From: %s  ->  To: %s', bode_in_names{in_idx}, bode_out_names{out_idx}), 'FontSize', 10);
        ylabel('Magnitude (dB)', 'FontSize', 9);
        
        if in_idx == n_inputs
            xlabel('Frequency (rad/s)', 'FontSize', 10);
        else
            set(gca, 'XTickLabel', []);
        end
        
        if in_idx == 1
            legend('Open-Loop (OL)', 'LQR Regulator', 'PP Regulator', 'Location', 'best');
        end
    end
    hold off;
    
    % Save Bode Figures
    clean_name = strrep(strrep(bode_out_names{out_idx}, ' ', '_'), '(', '');
    clean_name = strrep(clean_name, ')', '');
    saveas(fig_bode, fullfile(folder_path, ['Comparison_Bode_Magnitude_' clean_name '.png']));
end

fprintf('=> Comparison plots for regulators successfully saved to "%s"\n\n', folder_path);
end