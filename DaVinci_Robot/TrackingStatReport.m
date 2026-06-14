function TrackingStatReport(poles_cl, poles_obs_red, u_seg, t_sim, y_ref, y_real)
% TRACKINGSTATREPORT - Compiles and displays a comprehensive statistical report
% of tracking performance, actuator control efforts, and closed-loop poles.
%
% Inputs: poles_cl      - Closed-loop controller poles vector (5x1)
%         poles_obs_red - Reduced-order observer poles vector (3x1)
%         u_seg         - Control efforts time history matrix (N_points x 3)
%         t_sim         - Simulation time vector (N_points x 1)
%         y_ref         - Target reference positions matrix (N_points x 2) -> [q2_ref, q3_ref]
%         y_real        - Real plant output positions matrix (N_points x 2)  -> [q2_real, q3_real]

disp('==================================================');
disp('   COMPILING TRACKING STATISTICAL REPORT...       ');
disp('==================================================');

% 1. Calculate Actuator Effort Metrics (Max, Mean, and RMS values)
max_u  = max(abs(u_seg));
mean_u = mean(abs(u_seg));
rms_u  = rms(u_seg);

% 2. Isolate Steady-State Tracking Errors (Evaluating strictly for t >= 0.5s)
idx_steady = (t_sim >= 0.5);
if ~any(idx_steady)
    idx_steady = true(size(t_sim)); % Fallback to avoid empty indexing bugs if horizon is too short
end

% Conversion factors for intuitive clinical analysis
rad_to_deg = 180 / pi;
m_to_mm    = 1000;

% Compute tracking errors scaled to analytical units
error_q2 = (y_ref(idx_steady, 1) - y_real(idx_steady, 1)) * rad_to_deg;
error_q3 = (y_ref(idx_steady, 2) - y_real(idx_steady, 2)) * m_to_mm;

max_error_steady = [max(abs(error_q2)), max(abs(error_q3))];
rms_error_steady = [rms(error_q2), rms(error_q3)];

% =====================================================================
%  CONSOLE REPORT PRINTING
% =====================================================================

% --- TABLE 1: CLOSED-LOOP POLES (SEPARATION PRINCIPLE SUMMARY) ---
fprintf('\n===================================================================\n');
fprintf('        UPDATED SYSTEM POLES SUMMARY (FINE-TUNING PROFILE) \n');
fprintf('===================================================================\n');
fprintf(' Mode |     Controller Closed-Loop     |     Reduced-Order Observer   \n');
fprintf('-------------------------------------------------------------------\n');
for i = 1:5
    str_ctrl = sprintf('%+6.2f %+6.2fi', real(poles_cl(i)), imag(poles_cl(i)));
    if i <= length(poles_obs_red)
        str_red = sprintf('%+6.2f %+6.2fi', real(poles_obs_red(i)), imag(poles_obs_red(i)));
        fprintf('  %d   | %-30s | %-20s \n', i, str_ctrl, str_red);
    else
        fprintf('  %d   | %-30s | %-20s \n', i, str_ctrl, '         -         ');
    end
end
fprintf('===================================================================\n');

% --- TABLE 2: ACTUATOR CONTROL EFFORT METRICS ---
fprintf('\n===================================================================\n');
fprintf('            ACTUATOR METRICS REPORT (CONTROL EFFORT) \n');
fprintf('===================================================================\n');
fprintf(' Actuator       |   Maximum Value  |    Mean Value    |    RMS Value     \n');
fprintf('-------------------------------------------------------------------\n');
fprintf(' \\tau_1 (Base)   | %12.4f Nm | %11.4f Nm | %11.4f Nm \n', max_u(1), mean_u(1), rms_u(1));
fprintf(' \\tau_2 (Elev.)  | %12.4f Nm | %11.4f Nm | %11.4f Nm \n', max_u(2), mean_u(2), rms_u(2));
fprintf(' F_3 (Insertion) | %12.4f N  | %11.4f N  | %11.4f N  \n', max_u(3), mean_u(3), rms_u(3));
fprintf('===================================================================\n');

% --- TABLE 3: STEADY-STATE TRACKING ERRORS ---
fprintf('\n===================================================================\n');
fprintf('          STEADY-STATE TRACKING ERRORS SUMMARY (t >= 0.5 s) \n');
fprintf('===================================================================\n');
fprintf(' Joint         |     Maximum Absolute Error  |          RMS Error        \n');
fprintf('-------------------------------------------------------------------\n');
fprintf(' q2 (Elevation) | %18.6f deg       | %16.6f deg  \n', max_error_steady(1), rms_error_steady(1));
fprintf(' q3 (Insertion) | %18.6f mm        | %16.6f mm   \n', max_error_steady(2), rms_error_steady(2));
fprintf('===================================================================\n\n');
end