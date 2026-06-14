function PlotSimControl(t_sim, vec_ref, q_ref, X1, X2, u_lqt, u_ff, scenario_id, pert_id)
% PLOTSIMCONTROL - Performs a comparative tracking analysis between 
% Finite-Horizon LQT and Assumed Model state feedback controllers.
%
% Saves the consolidated performance figure automatically into 'Figures/Tracking'.

a2_val = 0.4;

% Conversion multipliers
rad_to_deg = 180 / pi;
m_to_mm    = 1000;

% Define and verify destination folder structure
folder_path = fullfile('Figures', 'Tracking');
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

fig_name = sprintf('Control Analysis - Scenario %s (%s)', scenario_id, pert_id);
fig_handle = figure('Name', fig_name, 'Color', 'w', 'Position', [100 100 1200 850]);
sgtitle(sprintf('Simulation Scenario %s | %s', scenario_id, pert_id), 'FontSize', 14, 'FontWeight', 'bold');


% COLUMNS 1: REAL JOINT SPACE TRACKING


% Subplot 1: Joint 1 (Base - Rotation) -> Scaled to Degrees
subplot(3,2,1);
plot(t_sim, q_ref(1,:) * rad_to_deg, 'Color', [0.6 0.6 0.6], 'LineWidth', 4); hold on;
plot(t_sim, X1(:,1) * rad_to_deg, 'b-', 'LineWidth', 1.5);
plot(t_sim, X2(:,1) * rad_to_deg, 'g-.', 'LineWidth', 1.5);
ylabel('Joint 1 (deg)'); title('Joint 1 Tracking Performance'); 
grid on; legend('Reference', 'Finite-Horizon LQT', 'Assumed Model', 'Location', 'best');

% Subplot 2: Joint 2 (Shoulder - Elevation) -> Scaled to Degrees
subplot(3,2,3);
plot(t_sim, q_ref(2,:) * rad_to_deg, 'Color', [0.6 0.6 0.6], 'LineWidth', 4); hold on;
plot(t_sim, X1(:,2) * rad_to_deg, 'b-', 'LineWidth', 1.5);
plot(t_sim, X2(:,2) * rad_to_deg, 'g-.', 'LineWidth', 1.5);
ylabel('Joint 2 (deg)'); title('Joint 2 Tracking Performance'); 
grid on;

% Subplot 3: Joint 3 (Prismatic - Insertion) -> Scaled to Millimeters
subplot(3,2,5);
plot(t_sim, q_ref(3,:) * m_to_mm, 'Color', [0.6 0.6 0.6], 'LineWidth', 4); hold on;
plot(t_sim, X1(:,3) * m_to_mm, 'b-', 'LineWidth', 1.5);
plot(t_sim, X2(:,3) * m_to_mm, 'g-.', 'LineWidth', 1.5);
ylabel('Joint 3 (mm)'); title('Joint 3 Tracking Performance'); xlabel('Time (s)');
grid on;


% COLUMN 2: 3D CARTESIAN PATH RECONSTRUCTION + ACTUATOR EFFORTS


% Subplot 4: Spatial End-Effector Trajectory (3D Isometric View)
subplot(3,2,[2,4]);

% Kinematic mapping for Controller 1 (LQT)
x_lqt_cart = zeros(size(t_sim)); y_lqt_cart = zeros(size(t_sim)); z_lqt_cart = zeros(size(t_sim));
for i = 1:length(t_sim)
    q1_i = X1(i,1); q2_i = X1(i,2); q3_i = X1(i,3);
    r_i = a2_val * cos(q2_i) + q3_i * sin(q2_i); 
    x_lqt_cart(i) = cos(q1_i) * r_i; 
    y_lqt_cart(i) = sin(q1_i) * r_i; 
    z_lqt_cart(i) = a2_val * sin(q2_i) - q3_i * cos(q2_i);
end

% Kinematic mapping for Controller 2 (Assumed Model)
x_ma_cart = zeros(size(t_sim)); y_ma_cart = zeros(size(t_sim)); z_ma_cart = zeros(size(t_sim));
for i = 1:length(t_sim)
    q1_i = X2(i,1); q2_i = X2(i,2); q3_i = X2(i,3);
    r_i = a2_val * cos(q2_i) + q3_i * sin(q2_i); 
    x_ma_cart(i) = cos(q1_i) * r_i; 
    y_ma_cart(i) = sin(q1_i) * r_i; 
    z_ma_cart(i) = a2_val * sin(q2_i) - q3_i * cos(q2_i);
end

% Render 3D Path comparison
plot3(vec_ref(1,:), vec_ref(2,:), vec_ref(3,:), 'Color', [0.6 0.6 0.6], 'LineWidth', 4); hold on;
plot3(x_lqt_cart, y_lqt_cart, z_lqt_cart, 'b-', 'LineWidth', 1.5);
plot3(x_ma_cart, y_ma_cart, z_ma_cart, 'g-.', 'LineWidth', 1.5);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(sprintf('Reconstructed Cartesian Space Trajectory | Scenario %s', scenario_id));
axis equal; grid on; view(3);
legend('Reference', 'Finite-Horizon LQT', 'Assumed Model', 'Location', 'best');

% Subplot 5: Real-Time Active Actuator Control Effort (Dual Vertical Axes)
subplot(3,2,6);

% Left Axis: Rotational Motor Torques (Solid lines for LQT, Dashed for Assumed Model)
yyaxis left
plot(t_sim, u_lqt(:,1), 'b-', 'LineWidth', 1.2); hold on;
plot(t_sim, u_lqt(:,2), 'r-', 'LineWidth', 1.2);
plot(t_sim, u_ff(:,1), 'b--', 'LineWidth', 1.2);
plot(t_sim, u_ff(:,2), 'r--', 'LineWidth', 1.2);
ylabel('Motor Torques (Nm)'); set(gca, 'YColor', 'k');

% Right Axis: Prismatic Joint Linear Insertion Force
yyaxis right
plot(t_sim, u_lqt(:,3), 'g-', 'LineWidth', 1.5);
plot(t_sim, u_ff(:,3), 'g--', 'LineWidth', 1.5);
ylabel('Insertion Force F_3 (N)'); set(gca, 'YColor', 'k');

title('Actuator Control Effort'); 
xlabel('Time (s)'); grid on;
hold off;

% Export graphic asset automatically
file_output = sprintf('Tracking_Analysis_Scenario_%s_%s.png', char(scenario_id), char(pert_id));
saveas(fig_handle, fullfile(folder_path, file_output));
end