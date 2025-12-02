function [] = PlotSim(t_traj, vec_ref, q_ref, tau_feedforward, x_lin, t_non, x_non)
%PLOTSIM Plot time domain simulations for a trajectory
 
a2_val = 0.4;

figure('Name', 'Simulation Results Analysis', 'Color', 'w', 'Position', [100 100 1200 800]);

% COLUMN 1: JOINT TRACKING + ERROR ---

% Subplot 1: Joint 1 (Base)
subplot(3,2,1);
plot(t_traj, q_ref(1,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 4); hold on;
plot(t_non, x_non(:,1), 'r-', 'LineWidth', 1.5);
plot(t_traj, x_lin(:,1), 'b--', 'LineWidth', 1.5);
ylabel('Joint 1 (rad)'); title('Joint 1 Tracking'); 
grid on; legend('Ref', 'Linear', 'Nonlinear', 'Location', 'best');

% Subplot 2: Joint 2 (Arm)
subplot(3,2,3);
plot(t_traj, q_ref(2,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 4); hold on;
plot(t_non, x_non(:,2), 'r-', 'LineWidth', 1.5);
plot(t_traj, x_lin(:,2), 'b--', 'LineWidth', 1.5);
ylabel('Joint 2 (rad)'); title('Joint 2 Tracking'); 
grid on;

% Subplot 3: Joint 3 (Prismatic)
subplot(3,2,5);
plot(t_traj, q_ref(3,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 4); hold on;
plot(t_non, x_non(:,3), 'r-', 'LineWidth', 1.5);
plot(t_traj, x_lin(:,3), 'b--', 'LineWidth', 1.5);
ylabel('Joint 3 (m)'); title('Joint 3 Tracking'); xlabel('Time (s)');
grid on;


% COLUMN 2: TRAJECTORY + TORQUES

% Subplot 4: Cartesian Path (3D View)
subplot(3,2,[2,4]);

% Linear model
x_lin_cart = zeros(size(t_traj));
y_lin_cart = zeros(size(t_traj)); 
z_lin_cart = zeros(size(t_traj));
for i=1:length(t_traj)
    q1_i = x_lin(i,1); q2_i = x_lin(i,2); q3_i = x_lin(i,3);
    r_i = a2_val*cos(q2_i) + q3_i*sin(q2_i); 
    x_lin_cart(i) = cos(q1_i) * r_i;
    y_lin_cart(i) = sin(q1_i) * r_i;
    z_lin_cart(i) = a2_val*sin(q2_i) - q3_i*cos(q2_i);
end


% Nonlinear model
x_non_cart = zeros(size(t_non));
y_non_cart = zeros(size(t_non));
z_non_cart = zeros(size(t_non));
for i=1:length(t_non)
    q1_i = x_non(i,1); q2_i = x_non(i,2); q3_i = x_non(i,3);
    r_i = a2_val*cos(q2_i) + q3_i*sin(q2_i);
    x_non_cart(i) = cos(q1_i) * r_i;
    y_non_cart(i) = sin(q1_i) * r_i;
    z_non_cart(i) = a2_val*sin(q2_i) - q3_i*cos(q2_i);
end


x_ref = vec_ref(1,:);
y_ref = vec_ref(2,:);
z_ref = vec_ref(3,:);

plot3(x_ref, y_ref, z_ref, 'Color', [0.7 0.7 0.7], 'LineWidth', 4); hold on;
plot3(x_lin_cart, y_lin_cart, z_lin_cart, 'g-', 'LineWidth', 1.5);
plot3(x_non_cart, y_non_cart, z_non_cart, 'm--', 'LineWidth', 1.5);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Cartesian Trajectory (3D)');
axis equal; grid on; view(3); legend('Ref', 'Linear', 'Nonlinear');


% Subplot 5: Computed Torques
subplot(3,2,6);
plot(t_traj, tau_feedforward(1,:), 'r', 'LineWidth', 1.5); hold on;
plot(t_traj, tau_feedforward(2,:), 'g', 'LineWidth', 1.5);
plot(t_traj, tau_feedforward(3,:), 'b', 'LineWidth', 1.5);
title('Computed Torques (Inverse Dynamics)'); ylabel('Input'); xlabel('Time (s)');
legend('Tau 1', 'Tau 2', 'Force 3', 'Location', 'best'); grid on;

end

