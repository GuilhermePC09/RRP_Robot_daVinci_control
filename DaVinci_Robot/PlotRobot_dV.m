function h = PlotRobot_dV(q, a2)
% PlotRobot Calculates robot kinematics and plots its 3D skeleton.
%   This is a self-contained, purely numeric function.
%
%   Inputs:
%       q   - [3x1] NUMERIC joint vector [q1_val; q2_val; q3_val]
%
%   Outputs:
%       h      - Cell array of plot handles for animation

% 1. Get Joint Positions

% Extract numeric joint variables
q1 = q(1);
q2 = q(2);
q3 = q(3);

% Numeric DH Matrix
auto_dh = @(d, theta, a, alpha) ...
         [cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta);
          sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
          0,           sin(alpha),             cos(alpha),            d;
          0,           0,                      0,                     1];

% Calculate individual transformation matrices
T_1_0 = auto_dh(0, q1, 0, pi/2);      % Frame {1} relative to {0}
T_2_1 = auto_dh(0, q2, a2, pi/2); % Frame {2} relative to {1}
T_3_2 = auto_dh(q3, 0, 0, 0);         % Frame {3} relative to {2}

% Calculate global (cumulative) transformation matrices
T_2_0 = T_1_0 * T_2_1;
T_3_0 = T_2_0 * T_3_2;

% Extract the 3 key points that define the robot
P_base  = [0; 0; 0];           % Origin of frame {0}
P_elbow = T_2_0(1:3, 4);       % Origin of frame {2} (end of link a2)
P_tip   = T_3_0(1:3, 4);       % Origin of frame {3} (end-effector)


% 2. Plot the Robot Skeleton

% Assemble X, Y, Z data for the links
Link1_X = [P_base(1), P_elbow(1)];
Link1_Y = [P_base(2), P_elbow(2)];
Link1_Z = [P_base(3), P_elbow(3)];

Link2_X = [P_elbow(1), P_tip(1)];
Link2_Y = [P_elbow(2), P_tip(2)];
Link2_Z = [P_elbow(3), P_tip(3)];

% 3. Plotting Logic

% Find or create a dedicated figure for the robot
fig_handle = findobj('Type', 'figure', 'Name', 'RRP Robot Visualization', 'Color', 'w');
if isempty(fig_handle)
    fig_handle = figure('Name', 'RRP Robot Visualization',  'Color', 'w');
end
figure(fig_handle);

clf; 

% Plot the first 3D object
h_link1 = plot3(Link1_X, Link1_Y, Link1_Z, 'Color', [0.2 0.2 0.8], 'LineWidth', 5);

hold on; %

% Plot the "Patient Bed" Plane
bed_z = -0.2;
bed_x_inf = 0.2;
bed_x_sup = 0.6;
bed_y_lim = 0.5;
bed_X = [bed_x_sup, bed_x_sup, bed_x_inf, bed_x_inf];
bed_Y = [bed_y_lim, -bed_y_lim, -bed_y_lim, bed_y_lim];
bed_Z = [bed_z, bed_z, bed_z, bed_z];
h_bed = patch(bed_X, bed_Y, bed_Z, [0.7 0.8 1.0], 'FaceAlpha', 0.3, 'EdgeColor', 'blue');

% Plot the rest of the robot
h_link2 = plot3(Link2_X, Link2_Y, Link2_Z, 'Color', [0.8 0.2 0.2], 'LineWidth', 5);
h_junta0 = scatter3(P_base(1), P_base(2), P_base(3), 100, 'k', 'filled', 'Marker', 's');
h_junta1 = scatter3(P_elbow(1), P_elbow(2), P_elbow(3), 100, 'k', 'filled', 'Marker', 'o');
h_junta_tip = scatter3(P_tip(1), P_tip(2), P_tip(3), 150, 'r', 'filled', 'Marker', 'h');

% Format the plot
xlabel('X Axis (m)');
ylabel('Y Axis (m)');
zlabel('Z Axis (m)');
title('RRP Robot 3D Visualization');
axis equal;
grid on;
view(3);

xlim([-0.2, 0.8]);
ylim([-0.5, 0.5]);
zlim([-0.4, 0.2]);

hold off;

% Return plot handles
h = {h_link1, h_link2, h_junta0, h_junta1, h_junta_tip};


end