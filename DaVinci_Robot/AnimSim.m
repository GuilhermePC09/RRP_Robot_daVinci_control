function [] = AnimSim(t_non, x_non, vec_ref)
%ANIMSIM Animate non-linear simulations

a2_val = 0.4;
x_ref = vec_ref(1,:);
y_ref = vec_ref(2,:);
z_ref = vec_ref(3,:);

% Setup Video Writer
video_filename = 'Robot_Animation.mp4';
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = 30;
open(v);
disp(['Recording video to: ' video_filename]);

% Animation Time Vector
fps = 30;
duration = t_non(end); 
t_anim = linspace(0, duration, round(duration * fps));

% Interpolate Simulation Data to Uniform Time Steps
q1_anim = interp1(t_non, x_non(:,1), t_anim, 'linear', 'extrap');
q2_anim = interp1(t_non, x_non(:,2), t_anim, 'linear', 'extrap');
q3_anim = interp1(t_non, x_non(:,3), t_anim, 'linear', 'extrap');

% Pre-calculate the "Trace" (The path drawn by the tip)
% Forward Kinematics to find X, Y, Z
trace_x = zeros(size(t_anim));
trace_y = zeros(size(t_anim));
trace_z = zeros(size(t_anim));

for i = 1:length(t_anim)
    r_i = a2_val*cos(q2_anim(i)) + q3_anim(i)*sin(q2_anim(i));
    
    trace_x(i) = cos(q1_anim(i)) * r_i;
    trace_y(i) = sin(q1_anim(i)) * r_i;
    trace_z(i) = a2_val*sin(q2_anim(i)) - q3_anim(i)*cos(q2_anim(i));
end

% Animation Loop
disp('Playing animation...');

for i = 1:length(t_anim)
    % Current Joint Configuration
    q_curr = [q1_anim(i); q2_anim(i); q3_anim(i)];
    PlotRobot_dV(q_curr, a2_val);
    hold on;
    
    % Plot Reference Trajectory (if it exists in workspace)
    if exist('x_ref', 'var') && exist('y_ref', 'var') && exist('z_ref', 'var')
        plot3(x_ref, y_ref, z_ref, 'Color', [0.8 0.8 0.8], 'LineWidth', 1, 'DisplayName', 'Reference');
    end
    
    % Plot the path taken SO FAR in Magenta
    plot3(trace_x(1:i), trace_y(1:i), trace_z(1:i), 'm-', 'LineWidth', 2, 'DisplayName', 'Trace');
    
    % Adjust View 
    %    Auto-adjust view based on scenario or keep standard isometric
    if max(abs(trace_y)) < 0.01 && max(abs(q1_anim)) < 0.01 
        % Mostly planar motion
        view(100, 10); 
    else
        % 3D motion
        view(135, 30);
    end
    
    title(sprintf('Simulation Time: %.2f s', t_anim(i)));
    grid on; axis equal;
    
    % Force draw and wait
    drawnow;

    % Capture frames for video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

disp('Animation finished.');


end

