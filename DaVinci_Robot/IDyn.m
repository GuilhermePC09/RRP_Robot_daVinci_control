function [q_ref, dq_ref, ddq_ref, tau_feedforward] = IDyn(num_points, vec_ref, dt)
%IDYN Inverse dynamics function

x_ref = vec_ref(1,:);
y_ref = vec_ref(2,:);
z_ref = vec_ref(3,:);

% Inverse Kinematics Loop
q_ref = zeros(3, num_points);

for k = 1:num_points
    % Get all possible IK solutions for current point
    q_sols = GetIKnumeric(x_ref(k), y_ref(k), z_ref(k));
    
    % Select configuration (Row 1 for d3 > 0)
    q_ref(:, k) = q_sols(1, :)'; 
end

% Kinematics Derivatives (Numerical)
dq_ref = gradient(q_ref, dt);
ddq_ref = gradient(dq_ref, dt);

% Inverse Dynamics (Compute Feedforward Torque)
tau_feedforward = zeros(3, num_points);

for k = 1:num_points
    x_current = [q_ref(:,k); dq_ref(:,k)];
    
    B_now = GetBmatrix(x_current);
    n_now = Getnvector(x_current);
    
    % Inverse Dynamics: Tau = B*ddq + n
    tau_feedforward(:, k) = B_now * ddq_ref(:,k) + n_now;
end


end

