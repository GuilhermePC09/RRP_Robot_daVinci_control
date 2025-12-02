function [ tau ] = DDyn_NE(Robot, M, CoM, I, g_vec)
% DDYN Direct Dynamics following forward Newton-Euler approach
%   * Forward recursion relative to the propagation of velocities and accelerations
%   * Backward recursion for the propagation of forces and moments
%
%   Inputs:
%       Robot - Robot D-H table (symbolic)
%       M     - [n x 1] vector of link masses [m1; m2; ...]
%       CoM   - [n x 3] matrix of link CoM vectors [rc1; rc2; ...]
%       I     - {n x 1} cell array of 3x3 inertia tensors {I1; I2; ...}
%       g_vec - [3 x 1] gravity vector (e.g., [0; 0; -g])
%
%
%   Output:
%       tau   - [n x 1] symbolic vector of joint torques/forces



% Extract number of links from the robot DH table
n = size(Robot, 1);

% Initialize symbolic variables
syms dq [n, 1] real;
syms ddq [n, 1] real;

% Initialize variables
z0 = sym([0; 0; 1]); % z axis of the base frame

% Preallocate arrays
w = sym(zeros(3, n+1)); % Angular velocity
wdot = sym(zeros(3, n+1)); % Angular acceleration
vdot = sym(zeros(3, n+1)); % Linear acceleration
vcdot = sym(zeros(3, n)); % Linear acceleration of CoM
vdot(:,1) = -g_vec;

f = sym(zeros(3, n+1)); % Force
mu = sym(zeros(3, n+1)); % Force
tau = sym(zeros(n, 1)); % Torque (adjusted to be 1D array)


% Forward recursion: Compute velocities and accelerations
for i = 1:n

    % Compute transformation matrix from frame i-1 to i
    p = Robot(i,:);
    T = DHTransf(p);
    
    % Extract rotation matrix2
    R = T(1:3, 1:3);
    r_i = [p(3); 0; p(1)];

    if ~isempty(symvar(p(2))) % Revolute joint 
        w(:, i+1) = R' * (w(:, i) + dq(i) * z0);
        wdot(:, i+1) = R' * (wdot(:, i) + ddq(i) * z0 + dq(i) * cross(w(:, i), z0));
        vdot(:, i+1) = R' * vdot(:, i) + cross(wdot(:, i+1), r_i) + cross(w(:, i+1), cross(w(:, i+1), r_i));
        % vdot(:, i+1) = R' * vdot(:, i) + cross(wdot(:, i), r_i) + cross(w(:, i), cross(w(:, i), r_i));


    else % Prismatic joint 
        w(:, i+1) = R' * w(:, i);
        wdot(:, i+1) = R' * wdot(:, i);
        vdot(:, i+1) = R' * (vdot(:, i) + ddq(i) * z0) + cross(2 * dq(i) * w(:, i+1), R' * z0) + cross(wdot(:, i+1), r_i) + cross(w(:, i+1), cross(w(:, i+1), r_i));
        
  
    end

    vcdot(:, i) = vdot(:, i+1) + cross(wdot(:, i+1), CoM(i, :)') + cross(w(:, i+1), cross(w(:, i+1), CoM(i, :)'));

end


% Backward recursion: Compute forces and torques
for i = n:-1:1

    % Compute position matrix from frame i-1 to i
    p1 = Robot(i,:);
    r_i = [p1(3); 0; p1(1)];

    % Compute transformation matrix from frame i to i+1 and extract R
    if i < n 
        p2 = Robot(i+1,:);
        T1 = DHTransf(p2);
        R = T1(1:3, 1:3);
    else
        R = eye(3);
    end

    % Compute forces and torques
    f(:, i) = R * f(:, i+1) + M(i) * vcdot(:, i);
    % mu(:, i) = - cross(f(:, i), r_i + CoM(i, :)') + R * mu(:, i+1) + cross(R * f(:, i+1), CoM(i, :)') + ...
    % I{i} * wdot(:, i+1) + cross(w(:, i+1), I{i} * w(:, i+1));

    mu(:, i) = I{i} * wdot(:, i+1) + cross(w(:, i+1), I{i} * w(:, i+1)) + ...
        R * mu(:, i+1) + ...
        cross(CoM(i, :)', M(i) * vcdot(:, i)) + ...
        cross(r_i, R * f(:, i+1));

    p = Robot(i,:);
    T = DHTransf(p);
    R = T(1:3, 1:3);

    if ~isempty(symvar(p(2))) % Revolute joint 
        tau(i) = mu(:, i)' * R' * z0;

    else % Prismatic joint
        tau(i) = f(:, i)' * R' * z0;
    end

end


end
