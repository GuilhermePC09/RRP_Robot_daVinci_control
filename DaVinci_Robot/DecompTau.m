function [B, phi, g_mat] = DecompTau(tau)
% Inputs:
% tau: symbolic expression for tau

% Initialize symbolic variables
syms g real;

% Determine the number of joints from tau
n = length(tau);

% Initialize B, phi, g
B = sym(zeros(n));
phi = sym(zeros(n, 1));
g_mat = sym(zeros(n, 1));

% Create symbolic vectors for velocities and accelerations
ddq = sym('ddq', [n, 1], 'real');

% Iterate over each joint to decompose tau
tau = expand(tau);
for i = 1:n
    tau_i = tau(i);
    
    % Extract the terms involving ddq
    for j = 1:n
        [coeffs_ddq, terms_ddq] = coeffs(tau_i, ddq(j));
        if terms_ddq(1) ~= 1
            B(i, j) = coeffs_ddq(1);
            tau_i = expand(tau_i - B(i, j) * ddq(j));
        end
    end
    
    % Extract the terms involving the gravitational constant g
    [coeffs_g, terms_g] = coeffs(tau_i, g);
    if terms_g(1) ~= 1
        g_mat(i) = coeffs_g(1) * g;
        tau_i = expand(tau_i - g_mat(i));
    end
    
    % The remaining terms are the Coriolis and centrifugal terms
    phi(i) = tau_i;
end

end
