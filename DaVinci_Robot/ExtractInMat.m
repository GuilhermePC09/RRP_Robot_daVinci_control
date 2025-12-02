function [ B ] = ExtractInMat( tau )
% Inputs:
% tau: symbolic expression for tau
% ddq: symbolic vector of joint accelerations [ddq1, ddq2, ..., ddqn]

% Determine the number of joints from tau
n = length(tau);

% Initialize B
B = sym(zeros(n));

% Create symbolic vectors for accelerations
ddq = sym('ddq', [n, 1], 'real');

% Iterate over each joint to extract B
for i = 1:n
    tau_i = tau(i);
    
    % Extract the terms involving ddq
    for j = 1:n
        [coeffs_ddq, terms_ddq] = coeffs(tau_i, ddq(j));
        if length(terms_ddq) == 2
            B(i, j) = coeffs_ddq(1);
            tau_i = tau_i - B(i, j) * ddq(j);
        end
    end
end

end

