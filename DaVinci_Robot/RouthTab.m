function rh_table = RouthTab(coeffs)
% ROUTHTAB Generates the Routh array for a given polynomial
%
% Input:
%   coeffs - Vector of polynomial coefficients [an, an-1, ..., a0]
%            (Output of poly(A))
%
% Output:
%   rh_table - The Routh matrix

n = length(coeffs); % Number of coefficients
l = ceil(n/2);      % Number of columns

rh_table = zeros(n, l);

% Initialize first two rows
rh_table(1, :) = [coeffs(1:2:n), zeros(1, l - length(coeffs(1:2:n)))];
rh_table(2, :) = [coeffs(2:2:n), zeros(1, l - length(coeffs(2:2:n)))];

% Calculate remaining rows
for i = 3:n
    for j = 1:l-1
        % Determinant logic for Routh array
        %      | a  b |
        %  x = | c  d |  -> (c*b - a*d) / c
        
        a = rh_table(i-2, 1);
        b = rh_table(i-2, j+1);
        c = rh_table(i-1, 1);
        d = rh_table(i-1, j+1);
        
        if c == 0
            c = 1e-6; % Avoid division by zero (epsilon method)
        end
        
        rh_table(i, j) = (c*b - a*d) / c;
    end
end


end
