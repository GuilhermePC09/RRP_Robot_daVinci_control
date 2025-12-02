function [f, x, u] = StateSpaceFunc(B, n)
% StateSpaceModel Generates the symbolic nonlinear state-space model f(x,u).
%
% Inputs:
%   B - Symbolic Mass Matrix B(q)
%   n - Symbolic vector C(q,dq)dq + G(q)
%
% Outputs:
%   f - The symbolic state function vector f(x, u) [6x1]
%   x - The symbolic state vector x [6x1]
%   u - The symbolic input vector u [3x1]
%

n_joints = size(B, 1);

syms q1 q2 q3 real
syms dq1 dq2 dq3 real

vars_q = [q1; q2; q3];
vars_dq = [dq1; dq2; dq3];

% Define the state (x) and input (u) vectors
x = sym('x', [2 * n_joints, 1], 'real'); % [x1..x6]
u = sym('u', [n_joints, 1], 'real');   % [u1..u3]

% Map the states to the (q, dq) variables
x_q = x(1:n_joints);
x_dq = x(n_joints+1:end);

% Substitute the (q, dq) variables with (x) states in B and n
B_x = subs(B, vars_q, x_q);
n_x = subs(n, [vars_q; vars_dq], [x_q; x_dq]);

% 6. Isolate the acceleration (ddq)
% The dynamics equation is: u = B*ddq + n
% Therefore: ddq = B^-1 * (u - n)

% (We use inv(B_x) as it is symbolic)
ddq_x = inv(B_x) * (u - n_x);

% Build the final state-space vector f(x,u)
% f = x_dot = [dq; ddq]

f_top = x_dq;       % The derivative of position is velocity
f_bottom = ddq_x;   % The derivative of velocity is acceleration

f = [f_top; f_bottom];

end