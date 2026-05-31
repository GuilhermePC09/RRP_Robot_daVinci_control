function [A, B, C, D, E] = LinModel(f_ss, x_ss, u_ss, w_ss, G, n)
% linearizeModel Linearizes the nonlinear state-space model f(x,u)
%                and h(x,u) around a defined operating point.
%
% Inputs:
%   f_ss  - Symbolic state function f(x,u) [6x1]
%   x_ss  - Symbolic state vector x [6x1]
%   u_ss  - Symbolic input vector u [3x1]
%   w_ss - Symbolic perturbation scalar w [1x1]
%   g_mat - Symbolic gravity vector G(q) [3x1]
%
% Outputs:
%   A, B, C, D, E - The symbolic linearized state-space matrices


% Output equation y = h(x,u)
% y = [x1; x2; x3]
h_ss = x_ss(1:3); % y = q

% Symbolic equilibrium point
% q_bar = [0, 0, q3_bar] (non-singular, horizontal pose)
% dq_bar = [0, 0, 0] (at rest)
syms q3_bar real
q_bar = [0; 0; q3_bar];
dq_bar = zeros(3, 1);

x_bar = [q_bar; dq_bar];

% Equilibrium input u_bar = G(q_bar)
syms q1 q2 q3 real % The original q variables
vars_q_orig = [q1; q2; q3];
u_bar = subs(G, vars_q_orig, q_bar);


% 3. Calculate Symbolic Jacobians (A, B, C, D)
% Taylor expansion

A_sym = jacobian(f_ss, x_ss);

B_sym = jacobian(f_ss, u_ss);

C_sym = jacobian(h_ss, x_ss); % C = d(h)/d(x)

D_sym = jacobian(h_ss, u_ss); % D = d(h)/d(u)

E_sym = jacobian(f_ss, w_ss);

% 4. Evaluate Jacobians at the Equilibrium Point
A = subs(A_sym, [x_ss; u_ss; w_ss], [x_bar; u_bar; 0]);
B = subs(B_sym, [x_ss; u_ss; w_ss], [x_bar; u_bar; 0]);
C = subs(C_sym, [x_ss; u_ss; w_ss], [x_bar; u_bar; 0]);
D = subs(D_sym, [x_ss; u_ss; w_ss], [x_bar; u_bar; 0]);
E = subs(E_sym, [x_ss; u_ss; w_ss], [x_bar; u_bar; 0]);

% 5. Simplify and finalize
A = simplify(A);
B = simplify(B);
C = simplify(C);
D = simplify(D); 
E = simplify(E);

end