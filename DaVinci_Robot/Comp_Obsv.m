function Comp_Obsv(A, C, L, obs_red, x0_real, t_sim, state_names, x0_hat)
%Compares error dynamics and poles of two observers.
% ---------------------------------------------------------------------------
% INPUTS
%   A           : plant state matrix                      (n x n)
%   C           : output matrix                           (m x n)
%   L           : full-order observer gain                (n x m)
%   obs_red     : struct with reduced-order observer (Friedland), fields:
%                   .F_red        (n-m) x (n-m)    reduced error dynamics
%                   .A_aug_red    (2n-m) x (2n-m)  augmented system [x; z]
%                   .recover_xhat n x (2n-m)       reconstructs x_hat = [S*C, N]*[x;z]
%                   .J            (n-m) x m         reduced-order observer gain
%   x0_real     : REAL initial state of the plant         (n x 1)
%   t_sim       : simulation time vector
%   state_names : (optional) cell array with state names for labels.
%                 Default = {'x_1', ..., 'x_n'}
%   x0_hat      : (optional) initial state estimate (n x 1).
%                 Default = zeros(n,1)  -> "uninformed" observer
%
% ---------------------------------------------------------------------------

    % Dimensions and defaults
    n_states  = size(A, 1);
    n_outputs = size(C, 1);

    if nargin < 7 || isempty(state_names)
        state_names = arrayfun(@(k) sprintf('x_%d', k), 1:n_states, ...
                               'UniformOutput', false);
    end
    if nargin < 8 || isempty(x0_hat)
        x0_hat = zeros(n_states, 1);
    end

    % Unpack the reduced-order observer
    F_red        = obs_red.F_red;
    A_aug_red    = obs_red.A_aug_red;
    recover_xhat = obs_red.recover_xhat;
    J            = obs_red.J;

    n_aug = 2*n_states - n_outputs;   % order of the augmented system [x; z]

    % Basic consistency checks
    assert(isequal(size(F_red), [n_states-n_outputs, n_states-n_outputs]), ...
        'obs_red.F_red deve ser (n-m)x(n-m) = %dx%d.', ...
        n_states-n_outputs, n_states-n_outputs);
    assert(isequal(size(A_aug_red), [n_aug, n_aug]), ...
        'obs_red.A_aug_red deve ser (2n-m)x(2n-m) = %dx%d.', n_aug, n_aug);
    assert(size(recover_xhat,1) == n_states && size(recover_xhat,2) == n_aug, ...
        'obs_red.recover_xhat deve ser n x (2n-m) = %dx%d.', n_states, n_aug);

    %  Open-loop error convergence

    % Full-order observer:  de/dt = (A - L*C)e
    A_obs = A - L*C;
    e0    = x0_real - x0_hat;

    sys_err_pp = ss(A_obs, zeros(n_states,1), eye(n_states), zeros(n_states,1));
    [e_pp, ~]  = initial(sys_err_pp, e0, t_sim);

    % Reduced-order observer: plant + observer in parallel
    %   z0 = -J*C*x0_real
    z0 = -J * C * x0_real;

    sys_aug_open = ss(A_aug_red, zeros(n_aug,1), eye(n_aug), zeros(n_aug,1));
    [X_aug, ~]   = initial(sys_aug_open, [x0_real; z0], t_sim);

    x_real   = X_aug(:, 1:n_states);
    xhat_red = (recover_xhat * X_aug.').';
    e_red    = x_real - xhat_red;

    % Plots per state
    for k = 1:n_states
        figure('Name', ['Erro de estimacao - ' state_names{k}], 'Color', 'w');
        plot(t_sim, e_pp(:, k),  'b-',  'LineWidth', 1.6); hold on;
        plot(t_sim, e_red(:, k), 'g-.', 'LineWidth', 1.6);
        grid on;
        xlabel('Tempo [s]');
        ylabel(['e_{' state_names{k} '}']);
        title(['Erro de estimacao: ' state_names{k} ' (malha aberta, sem entrada)']);
        legend('PP (ordem completa)', 'Reduzido', 'Location', 'best');
    end

    %  Observer pole map

    figure('Name', 'Polos dos observadores', 'Color', 'w');
    hold on; grid on;
    plot(real(eig(A_obs)), imag(eig(A_obs)), 'bx', ...
         'MarkerSize', 12, 'LineWidth', 2);
    plot(real(eig(F_red)), imag(eig(F_red)), 'gs', ...
         'MarkerSize', 12, 'LineWidth', 2);
    xline(0, 'k--'); yline(0, 'k--');
    xlabel('Re'); ylabel('Im');
    title('Polos: observador PP vs observador reduzido');
    legend(sprintf('Obs identidade (n=%d)', n_states), ...
           sprintf('Obs Reduzido (n=%d)', n_states - n_outputs), ...
           'Location', 'best');
    axis equal;
end