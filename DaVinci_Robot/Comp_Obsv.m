function Comp_Obsv(A, C, L, obs_red, x0_real, t_sim, state_names, x0_hat)
%Compara dinamica de erro e polos de dois observadores.
% ---------------------------------------------------------------------------
% ENTRADAS
%   A           : matriz de estados da planta            (n x n)
%   C           : matriz de saida                        (m x n)
%   L           : ganho do observador de ordem completa  (n x m)
%   obs_red     : struct com o observador reduzido (Friedland), campos:
%                   .F_red        (n-m) x (n-m)    dinamica do erro reduzido
%                   .A_aug_red    (2n-m) x (2n-m)  sistema aumentado [x; z]
%                   .recover_xhat n x (2n-m)       reconstroi x_hat = [S*C, N]*[x;z]
%                   .J            (n-m) x m         ganho do observador reduzido
%   x0_real     : estado inicial REAL da planta          (n x 1)
%   t_sim       : vetor de tempo da simulacao
%   state_names : (opcional) cell com nomes dos estados p/ rotulos.
%                 Default = {'x_1', ..., 'x_n'}
%   x0_hat      : (opcional) estimativa inicial do estado (n x 1).
%                 Default = zeros(n,1)  -> observador "ignorante"
%
% ---------------------------------------------------------------------------

    % Dimensoes e defaults 
    n_states  = size(A, 1);
    n_outputs = size(C, 1);

    if nargin < 7 || isempty(state_names)
        state_names = arrayfun(@(k) sprintf('x_%d', k), 1:n_states, ...
                               'UniformOutput', false);
    end
    if nargin < 8 || isempty(x0_hat)
        x0_hat = zeros(n_states, 1);
    end

    % Desempacota o observador reduzido
    F_red        = obs_red.F_red;
    A_aug_red    = obs_red.A_aug_red;
    recover_xhat = obs_red.recover_xhat;
    J            = obs_red.J;

    n_aug = 2*n_states - n_outputs;   % ordem do sistema aumentado [x; z]

    % Checagens basicas de consistencia
    assert(isequal(size(F_red), [n_states-n_outputs, n_states-n_outputs]), ...
        'obs_red.F_red deve ser (n-m)x(n-m) = %dx%d.', ...
        n_states-n_outputs, n_states-n_outputs);
    assert(isequal(size(A_aug_red), [n_aug, n_aug]), ...
        'obs_red.A_aug_red deve ser (2n-m)x(2n-m) = %dx%d.', n_aug, n_aug);
    assert(size(recover_xhat,1) == n_states && size(recover_xhat,2) == n_aug, ...
        'obs_red.recover_xhat deve ser n x (2n-m) = %dx%d.', n_states, n_aug);

    %  Convergencia do erro em malha aberta
    
    % Observador de ordem completa:  de/dt = (A - L*C)e 
    A_obs = A - L*C;
    e0    = x0_real - x0_hat;

    sys_err_pp = ss(A_obs, zeros(n_states,1), eye(n_states), zeros(n_states,1));
    [e_pp, ~]  = initial(sys_err_pp, e0, t_sim);

    % Observador reduzido: planta + observador em paralelo 
    %   z0 = -J*C*x0_real
    z0 = -J * C * x0_real;

    sys_aug_open = ss(A_aug_red, zeros(n_aug,1), eye(n_aug), zeros(n_aug,1));
    [X_aug, ~]   = initial(sys_aug_open, [x0_real; z0], t_sim);

    x_real   = X_aug(:, 1:n_states);
    xhat_red = (recover_xhat * X_aug.').';
    e_red    = x_real - xhat_red;

    % Graficos por estado
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

    %  Mapa de polos dos observadores

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