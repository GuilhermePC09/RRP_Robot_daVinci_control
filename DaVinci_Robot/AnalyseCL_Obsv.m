function AnalyseCL_Obsv(A, C, L, obs_red, poles_control, x0_real, t_sim, state_labels, x0_hat)
% Tabela de polos, convergencia de estados e mapa global.
%
% ---------------------------------------------------------------------------
% ENTRADAS
%   A             : matriz de estados da planta            (n x n)
%   C             : matriz de saida                        (m x n)
%   L             : ganho do observador de ordem completa  (n x m)
%   obs_red       : struct com o observador reduzido (Friedland), campos:
%                     .F_red        (n-m) x (n-m)    dinamica do erro reduzido
%                     .A_aug_red    (2n-m) x (2n-m)  sistema aumentado [x; z]
%                     .recover_xhat n x (2n-m)       reconstroi x_hat
%                     .J            (n-m) x m         ganho do observador reduzido
%   poles_control : polos de malha fechada do controlador (LQR/PP), vetor n x 1
%                   (no script original: poles_cl)
%   x0_real       : estado inicial REAL da planta          (n x 1)
%   t_sim         : vetor de tempo da simulacao
%   state_labels  : (opcional) cell com rotulos dos estados. Default {'x_1'..}
%   x0_hat        : (opcional) estimativa inicial do estado (n x 1).
%                   Default = zeros(n,1).
%
% ---------------------------------------------------------------------------


    % Dimensoes e defaults
    n_states  = size(A, 1);
    n_outputs = size(C, 1);
    n_red     = n_states - n_outputs;
    n_aug     = 2*n_states - n_outputs;

    if nargin < 8 || isempty(state_labels)
        state_labels = arrayfun(@(k) sprintf('x_%d', k), 1:n_states, ...
                                'UniformOutput', false);
    end
    if nargin < 9 || isempty(x0_hat)
        x0_hat = zeros(n_states, 1);
    end

    % Desempacota o observador reduzido 
    F_red        = obs_red.F_red;
    A_aug_red    = obs_red.A_aug_red;
    recover_xhat = obs_red.recover_xhat;
    J            = obs_red.J;

    % Checagens de consistencia 
    assert(isequal(size(F_red), [n_red, n_red]), ...
        'obs_red.F_red deve ser (n-m)x(n-m) = %dx%d.', n_red, n_red);
    assert(isequal(size(A_aug_red), [n_aug, n_aug]), ...
        'obs_red.A_aug_red deve ser (2n-m)x(2n-m) = %dx%d.', n_aug, n_aug);
    assert(size(recover_xhat,1) == n_states && size(recover_xhat,2) == n_aug, ...
        'obs_red.recover_xhat deve ser n x (2n-m) = %dx%d.', n_states, n_aug);

    A_obs          = A - L*C;
    poles_obs_full = eig(A_obs);
    poles_obs_red  = eig(F_red);
    poles_control  = poles_control(:);


    % TABELA COMPARATIVA DE POLOS (PRINCIPIO DA SEPARACAO)

    fmt = @(p) sprintf('%+6.2f %+6.2fi', real(p), imag(p));
    dash = '       -       ';

    fprintf('\n===================================================================\n');
    fprintf('        TABELA COMPARATIVA DE POLOS DO SISTEMA GLOBAL \n');
    fprintf('===================================================================\n');
    fprintf(' Modo |   Controlador (MF)    |   Obs. Completo   |   Obs. Reduzido   \n');
    fprintf('-------------------------------------------------------------------\n');
    for i = 1:n_states
        if i <= numel(poles_control),  s_ctrl = fmt(poles_control(i));  else, s_ctrl = dash; end
        if i <= numel(poles_obs_full), s_full = fmt(poles_obs_full(i)); else, s_full = dash; end
        if i <= numel(poles_obs_red),  s_red  = fmt(poles_obs_red(i));  else, s_red  = dash; end
        fprintf('  %d   | %-21s | %-17s | %-17s \n', i, s_ctrl, s_full, s_red);
    end
    fprintf('===================================================================\n\n');


    % CONVERGENCIA DOS ESTADOS: PLANTA VS ESTIMADO (MALHA ABERTA)
    e0 = x0_real - x0_hat;

    % Erro do observador completo -> reconstroi x_hat_pp = x_real - e_pp
    sys_err_pp = ss(A_obs, zeros(n_states,1), eye(n_states), zeros(n_states,1));
    [e_pp, ~]  = initial(sys_err_pp, e0, t_sim);

    % Observador reduzido (planta + observador em paralelo)
    z0 = -J * C * x0_real;
    sys_aug_open = ss(A_aug_red, zeros(n_aug,1), eye(n_aug), zeros(n_aug,1));
    [X_aug, ~]   = initial(sys_aug_open, [x0_real; z0], t_sim);

    x_real   = X_aug(:, 1:n_states);
    xhat_red = (recover_xhat * X_aug.').';
    xhat_pp  = x_real - e_pp;

    for k = 1:n_states
        figure('Name', ['Convergencia de Estado - ' state_labels{k}], 'Color', 'w');
        plot(t_sim, x_real(:, k),   'r-.', 'LineWidth', 2.5); hold on;
        plot(t_sim, xhat_pp(:, k),  'b--', 'LineWidth', 1.5);
        plot(t_sim, xhat_red(:, k), 'g-',  'LineWidth', 1.5);
        grid on;
        xlabel('Tempo (s)'); ylabel(state_labels{k});
        title(['Rastreamento e Convergencia do Estado: ' state_labels{k}]);
        legend('Planta Real', 'Obs Ordem Completa', 'Obs Reduzido', 'Location', 'best');
    end

    % MAPA GLOBAL DE POLOS EXPANDIDO

    figure('Name', 'Mapa Global de Polos do Sistema', 'Color', 'w');
    hold on; grid on;
    plot(real(poles_control),  imag(poles_control),  'rx', 'MarkerSize', 10, 'LineWidth', 2);
    plot(real(poles_obs_full), imag(poles_obs_full), 'bs', 'MarkerSize', 8,  'LineWidth', 1.5);
    plot(real(poles_obs_red),  imag(poles_obs_red),  'g+', 'MarkerSize', 8,  'LineWidth', 1.5);
    xline(0, 'k--'); yline(0, 'k--');
    xlabel('Eixo Real (1/s)'); ylabel('Eixo Imaginario (1/s)');
    title('Mapa de Polos Global: Principio da Separacao');
    legend('Controlador (Malha Fechada)', ...
           sprintf('Obs Ordem Completa (n=%d)', n_states), ...
           sprintf('Obs Reduzido (n=%d)', n_red), 'Location', 'best');
end