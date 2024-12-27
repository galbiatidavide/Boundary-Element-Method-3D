
% Subplot configuration
figure;
colors = ['r', 'b']; % Colors for errors and errors_quad

name = [{'u_n boundary'}, {'u boundary'}, {'u domain'}];


for j = 1:3
    % Create subplot for this iteration
    subplot(1, 3, j);
    
    % Extract data for this iteration
    errors = squeeze(errors_linear_to_squeeze(:,:,j));
    errors_quadratic = squeeze(errors_quadratic_to_squeeze(:,:,j));
    
    % Plot errors and errors_quad with log-log scaling
    loglog(sizes, errors/errors(1), 'o-', 'Color', colors(1), 'DisplayName', 'Linear', 'LineWidth', 2);
    hold on;
    loglog(sizes_quadratic, errors_quadratic/errors_quadratic(1), 's--', 'Color', colors(2), 'DisplayName', 'Quadratic', 'LineWidth', 2);

    % Plot theoretical lines with log-log scaling
    sizes_fine = logspace(log10(min(sizes)), log10(max(sizes)), 100); % Finer grid for theoretical lines
    orders = [0.5, 1, 2, 3];
    lineStyles = {'--', '-.', ':', '-'}; % Different line styles for each order
    for i = 1:length(orders)
        order = orders(i);
        loglog(sizes_fine, (sizes_fine/sizes_fine(end)).^order, ...
               'k', 'LineStyle', lineStyles{i}, 'DisplayName', [num2str(order), '-order']);
    end
    % Annotations
    xlabel('Mesh Size');
    ylabel('Relative Error Norm 2');
    title(name(j));
    legend('show', 'Location', 'Best');
    grid on;
    hold off;
end


% Overall figure title
sgtitle('Convergence Analysis in Log-Log Scale');

warning('The reason why some of the relative errors are null may be the fact that you imposed the values via bcs')
