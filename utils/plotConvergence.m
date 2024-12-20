function plotConvergence(problem, algorithms, convergence_data)
%PLOTCONVERGENCE Plot convergence curves for all algorithms
%   Inputs:
%       problem: name of the optimization problem
%       algorithms: cell array of algorithm names
%       convergence_data: matrix of convergence curves (algorithms x iterations)

figure('Name', sprintf('%s Convergence Comparison', problem));
hold on;

colors = {'b', 'r', 'g', 'm'};
styles = {'-', '--', ':', '-.'};

for i = 1:length(algorithms)
    plot(convergence_data(i,:), [colors{i} styles{i}], 'LineWidth', 1.5, ...
         'DisplayName', algorithms{i});
end

xlabel('Iteration');
ylabel('Objective Value');
title(sprintf('%s Problem - Convergence Comparison', problem));
legend('Location', 'best');
grid on;

% Save figure
saveas(gcf, sprintf('results/%s_convergence.png', lower(problem)));
end
