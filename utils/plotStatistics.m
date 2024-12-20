function plotStatistics(problem, algorithms, results)
%PLOTSTATISTICS Plot statistical results for all algorithms
%   Inputs:
%       problem: name of the optimization problem
%       algorithms: cell array of algorithm names
%       results: matrix of results (algorithms x runs)

% Calculate statistics
means = mean(results, 2);
stds = std(results, 0, 2);

% Create bar plot
figure('Name', sprintf('%s Statistical Comparison', problem));

% Plot means with error bars
bar_h = bar(means);
hold on;
errorbar(1:length(algorithms), means, stds, 'k.', 'LineWidth', 1.5);

% Customize plot
set(gca, 'XTick', 1:length(algorithms));
set(gca, 'XTickLabel', algorithms);
xlabel('Algorithm');
ylabel('Objective Value');
title(sprintf('%s Problem - Statistical Comparison', problem));
grid on;

% Add value labels on top of bars
for i = 1:length(means)
    text(i, means(i), sprintf('%.2f±%.2f', means(i), stds(i)), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom');
end

% Save figure
saveas(gcf, sprintf('results/%s_statistics.png', lower(problem)));
end
