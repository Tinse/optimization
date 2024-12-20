function displayResults(problems, algorithms, best_values, execution_time)
%DISPLAYRESULTS Display final optimization results
%   Inputs:
%       problems: cell array of problem names
%       algorithms: cell array of algorithm names
%       best_values: 3D matrix of best values (problems x algorithms x runs)
%       execution_time: matrix of execution times (problems x algorithms)

fprintf('\nOptimization Results Summary\n');
fprintf('==========================\n\n');

for p = 1:length(problems)
    fprintf('%s Problem:\n', problems{p});
    fprintf('-------------\n');
    
    % Calculate statistics
    means = mean(squeeze(best_values(p,:,:)), 2);
    stds = std(squeeze(best_values(p,:,:)), 0, 2);
    best_runs = min(squeeze(best_values(p,:,:)), [], 2);
    
    % Display results for each algorithm
    for a = 1:length(algorithms)
        fprintf('%s:\n', algorithms{a});
        fprintf('  Best Value: %.6f\n', best_runs(a));
        fprintf('  Mean ± Std: %.6f ± %.6f\n', means(a), stds(a));
        fprintf('  Time: %.2f seconds\n', execution_time(p,a));
    end
    fprintf('\n');
end

% Save results to file
fid = fopen('results/optimization_summary.txt', 'w');
fprintf(fid, 'Optimization Results Summary\n');
fprintf(fid, '==========================\n\n');

for p = 1:length(problems)
    fprintf(fid, '%s Problem:\n', problems{p});
    fprintf(fid, '-------------\n');
    
    means = mean(squeeze(best_values(p,:,:)), 2);
    stds = std(squeeze(best_values(p,:,:)), 0, 2);
    best_runs = min(squeeze(best_values(p,:,:)), [], 2);
    
    for a = 1:length(algorithms)
        fprintf(fid, '%s:\n', algorithms{a});
        fprintf(fid, '  Best Value: %.6f\n', best_runs(a));
        fprintf(fid, '  Mean ± Std: %.6f ± %.6f\n', means(a), stds(a));
        fprintf(fid, '  Time: %.2f seconds\n', execution_time(p,a));
    end
    fprintf(fid, '\n');
end

fclose(fid);
end
