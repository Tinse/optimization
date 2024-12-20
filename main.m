%% Main script for optimization algorithm comparison
% Author: Tinse
% Date: 2024-12-19

%% Clear workspace and set random seed
clear;
clc;
rng('default');
rng(42);

% Add paths
addpath('algorithms', 'problems', 'utils');

%% Initialize global variables for Knapsack problem
global weights values capacity

% 背包问题参数
weights = [871,827,369,306,435,349,843,270,540,504,838,545,191,771,701,258,946,888,838,924,...
    326,708,873,781,380,5,528,548,838,909,538,954,165,791,469,30,664,429,795,605,...
    52,401,788,342,266,843,239,636,323,902,839,542,11,356,577,901,60,571,695,483,...
    809,523,978,606,550,416,665,335,588,193,239,154,522,613,425,182,728,834,896,...
    933,926,948,584,369,689,306,712,788,863,538,978,556,526,675,980,897,835,77,...
    646,396];

values = [845,758,421,259,512,405,784,304,477,584,909,505,282,756,619,251,910,983,811,903,...
    311,730,899,684,473,101,435,611,914,967,478,866,261,806,549,15,720,399,825,669,...
    2,494,868,244,326,871,192,568,239,968,804,448,81,321,508,933,110,552,707,548,...
    815,541,964,604,588,445,597,385,576,291,190,187,613,657,477,90,758,877,924,...
    843,899,924,541,392,706,276,812,850,896,590,950,580,451,661,997,917,794,83,...
    613,487];

total_weights = sum(weights);
capacity = total_weights * 0.1;  % 背包容量为总体积的10%

%% Parameters
max_iter = 200;      % Increased maximum iterations
pop_size = 50;       % Increased population size
runs = 30;            % Number of independent runs

% Problem dimensions
dim_ackley = 30;      % Dimension for Ackley function
dim_g06 = 2;         % Dimension for G06 problem
dim_knapsack = 100;   % Dimension for Knapsack problem (number of items)

%% Initialize result storage
problems = { 'Ackley', 'G06', 'Knapsack'};
% problems = { 'Knapsack'};
algorithms = {'PSO', 'SecVibratPSO', 'SA', 'GA'};
best_values = zeros(length(problems), length(algorithms), runs);
convergence = zeros(length(problems), length(algorithms), max_iter);
execution_time = zeros(length(problems), length(algorithms));

%% Create results directory if it doesn't exist
if ~exist('results', 'dir')
    mkdir('results');
end

%% Run experiments
for p = 1:length(problems)
    problem = problems{p};
    
    % Set problem dimension
    switch problem
        case 'Ackley'
            dim = dim_ackley;
        case 'G06'
            dim = dim_g06;
        case 'Knapsack'
            dim = dim_knapsack;
    end
    
    % Run each algorithm multiple times
    for a = 1:length(algorithms)
        algorithm = algorithms{a};
        
        % Initialize timer
        tic;
        
        % Run multiple times
        for r = 1:runs
            % Run the selected algorithm
            switch algorithm
                case 'PSO'
                    [best_val, conv_curve] = PSO(problem, dim, pop_size, max_iter);
                case 'SecVibratPSO'
                    [best_val, conv_curve] = SecVibratPSO(problem, dim, pop_size, max_iter);
                case 'SA'
                    [best_val, conv_curve] = SA(problem, dim, pop_size, max_iter);
                case 'GA'
                    [best_val, conv_curve] = GA(problem, dim, pop_size, max_iter);
            end
            
            % Store results
            best_values(p, a, r) = best_val;
            if r == 1  % Only store convergence curve for first run
                convergence(p, a, :) = conv_curve;
            end
        end
        
        execution_time(p, a) = toc; % End timing
    end
    
    % Plot convergence curves for current problem
    plotConvergence(problem, algorithms, squeeze(convergence(p,:,:)));
    
    % Plot statistical results
    plotStatistics(problem, algorithms, squeeze(best_values(p,:,:)));
end

%% Save results
save('results/optimization_results.mat', 'best_values', 'convergence', 'execution_time');

%% Display final results
displayResults(problems, algorithms, best_values, execution_time);
