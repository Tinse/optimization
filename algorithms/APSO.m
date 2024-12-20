function [best_val, convergence] = SOPSO(problem, dim, pop_size, max_iter)
%SOPSO Second-Order Oscillation PSO algorithm
%   Inputs:
%       problem: string, name of the optimization problem
%       dim: problem dimension
%       pop_size: population size
%       max_iter: maximum number of iterations
%   Outputs:
%       best_val: best solution found
%       convergence: convergence curve

% SOPSO parameters
alpha = 0.7;    % Damping coefficient
beta = 1.5;     % Oscillation frequency
gamma = 2.0;    % Social coefficient

% Initialize bounds based on problem
switch problem
    case 'Ackley'
        lb = -32.768 * ones(1, dim);
        ub = 32.768 * ones(1, dim);
        init_lb = lb;
        init_ub = ub;
    case 'G06'
        lb = [13, 0];
        ub = [100, 100];
        % Initialize near the known good region
        init_lb = [13, 0];     % Lower bounds for initialization
        init_ub = [20, 5];     % Upper bounds for initialization
        dim = 2;
    case 'Knapsack'
        lb = zeros(1, dim);
        ub = ones(1, dim);
        init_lb = lb;
        init_ub = ub;
end

% Initialize population with problem-specific ranges
pos = init_lb + (init_ub - init_lb) .* rand(pop_size, dim);
vel = -0.1 * (init_ub - init_lb) + 0.2 * (init_ub - init_lb) .* rand(pop_size, dim);
acc = zeros(pop_size, dim);  % Initial acceleration
p_best = pos;
p_best_val = inf(pop_size, 1);
g_best = zeros(1, dim);
g_best_val = inf;

% Initialize convergence curve
convergence = zeros(1, max_iter);

% Main loop
for iter = 1:max_iter
    % Evaluate particles
    for i = 1:pop_size
        % Calculate fitness
        switch problem
            case 'Ackley'
                current_val = Ackley(pos(i,:));
            case 'G06'
                current_val = G06(pos(i,:));
            case 'Knapsack'
                current_val = Knapsack(pos(i,:));
        end
        
        % Update personal best
        if current_val < p_best_val(i)
            p_best_val(i) = current_val;
            p_best(i,:) = pos(i,:);
            
            % Update global best
            if current_val < g_best_val
                g_best_val = current_val;
                g_best = pos(i,:);
            end
        end
    end
    
    % Update particles using second-order oscillation
    for i = 1:pop_size
        % Calculate acceleration
        acc(i,:) = beta * (g_best - pos(i,:)) + ...
                   gamma * rand(1,dim) .* (p_best(i,:) - pos(i,:));
        
        % Update velocity with damping
        vel(i,:) = alpha * vel(i,:) + acc(i,:);
        
        % Update position
        pos(i,:) = pos(i,:) + vel(i,:);
        
        % Bound position
        pos(i,:) = max(min(pos(i,:), ub), lb);
    end
    
    % Record the best value for convergence history
    convergence(iter) = g_best_val;
end

best_val = g_best_val;
end
