function [best_val, convergence] = SA(problem, dim, pop_size, max_iter)
%SA Simulated Annealing algorithm
%   Inputs:
%       problem: string, name of the optimization problem
%       dim: dimension of the problem
%       pop_size: not used in SA, kept for interface consistency
%       max_iter: maximum number of iterations
%   Outputs:
%       best_val: best solution found
%       convergence: convergence curve

% For Knapsack problem
global weights values capacity

% SA parameters
T0 = 100;          % Initial temperature
alpha = 0.95;      % Cooling rate
inner_iter = 20;   % Number of iterations at each temperature

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

% Initialize solution near the promising region
current_pos = init_lb + (init_ub - init_lb) .* rand(1, dim);

switch problem
    case 'Ackley'
        current_val = Ackley(current_pos);
    case 'G06'
        current_val = G06(current_pos);
    case 'Knapsack'
        current_val = Knapsack(current_pos);
end

best_pos = current_pos;
best_val = current_val;

% Initialize convergence curve
convergence = zeros(1, max_iter);

% Initialize temperature
T = T0;
iter = 1;

% Main loop
while iter <= max_iter
    % Inner loop
    for j = 1:inner_iter
        if iter > max_iter
            break;
        end
        
        % Generate new solution
        if strcmp(problem, 'G06')
            % Special handling for G06 problem
            step = zeros(1, dim);
            for d = 1:dim
                if d == 1
                    % For x1: 13 <= x1 <= 100
                    range = min(100 - current_pos(d), current_pos(d) - 13) * 0.1;
                    step(d) = range * randn;
                    new_pos(d) = current_pos(d) + step(d);
                    new_pos(d) = max(min(new_pos(d), 100), 13);
                else
                    % For x2: 0 <= x2 <= 100
                    range = min(100 - current_pos(d), current_pos(d)) * 0.1;
                    step(d) = range * randn;
                    new_pos(d) = current_pos(d) + step(d);
                    new_pos(d) = max(min(new_pos(d), 100), 0);
                end
            end
        else
            % Normal handling for other problems
            step = (ub - lb) .* randn(1, dim) * 0.1;
            new_pos = current_pos + step;
            new_pos = max(min(new_pos, ub), lb);
        end
        
        % Evaluate new solution
        switch problem
            case 'Ackley'
                new_val = Ackley(new_pos);
            case 'G06'
                new_val = G06(new_pos);
            case 'Knapsack'
                new_val = Knapsack(new_pos);
        end
        
        % Accept or reject new solution
        delta = new_val - current_val;
        if delta < 0 || rand < exp(-delta/T)
            current_pos = new_pos;
            current_val = new_val;
            
            % Update best solution
            if new_val < best_val
                best_pos = new_pos;
                best_val = new_val;
            end
        end
        
        % Record the best value for convergence history
        convergence(iter) = best_val;
        iter = iter + 1;
    end
    
    % Cool down
    T = alpha * T;
end

end
