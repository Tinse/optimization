function [best_val, convergence] = PSO(problem, dim, pop_size, max_iter)
%PSO Particle Swarm Optimization algorithm
%   Inputs:
%       problem: string, name of the optimization problem
%       dim: dimension of the problem
%       pop_size: population size
%       max_iter: maximum number of iterations
%   Outputs:
%       best_val: best solution found
%       convergence: convergence curve

% Algorithm parameters
w = 0.9;             % Inertia weight
c1 = 1.4;            % Personal learning coefficient
c2 = 1.4;            % Global learning coefficient

% For Knapsack problem
global weights values capacity

% Initialize bounds based on problem
switch problem
    case 'Ackley'
        lb = -32.768 * ones(1, dim);
        ub = 32.768 * ones(1, dim);
        init_lb = lb;
        init_ub = ub;
        vmax = 0.2 * (ub - lb);  % Maximum velocity for Ackley
    case 'G06'
        lb = [13, 0];
        ub = [100, 100];
        init_lb = [13, 0];
        init_ub = [20, 5];
        dim = 2;
        vmax = 0.1 * (ub - lb);  % More restricted velocity for G06
    case 'Knapsack'
        lb = zeros(1, dim);
        ub = ones(1, dim);
        init_lb = lb;
        init_ub = ub;
        vmax = 0.15 * ones(1, dim);  % Restricted velocity for Knapsack
end

% Initialize population
pos = init_lb + (init_ub - init_lb) .* rand(pop_size, dim);
vel = -vmax + 2 * vmax .* rand(pop_size, dim);

% Initialize personal best
p_best = pos;
p_best_val = inf(pop_size, 1);

% Initialize global best
g_best = zeros(1, dim);
g_best_val = inf;

% Evaluate initial population
for i = 1:pop_size
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
    end
    
    % Update global best
    if current_val < g_best_val
        g_best_val = current_val;
        g_best = pos(i,:);
    end
end

% Initialize convergence history
convergence = zeros(1, max_iter);
convergence(1) = g_best_val;

% Main loop
for iter = 2:max_iter
    % Update each particle
    for i = 1:pop_size
        % Update velocity
        vel(i,:) = w * vel(i,:) + ...
                   c1 * rand(1, dim) .* (p_best(i,:) - pos(i,:)) + ...
                   c2 * rand(1, dim) .* (g_best - pos(i,:));
        
        % Apply velocity clamping
        % vel(i,:) = min(max(vel(i,:), -vmax), vmax);
        
        % Update position
        pos(i,:) = pos(i,:) + vel(i,:);
        
        % Bound position and handle special cases
        if strcmp(problem, 'G06')
            % Strict boundary handling for G06
            pos(i,1) = max(min(pos(i,1), 100), 13);  % 13 <= x1 <= 100
            pos(i,2) = max(min(pos(i,2), 100), 0);   % 0 <= x2 <= 100
            
            % Adjust velocity if position hits boundary
            if pos(i,1) == 13 || pos(i,1) == 100
                vel(i,1) = 0;
            end
            if pos(i,2) == 0 || pos(i,2) == 100
                vel(i,2) = 0;
            end
        elseif strcmp(problem, 'Knapsack')
            % Special handling for Knapsack problem
            % Ensure values stay in [0,1]
            pos(i,:) = max(min(pos(i,:), 1), 0);
            
            % Reset velocity if position hits boundary
            vel(i,:) = vel(i,:) .* (pos(i,:) > 0 & pos(i,:) < 1);
        else
            % Normal boundary handling for other problems
            pos(i,:) = max(min(pos(i,:), ub), lb);
        end
        
        % Evaluate new position
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
    
    % Record convergence
    convergence(iter) = g_best_val;
    
    % Optional: Adaptive parameter adjustment
    % w = max(0.4, w * 0.999);  % Gradually decrease inertia weight
end

best_val = g_best_val;
end
