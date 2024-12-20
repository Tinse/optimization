function [best_val, convergence] = SecVibratPSO(problem, dim, pop_size, max_iter, weights, capacity, values)
%SecVibratPSO Second-Order Vibrating Particle Swarm Optimization
%   Inputs:
%       problem: string, name of the optimization problem
%       dim: problem dimension
%       pop_size: population size
%       max_iter: maximum number of iterations
%       weights: weights of items in the knapsack problem
%       capacity: capacity of the knapsack
%       values: values of items in the knapsack problem
%   Outputs:
%       best_val: best fitness value found
%       convergence: convergence curve

% For Knapsack problem
global weights values capacity

% Algorithm parameters
w_max = 0.9;         % Maximum inertia weight
w_min = 0.2;         % Minimum inertia weight (reduced for better convergence)
phi1 = 2.0;          % Personal learning coefficient
phi2 = 2.0;          % Global learning coefficient
xi1_max = 0.5;       % Maximum first vibration coefficient
xi1_min = 0.05;      % Minimum first vibration coefficient (reduced)
xi2_max = 0.5;       % Maximum second vibration coefficient
xi2_min = 0.05;      % Minimum second vibration coefficient (reduced)

% Initialize adaptive parameters
w = w_max;           % Start with maximum exploration
xi1 = xi1_max;       % Start with strong vibration
xi2 = xi2_max;       % Start with strong vibration

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
        vmax = 1.5 * ones(1, dim);  % Restricted velocity for Knapsack
end

% Initialize population
pos = init_lb + (init_ub - init_lb) .* rand(pop_size, dim);
vel = -0.1 * (init_ub - init_lb) + 0.2 * (init_ub - init_lb) .* rand(pop_size, dim);
pos_prev = pos;  % Previous position for second-order terms

% Initialize personal and global best
p_best = pos;
p_best_val = inf(pop_size, 1);
g_best = zeros(1, dim);
g_best_val = inf;

% Initialize convergence curve
convergence = zeros(1, max_iter);

% Evaluate initial population
for i = 1:pop_size
    switch problem
        case 'Ackley'
            val = Ackley(pos(i,:));
        case 'G06'
            val = G06(pos(i,:));
        case 'Knapsack'
            val = Knapsack(pos(i,:));
    end
    
    % Update personal best
    if val < p_best_val(i)
        p_best_val(i) = val;
        p_best(i,:) = pos(i,:);
        
        % Update global best
        if val < g_best_val
            g_best_val = val;
            g_best = pos(i,:);
        end
    end
end

% Initialize history for stagnation detection
best_val_history = inf(max_iter, 1);
best_val_history(1) = g_best_val;
last_improve = 1;
stagnation_count = 0;

% Variables for tracking improvement
% last_improve = 1;
% stagnation_count = 0;
% best_val_history = inf(1, max_iter);
% best_val_history(1) = g_best_val;

% Main loop
for iter = 1:max_iter
    % Calculate evolutionary state
    progress_ratio = iter / max_iter;
    
    % Adaptive parameter adjustment based on iteration progress and improvement
    if iter > 1
        if best_val_history(iter-1) < best_val_history(last_improve)
            last_improve = iter-1;
            stagnation_count = 0;
        else
            stagnation_count = stagnation_count + 1;
        end
    end
    
    % Early stage: Strong exploration but with controlled vibration
    if progress_ratio < 0.1
        w = w_max - (w_max - 0.6) * (progress_ratio/0.1);  % Gradual decrease
        xi1 = xi1_max * (1 - progress_ratio/0.1 * 0.3);    % Slightly reduce vibration
        xi2 = xi2_max * (1 - progress_ratio/0.1 * 0.3);
    % Middle stage: Balanced search
    elseif progress_ratio < 0.5
        w = 0.6 - (0.6 - 0.4) * ((progress_ratio - 0.1)/0.4);  % Smoother transition
        xi1 = xi1_max * 0.7 - (xi1_max * 0.7 - xi1_min * 2) * ((progress_ratio - 0.1)/0.4);
        xi2 = xi2_max * 0.7 - (xi2_max * 0.7 - xi2_min * 2) * ((progress_ratio - 0.1)/0.4);
    % Late stage: Refined local search
    else
        w = 0.4 - (0.4 - w_min) * ((progress_ratio - 0.5)/0.5)^2;  % Quadratic decrease
        xi1 = xi1_min * 2 - xi1_min * ((progress_ratio - 0.5)/0.5);
        xi2 = xi2_min * 2 - xi2_min * ((progress_ratio - 0.5)/0.5);
    end
    
    % Adjust parameters if stagnation is detected
    if stagnation_count > 10
        w = min(w_max, w * 1.5);  % Increase w but with limit
        xi1 = min(xi1_max, xi1 * 1.2);  % Moderate increase in vibration
        xi2 = min(xi2_max, xi2 * 1.2);
        stagnation_count = 0;
    end
    
    for i = 1:pop_size
        % Calculate velocity according to the second-order vibrating PSO formula
        vel(i,:) = w * vel(i,:) + ...
                   phi1 * rand(1, dim) .* (p_best(i,:) - (1 + xi1) * pos(i,:) + xi1 * pos_prev(i,:)) + ...
                   phi2 * rand(1, dim) .* (g_best - (1 + xi2) * pos(i,:) + xi2 * pos_prev(i,:));
        
        % Apply velocity clamping
        vel(i,:) = min(max(vel(i,:), -vmax), vmax);
        
        % Store current position before update
        pos_prev(i,:) = pos(i,:);
        
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
                % Also adjust previous position to avoid oscillation
                pos_prev(i,1) = pos(i,1);
            end
            if pos(i,2) == 0 || pos(i,2) == 100
                vel(i,2) = 0;
                % Also adjust previous position to avoid oscillation
                pos_prev(i,2) = pos(i,2);
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
                val = Ackley(pos(i,:));
            case 'G06'
                val = G06(pos(i,:));
            case 'Knapsack'
                val = Knapsack(pos(i,:));
        end
        
        % Update personal best
        if val < p_best_val(i)
            p_best(i,:) = pos(i,:);
            p_best_val(i) = val;
            
            % Update global best
            if val < g_best_val
                g_best = pos(i,:);
                g_best_val = val;
            end
        end
    end
    
    % Record the best value for convergence history
    convergence(iter) = g_best_val;
    best_val_history(iter) = g_best_val;
end

best_val = g_best_val;
end
