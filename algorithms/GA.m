function [best_val, convergence] = GA(problem, dim, pop_size, max_iter)
%GA Genetic Algorithm with real-coded operators
%   Inputs:
%       problem: string, name of the optimization problem
%       dim: problem dimension
%       pop_size: population size
%       max_iter: maximum number of iterations
%   Outputs:
%       best_val: best solution found
%       convergence: convergence curve

% GA parameters
p_crossover = 0.8;    % Crossover probability
p_mutation = 0.1;     % Mutation probability
elite_size = 2;       % Number of elite solutions
eta_c = 20;          % Distribution index for SBX crossover
eta_m = 20;          % Distribution index for polynomial mutation

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
        % Adjust mutation parameters for G06
        eta_m = 10;  % Smaller value for more exploration
        p_mutation = 0.2;  % Higher mutation rate
    case 'Knapsack'
        lb = zeros(1, dim);
        ub = ones(1, dim);
        init_lb = lb;
        init_ub = ub;
end

% Initialize convergence curve
convergence = zeros(1, max_iter);

% Initialize population with problem-specific ranges
pop = init_lb + (init_ub - init_lb) .* rand(pop_size, dim);

% Add some diversity to the initial population
if strcmp(problem, 'G06')
    % Add some individuals exploring wider range
    num_diverse = round(pop_size * 0.2);  % 20% of population
    pop(end-num_diverse+1:end, :) = lb + (ub - lb) .* rand(num_diverse, dim);
end

% Initialize fitness
fitness = zeros(pop_size, 1);
best_val = inf;

% Main loop
for iter = 1:max_iter
    % Evaluate fitness
    for i = 1:pop_size
        if strcmp(problem, 'Knapsack')
            eval_pop = round(pop(i,:));  % Round for discrete problems
        else
            eval_pop = pop(i,:);
        end
        
        switch problem
            case 'Ackley'
                fitness(i) = Ackley(eval_pop);
            case 'G06'
                fitness(i) = G06(eval_pop);
            case 'Knapsack'
                fitness(i) = Knapsack(eval_pop);
        end
    end
    
    % Update best solution
    [min_fit, min_idx] = min(fitness);
    if min_fit < best_val
        best_val = min_fit;
        best_solution = pop(min_idx,:);
    end
    
    % Record the best value for convergence history
    convergence(iter) = best_val;
    
    % Selection (Tournament selection)
    new_pop = zeros(pop_size, dim);
    [~, sorted_idx] = sort(fitness);
    
    % Elitism
    new_pop(1:elite_size,:) = pop(sorted_idx(1:elite_size),:);
    
    % Generate rest of new population
    for i = (elite_size+1):2:pop_size
        % Select parents using tournament selection
        p1_idx = tournament_select(fitness);
        p2_idx = tournament_select(fitness);
        
        % Simulated Binary Crossover (SBX)
        if rand < p_crossover
            [child1, child2] = sbx_crossover(pop(p1_idx,:), pop(p2_idx,:), lb, ub, eta_c);
        else
            child1 = pop(p1_idx,:);
            child2 = pop(p2_idx,:);
        end
        
        % Mutation
        for j = 1:2
            if rand < p_mutation
                if strcmp(problem, 'G06')
                    % Special mutation for G06
                    for d = 1:dim
                        if rand < p_mutation
                            if d == 1
                                % For x1: 13 <= x1 <= 100
                                range = min(100 - child1(d), child1(d) - 13) * 0.1;
                                delta = range * randn;
                                child1(d) = child1(d) + delta;
                                child1(d) = max(min(child1(d), 100), 13);
                            else
                                % For x2: 0 <= x2 <= 100
                                range = min(100 - child1(d), child1(d)) * 0.1;
                                delta = range * randn;
                                child1(d) = child1(d) + delta;
                                child1(d) = max(min(child1(d), 100), 0);
                            end
                        end
                    end
                else
                    % Normal mutation for other problems
                    r = rand(1, dim);
                    pos = r < p_mutation;
                    if any(pos)
                        child1(pos) = child1(pos) + 0.1 * (ub(pos) - lb(pos)) .* randn(1,sum(pos));
                        child1 = max(min(child1, ub), lb);
                    end
                end
            end
            if i+1 <= pop_size
                if j == 1
                    new_pop(i,:) = child1;
                else
                    new_pop(i+1,:) = child2;
                end
            else
                new_pop(i,:) = child1;
            end
        end
    end
    
    pop = new_pop;
    
    % Adaptive parameter adjustment
    if mod(iter, 100) == 0
        eta_c = min(eta_c * 1.05, 30);  % Gradually increase distribution index
        p_mutation = max(p_mutation * 0.95, 0.01);  % Gradually decrease mutation rate
    end
end

end

function idx = tournament_select(fitness)
    k = 3;  % Tournament size
    pop_size = length(fitness);
    tournament_idx = randi(pop_size, 1, k);
    [~, min_idx] = min(fitness(tournament_idx));
    idx = tournament_idx(min_idx);
end

function [child1, child2] = sbx_crossover(parent1, parent2, lb, ub, eta_c)
    % Simulated Binary Crossover (SBX)
    dim = length(parent1);
    child1 = zeros(1, dim);
    child2 = zeros(1, dim);
    
    for i = 1:dim
        if rand < 0.5
            if abs(parent1(i) - parent2(i)) > 1e-14
                if parent1(i) < parent2(i)
                    y1 = parent1(i);
                    y2 = parent2(i);
                else
                    y1 = parent2(i);
                    y2 = parent1(i);
                end
                
                % Calculate beta_q
                rand_num = rand();
                beta = 1.0 + (2.0 * (y1 - lb(i)) / (y2 - y1));
                alpha = 2.0 - beta^(-(eta_c + 1.0));
                if rand_num <= 1.0/alpha
                    beta_q = (rand_num * alpha)^(1.0/(eta_c + 1.0));
                else
                    beta_q = (1.0/(2.0 - rand_num * alpha))^(1.0/(eta_c + 1.0));
                end
                
                % Generate children
                c1 = 0.5 * ((y1 + y2) - beta_q * (y2 - y1));
                c2 = 0.5 * ((y1 + y2) + beta_q * (y2 - y1));
                
                % Bound check
                c1 = max(lb(i), min(ub(i), c1));
                c2 = max(lb(i), min(ub(i), c2));
                
                if rand < 0.5
                    child1(i) = c2;
                    child2(i) = c1;
                else
                    child1(i) = c1;
                    child2(i) = c2;
                end
            else
                child1(i) = parent1(i);
                child2(i) = parent2(i);
            end
        else
            child1(i) = parent1(i);
            child2(i) = parent2(i);
        end
    end
end

function child = polynomial_mutation(parent, lb, ub, p_mutation, eta_m)
    % Polynomial Mutation
    dim = length(parent);
    child = parent;
    
    for i = 1:dim
        if rand < p_mutation
            y = parent(i);
            delta1 = (y - lb(i)) / (ub(i) - lb(i));
            delta2 = (ub(i) - y) / (ub(i) - lb(i));
            
            rand_num = rand();
            mut_pow = 1.0/(eta_m + 1.0);
            
            if rand_num <= 0.5
                xy = 1.0 - delta1;
                val = 2.0 * rand_num + (1.0 - 2.0 * rand_num) * (xy^(eta_m + 1.0));
                deltaq = val^mut_pow - 1.0;
            else
                xy = 1.0 - delta2;
                val = 2.0 * (1.0 - rand_num) + 2.0 * (rand_num - 0.5) * (xy^(eta_m + 1.0));
                deltaq = 1.0 - val^mut_pow;
            end
            
            y = y + deltaq * (ub(i) - lb(i));
            child(i) = max(lb(i), min(ub(i), y));
        end
    end
end
