function f = G06(x)
%G06 G06 constrained optimization problem
%   Minimize: f(x) = (x(1) - 10)^3 + (x(2) - 20)^3
%   Subject to:
%   g1(x) = -(x(1) - 5)^2 - (x(2) - 5)^2 + 100 <= 0
%   g2(x) = (x(1) - 6)^2 + (x(2) - 5)^2 - 82.81 <= 0
%   13 <= x(1) <= 100
%   0 <= x(2) <= 100
%   Global minimum: f(14.095, 0.84296) = -6961.81388

% Objective function
f_obj = (x(1) - 10)^3 + (x(2) - 20)^3;

% Constraints evaluation
g1 = -(x(1) - 5)^2 - (x(2) - 5)^2 + 100;
g2 = (x(1) - 6)^2 + (x(2) - 5)^2 - 82.81;

% Calculate constraint violations
v_bound1 = max(0, 13 - x(1)) + max(0, x(1) - 100);
v_bound2 = max(0, -x(2)) + max(0, x(2) - 100);
v_g1 = max(0, g1);
v_g2 = max(0, g2);

% Total violation
total_violation = 100*v_bound1 + 100*v_bound2 + v_g1 + v_g2;

% Check feasibility
if total_violation == 0
    f = f_obj;
else
    % Use static penalty with reasonable coefficient
    penalty = 1e4;
    f = f_obj + penalty * total_violation;
end

end
