function f = Knapsack(x)
%KNAPSACK 0-1 Knapsack problem
%   Maximize the value of items in the knapsack while respecting weight constraint
%   For continuous optimization algorithms, we use sigmoid function to convert
%   continuous values to binary decisions

% Problem parameters (example instance)
global weights values capacity  % Make these global so they're consistent across calls

% Initialize if not already defined
if isempty(weights) || isempty(values) || isempty(capacity)
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
end

% Ensure x has the correct dimension
if length(x) ~= length(weights)
    x = x(1:length(weights));  % Truncate if too long
    % Pad with zeros if too short
    if length(x) < length(weights)
        x = [x, zeros(1, length(weights) - length(x))];
    end
end

% Convert continuous values to binary using sigmoid function
x_binary = 1./(1 + exp(-10*(x - 0.5)));

% Calculate total weight and value using element-wise multiplication
total_weight = sum(weights .* x_binary);
total_value = sum(values .* x_binary);

% Penalty for exceeding capacity
penalty = 1e1;
if total_weight > capacity
    f = -total_value + penalty * (total_weight - capacity);
else
    f = -total_value;  % Negative because we're minimizing
end

end
