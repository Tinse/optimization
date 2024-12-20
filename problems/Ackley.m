function f = Ackley(x)
%ACKLEY Ackley function
%   f(x) = -20*exp(-0.2*sqrt(1/d * sum(x.^2))) - exp(1/d * sum(cos(2*pi*x))) + 20 + e
%   Global minimum: f(0,...,0) = 0

d = length(x);
sum1 = sum(x.^2);
sum2 = sum(cos(2*pi*x));

term1 = -20 * exp(-0.2 * sqrt(sum1/d));
term2 = -exp(sum2/d);

f = term1 + term2 + 20 + exp(1);
end
