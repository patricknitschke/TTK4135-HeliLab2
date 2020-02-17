function [c, ceq] = nonlcon_e(x)
%NONLCON Summary of this function goes here
%   Detailed explanation goes here
N = 40;
mx = 6;

alpha = 0.2;
beta = 20;
lambda_t = 2*pi/3;

c = alpha * exp(-beta*(x(1:mx:N*mx)-lambda_t).^2) - x(5:mx:N*mx);
ceq = [];
end

