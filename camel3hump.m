function [y] = camel3hump(x)
% [y] = camel3(x)
%
% 3 Hump Camel Function
%
% Lachlan Moore
% 2020 December

x1 = x(1);
x2 = x(2);

y = 2*x1^2 + -1.05*x1^4 + x1^6 / 6 + x1*x2 + x2^2;

end
