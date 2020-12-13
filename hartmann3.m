function [y] = hartmann3(x)
% [y] = hartmann3(x)
%
% 3 Hump Camel Function
%
% Lachlan Moore
% 2020 December

alpha = [1.0; 1.2; 3.0; 3.2];

A = [3, 10, 30;
     0.1, 10, 35;
     3, 10, 30;
     0.1, 10, 35];
 
P = 1e-4 * [3689, 1170, 2673;
               4699, 4387, 7470;
               1091, 8732, 5547;
               381, 5743, 8828];

y = 0;
for i = 1:4
   si = 0;
   for j = 1:3
       si = si + A(i,j)*(x(j) - P(i,j))^2;  
   end
   si = -1 * si;
   y = y + alpha(i) * exp(si);
end
y = -1 * y;

end
