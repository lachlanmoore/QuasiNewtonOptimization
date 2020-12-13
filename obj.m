function [f, g] = obj(x)
% [f, g] = obj(x)
%
% Objective function and gradient evaluation
% 
% -----------------
% Inputs:
% x: Evaluation location
% 
% Outputs:
% f: function eval
% g: function gradient, complex step approximation
% ------------------
%
% Lachlan Moore
% 2020 December




h = 1e-60; % complex step 
f = sub(x);
g = zeros(length(x), 1);

for i = 1:length(x)
    
   xc = x;
   xc(i) = complex(xc(i), h);  % complex step
   g(i) = imag(sub(xc)/h);     % complex step
   
end


    function [val] = sub(xc)
     % [val] = sub(xc)
     % Objective Function Definition
        
        
%% Line Search Test Cases
%       val = xc(1)^2;                                  % 1D Test
%       val = xc(1)*sin(xc(1)) + xc(1)*cos(2*xc(1));    % 1D Test
%       val = sin(xc(1)) + sin(10/3*xc(1));             % 1D Test

%% Final Project Test Cases
      val = -cos(xc(1))*cos(xc(2))*exp(-(xc(1)-pi)^2-(xc(2)-pi)^2); %Easom Function 2D
%       val = hartmann3(xc);                                          %Hartman 3D
%       val = camel3hump(xc);                                         %3 Hump Camel Function 2D

    end
end