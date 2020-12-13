function bfgs(xk, n, fun_input)
% bfgs(xk, n, fun_input)
%
% A BFGS quasi-Newton optimization algorithm for unconstrained minimization
%
% -----------------
% Inputs:
% xk: Initial guess
% n: Number of design variables
% fun_input: Objective function 
% ------------------
%
% Taken from "Numerical Optimization, Nocedal and Wright, 2006"
% Algorithm 6.1
% 
% Lachlan Moore
% 2020 December

close all

H  = eye(n);

func      = @(x) fun_input(x);

[func_eval, func_deriv] = func(xk);

i = 1;
iter = 300;
eps = 1e-6; 


while norm(func_deriv) > eps
 
    pk     = -H * func_deriv; % search direction
    pk     = pk/norm(pk);
    
    alpha  = 1; %compute from search direction, will be line search
    alpha  = line_search(xk, pk, alpha);
    xk_new = xk + alpha * pk;
    
    [func_eval_new, func_deriv_new] = func(xk_new);
    
    sk     = xk_new - xk;
    yk     = func_deriv_new - func_deriv;
    
    rk     = 1/(yk'*sk);
    H      = (eye(n) - rk*(sk*yk'))*H*(eye(n) - rk*(yk*sk')) + rk*(sk*sk');
    
    xk     = xk_new;
    func_eval  = func_eval_new;
    func_deriv = func_deriv_new;
    
    
    i      = i + 1;
    if i >= iter
        break
    end

end

fprintf('Final Evaluation %d\n', func_eval)
fprintf('Iterations: %d\n', i-1)
fprintf('Gradient %d\n',norm(func_deriv))
fprintf('Local Minimum at \n')
fprintf('%c \n', xk)
disp(H)

end