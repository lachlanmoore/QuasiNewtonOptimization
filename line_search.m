function [alphastar] = line_search(x, dir, alpha1)
% [alphastar] = line_search(x, dir, alpha1)
%
% A line search algorithm that satisfies the strong wolfe conditions
%
% -----------------
% Inputs:
% x: Initial starting location
% dir: Search direction
% alpha1: Initial step length
%
% Outputs: 
% alphastar: Step length meeting the strong wolfe conditions
% ------------------
%
% Taken from "Numerical Optimization, Nocedal and Wright, 2006"
% Algorithms 3.5, 3.6
% 
% Lachlan Moore
% 2020 December


    a0 = 0;
    a1 = alpha1;
    amax = 10*a1;
    
    c1 = 1e-4;
    c2 = .9;
    i = 1;
    
    [fun0, grad0]     = obj(x);
    [fold, ~] = obj(x+a0*dir);
    
    d = 1;
    while 1
        
        if d >= 1000
                error('Runtime too long')
        end
        
        [fun1, grad1]   = obj(x+a1*dir);
        
        if (fun1 > (fun0 + c1*a1*grad0'*dir)) || ((fun1 > fold) && (i > 1)) %Violates sufficient decrease
           alphastar = zoom(a0, a1); 
           break
        end

        if abs(grad1'*dir) <= -c2*grad0'*dir %Sufficient Decrease
            alphastar = a1;
            break
        end
        if grad1'*dir >= 0 %Curvature condition
            alphastar = zoom(a1, a0);
            break
        end
        
        a0 = a1;
        a1 = (a0 + amax )/2;
        i = i+1;
        fold = fun1;
        d = d+1;
        
    end
    
    
    function [alphastar] = zoom(alphalo, alphahi)
        % [alphastar] = zoom(alphalo, alphahi)
        % 
        % Inputs:
        % alphalo, alphahi: Bounds of search location
        %
        % Output: 
        % alphastar: Step length meeting stong wolfe conditions
        
        k = 1;
        while 1
            if k >= 5000
                error('Runtime too long')
            end
            
            alphaj = (alphalo+alphahi)/2;
            [funj, gradj] = obj(x+alphaj*dir);
            [funlo, ~] = obj(x+alphalo*dir);
            
            if (funj > fun0 + c1*alphaj*grad0'*dir) || (funj >= funlo)
               alphahi = alphaj; 
            else
                if abs(gradj'*dir) <= -c2 * grad0'*dir
                    alphastar = alphaj;
                    return
                end
                
                if gradj'*dir*(alphahi-alphalo) >= 0
                    alphahi = alphalo;
                end
                alphalo = alphaj; 
            end 
            k = k+1;
        end 
    end

end