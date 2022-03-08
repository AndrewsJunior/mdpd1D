%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     (C) Andrews T. Anum and Michael Pokojovy (2022)       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Hn, mu, sigma] = ddiv_estimator3(mu0, sigma0, X, alpha, plot_flag)
    if (~exist('mu0'))
       mu0 = median(X);
    end
    if (~exist('sigma0'))
       sigma0 = median(abs(diff(x)))/(2*erfinv(0.5));
    end
    if (nargin < 4) 
        error('This function requires at least four input arguments.');
    end
    if (nargin == 4)
        plot_flag = 0;
    end
    if (nargin == 5)
        plot_flag = 1;
    end
    
    %min_crit_val = Inf;
    Hnc = sum(mu0^2 + sigma0^2); 
    [Hn, mu, sigma] = NM_GD_descent([mu0; sigma0], alpha, X, plot_flag);
    
    function [Hn, mu, sigma] = NM_GD_descent(theta, alpha, X, plot_flag)
        num_flag = 0;
        step_prev      = 1;
        dtheta_norm_sq = Inf; 
    
        max_glob_it = 1E6;
        glob_it = 1;
        
        while (glob_it <= max_glob_it)
            [dir, step] = armijo_rule(dtheta_norm_sq, step_prev, theta, alpha, X);
            step_prev = step;
            %stopping criterion value
            dtheta_norm_sq = sum((project(theta + step*dir) - theta).^2);

            if (dtheta_norm_sq < Hnc*(1E-12) )
                [Hn, ~, ~] = DenPow_Divergence(theta, alpha, X, num_flag);
                break;
            else
                [Hn, ~, ~] = DenPow_Divergence(theta, alpha, X, num_flag);
                
                if (plot_flag == 1)
                   hold on 
                   plot3(theta(1), theta(2), Hn,'r*')
                end 
                %theta update
                theta = project(theta + step*dir);
            end
            glob_it = glob_it + 1;
        end
        mu    = theta(1);
        sigma = theta(2);
    end
end
  
function y = project(x)
    y = [x(1); max(1E-20, x(2))];
end

function[dir, step] = NM(theta, alpha, X, Hn_Hessian, Hn_grad, vartheta)
    global Hn_prev;
    num_flag = 0;
    dir  = -Hn_Hessian\Hn_grad;
    step = 1.0;

    Hn = DenPow_Divergence(project(theta + step*dir), alpha, X, num_flag);

    if((isequal(project(theta + step*dir),(theta + step*dir)  )) && (Hn - Hn_prev <= -vartheta*sum((Hn_grad).^2))) 
        return;
    else
        [dir, step] = GD(theta,alpha, X, step);
    end
end

function[dir, step] = GD(theta, alpha, X, step )
    num_flag = 0;
    global Hn_prev;
    global Hn_grad;
    dir = -Hn_grad;
    max_it = 50000;
    it = 1;
    
    while (it <= max_it)
        Hn = DenPow_Divergence(project(theta + step*dir), alpha, X, num_flag);
        gamma = 1E-4;
        r = -gamma/step*sum((project(theta + step*dir) - theta).^2);
        %choosing maximum step size for which Hn - Hn_prev <= r in {1, 1/2, 1/4, ..}
        if (Hn - Hn_prev > r )
            step = 0.5*step;
        else
            break;
        end
        it = it + 1;
    end
end 

function [dir, step] = armijo_rule(dtheta_norm_sq, step_prev, theta, alpha, X)
    num_flag = 0;
    gamma = 1E-4;
    vartheta = 1E-6;
    global Hn_prev;
    global Hn_grad;
    
    [Hn_prev, Hn_grad, Hn_Hessian] = DenPow_Divergence(theta, alpha, X, num_flag);
    % Check for positive definiteness of Hessian
    eig_Hessian = eig(Hn_Hessian);
    Hessian_PD_flag = ((min(eig_Hessian) >= vartheta) && (max(eig_Hessian) <= (1/vartheta)));
    step = step_prev;
    
    % Check if Newton's method can be employed (both positive definitness
    % and convergence check
    if(Hessian_PD_flag)
        [dir,step] = NM(theta, alpha, X, Hn_Hessian, Hn_grad, vartheta);
         return; 
    else
        [dir,step] = GD(theta,alpha, X, step );
    end
end 
