%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     (C) Andrews T. Anum and Michael Pokojovy (2022)       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Hn, mu, sigma] = newddiv_estimator(mu_init, sigma_init, X, alpha)

    if (nargin < 4) 
        error('This function requires four input arguments.');
    end
    min_crit_val = Inf;
    hnc = sum(mu_init^2 + sigma_init^2);
    [Hn, theta_hat.loc, theta_hat.cov] = gradient_decent(mu_init, sigma_init, alpha, X);
    mu = theta_hat.loc;
    sigma = theta_hat.cov;
    
    function [Hn, mu_init, sigma_init] = gradient_decent(mu_init, sigma_init, alpha, X)
        num_flag = 0;
        step_prev = 1;
        max_glob_it = 1E6;
        glob_it = 1;
        
        while (glob_it <= max_glob_it)
            [dir_mu, dir_sigma, step] = armijo_rule(step_prev, mu_init, sigma_init, alpha, X);
            step_prev = step;
            delta = (step*dir_mu)^2 + (max(1E-100, sigma_init + step*dir_sigma) - sigma_init)^2;

            if (delta < hnc*(1E-12) )
                [Hn,~,~] = DenPown_Divergence(mu_init, sigma_init, alpha, X, num_flag);                
                break;
            else
                [Hn,~,~] = DenPown_Divergence(mu_init, sigma_init, alpha, X, num_flag);
                mu_init    = mu_init    + step*dir_mu;
                sigma_init = max(1E-100, sigma_init + step*dir_sigma);
            end

            glob_it = glob_it + 1;
        end        
    end
end
      
function [dir_mu, dir_sigma, step] = armijo_rule(step_prev, mu_init, sigma_init, alpha, X)
    num_flag = 0;
    gamma = 1E-4;
    [Hn_prev, Hn_grad_prev] = DenPown_Divergence(mu_init, sigma_init, alpha, X, num_flag);
    dir_mu    = -Hn_grad_prev(1);
    dir_sigma = -Hn_grad_prev(2);
    step = step_prev;
    
    if(dir_mu^2 + dir_sigma^2 < 1E-18)
        dir_mu = 0;
        dir_sigma = 0;
        return;
    end
    max_it = 500;
    it = 1;

    while (it <= max_it)
        Hn = DenPown_Divergence(mu_init + step*dir_mu, max(1E-100, sigma_init + step*dir_sigma), alpha, X, num_flag);        
        r = -gamma/step*((step*dir_mu)^2 + (max(1E-100, sigma_init + step*dir_sigma) - sigma_init)^2);
        %choosing maximum step size for which Hn - Hn_prev <= r in {1, 1/2, 1/4, ..}
            if (Hn - Hn_prev > r )
                step = 0.5*step;
            else
                break;
            end

        it = it + 1;
    end
end
