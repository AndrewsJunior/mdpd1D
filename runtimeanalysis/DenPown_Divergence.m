%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     (C) Andrews T. Anum and Michael Pokojovy (2022)       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H, H_grad, H_Hessian] = DenPown_Divergence(mu_init, sigma_init, alpha, X, num_flag)

    if (nargin < 4) 
        error('This function requires at least four input arguments.');
    end
    if (nargin == 4)
        num_flag = 0;
    end
    if (nargin == 5)
        num_flag =1;
    end 
    if(size(X, 2) > size(X, 1))
        X = X';
    end
    
    n = size(X, 1);
    isig = inv(sigma_init);
    sigmasq = sigma_init * sigma_init;
    tp_aot = (2 * pi)^(alpha/2);
    expon = exp((-alpha * ((X-mu_init).^2)) / (2 * sigmasq));
    phi = (1 / sqrt(2 * pi))^(1 + alpha) * (isig).^(alpha) * ((2 * pi) / (1 + alpha))^0.5;
    
   %%% ANALYTIC SOLUTION     
    if (num_flag == 0)
        if (nargout == 1)
            H = phi -(1 + 1/alpha)/n *sum( ((isig)^alpha)/tp_aot * expon) ;
        end

        if (nargout == 2)
            H = phi -(1 + 1/alpha)/n *sum( ((isig)^alpha)/tp_aot * expon) ;            
            dHn_dmu = -(1 + 1/alpha)/n * sum((alpha * ((isig).^(alpha + 2))/ tp_aot) * (X-mu_init) .* expon);
            phi_Sigma = - alpha * (isig).^(alpha + 1) * (1 / sqrt(2 * pi))^(1 + alpha) * ((2 * pi) / (1 + alpha))^0.5 ;
            dHn_dsigma = phi_Sigma - (1 + 1/alpha)/n * sum((- alpha / (tp_aot) * ((isig).^(alpha + 1)) .* expon ...
           +  alpha / (tp_aot) * ((isig).^(alpha + 3)) .* (X-mu_init).^2 .* expon )) ;
            H_grad = [dHn_dmu ; dHn_dsigma];
        end
        
        if (nargout == 3)
            H = phi -(1 + 1/alpha)/n *sum( ((isig)^alpha)/tp_aot * expon) ;
            dHn_dmu = -(1 + 1/alpha)/n * sum((alpha * ((isig).^(alpha + 2))/ tp_aot) * (X-mu_init) .* expon);
            phi_Sigma = - alpha * (isig).^(alpha + 1) * (1 / sqrt(2 * pi))^(1 + alpha) * ((2 * pi) / (1 + alpha))^0.5 ;
            dHn_dsigma = phi_Sigma - (1 + 1/alpha)/n * sum((- alpha / (tp_aot) * ((isig).^(alpha + 1)) .* expon ...
           +  alpha / (tp_aot) * ((isig).^(alpha + 3)) .* (X-mu_init).^2 .* expon )) ;
            H_grad = [dHn_dmu ; dHn_dsigma];
            
            d2Hn_dmu2 = -(1 + 1/alpha)/n * sum(((alpha * (isig).^(alpha + 2) / tp_aot) .* (- expon + (alpha ...
                 .* ((X - mu_init).^2)/ sigmasq) .* expon ) ));
            d2Hn_dmudsigma = -(1 + 1/alpha)/n * sum((X-mu_init) .* ( (( alpha^2 * (isig).^(alpha + 5)) / tp_aot ) ...
            .* (X-mu_init).^2 .* expon - (alpha/tp_aot  * (alpha + 2) * (isig).^(alpha + 3)) .* expon  )) ;
             
            phi_Sigmasq = alpha * (alpha  + 1) * ((isig).^(alpha + 2)) ... 
            * (1 / sqrt(2 * pi))^(1 + alpha) * ((2 * pi) / (1 + alpha))^0.5;
        
            d2Hn_dsigma2 = phi_Sigmasq  - (1 + 1/alpha)/n * sum((  alpha * (alpha + 1)/ tp_aot  ...
            * (isig).^(alpha + 2) .* expon  - alpha * alpha / tp_aot  * (isig).^(alpha + 4) ...
            .* (X-mu_init).^2 .* expon + alpha * (-alpha - 3)/ tp_aot  * (isig).^(alpha + 4) ...
            .* (X - mu_init).^2 .* expon +  alpha * alpha / tp_aot  * (isig).^(alpha + 6) .* (X - mu_init).^4 ...
            .* expon ));
             
            H_Hessian = [d2Hn_dmu2 d2Hn_dmudsigma ; d2Hn_dmudsigma d2Hn_dsigma2];
        end
    end
    
    %%% NUMERICAL  SOLUTION    
    if (num_flag == 1)
        h = 1E-6;
        H = DenPown_Divergence(mu_init, sigma_init, alpha, X);
        
        if (nargout == 1)
            [H] = DenPown_Divergence(mu_init, sigma_init, alpha, X);
        end
        if (nargout == 2)
            H = DenPown_Divergence(mu_init, sigma_init, alpha, X);
            
            [H_mu, ~] = DenPown_Divergence(mu_init + h,   sigma_init,      alpha,   X);
            [H_sigma, ~] = DenPown_Divergence(mu_init,    sigma_init + h,  alpha,   X);
            
            dH_dmu_num = (H_mu - H)/h;
            dH_dsigma_num = (H_sigma - H)/h;
            
            H_grad = [dH_dmu_num; dH_dsigma_num];
        end
        
        if (nargout == 3)
            H = DenPown_Divergence(mu_init, sigma_init, alpha, X);
            
            [H_mu, ~] = DenPown_Divergence(mu_init + h,   sigma_init,      alpha,   X);
            [H_sigma, ~] = DenPown_Divergence(mu_init,    sigma_init + h,  alpha,   X);
            
            dH_dmu_num = (H_mu - H)/h;
            dH_dsigma_num = (H_sigma - H)/h;
            
            H_grad = [dH_dmu_num; dH_dsigma_num];
            
            [H_mu_m, ~]    = DenPown_Divergence(mu_init - h, sigma_init,     alpha,   X);
            [H_musigma, ~] = DenPown_Divergence(mu_init + h, sigma_init + h, alpha,   X);
            [H_sigma_m, ~] = DenPown_Divergence(mu_init,     sigma_init - h, alpha,   X);

            d2H_dmu2_num = (H_mu_m -2*H + H_mu)/(h^2);
            d2H_dsigma2_num = (H_sigma_m -2*H + H_sigma)/(h^2);
            d2H_dmudsigma_num= (H_musigma - H_mu - H_sigma + H)/(h^2);
            
            H_Hessian = [d2H_dmu2_num d2H_dmudsigma_num; d2H_dmudsigma_num d2H_dsigma2_num];
        end 
    end
end
