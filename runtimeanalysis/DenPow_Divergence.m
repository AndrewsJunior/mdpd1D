%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     (C) Andrews T. Anum and Michael Pokojovy (2022)       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H, H_grad, H_Hessian] = DenPow_Divergence(theta, alpha, X, num_flag)
    mu    = theta(1);
    sigma = theta(2);

    if (nargin < 4) 
        error('This function requires at least four input arguments.');
    end
    if (nargin == 4)
        num_flag = 0;
    end
    if (nargin == 5)
        num_flag =1;
    end 
    if (size(X, 1) > size(X, 2))
        X = X';
    end

    n = size(X, 2);    
    isig = inv(sigma);
    sigmasq = sigma * sigma;
    ttp_php = (2 * pi)^(alpha/2);
    expon = exp((-alpha * ((X-mu).^2)) / (2 * sigmasq));
    phi = (1 / sqrt(2 * pi))^(1 + alpha) * (isig).^(alpha) * ((2 * pi) / (1 + alpha))^0.5;

    if (num_flag == 0)
        
%         disp('Analytic Solutions:');
        
        if (nargout == 1)
            H = phi -(1 + 1/alpha)/n *sum( ((isig)^alpha)/ttp_php * expon) ;
        end
        
        if (nargout == 2)
            H = phi -(1 + 1/alpha)/n *sum( ((isig)^alpha)/ttp_php * expon) ;
            dHn_dmu = -(1 + 1/alpha)/n * sum((alpha * ((isig).^(alpha + 2))/ ttp_php) * (X-mu) .* expon);
            phi_Sigma = - alpha * (isig).^(alpha + 1) * (1 / sqrt(2 * pi))^(1 + alpha) * ((2 * pi) / (1 + alpha))^0.5 ;
            dHn_dsigma = phi_Sigma - (1 + 1/alpha)/n * sum((- alpha / (ttp_php) * ((isig).^(alpha + 1)) .* expon ...
           +  alpha / (ttp_php) * ((isig).^(alpha + 3)) .* (X-mu).^2 .* expon )) ;
            H_grad = [dHn_dmu ; dHn_dsigma];
        end
        
        if (nargout == 3)
            H = phi -(1 + 1/alpha)/n *sum( ((isig)^alpha)/ttp_php * expon) ;
            dHn_dmu = -(1 + 1/alpha)/n * sum((alpha * ((isig).^(alpha + 2))/ ttp_php) * (X-mu) .* expon);
            phi_Sigma = - alpha * (isig).^(alpha + 1) * (1 / sqrt(2 * pi))^(1 + alpha) * ((2 * pi) / (1 + alpha))^0.5 ;
            dHn_dsigma = phi_Sigma - (1 + 1/alpha)/n * sum((- alpha / (ttp_php) * ((isig).^(alpha + 1)) .* expon ...
           +  alpha / (ttp_php) * ((isig).^(alpha + 3)) .* (X-mu).^2 .* expon )) ;
            H_grad = [dHn_dmu ; dHn_dsigma];
            
            d2Hn_dmu2 = -(1 + 1/alpha)/n * sum(((alpha * (isig).^(alpha + 2) / ttp_php) .* (- expon + (alpha ...
                 .* ((X - mu).^2)/ sigmasq) .* expon ) ));
            d2Hn_dmudsigma = -(1 + 1/alpha)/n * sum((X-mu) .* ( (( alpha^2 * (isig).^(alpha + 5)) / ttp_php ) ...
            .* (X-mu).^2 .* expon - (alpha/ttp_php  * (alpha + 2) * (isig).^(alpha + 3)) .* expon  )) ;
             
            phi_Sigmasq = alpha * (alpha  + 1) * ((isig).^(alpha + 2)) ... 
            * (1 / sqrt(2 * pi))^(1 + alpha) * ((2 * pi) / (1 + alpha))^0.5;
        
            d2Hn_dsigma2 = phi_Sigmasq  - (1 + 1/alpha)/n * sum((  alpha * (alpha + 1)/ ttp_php  ...
            * (isig).^(alpha + 2) .* expon  - alpha * alpha / ttp_php  * (isig).^(alpha + 4) ...
            .* (X-mu).^2 .* expon + alpha * (-alpha - 3)/ ttp_php  * (isig).^(alpha + 4) ...
            .* (X - mu).^2 .* expon +  alpha * alpha / ttp_php  * (isig).^(alpha + 6) .* (X - mu).^4 ...
            .* expon ));
             
            H_Hessian = [d2Hn_dmu2 d2Hn_dmudsigma ; d2Hn_dmudsigma d2Hn_dsigma2];
        end
    end
    
    if (num_flag == 1)
%         disp('Numeric Solution');
        h = 1E-6;
        H = DenPow_Divergence(mu, sigma, alpha, X);
        
        if (nargout == 1)
            [H] = DenPow_Divergence(mu, sigma, alpha, X);
        end
        
        if (nargout == 2)
            H = DenPow_Divergence(mu, sigma, alpha, X);
            
            [H_mu, ~] = DenPow_Divergence(mu + h,   sigma,      alpha,   X);
            [H_sigma, ~] = DenPow_Divergence(mu,    sigma + h,  alpha,   X);
            
            dH_dmu_num = (H_mu - H)/h;
            dH_dsigma_num = (H_sigma - H)/h;
            
            H_grad = [dH_dmu_num; dH_dsigma_num];
        end
            
        if (nargout == 3)
            H = DenPow_Divergence(mu, sigma, alpha, X);
            
            [H_mu, ~] = DenPow_Divergence(mu + h,   sigma,      alpha,   X);
            [H_sigma, ~] = DenPow_Divergence(mu,    sigma + h,  alpha,   X);
            
            dH_dmu_num = (H_mu - H)/h;
            dH_dsigma_num = (H_sigma - H)/h;
            
            H_grad = [dH_dmu_num; dH_dsigma_num];
            
            [H_mu_m, ~]    = DenPow_Divergence(mu - h, sigma,     alpha,   X);
            [H_musigma, ~] = DenPow_Divergence(mu + h, sigma + h, alpha,   X);
            [H_sigma_m, ~] = DenPow_Divergence(mu,     sigma - h, alpha,   X);

            d2H_dmu2_num = (H_mu_m -2*H + H_mu)/(h^2);
            d2H_dsigma2_num = (H_sigma_m -2*H + H_sigma)/(h^2);
            d2H_dmudsigma_num= (H_musigma - H_mu - H_sigma + H)/(h^2);
            
            H_Hessian = [d2H_dmu2_num d2H_dmudsigma_num; d2H_dmudsigma_num d2H_dsigma2_num];
        end  
    end
end
