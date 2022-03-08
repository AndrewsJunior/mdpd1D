%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     (C) Andrews T. Anum and Michael Pokojovy (2022)       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[raw] = mcd1D(vecX, bdp)

   [n, v]=size(vecX);
   h=floor(2*floor((n+v+1)/2)-n+2*(n-floor((n+v+1)/2))*(1-bdp));

   if(h==n)
      [~,rank_vecX] = sort(vecX);
      raw.loc = mean(vecX);
      raw.best = rank_vecX;
      raw.cov = var(vecX);
   else
      [~,rank_vecX] = sort(vecX);
      y = vecX(rank_vecX);
      s_sq  = zeros(1, (n-h+1));
      y_bar = zeros(1, (n-h+1));
      y_bar(1) = mean(y(1:h)) ;
      s_sq(1) = var(y(1:h))*(h-1);

      for j = 2:(n-h+1)
         y_bar(j) = (h*y_bar(j-1) - y(j-1) + y(j+h-1) )/h;
         s_sq(j)  = s_sq(j - 1) - y(j - 1)*y(j - 1) + y(j + h - 1)*y(j + h - 1) - h*(y_bar(j))*(y_bar(j)) + h*(y_bar(j - 1))*(y_bar(j - 1));            
      end
        
      [~, j_ast] = min(s_sq);
      raw.loc = mean(y_bar(j_ast));
      raw.best = rank_vecX(j_ast:j_ast + (h-1)) + 1;
      alpha = 1- bdp;
      smallsamplecorfactor =1;
        
      if smallsamplecorfactor==1
         c1factor = corfactorRAW(1,n,alpha);
      end
        
      c1factor = c1factor * consistencyfactor(h,n,1);
      raw.cov =c1factor*(s_sq(j_ast)/(h-1));
      raw.h = h;
   end
end

function rawconsfac = consistencyfactor(h,n,v)
   a = chi2inv(h/n,v);
   rawconsfac = (h/n)/(chi2cdf(a,v+2));
end
  
% corfactorRAW function
function rawcorfac = corfactorRAW(p,n,alpha)        
   if p == 1
      fp_500_n = 1-(exp(0.262024211897096)*1)/n^0.604756680630497;
      fp_875_n = 1-(exp(-0.351584646688712)*1)/n^1.01646567502486;
   end
        
   if 0.5 <= alpha && alpha <= 0.875
      fp_alpha_n=fp_500_n+(fp_875_n-fp_500_n)/0.375*(alpha-0.5);
   end
   if 0.875 < alpha && alpha < 1
      fp_alpha_n=fp_875_n+(1-fp_875_n)/0.125*(alpha-0.875);
   end
    rawcorfac=1/fp_alpha_n;
   if rawcorfac <=0 || rawcorfac>50
      rawcorfac=1;
   end
end
