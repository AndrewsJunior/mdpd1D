%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     (C) Andrews T. Anum and Michael Pokojovy (2022)       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
k = 10;
alpha = 0.5;
%nsamp = [30, 50, 100, 200, 300, 400];
nsamp = 200;
N = 5000; 

for i= 1: size(nsamp, 2)    
    testdata = randn(N, nsamp(i));
    dpd_mean = zeros(1, N);
    dpd_sigma = zeros(1,N);
    mcd_mean = zeros(1, N);
    mcd_sigma = zeros(1,N);
    
    for l = 1:N 
        cdata = testdata(l,:);
        cdata(1:61) = 1000; %increase until breakdown occurs
        [~, mu, sigma] = ddiv_estimator2(mean(cdata), std(cdata), cdata, alpha);
        dpd_mean(l) =  mu;
        dpd_sigma(l) = sigma;
        [raw] = mcd1D(cdata', 0.272);
        mcd_mean(l) = raw.loc;
        mcd_sigma(l) = raw.cov;
    end
    pos_dpd_mean = size(find(abs(dpd_mean) <= k),2);
    pos_mcd_mean = size(find(abs(mcd_mean) <= k),2);
    pos_dpd_sigma = size(intersect(find(1/k <= abs(dpd_sigma)), find(abs(dpd_sigma)<= k) ),2);
    pos_mcd_sigma = size(intersect(find(1/k <= abs(mcd_sigma)), find(abs(mcd_sigma)<= k) ),2);
end
