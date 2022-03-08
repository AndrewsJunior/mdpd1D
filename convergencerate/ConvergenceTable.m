%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     (C) Andrews T. Anum and Michael Pokojovy (2022)       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
nsamp = [30 50 100 200 300 400];
N = 1000000;
alpha = [0.1, 0.25, 0.5, 0.75, 1];
fh = @(x) x./((1+x).^(3/2));
alphamcd = round(feval(fh, alpha),3);
mean_dpd   = zeros(size(nsamp, 2)*size(alpha, 2), N);
sigma_dpd  = zeros(size(nsamp, 2)*size(alpha, 2), N);
mean_mcd   = zeros(size(nsamp, 2)*size(alpha, 2), N);
sigma_mcd  = zeros(size(nsamp, 2)*size(alpha, 2), N);
m = 0;
for k = 1:size(nsamp, 2)%controls sample size nsapm
    for j = 1:size(alpha, 2)%controls alpha
        testdata = randn(N, nsamp(k) );
        dpd_mean = zeros(1, N);
        dpd_sigma = zeros(1,N);
        mcd_mean = zeros(1, N);
        mcd_sigma = zeros(1,N);

        for l = 1:N %controls the replicatins
            %with non robust warmstart
            [~, mu, sigma] = ddiv_estimator3(mean(testdata(l,:)), std(testdata(l,:) ), testdata(l,:), alpha(j));
            dpd_mean(l) =  mu;
            dpd_sigma(l) = sigma;
            %[raw, ~] = mcd(testdata(l,:)', 'bdp', alphamcd(j), 'msg',0);
            [raw]    = mcd1D(testdata(l,:)', alphamcd(j));
            mcd_mean(l) = raw.loc;
            mcd_sigma(l) = raw.cov;
        end
        m = m+1;
        mean_dpd(m, :)  = dpd_mean;
        sigma_dpd(m, :) = dpd_sigma;
        mean_mcd(m, :)  = mcd_mean;
        sigma_mcd(m, :)  = mcd_sigma;
    end
end
% COMPUTING THE CONVERGENCE RATES FOR MDPD AND MCD
% MEAN

d = 5;
MSE_mean_dpd  = zeros(size(alpha,2), size(nsamp,2));
MSE_sigma_dpd = zeros(size(alpha,2), size(nsamp,2));

MSE_mean_mcd  = zeros(size(alpha,2), size(nsamp,2));
MSE_sigma_mcd = zeros(size(alpha,2), size(nsamp,2));

line_coeffs_dpd = zeros(size(alpha,2), 2);
line_coeffs_mcd = zeros(size(alpha,2), 2);

line_coeffs_dpd_sigma = zeros(size(alpha,2), 2);
line_coeffs_mcd_sigma = zeros(size(alpha,2), 2);

for j = 1:size(alpha,2)%controls alpha
    for i = 1: size(nsamp,2)%controls nsamp
        MSE_mean_dpd(j,i)  = (1/N)*sum( mean_dpd(j + (i-1)*d, :).^2);
        MSE_sigma_dpd(j,i) = (1/N)*sum((log(sigma_dpd(j + (i-1)*d, :))).^2 );
        MSE_mean_mcd(j,i)  = (1/N)*sum( mean_mcd(j + (i-1)*d, :)'.^2);
        MSE_sigma_mcd(j,i) = (1/N)*sum((log(sigma_mcd(j + (i-1)*d, :)))'.^2 );
    end

    line_coeffs_dpd(j,:) = polyfit(log(nsamp), log(MSE_mean_dpd(j,:)), 1);
    line_coeffs_mcd(j,:) = polyfit(log(nsamp), log(MSE_mean_mcd(j,:)), 1);    
    line_coeffs_dpd_sigma(j,:) = polyfit(log(nsamp), log(MSE_sigma_dpd(j,:)), 1);
    line_coeffs_mcd_sigma(j,:) = polyfit(log(nsamp), log(MSE_sigma_mcd(j,:)), 1);
end

A = [alpha; alphamcd; line_coeffs_dpd'; line_coeffs_mcd'];
B = [alpha; alphamcd; line_coeffs_dpd_sigma'; line_coeffs_mcd_sigma'];

fileID = fopen('conv_table_mean.txt','w');
fprintf(fileID,'%6s %18s %18s \n','alpha','MDPD', 'MCD');
fprintf(fileID,'%6.5f(%6.5f) %6.5f %6.5f  %6.5f %6.5f\n',A);
fclose(fileID);

fileID = fopen('conv_table_sigma.txt','w');
fprintf(fileID,'%6s %18s %18s \n','alpha','MDPD', 'MCD');
fprintf(fileID,'%6.5f(%6.5f) %6.5f %6.5f  %6.5f %6.5f\n', B);
fclose(fileID);
