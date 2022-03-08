%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     (C) Andrews T. Anum and Michael Pokojovy (2022)       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

FILE = readtable('EP.STDS.csv');
FILEDATE = table2array(FILE(:,1));
MONTHS = datenum(FILEDATE, "yyyy-mm");
new_cases = table2array(FILE(:,2));
%n = size(closedP,1);
difflogx = diff(log(new_cases));
n = size(difflogx,1);

rem = 24;
newn = n - rem;
delta_t  = 1;
total_time = rem;
nsim = 50000;
alpha = 0.5;
[rate, vola, mu, sigma]     = fit_GBM(new_cases(1:newn), delta_t, alpha);
[rateu, volau, muu, sigmau] = fit_GBM_usual(new_cases(1:newn), delta_t);
[ratem, volam, mum, sigmam] = fit_GBM_mcd(new_cases(1:newn), delta_t, alpha);
S0 = new_cases(newn);

[x_pred]  = pred_GBM(rate,  vola,  delta_t, total_time, S0, nsim);
[x_predu] = pred_GBM(rateu, volau, delta_t, total_time, S0, nsim);
[x_predm] = pred_GBM(ratem, volam, delta_t, total_time, S0, nsim);

mse  = zeros(1, size(x_pred,2));
umse = zeros(1, size(x_pred,2));
mmse = zeros(1, size(x_pred,2));
% % 
for l = 1:size(x_pred,2)
    mse(l)  = (1/nsim)*sum((x_pred(:,l) -  new_cases(newn+l) ).^2);
    umse(l) = (1/nsim)*sum((x_predu(:,l) - new_cases(newn+l) ).^2);
    mmse(l) = (1/nsim)*sum((x_predm(:,l) - new_cases(newn+l) ).^2);
end
 
figure()
plot(1:rem, mse,'k','LineWidth', 3.5)
hold on 
plot(1:rem, mmse,'b--', 'LineWidth', 3.5)
hold on 
plot(1:rem, umse,'r:', 'LineWidth', 3.5)
 
ylabel('MSE ','fontsize',28,'interpreter','latex');
xlabel('Future month index','fontsize',28,'interpreter','latex');
%title('MSE comparison among the estimators','fontsize', 45,'interpreter','latex');
legend('MDPD($\alpha = 0.5$)','MCD(bdp = 0.272)','Usual','MCD','fontsize',30,'interpreter','latex','Location','Best');
ax.FontSize = 25;
ax.LabelFontSizeMultiplier = 1.2; 
hold off
 
function[rate, vola, mu, sigma] = fit_GBM(x_data, dt, alpha)
    dlogx = diff(log(x_data));
    [~, mu, sigma] = ddiv_estimator3(mean(dlogx), std(dlogx), dlogx, alpha);
    %[~, mu, sigma] = ddiv_estimator3(median(dlogx), iqr(dlogx)/1.349, dlogx, alpha);
    vola = sigma/sqrt(dt);
    rate = mu/dt + 0.5*vola^2;
end

function[rate, vola, mu, sigma] = fit_GBM_mcd(x_data, dt, alpha)
    dlogx = diff(log(x_data));
    bdp = round(alpha/((1+alpha)^(3/2)),3);
    [raw] = mcd1D(dlogx, bdp);
    mu = raw.loc;
    sigma = raw.cov;
    vola = sigma/sqrt(dt);
    rate = mu/dt + 0.5*vola^2;
end

function[rate, vola, mu, sigma] = fit_GBM_usual(x_data, dt)
    dlogx = diff(log(x_data));
    mu = mean(dlogx);
    sigma = std(dlogx);
    vola = sigma/sqrt(dt);
    rate = mu/dt + 0.5*vola^2;
end

function[x_pred] = pred_GBM(rate, vola, delta_t, total_time, S0, nsim)
    m = rate - 0.5*(vola^2);
    x_pred = zeros(nsim, total_time);
    x_pred(:,1) = S0;
    for rep = 1:nsim
        if(total_time > 1)
            for i = 2:(total_time)
                x_pred(rep, i) = x_pred(rep,i-1)*exp(delta_t*m + sqrt(delta_t)*vola*randn(1));
            end
        end
    end
end
