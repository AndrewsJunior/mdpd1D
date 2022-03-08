%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     (C) Andrews T. Anum and Michael Pokojovy (2022)       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

FILE = readtable('EP.STDS.csv');
FILEDATE = table2array(FILE(:,1));
MONTHS = datenum(FILEDATE, "yyyy-mm");
new_cases = table2array(FILE(:,2));
difflogx = diff(log(new_cases));

figure(1)
hold on
ax = gca;
plot(MONTHS, new_cases, 'LineWidth', 2)
datetick('x', 'mmm-yy', 'keeplimits')
ylabel('New chlamydia cases','fontsize',34,'interpreter','latex');
xlabel('Date','fontsize',34,'interpreter','latex');%change to months
title('Number of new chlamydia cases','fontsize',48,'interpreter','latex')
ax.FontSize = 25;
ax.LabelFontSizeMultiplier = 1.5;
hold on
plot(MONTHS([100 111 146 164]),new_cases([100 111 146 164]), 'O', 'MarkerSize', 12)
hold off

n = size(difflogx, 1);
rem = 24;
newn = n - rem;
delta_t = 1;
total_time = rem;
nsim = 50000;
alpha = 0.5;
[rate, vola, mu, sigma] = fit_GBM(new_cases(1:newn), delta_t, alpha);
S0 = new_cases(newn);
[x_pred] = pred_GBM(rate, vola, delta_t, total_time, S0, nsim);
lb = zeros(1, size(x_pred,2));
ub = zeros(1, size(x_pred,2));

for l = 1:size(x_pred,2)
    lb(l) = quantile(x_pred(:,l), 0.05);
    
    ub(l) = quantile(x_pred(:,l), 0.95);
end
 
color_order = get(gca, 'colororder');
figure(2)
%color_order = get(gca, 'colororder');
ax = gca;
hold on 
%ylim([-50, 500+max(lb)]);
plot(MONTHS(1:newn+1),new_cases(1:newn+1),'LineWidth',2)
hold on
plot(MONTHS(newn:(newn+size(x_pred,2)-1) ),lb,'Color', color_order(7, :),'LineWidth',3)
hold on
plot(MONTHS(newn:(newn+size(x_pred,2)-1)), mean(x_pred),'b', 'LineWidth',3)
hold on
plot(MONTHS(newn:(newn+size(x_pred,2)-1) ),ub,'Color', color_order(7, :),'LineWidth',3)
hold on
plot(MONTHS(newn:(newn+size(x_pred,2)-1)),new_cases(newn:newn+total_time-1),	'k','LineWidth',3)
datetick('x', 'mmm-yy') 
ylabel('New chlamydia cases','fontsize',28,'interpreter','latex');
xlabel('Date','fontsize',28,'interpreter','latex');
ax.FontSize = 25;
ax.LabelFontSizeMultiplier = 1.2;
shade(MONTHS(newn:(newn+size(x_pred,2)-1)), lb, MONTHS(newn:(newn+size(x_pred,2)-1) ), ub, 'FillType', [1 2;2 1], 'FillAlpha',0.3);
legend('Historic data','0.05 and 0.95 quantiles','Predicted mean','fontsize',30,'interpreter','latex','Location','Best');
hold off
 
figure(3)
hold on
ax = gca;
%ylim([-50, 500+max(lb)]);
plot(MONTHS(1:newn+1),new_cases(1:newn+1),'LineWidth',2)
hold on
plot(MONTHS(newn:(newn+size(x_pred,2)-1) ), new_cases(newn:newn+total_time-1),'k','LineWidth', 3)%future
hold on 
%rsim = randperm(nsim,25);
[a1, a2]= sort(x_pred(:,rem));
%rsim = randperm(50,25);
rsim =[11, 21, 29, 1, 26, 19, 22, 15, 6, 12, 4, 7, 8, 49, 35, 41, 18, 43, 16, 48, 13, 27, 34, 47, 38];
selfive = randperm(size(rsim,2), 5);
for i = rsim(selfive)
    plot(MONTHS(newn:(newn+size(x_pred,2)-1) ),x_pred(i,:),'k--', 'LineWidth',1)
    hold on
end 
datetick('x', 'mmm-yy') 
ylabel('New chlamydia cases','fontsize',28,'interpreter','latex');
xlabel('Date','fontsize',28,'interpreter','latex');
%title('Number of new gonorrhea cases','fontsize', 45,'interpreter','latex');
legend('Historic data','Observed ``future" data','Possible future paths','fontsize',30,'interpreter','latex','Location','Best');
ax.FontSize = 25;
ax.LabelFontSizeMultiplier = 1.2;
hold off 

function[rate, vola, mu, sigma] = fit_GBM(x_data, dt, alpha)
    dlogx = diff(log(x_data));
    [~, mu, sigma] = ddiv_estimator3(mean(dlogx), std(dlogx), dlogx, alpha);
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
