%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     (C) Andrews T. Anum and Michael Pokojovy (2022)       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

FILE = readtable('EP.STDS.csv');
new_cases = table2array(FILE(:,2));

logdiff = diff(log(new_cases));
m = size(logdiff,1);
rem = 24;
newn = m - rem;

alpha = 0.5;
[~, mu, sigma] = ddiv_estimator3(mean(logdiff(1:newn)), std(logdiff(1:newn)), logdiff(1:newn), alpha);

z = (logdiff - mu)/sigma;

testalpha = 0.05;
n = newn;
testalphastar = (1- (1-testalpha)^(1/n));%^(1/n);

%2P(z>|z_i|)
pvals1 = 2*(erfc(abs(z)/sqrt(2))/2);

issignif = find(pvals1 < testalphastar);
% significant if 2P < testalphastar
% no outliers if none of them is significant
outnumb = size(issignif,1);

issignif2 = find(pvals1 < testalpha);
% significant if 2P < testalphastar
% no outliers if none of them is significant
outnumb2 = size(issignif2,1);

color_order = get(gca, 'colororder');
figure(1)
hold on
ax = gca;
plot(1:size(z,1), z, 'ro', 'MarkerSize', 9);
ylim([-8 10]);
xlim([1 size(z,1)]);
hold on 
plot(1:size(z,1),  icdf('Normal', (1-testalpha/2),0 ,1)*ones(1,size(z,1)), 'k--', 'LineWidth',2);
hold on
plot(1:size(z,1),  icdf('Normal', (1-testalphastar/2),0 ,1)*ones(1,size(z,1)), 'LineWidth',2, 'Color', color_order(1, :));
hold on
plot(1:size(z,1), -icdf('Normal', (1-testalpha/2),0 ,1)*ones(1,size(z,1)), 'k--', 'LineWidth',2);
hold on
plot(1:size(z,1), -icdf('Normal', (1-testalphastar/2),0 ,1)*ones(1,size(z,1)), 'LineWidth',2, 'Color', color_order(1, :));
hold on
%plot(1:n, -quantile(z, 1-testalphastar/2)*ones(1,n), 'b', 'LineWidth',2);
%plot(1:n, -icdf('Normal', (1-testalphastar/2),0 ,1)*ones(1,n),  'LineWidth',2);
plot(newn*ones(1,2), [-8 10],'b:', 'LineWidth',2)
hold on
plot(issignif2, z(issignif2), 'r*', 'MarkerSize', 15); 
xlabel('Index', 'interpreter','latex')
ylabel('Robust $z$-scores', 'interpreter','latex');
legend('$z$-scores','per-comparison limit','family-wise limit','fontsize',30,'interpreter','latex','Location','southwest');
ax.FontSize = 20;
ax.LabelFontSizeMultiplier = 1.5;
hold off

alpha = 0.5;
[~, mu, sigma] = ddiv_estimator3(median(logdiff(1:newn)), iqr(logdiff(1:newn))/1.349, logdiff(1:newn), alpha);
yval    = logdiff(1:newn);
z_score = (yval - mu)./sigma;
y       = sort(z_score);
probs = linspace(0, 1, size(yval,1)+2);
probs = probs(2:(size(probs,2)-1) );
xval    = norminv(probs);
line = xval;

% robust qqplot 
figure(2)
hold on
plot(xval, y,'ro')
hold on
plot(xval, line, 'LineWidth', 2)
ylabel('Empirical quantiles','interpreter','latex', 'fontsize',30)
xlabel('Theoretical quantiles', 'interpreter','latex', 'fontsize',30)
title('Robust $z$-scores', 'interpreter','latex', 'fontsize', 25)
hold off

figure(3)
hold on
ax = gca;
[f, xi] = ksdensity(logdiff(1:newn));
plot(xi, f, 'LineWidth', 3, 'Color', color_order(1, :))
hold on
plot(xi, normpdf(xi, mu, sigma),'k--', 'LineWidth', 3)
ylabel('Kernel density estimate')
xlabel('sample data')
legend('KDE','Guassian pdf','fontsize',30,'interpreter','latex','Location','best');
ax.FontSize = 20;
ax.LabelFontSizeMultiplier = 1.5;
hold off
