%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     (C) Andrews T. Anum and Michael Pokojovy (2022)       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sam_sizes = [30 50 100 200 300 400];
mean_mdpd_bdp = [10 16 31 61 91 121]; 
sigma_mdpd_bdp = [10 16 31 61 91 121]; 
mean_mcd_bdp = [9 15 28 56 85 112]; 
sigma_mcd_bdp = [9 15 28 55 83 110]; 

vals = [mean_mdpd_bdp; mean_mcd_bdp]';
vals2 = [sigma_mdpd_bdp; sigma_mcd_bdp]';
subplot(2,1,1)
ax = gca;
bar(sam_sizes, vals, 2);
xlabel('Sample size','fontsize',22,'interpreter','latex');
ylabel('Breakdown point','fontsize',22,'interpreter','latex');
legend('MDPD','MCD','fontsize',10, 'Location', 'Best');
ax.FontSize = 25;
ax.LabelFontSizeMultiplier = 1.2;
hold on 
subplot(2,1,2)
ax = gca;
bar(sam_sizes, vals2, 2,'grouped')
xlabel('Sample size','fontsize',22,'interpreter','latex');
ylabel('Breakdown point','fontsize',22,'interpreter','latex');
legend('MDPD','MCD','fontsize',10, 'Location', 'Best');
ax.FontSize = 25;
ax.LabelFontSizeMultiplier = 1.2;
hold off
