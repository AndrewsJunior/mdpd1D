%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              (C) Andrews T. Anum and Michael Pokojovy (2022)                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

load('boxplots.mat')
figure(1)
set(gcf, 'PaperUnits', 'centimeters');
xSize = 24; ySize = 15;
xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
set(gcf, 'PaperPosition', [xLeft yTop xSize ySize]);
set(gcf, 'Position',[0 0 xSize*50 ySize*50]);
        
for i=1:length(sam_sizes)
    subplot_tight(3,2,i,[0.08 0.095]);
    begin = 1 + (i-1)*5;
    ending = 5*i;
    h1= boxplot(mean_dpd(begin:ending,:)','Widths',0.4,'Labels',...
            {' 0.1','0.25','0.50','0.75','1.00'},'Color', 'r', 'symbol', '', 'whisker', inf);
    set(h1,{'linew'},{1.0},'Linestyle','-');
    hold on
    h = boxplot(mean_mcd(begin:ending,:)','Widths',0.4,'Labels',...
            {' 0.1','0.25','0.50','0.75','1.00'},'Color', 'b', 'symbol', '', 'whisker', inf);
    set(h,{'linew'},{1.0},'Linestyle','--');
    n = num2str(sam_sizes(i));
    title({['Location: MDPD vs MCD, $n=$ ', n]}, 'fontsize',24,'interpreter','latex');
    ylabel('${\mathrm{\hat{\mu}}}$','fontsize',22,'interpreter','latex');
    xlabel('${\mathrm{\alpha}}$','fontsize',22,'interpreter','latex');
    legend([h1(1);h(1)],'MDPD','MCD','FontSize',0.01, 'Position', [0.5 0.98 0.001 0.002],'Orientation', 'horizontal');
    hold off
end

%%%for sigma
figure(2)
set(gcf, 'PaperUnits', 'centimeters');
xSize = 24; ySize = 15;
xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
set(gcf, 'PaperPosition', [xLeft yTop xSize ySize]);
set(gcf, 'Position',[0 0 xSize*50 ySize*50]);

for i=1:length(sam_sizes)
    subplot_tight(3,2,i, [0.08 0.095]);    
    begin = 1 + (i-1)*5;
    ending = 5*i;
    h1= boxplot(sigma_dpd(begin:ending,:)','Widths',0.5,'Labels',...
            {' 0.1','0.25','0.50','0.75','1.00'},'Color', 'r', 'symbol', '', 'whisker', inf);
    set(h1,{'linew'},{1.0},'Linestyle','-');
    hold on

    h = boxplot(sigma_mcd(begin:ending,:)','Widths',0.5,'Labels',...
            {' 0.1','0.25','0.50','0.75','1.00'},'Color', 'b', 'symbol', '', 'whisker', inf);
    set(h,{'linew'},{1.0},'Linestyle','--');
        
    n = num2str(sam_sizes(i));
    title({['Scale: MDPD vs MCD, $n=$ ', n]}, 'fontsize',24,'interpreter','latex');
    ylabel('${\mathrm{\hat{\sigma}}}$','fontsize',22,'interpreter','latex');
    xlabel('${\mathrm{\alpha}}$','fontsize',22,'interpreter','latex');
    legend([h1(1);h(1)],'MDPD','MCD','fontsize',0.01, 'Position', [0.5 0.98 0.001 0.002],'Orientation', 'horizontal');
    hold off
end
