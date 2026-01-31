%% Figures
% Description: Plots official sov. debt
% Author: Leonardo Barreto
% This version: 6/1/2023
% Modified:

% Housekeeping
clear
% close all
clc

cd("C:\Users\Leonardo Barreto\Google Drive\GoodNotes\Fourth Year\Aiyagari RPI\AiyagariOpen\src\Matlab")

[year, ~, ~]= xlsread('Data_intro.xlsx',3,'A2:A143');
[quarter, ~, ~]= xlsread('Data_intro.xlsx',3,'B2:B143');
[unemp, ~, ~]= xlsread('Data_intro.xlsx',3,'C2:C143');
[inf, ~, ~]= xlsread('Data_intro.xlsx',3,'D2:D143');
[misery, ~, ~]= xlsread('Data_intro.xlsx',3,'E2:E143');
[csi, ~, ~]= xlsread('Data_intro.xlsx',3,'F2:F143');
[debt_gdp, ~, ~]= xlsread('Data_intro.xlsx',3,'G2:G143');

t_init = 61; %2005
xlabels1 = cell(1,length(debt_gdp));
init_y = year(t_init);
init_q = quarter(1);
end_y = year(end);
end_q = quarter(end);

ticks = init_y:5:end_y;
xticks_every5yrs = 1:5*4:length(unemp(t_init:end));
T = length(unemp(t_init:end));   % final x value
x0 = T - 7;                     % start 20 periods before the end
y0 = 2;                         % bottom of your y-axis
w = 7;                          % width: 20 periods
h = 3;  

figure('WindowState','maximized','Color',[1 1 1]);
    fig = subplot(1,2,1)
    ax = gca;
    %bar(parr_g,1, 'FaceColor',[217/255 217/255 217/255], "EdgeColor",[217/255 217/255 217/255])
    %text(62,12.5,'Corr. = 0.76','interpreter','latex','fontsize',20)
    f1 = plot(unemp(t_init:end), 'b-.', 'LineWidth', 3);
    hold on;
    f2 = plot(inf(t_init:end), 'r-', 'LineWidth', 3);
    rectangle('Position',[x0, y0, w, h], ...
          'EdgeColor',[0 0 0], ...   % black
          'LineWidth',4, ...
          'FaceColor','none');
    hold off    
    grid on
    ax.TickLabelInterpreter='latex';
    ax.FontSize = 22;
    xlim([1 length(unemp(t_init:end))]);
    xticks(xticks_every5yrs)
    xticklabels(ticks)
    ylim([-2 14])
    ylabel('\%','fontsize',23,'interpreter','latex');
    legend([f1 f2],'Unemployment','Inflation',...
        'location', 'northwest');
    legend('boxoff')
    set(legend,'FontSize',20,'interpreter','latex')  
% after you call title(...)
h = title('{(a). Unemployment \& Inflation}','fontsize',26,'interpreter','latex');
set(h,'Units','normalized');
h.Position(2) = 1.04;   % >1 pushes the title further above the axes (try 1.04–1.10)
    axis square

%     fig = subplot(1,3,2)
%     ax = gca;
%     %bar(parr_g,1, 'FaceColor',[217/255 217/255 217/255], "EdgeColor",[217/255 217/255 217/255])
%     %text(62,12.5,'Corr. = 0.76','interpreter','latex','fontsize',20)
%     f1 = plot(inf, 'b-', 'LineWidth', 3);
%     grid on
%     ax.TickLabelInterpreter='latex';
%     ax.FontSize = 18;
%     xlim([1 length(debt_gdp)]);
%     xticks(xticks_every5yrs)
%     xticklabels(ticks)
%     %ylim([-2 16])
%     ylabel('\%','fontsize',22,'interpreter','latex');
%     title('{(b). Inflation}','fontsize',24,'interpreter','latex')
%     axis square

fig = subplot(1,2,2)
    ax = gca;
    %bar(parr_g,1, 'FaceColor',[217/255 217/255 217/255], "EdgeColor",[217/255 217/255 217/255])
    %text(62,12.5,'Corr. = 0.76','interpreter','latex','fontsize',20)
%     bar([13:18], 200, 10, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
%     bar([13:18], -200, 10, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
    f1 = plot((misery(t_init:end) - mean(misery(t_init:end)))/std(misery(t_init:end)), 'b-', 'LineWidth', 3, 'LineStyle','-.');
    hold on;
    f2 = plot(-(csi(t_init:end) - mean(csi(t_init:end)))/std(csi(t_init:end)), 'r', 'LineWidth', 3, 'LineStyle','-');
    f3 = plot(zeros(100,1), 'k', 'LineWidth', 2.0);
    hold off
    grid on
    ax.TickLabelInterpreter='latex';
    ax.FontSize = 22;
    xlim([1 length(unemp(t_init:end))]);
    xticks(xticks_every5yrs)
    xticklabels(ticks)
    ylim([-2 3])
    %ylabel('\%','fontsize',22,'interpreter','latex');
    ylabel('Standardized units','fontsize',23,'interpreter','latex');
    legend([f1 f2],'Unempl. + Inflation','(-) CSI',...
        'location', 'northwest');
    legend('boxoff')
    set(legend,'FontSize',20,'interpreter','latex')  
h = title('{(b). Consumer Sentiment Index}','fontsize',26,'interpreter','latex');
set(h,'Units','normalized');
h.Position(2) = 1.04;
    axis square
set(gcf,'PaperPositionMode','auto')
print(gcf, '-dpng', '../Figures/Figure_Intro')
