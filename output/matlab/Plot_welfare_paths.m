%% Figures
% Description: Plots official sov. debt
% Author: Leonardo Barreto
% This version: 6/1/2023
% Modified:

% Housekeeping
clear
close all
clc
% Load the simulation data if not already loaded
load('welfare_time_paths.mat');  % contains u_sim

group_1 = 1;
group_2 = 2:4;
group_3 = 5:6;
group_4 = 7:9;
group_5 = 10;

table_means = nan(5,5);
for i=1:5
    eval(['table_means(', int2str(i) ,',1) = mean(mean(welfare_mat(group_', int2str(i) ,',1:3)))'])
    eval(['table_means(', int2str(i) ,',2) = mean(mean(welfare_mat(group_', int2str(i) ,',4:7)))'])
    eval(['table_means(', int2str(i) ,',3) = mean(mean(welfare_mat(group_', int2str(i) ,',8:11)))'])
    eval(['table_means(', int2str(i) ,',4) = mean(mean(welfare_mat(group_', int2str(i) ,',12:14)))'])
    eval(['table_means(', int2str(i) ,',5) = mean(mean(welfare_mat(group_', int2str(i) ,',15:18)))'])
end
%%
xtick_labels = {'2020','2021','2022','2023','2024','2025'};
T = 20;               % original time series length
T_extra = 0;          % quarters before 2020Q2 (2019Q1–2020Q1)
T_total = T + T_extra;
x = 1:T_total;

figure('WindowState','maximized','Color',[1 1 1]);
%subplot(2,3,1)
ax = gca;
f1 = plot([NaN, welfare_mat(1,1:T_total)], 'r-', 'LineWidth', 3);
hold on;
f2 = plot([NaN, welfare_mat(3,1:T_total)], 'k-', 'LineWidth', 3);
f3 = plot([NaN, welfare_mat(5,1:T_total)], 'b-', 'LineWidth', 3);
f4 = plot([NaN, welfare_mat(7,1:T_total)], 'r-', 'LineWidth', 3);
f5 = plot([NaN, welfare_mat(10,1:T_total)], 'r-', 'LineWidth', 3);
grid on
ax.TickLabelInterpreter='latex';
ax.FontSize = 18;
xlim([1 T_total]);
yline(0.0, '-k', 'LineWidth', 1.5);
%ylim([-2 6]);
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('Change in welfare (\%)','fontsize',22,'interpreter','latex');
%title('{Baseline}','fontsize',25,'interpreter','latex')
%numerator   = sum(sum(u_gain(t0:t0+3,:) >= 1, 'omitnan')); 
%denominator = sum(sum(~isnan(u_gain(t0:t0+3,:))));
%ratio_better = round(100 * numerator / denominator);
%label_str = sprintf('Ratio $>1$ = %d\\%%', ratio_better);
%text(1.02,30,label_str,'fontsize',15,'interpreter','latex')

return
subplot(2,3,2)
ax = gca;
f1 = plot(summary_paths(3:end,2), 'r-', 'LineWidth', 3);
hold on;
grid on
ax.TickLabelInterpreter='latex';
ax.FontSize = 18;
xlim([1 10]);
ylim([-2 6]);
xticks(1:1:10);
xtickangle(ax,0)
xticklabels({'','P20','','P40','','P60','','P80','','P99'})
ylabel('Change in welfare (\%)','fontsize',22,'interpreter','latex');
xlabel('Wealth percentile','fontsize',22,'interpreter','latex')
title('{Constant $w$}','fontsize',25,'interpreter','latex')

subplot(2,3,3)
ax = gca;
f1 = plot(summary_paths(3:end,3), 'r-', 'LineWidth', 3);
hold on;
grid on
ax.TickLabelInterpreter='latex';
ax.FontSize = 18;
xlim([1 10]);
ylim([-2 6]);
xticks(1:1:10);
xtickangle(ax,0)
xticklabels({'','P20','','P40','','P60','','P80','','P99'})
ylabel('Change in welfare (\%)','fontsize',22,'interpreter','latex');
xlabel('Wealth percentile','fontsize',22,'interpreter','latex')
title('{Constant $r$}','fontsize',25,'interpreter','latex')

subplot(2,3,4)
ax = gca;
f1 = plot(summary_paths(3:end,4), 'r-', 'LineWidth', 3);
hold on;
grid on
ax.TickLabelInterpreter='latex';
ax.FontSize = 18;
xlim([1 10]);
ylim([-2 6]);
xticks(1:1:10);
xtickangle(ax,0)
xticklabels({'','P20','','P40','','P60','','P80','','P99'})
ylabel('Change in welfare (\%)','fontsize',22,'interpreter','latex');
xlabel('Wealth percentile','fontsize',22,'interpreter','latex')
title('{Constant $\tau^n$}','fontsize',25,'interpreter','latex')

subplot(2,3,5)
ax = gca;
f1 = plot(summary_paths(3:end,5), 'r-', 'LineWidth', 3);
hold on;
grid on
ax.TickLabelInterpreter='latex';
ax.FontSize = 18;
xlim([1 10]);
ylim([-2 6]);
xticks(1:1:10);
xtickangle(ax,0)
xticklabels({'','P20','','P40','','P60','','P80','','P99'})
ylabel('Change in welfare (\%)','fontsize',22,'interpreter','latex');
xlabel('Wealth percentile','fontsize',22,'interpreter','latex')
title('{Constant $T^u$}','fontsize',25,'interpreter','latex')

subplot(2,3,6)
ax = gca;
f1 = plot(summary_paths(3:end,6), 'r-', 'LineWidth', 3);
hold on;
grid on
ax.TickLabelInterpreter='latex';
ax.FontSize = 18;
xlim([1 10]);
ylim([-2 6]);
xticks(1:1:10);
xtickangle(ax,0)
xticklabels({'','P20','','P40','','P60','','P80','','P99'})
ylabel('Change in welfare (\%)','fontsize',22,'interpreter','latex');
xlabel('Wealth percentile','fontsize',22,'interpreter','latex')
title('{Constant $T$}','fontsize',25,'interpreter','latex')

%% FIGURE IN ONLY ONE PLOT

figure('WindowState','maximized','Color',[1 1 1]);
ax = gca;
f1 = plot(summary_paths(3:end,1), 'r-','LineWidth', 6);
hold on;
f2 = plot(summary_paths(3:end,2), 'b-','LineWidth', 4);
f3 = plot(summary_paths(3:end,3), 'Color',[110/255 110/255 110/255],'LineStyle','--','LineWidth', 4);
f4 = plot(summary_paths(3:end,4), 'Color',[0.2 0.5 0.0],'LineStyle','--','LineWidth', 4);
f5 = plot(summary_paths(3:end,5), 'Color',[115/255*0 43/255*0 245/255*0],'LineStyle','--','LineWidth', 4);
f6 = plot(summary_paths(3:end,6), 'Color',[142/255 64/255 58/255],'LineStyle',':','LineWidth', 4);
grid on
ax.TickLabelInterpreter='latex';
ax.FontSize = 18;
xlim([1 10]);
ylim([-2 6]);
xticks(1:1:10);
yline(0.0, '-k', 'LineWidth', 1.5);
xtickangle(ax,0)
xticklabels({'P10','P20','P30','P40','P50','P60','P70','P80','P90','P99'})
axis square
ylabel('Change in welfare (\%)','fontsize',22,'interpreter','latex');
xlabel('Wealth percentile','fontsize',22,'interpreter','latex')
legend([f1 f2 f3 f4 f5 f6],'Baseline','Constant $w$','Constant $r$','Constant $\tau^n$','Constant $T^u$','Constant $T$',...
        'location', 'northeast');
    legend('boxoff')
    set(legend,'FontSize',18,'interpreter','latex')
set(gcf,'PaperPositionMode','auto')
print(gcf, '-dpng', '../Figures/welfare_decomp')