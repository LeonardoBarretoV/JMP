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
load('transition_decomp.mat');  % contains u_sim

%% FIGURE IN ONLY ONE PLOT

figure('WindowState','maximized','Color',[1 1 1]);
ax = gca;
f1 = plot(summary_paths(3:end,1), 'r-','LineWidth', 6);
hold on;
% Constant w
f2 = plot(summary_paths(3:end,2), 'b-','LineWidth', 4);
% Constant r
f3 = plot(summary_paths(3:end,3), 'Color','k','LineStyle','-','LineWidth', 4);
% Constant tau_n
f4 = plot(summary_paths(3:end,4), 'Color',[142/255 64/255 58/255],'LineStyle',':','LineWidth', 4);
% Constant Tu
f5 = plot(summary_paths(3:end,5), 'Color',[115/255 43/255 245/255],'LineStyle','--','LineWidth', 4);
% Constant T
f6 = plot(summary_paths(3:end,6), 'Color',[0.2 0.5 0.0],'LineStyle','--','LineWidth', 4);
grid on
ax.TickLabelInterpreter='latex';
ax.FontSize = 18;
xlim([1 10]);
ylim([-2 4]);
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
summary_paths

%% FIGURE FOR THE SLIDES 
%% FIGURE IN ONLY ONE PLOT

figure('WindowState','maximized','Color',[1 1 1]);
ax = gca;
f1 = plot(summary_paths(3:end,1), 'r-','LineWidth', 6);
hold on;
% Constant r
f3 = plot(summary_paths(3:end,3), 'Color','k','LineStyle','-','LineWidth', 4);
% Constant T
f6 = plot(summary_paths(3:end,6), 'Color',[0.2 0.5 0.0],'LineStyle','--','LineWidth', 4);
grid on
ax.TickLabelInterpreter='latex';
ax.FontSize = 22;
xlim([1 10]);
ylim([-2 4]);
xticks(1:1:10);
yline(0.0, '-k', 'LineWidth', 1.5);
xtickangle(ax,0)
xticklabels({'P10','P20','P30','P40','P50','P60','P70','P80','P90','P99'})
axis square
ylabel('Change in welfare (\%)','fontsize',26,'interpreter','latex');
xlabel('Wealth percentile','fontsize',26,'interpreter','latex')
legend([f1 f3 f6],'Baseline','Constant $r$','Constant $T$',...
        'location', 'northeast');
    legend('boxoff')
    set(legend,'FontSize',20,'interpreter','latex')
set(gcf,'PaperPositionMode','auto')
print(gcf, '-dpng', '../Figures/welfare_decomp_slides')
return

%% FIGURE IN ONLY ONE PLOT

figure('WindowState','maximized','Color',[1 1 1]);
ax = gca;
f1 = plot(summary_paths_diff(3:end,1), 'r-','LineWidth', 6);
hold on;
% Constant w
f2 = plot(summary_paths_diff(3:end,2), 'b-','LineWidth', 4);
% Constant r
f3 = plot(summary_paths_diff(3:end,3), 'Color','k','LineStyle','-','LineWidth', 4);
% Constant tau_n
f4 = plot(summary_paths_diff(3:end,4), 'Color',[0.2 0.5 0.0],'LineStyle','--','LineWidth', 4);
% Constant Tu
f5 = plot(summary_paths_diff(3:end,5), 'Color',[115/255 43/255 245/255],'LineStyle','--','LineWidth', 4);
% Constant T
f6 = plot(summary_paths_diff(3:end,6), 'Color',[142/255 64/255 58/255],'LineStyle',':','LineWidth', 4);
grid on
ax.TickLabelInterpreter='latex';
ax.FontSize = 18;
xlim([1 10]);
%ylim([-2 4]);
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
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-dpng', '../Figures/welfare_decomp')