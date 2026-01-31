%% Figure 1
% Description: Plot Figure 1
% Author: Leonardo Barreto
% This version: 1/17/2023
% Modified:

% Housekeeping
clear
%close all
clc

welfare_1 = [   3.0308372387383593
  2.643223567871582
  2.5365741199166316
  2.49527015848334
  2.184437118857585
  1.8181813041162003
  1.473384102038734
  0.9876306945686555
  0.49033844052506964
 -0.19520245305899886
];

init_y = 2016;
init_m = 3;
end_y = 2023;
end_m = 2;

ticks = init_y:1:end_y;

%%
figure('WindowState','maximized','Color',[1 1 1]);
    fig = subplot(1,1,1)
    ax = gca;
    f1 = plot(welfare_1(1:end), 'r-', 'LineWidth', 6);
    hold on;
    f2 = plot(zeros(10,1), 'k', 'LineWidth', 2.5);
    hold off
    grid on
    ax.TickLabelInterpreter='latex';
    ax.FontSize = 22;
    xlim([1 10]);
    xticks([1:1:10]);
    xticklabels({'P10','P20','P30','P40','P50','P60','P70','P80','P90','P99'})
    ylim([-2 4])
    ylabel('Change in welfare (\%)','fontsize',26,'interpreter','latex');
    xlabel('Wealth percentile','fontsize',26,'interpreter','latex')
    axis square
    %text(6.5,3.2,'Mean = 1.82\%','fontsize',24,'interpreter','latex')
    %text(6.5,2.9,'Minimum = -1.77\%','fontsize',24,'interpreter','latex')
    %legend([f1 f2],'No inflation','Inflation','location','southeast','orientation','vertical');
    %set(legend,'FontSize',22,'interpreter','latex')
    %legend('boxoff');
    %print(gcf,'-dpng','Figure_Welfare_locks')    
set(gcf,'PaperPositionMode','auto')
print(gcf, '-dpng', '../Figures/welfare_gains')    
return
%%
 welfare_2 = [ -1.0663115828731946
 -1.0126635657845926
 -0.9452922667853009
 -0.8916077546892128
 -0.8385313716774565
 -0.7707758734005528
 -0.6860886267303974
 -0.5620238953519167
 -0.36508355457009867
  0.2700898669788725
];

welfare_3 = [
 -1.0504470895576736
 -1.0097065670807113
 -0.9454273208622858
 -0.8776118121478393
 -0.8139070787807904
 -0.7260834485824705
 -0.6213315670715246
 -0.4710587074972834
 -0.24001579386915636
  0.46666723018131073];

welfare_4 = [
 -1.3317331814615092
 -1.2814443135918356
 -1.2100263914975116
 -1.1310431761075501
 -1.0484554429945057
 -0.9398113978581965
 -0.811110945281749
 -0.6363324637957345
 -0.37665809387070276
  0.4121116599568084];

welfare_5 = [ -1.4922578779798592
 -1.4348316530967977
 -1.3612130043116406
 -1.2828401853535234
 -1.1926318607412947
 -1.0791763358658695
 -0.9422219207442984
 -0.7600798620840878
 -0.49314885899330946
  0.32525315151172335
    ]
figure('WindowState','maximized','Color',[1 1 1]);
    fig = subplot(1,1,1)
    ax = gca;
    %bar(parr_g,1, 'FaceColor',[217/255 217/255 217/255], "EdgeColor",[217/255 217/255 217/255])
    %text(62,12.5,'Corr. = 0.76','interpreter','latex','fontsize',20)
    f1 = plot(welfare_1(1:end), 'r-', 'LineWidth', 6);
    hold on;
    f2 = plot(welfare_2(1:end), 'b--', 'LineWidth', 4);
    f3 = plot(welfare_3(1:end), 'Color',[0.2 0.5 0.0],'LineStyle','--', 'LineWidth', 4);
    f4 = plot(welfare_4(1:end), 'Color',[115/255 43/255 245/255],'LineStyle',':', 'LineWidth', 4);
    f5 = plot(welfare_5(1:end), 'Color',[100/255 100/255 100/255],'LineStyle',':', 'LineWidth', 4);    
    f6 = plot(zeros(10,1), 'k', 'LineWidth', 2.5);
    hold off
    grid on
    ax.TickLabelInterpreter='latex';
    ax.FontSize = 18;
    xlim([1 10]);
    xticks([1:1:10]);
    xticklabels({'P10','P20','P30','P40','P50','P60','P70','P80','P90','P99'})
    ylim([-2 4])
    ylabel('Change in welfare (\%)','fontsize',22,'interpreter','latex');
    xlabel('Wealth percentile','fontsize',22,'interpreter','latex')
    axis square
    legend([f1 f2 f3 f4 f5],'2020Q2','2021Q2','2022Q2','2023Q2', '2024Q2',...
        'location', 'northeast');
    legend('boxoff')
    set(legend,'FontSize',20,'interpreter','latex')
    %legend([f1 f2],'No inflation','Inflation','location','southeast','orientation','vertical');
    %set(legend,'FontSize',22,'interpreter','latex')
    %legend('boxoff');
    %print(gcf,'-dpng','Figure_Welfare_locks')    
set(gcf,'PaperPositionMode','auto')
print(gcf, '-dpng', '../Figures/welfare_gains') 
return

    %% Figure try

% Example data: rows = wealth percentiles, columns = time
% Replace this with your actual matrix
W = rand(10, 6);  % 10 percentiles x 6 time periods

% X-axis: percentiles
percentiles = linspace(10, 100, size(W,1));

% Create colormap (e.g., blue fading to light blue)
cmap = flipud(gray(size(W,2)));

figure;
hold on;
for t = 1:size(W,2)
    plot(percentiles, W(:,t), 'LineWidth', 2, 'Color', cmap(t,:));
end
hold off;

xlabel('Wealth Percentile');
ylabel('Welfare (Consumption Equivalent)');
title('Welfare Distribution over Time');
colormap(cmap);
colorbar('Ticks', linspace(0,1,size(W,2)), ...
         'TickLabels', arrayfun(@(x) sprintf('t=%d', x), 1:size(W,2), 'UniformOutput', false));
    
