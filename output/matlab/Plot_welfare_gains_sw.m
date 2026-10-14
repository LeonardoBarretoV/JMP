%% [SW] WELFARE GAINS -- sticky-wage version
% The ORIGINAL Plot_welfare_gains.m hardcodes its numbers as literal arrays
% rather than loading a .mat, so this is not a path-only copy: the decile
% vector below is the kappa_w = 0.05 result, computed from
% baseline_stage3_kw05.jld2 with the same bucketing statistics.jl uses.
% welfare_paper is the post-[Z-NORM] Stage 0 column, shown for reference.
% NOTE: the paper draft's figures were generated from the PRE-[Z-NORM]
% Fourth Year vintage; use the Stage 0 column here, not those.

clear
clc

% CE welfare in composite units (%), by wealth decile P10..P100
welfare_paper = [ 3.0468; 2.6746; 2.6033; 2.5260; 2.2147; 1.8556; 1.4999; 1.0056; 0.4978; -0.2019];
welfare_sw    = [ 2.9897; 2.6237; 2.5515; 2.4540; 2.1418; 1.7741; 1.4151; 0.9181; 0.4010; -0.3097];

figure('WindowState','maximized','Color',[1 1 1]);
hold on; box on; grid on;
f1 = plot(welfare_paper, 'r-',  'LineWidth', 6);
f2 = plot(welfare_sw,    'b--', 'LineWidth', 6);
yline(0, 'k:', 'LineWidth', 2);
xlabel('Wealth decile','fontsize',26,'interpreter','latex');
ylabel('Change in welfare (\%)','fontsize',26,'interpreter','latex');
xticks(1:10); xticklabels({'P10','P20','P30','P40','P50','P60','P70','P80','P90','P100'});
legend([f1 f2], {'Flexible wages (Stage 0)','Sticky wages ($\kappa_w=0.05$)'}, ...
       'interpreter','latex','fontsize',22,'Location','best');
set(gca,'fontsize',20,'TickLabelInterpreter','latex')
title('Welfare gains by wealth decile','fontsize',28,'interpreter','latex')
print(gcf, '-dpng', 'C:/Users/Leonardo Barreto/Google Drive/GoodNotes/Fifth Year/JMP_inflation/output/figures_sw/welfare_gains_sw')
