%% Consumer Sentiment Index -- "better off than a year ago" share
% Model vs data, sticky-wage / endogenous-inflation baselines.
%
% PROVENANCE.  The original figure is
%     Fourth Year\Aiyagari RPI\AiyagariOpen\src\Figures\CSI_answers_3.png
% produced by Plot_CSI_answers.m, whose `model` matrix was HARD-CODED.  Those
% numbers did NOT come from any share_better_* series in simulation.jl or
% simulation_interp1.jl -- every one of those uses N or count(~=) as the
% denominator and reproduces nothing.  They came from THIS estimator, in
%     Fourth Year\...\src\Matlab\Plots_histograms.m :
%       u_gain(t,:) = u_sim(t,:) ./ u_sim(t-4,:);
%       numerator   = sum(sum(u_gain(t:t+3,:) >= u_high, 'omitnan'));
%       denominator = sum(sum(~isnan(u_gain(t:t+3,:))));
%       better      = round(100 * numerator / denominator);
%     windows t = 5,9,13,17,21,25 (each t:t+3) -> 2019(burn-in),2020,...,2024
% With u_high anywhere in [1.0000001, 1.0005] this reproduces 30 43 44 29 27 28
% EXACTLY off sim_data_2.mat.  Plots_histograms.m now uses u_high = 1.01, which
% gives 25 35 35 22 22 22 -- it was widened AFTER the figure was made.
%
% NOTHING IS HARD-CODED HERE.  All numbers come from csi_shares_sw.mat, written
% by scripts/main_CSI_sw.jl.  Re-run that and this figure updates.

clear; clc;
S = load('csi_shares_sw.mat');

data   = S.data(:,1)';          % survey "better off" share
fourth = S.fourth_year(:)';     % Fourth Year model, reproduced exactly
epsm_v = S.epsm_fed_v(:)';      % sticky wages, pi fed from data (~9%)
run3_v = S.run3_v(:)';          % sticky wages, pi endogenous (~2%)
epsm_x = S.epsm_fed_x(:)';      % same, utility = flow GHH composite
run3_x = S.run3_x(:)';

%% ---- headline figure: DATA vs NEW BASELINE only ----------------------------
% New baseline = reduced_epsm.jld2 : sticky wages, Taylor rule, epsm fed so
% inflation matches the data.  Utility = value function (like-for-like with the
% original CSI_answers_3).  Same styling as the original figure.
figure('WindowState','maximized','Color',[1 1 1]);
ax = gca;
f1 = plot(epsm_v, 'b-',  'LineWidth', 5); hold on;
f2 = plot(data,   'r-.', 'LineWidth', 5); hold off; grid on
ax.TickLabelInterpreter = 'latex';
xlim([1 6]); xticks(1:6); xticklabels([]);
ax.Clipping = 'off';
ylim([20 50]);                      % set BEFORE the labels: y0 must come from the FINAL
                                     % limits, else the label collides with the y-tick
offset0 = 7; y00 = ax.YLim(1);   % AT the axis bottom; VerticalAlignment top hangs it below
text(1, y00,                'Mean',       'Interpreter','latex', ...
     'HorizontalAlignment','center','VerticalAlignment','top','FontSize',24);
text(1, y00 - 0.2*offset0,  '2008--2019', 'Interpreter','latex', ...
     'HorizontalAlignment','center','VerticalAlignment','top','FontSize',24);
yrs0 = {'2020','2021','2022','2023','2024'};
for k = 2:6
    text(k, y00, yrs0{k-1}, 'HorizontalAlignment','center', ...
        'VerticalAlignment','top','FontSize',24,'Interpreter','latex');
end
ylabel('\%','fontsize',24,'interpreter','latex');
axis square; ax.FontSize = 24;
legend([f1 f2], 'Model', 'Data', 'Location','northeast','Orientation','vertical');
legend('boxoff'); set(legend,'FontSize',24,'interpreter','latex');
set(gcf,'PaperPositionMode','auto');
print(gcf, '-dpng', 'CSI_answers_3_sw_baseline');

%% ---- main figure: value-function utility, the measure the original used ----
figure('WindowState','maximized','Color',[1 1 1]);
ax = gca;
f1 = plot(epsm_v, 'b-',  'LineWidth', 5); hold on;
f2 = plot(run3_v, 'k-',  'LineWidth', 5, 'Marker','o','MarkerSize',10);
f3 = plot(fourth, 'g--', 'LineWidth', 3);
f4 = plot(data,   'r-.', 'LineWidth', 5); hold off; grid on
ax.TickLabelInterpreter = 'latex';
xlim([1 6]); xticks(1:6); xticklabels([]);
ax.Clipping = 'off';
ylim([20 50]);                      % set BEFORE the labels (see above)
offset = 7; y0 = ax.YLim(1);     % AT the axis bottom (see above)
text(1, y0,               'Mean',       'Interpreter','latex', ...
     'HorizontalAlignment','center','VerticalAlignment','top','FontSize',24);
text(1, y0 - 0.2*offset,  '2008--2019', 'Interpreter','latex', ...
     'HorizontalAlignment','center','VerticalAlignment','top','FontSize',24);
years = {'2020','2021','2022','2023','2024'};
for k = 2:6
    text(k, y0, years{k-1}, 'HorizontalAlignment','center', ...
        'VerticalAlignment','top','FontSize',24,'Interpreter','latex');
end
ylabel('\%','fontsize',24,'interpreter','latex');
axis square; ax.FontSize = 24;
legend([f1 f2 f3 f4], 'Sticky wages, $\pi$ fed', 'Sticky wages, $\pi$ endogenous', ...
       'Fourth Year model', 'Data', 'Location','northeast','Orientation','vertical');
legend('boxoff'); set(legend,'FontSize',22,'interpreter','latex');
set(gcf,'PaperPositionMode','auto');
print(gcf, '-dpng', 'CSI_answers_3_sw');

%% ---- companion: the two utility measures, since they disagree ----
figure('WindowState','maximized','Color',[1 1 1]);
ax = gca;
g1 = plot(epsm_v, 'b-',  'LineWidth', 4); hold on;
g2 = plot(epsm_x, 'b--', 'LineWidth', 4);
g3 = plot(run3_v, 'k-',  'LineWidth', 4);
g4 = plot(run3_x, 'k--', 'LineWidth', 4);
g5 = plot(data,   'r-.', 'LineWidth', 4); hold off; grid on
ax.TickLabelInterpreter = 'latex';
xlim([1 6]); xticks(1:6); xticklabels([]);
ax.Clipping = 'off';
y0 = ax.YLim(1);
text(1, y0, 'Mean', 'Interpreter','latex','HorizontalAlignment','center', ...
     'VerticalAlignment','top','FontSize',22);
text(1, y0 - 0.2*offset, '2008--2019','Interpreter','latex', ...
     'HorizontalAlignment','center','VerticalAlignment','top','FontSize',22);
for k = 2:6
    text(k, y0, years{k-1}, 'HorizontalAlignment','center', ...
        'VerticalAlignment','top','FontSize',22,'Interpreter','latex');
end
ylabel('\%','fontsize',24,'interpreter','latex'); axis square; ax.FontSize = 22;
legend([g1 g2 g3 g4 g5], '$\pi$ fed, $V$', '$\pi$ fed, flow $x$', ...
       '$\pi$ endog, $V$', '$\pi$ endog, flow $x$', 'Data', ...
       'Location','northeast'); legend('boxoff');
set(legend,'FontSize',20,'interpreter','latex');
title('{Sensitivity to the utility measure}','fontsize',24,'interpreter','latex');
set(gcf,'PaperPositionMode','auto');
print(gcf, '-dpng', 'CSI_answers_3_sw_umeasure');

%% ---- console table ----
fprintf('\n  BETTER-OFF SHARE (%%)   u_high = %.7f, N = %g, T0 = %g, seed = %g\n', ...
        S.u_high, S.N, S.T0, S.seed);
fprintf('  %-22s %8s %6s %6s %6s %6s %6s\n','','2008-19','2020','2021','2022','2023','2024');
rows = {'Data',data; 'Fourth Year',fourth; 'Sticky, pi fed (V)',epsm_v; ...
        'Sticky, pi endog (V)',run3_v; 'Sticky, pi fed (x)',epsm_x; 'Sticky, pi endog (x)',run3_x};
for i = 1:size(rows,1)
    fprintf('  %-22s %8.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n', rows{i,1}, rows{i,2});
end
