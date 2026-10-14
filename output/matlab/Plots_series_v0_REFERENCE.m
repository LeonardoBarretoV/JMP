%% Figures
% Description: Plots official sov. debt
% Author: Leonardo Barreto
% This version: 6/1/2023
% Modified:

% Housekeeping
clear
% close all
clc

[year, ~, ~]= xlsread('Data_plots.xlsx',3,'A2:A37');
[quarter, ~, ~]= xlsread('Data_plots.xlsx',3,'B2:B37');
[un_series_data, ~, ~]= xlsread('Data_plots.xlsx',3,'C14:C37');
[pi_series_data, ~, ~]= xlsread('Data_plots.xlsx',3,'D14:D37');
[misery, ~, ~]= xlsread('Data_plots.xlsx',3,'E2:E37');
[gdp_g, ~, ~]= xlsread('Data_plots.xlsx',3,'F14:F37');
[csi, ~, ~]= xlsread('Data_plots.xlsx',3,'G2:G37');
[b_series_data, ~, ~]= xlsread('Data_plots.xlsx',3,'H14:H37');
[bstar_series_data, ~, ~]= xlsread('Data_plots.xlsx',3,'I14:I37');
[cons_g, ~, ~]= xlsread('Data_plots.xlsx',3,'J14:J37');
[ui, ~, ~]= xlsread('Data_plots.xlsx',3,'K14:K37');
[gov_exp, ~, ~]= xlsread('Data_plots.xlsx',3,'L14:L37');
[int_rate, ~, ~]= xlsread('Data_plots.xlsx',3,'N14:N37');
[tips_5yr, ~, ~]= xlsread('Data_plots.xlsx',3,'O14:O37');
[tips_10yr, ~, ~]= xlsread('Data_plots.xlsx',3,'P14:P37');
[wage_data, ~, ~]= xlsread('Data_plots.xlsx',3,'Q2:Q37');
wage_data = 100.*(wage_data ./ wage_data(16));

T_w = length(wage_data);
time = (1:T_w)';                        % Create time trend if not given
% Regress real wages on a constant and linear time trend
X = [ones(T_w,1), time];
b = X \ wage_data;                       % OLS coefficients
trend = X * b;                        % Fitted trend
w_detrended = (wage_data - trend)./100;
w_level = ones(T_w,1).*(1+w_detrended);
w_level = w_level(13:end);
%Load model
cd("C:\Users\Leonardo Barreto\Google Drive\GoodNotes\Fourth Year\Aiyagari RPI\AiyagariOpen\src\Matlab")

% Load data from .mat file
data_model = load('model_time_series_ratio.mat');

% Extract each series
ts = data_model.time_series;

y              = ts(:, 1);
pi_            = ts(:, 2);
un             = ts(:, 3);
q              = ts(:, 4);
b              = ts(:, 5);
bd             = ts(:, 6);
bstar          = ts(:, 7);
l              = ts(:, 8);
lump_sum_T     = ts(:, 9);
lump_sum_Tu    = ts(:,10);
r              = ts(:,11);
r_k            = ts(:,12);
k              = ts(:,13);
n              = ts(:,14);
w              = ts(:,15);
C              = ts(:,16);
lab_tax        = ts(:,17);
tau_k          = ts(:,18);
tau_n          = ts(:,19);
X              = ts(:,20);
G              = ts(:,21);

% Steady states
ss_init  = data_model.ss_init;
ss_final = data_model.ss_final;


% Read full matrix from Sheet 1
%data_model = xlsread('model_series.xlsx', 1);  % or use readmatrix for newer MATLAB versions

T = 20;               % original time series length
T_extra = 5;          % quarters before 2020Q2 (2019Q1–2020Q1)
T_total = T + T_extra;
x = 1:T_total;

% Time axis labels (adjust if you want quarterly labels)
xtick_labels = {'2019','2020','2021','2022','2023','2024','2025'};

% Steady-state values from ss_init (sheet 2 of Excel)
pi_ss     = 100 * ((1 + ss_init(2))^4 - 1);
un_ss    = 100 * ss_init(3);
debt_ss  = 100 * ss_init(4) * ss_init(5) / (4 * ss_init(1));
ext_debt_ss  = 100 * ss_init(4) * ss_init(7) / (4 * ss_init(1));
l_ss     = 100 * ss_init(4) * ss_init(8) / (4 * ss_init(1));
T_ss     = 100 * ss_init(9) / (4 * ss_init(1));
Tu_ss    = 100 * ss_init(10) / (4 * ss_init(1));
r_ss = 100 * ((1 + ss_init(11))^4 - 1);  % still annualized (as this is a level)
w_ss = ss_init(15);
c_ss = ss_init(16);
y_ss = ss_init(1);
k_ss = ss_init(13);
n_ss = ss_init(14);
G_ss = 100 *ss_init(21)/( ss_init(1));
tau_n_ss = 100*(1 - ss_init(19));
tau_k_ss = 100*ss_init(18);

% Precompute shifted debt: q(t) * B(t+1)
value_debt = q(1:T) .* [b(2:T); b(T)];
value_ext_debt = q(1:T) .* [bstar(2:T); bstar(T)];

% Construct padded time series with 5 quarters of steady state
pi_series   = [repmat(pi_ss, T_extra, 1); 100 * ((1 + pi_(1:T)).^4 - 1)];
un_series  = [repmat(un_ss, T_extra, 1); 100 * un(1:T)];
debt_series = [repmat(debt_ss, T_extra, 1); 100 * value_debt ./ (4 * y(1:T))];
ext_debt_series = [repmat(ext_debt_ss, T_extra, 1); 100 * value_ext_debt ./ (4 * y(1:T))];
l_series   = [repmat(l_ss, T_extra, 1); 100 * q(1:T) .* l(1:T) ./ (4 * y(1:T))];
T_series   = [repmat(T_ss, T_extra, 1); 100 * lump_sum_T(1:T) ./ (4 * y(1:T))];
Tu_series  = [repmat(Tu_ss, T_extra, 1); 100 * lump_sum_Tu(1:T) ./ (4 * y(1:T))];
r_series = [repmat(r_ss, T_extra, 1); 100 * ((1 + r(1:T)).^4 - 1)]; % leave interest rate annualized
w_series = [ones(T_extra, 1); w(1:T) ./ w_ss];
lab_income_series = [ones(T_extra, 1); w(1:T).* n(1:T) ./ (w_ss * n_ss)];
k_series = [repmat(k_ss / (4*y_ss), T_extra, 1); k(1:T) ./(4 *y_ss)];
n_series = [ones(T_extra, 1); n(1:T) ./ n_ss];  % normalized so n_ss = 1
gdp_growth = [100*(y(1) / y_ss - 1); 100 * (y(2:T) - y(1:T-1)) ./ y(1:T-1)];
cons_growth = [100*(C(1) / c_ss - 1); 100 * (C(2:T) - C(1:T-1)) ./ C(1:T-1)];
G_series  = [repmat(G_ss, T_extra, 1); 100 * G(1:T) ./ y(1:T)];
tau_n_series = [repmat(tau_n_ss, T_extra, 1); 100.*(1-tau_n(1:T))];
tau_k_series = [repmat(tau_k_ss, T_extra, 1); 100.*(tau_k(1:T))];
q_series = [ones(T_extra, 1); q(1:T)];  % normalized so n_ss = 1

% after_tax_lab_income = [ones(T_extra, 1); (w(1:T).* n(1:T) - lab_tax(1:T) + lump_sum_T(1:T) + lump_sum_Tu(1:T) )./ (w_ss * n_ss - 0.0011459707922496613 + 0.00022531029000003922)];
% 
% lab_income_series - lab_tax_series + T_series + Tu_series;
gdp_growth_series = [zeros(T_extra + 1, 1); gdp_growth];
cons_growth_series = [zeros(T_extra + 1, 1); cons_growth];

% Define shading range for 2020Q2
shock_start = 6 - 0.5;  % half step before quarter
shock_end = 6 + 0.5;    % half step after quarter
y_shade = [-1e6, 1e6];  % Large range to cover y-axis

% X coordinates for vertical bar
x_fill = [shock_start, shock_end, shock_end, shock_start];


% -------- PLOT --------
figure('WindowState','maximized','Color',[1 1 1]);

% (a) Inflation
subplot(2,3,1)
ax = gca;
%yl = [0.0, max(pi_series)+1];  % preserve original y-axis range
bar(6, 200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, pi_series, 'b-', 'LineWidth', 4);
grid on
plot(x, [pi_series_data; NaN], 'r--', 'LineWidth', 4);
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
yticks(0:2:10);
ylim([0 10])  % restore original y-limits
xlim([1 length(x)])
ylabel('\%','fontsize',22,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
title('{(a). Inflation}','fontsize',24,'interpreter','latex')


% (b) Unemployment
subplot(2,3,2)
ax = gca;
bar(6, 100, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, un_series, 'b-', 'LineWidth', 4);
grid on
plot(x, [un_series_data; NaN], 'r--', 'LineWidth', 4);
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',22,'interpreter','latex')
ylim([2 14])
xlim([1 length(x)])
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
title('{(b). Unemployment}','fontsize',24,'interpreter','latex')

% (c) Debt
subplot(2,3,3)
ax = gca;
bar(6, 200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, debt_series, 'b-', 'LineWidth', 4); 
grid on
plot(x, [b_series_data; NaN], 'r--', 'LineWidth', 4);
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',22,'interpreter','latex')
ylim([100 135])
xlim([1 length(x)])
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
title('{(c). Debt-to-GDP}','fontsize',24,'interpreter','latex')

% (d)  External Debt
subplot(2,3,4)
ax = gca;
bar(6, 200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
p1  = plot(x, ext_debt_series, 'b-', 'LineWidth', 4); 
grid on
p2 = plot(x, [bstar_series_data; NaN], 'r--', 'LineWidth', 4);
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',22,'interpreter','latex')
ylim([26 36])
xlim([1 length(x)])
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
title('{(d). External Debt-to-GDP}','fontsize',24,'interpreter','latex')
legend([p1 p2],'Model $\hspace{0.5em}$','Data$\hspace{0.5em}$');
set(legend,...
    'Position',[0.440118532992161 0.0110994873080163 0.160908767349009 0.04607330017689],...
    'Orientation','horizontal',...
    'Interpreter','latex',...
    'FontSize',25);
legend('boxoff');
set(gcf,'PaperPositionMode','auto')

% (c) GDP Growth (non-annualized)
subplot(2,3,5)
ax= gca;
bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
bar(6, -200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
plot(x, gdp_growth_series(2:end), 'b-', 'LineWidth', 4);
plot(x, [gdp_g-0.5; NaN], 'r--', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',22,'interpreter','latex')
ylim([-10 8])
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
title('{(e). GDP Growth}','fontsize',24,'interpreter','latex')
xlim([1 length(x)])

% (d) Consumption Growth (non-annualized)
subplot(2,3,6)
ax = gca; 
bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
bar(6, -200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
p1 = plot(x, cons_growth_series(2:end), 'b-', 'LineWidth', 4);
p2 = plot(x, [cons_g-0.5; NaN], 'r--', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',22,'interpreter','latex')
ylim([-12 10])
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
title('{(f). Consumption Growth}','fontsize',24,'interpreter','latex')
xlim([1 length(x)])
set(gcf,'PaperPositionMode','auto')
print(gcf, '-dpng', '../Figures/f1_baseline')

%% SECOND PLOT
wide_param = -0.03;
figure('WindowState','maximized','Color',[1 1 1]);

% (b) Wage Relative to SS
subplot(2,2,1)
yyaxis left
ax = gca; 
bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, w_series, 'b-', 'LineWidth', 4);
ylim([0.98 1.05])

yyaxis right
p2 = plot(x, [w_level; NaN], 'r--', 'LineWidth', 4);
ylim([0.95 1.1])
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('','fontsize',22,'interpreter','latex')
%ylim([0.98 1.05])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
title('{(a). Wage Rate}','fontsize',24,'interpreter','latex')

subplot(2,2,2)
ax = gca; 
yyaxis left
bar(6, 200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
bar(6, -200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
p1 = plot(x, [r_series(1:5); NaN; r_series(7:end)], 'b-', 'LineWidth', 4);
plot(x, [r_series(1:5); r_series(1); r_series(7:end)], 'b--', 'LineWidth', 4);
ylim([0.6 1.2])
ylabel('\%','fontsize',22,'interpreter','latex')

yyaxis right
p2 = plot(x, [tips_10yr; NaN], 'r--', 'LineWidth', 4);
ylim([-1 3])

hold off
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',22,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
title('{(b). Interest Rate}','fontsize',24,'interpreter','latex')
legend([p1 p2],'Model (left-axis) $\hspace{0.5em}$','Data (right-axis)$\hspace{0.5em}$');
set(legend,...
    'Position',[0.440118532992161 0.0110994873080163 0.160908767349009 0.04607330017689],...
    'Orientation','horizontal',...
    'Interpreter','latex',...
    'FontSize',25);
legend('boxoff');
set(gcf,'PaperPositionMode','auto')


% (b) Labor income
subplot(2,2,3)
ax = gca; 
bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, lab_income_series, 'b-', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('','fontsize',22,'interpreter','latex')
ylim([0.90 1.02])
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
title('{(c). Labor Income}','fontsize',24,'interpreter','latex')


% Panel (a) - Capital / Output_ss
subplot(2,2,4)
ax = gca; 
bar(6, 200, 1, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none'); hold on
plot(x, k_series, 'b-', 'LineWidth', 4);
grid on
xticks(1:4:T+T_extra)
xticklabels(xtick_labels)
xtickangle(ax,0)
xlim([1 length(x)])
%ylabel('$k_t / y_{ss}$','fontsize',22,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
ylim([4.2 4.35])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
title('{(d). Capital-to-Initial Output}','fontsize',24,'interpreter','latex')

print(gcf, '-dpng', '../Figures/f2_baseline')
return
%% THIRD PLOT
figure('WindowState','maximized','Color',[1 1 1]);
% (a) Interest Rate

% (a) Government expenditure
subplot(2,3,1)
ax = gca;
bar(6, 200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
p1 = plot(x, G_series, 'b-', 'LineWidth', 4); 
p2 = plot(x, [gov_exp; NaN], 'r--', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',22,'interpreter','latex')
title('{(a). Government Expenditure}','fontsize',24,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
%ylim([0 8])
ylim([4.5 6.5])
legend([p1 p2],'Model $\hspace{0.5em}$','Data$\hspace{0.5em}$');
set(legend,...
    'Location','northeast',...
    'Orientation','vertical',...
    'Interpreter','latex',...
    'FontSize',19);
legend('boxoff');

% (e) Stimulus Checks
subplot(2,3,2)
ax = gca;
bar(6, 200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, T_series, 'b-', 'LineWidth', 4); 
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',22,'interpreter','latex')
title('{(b). Stimulus Checks}','fontsize',24,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
ylim([0 8])

subplot(2,3,3)
ax = gca;
bar(6, 200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, Tu_series, 'b-', 'LineWidth', 4); 
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',22,'interpreter','latex')
title('{(c). Unemployment Insurance}','fontsize',24,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
ylim([0 8])

subplot(2,3,4)
ax = gca;
bar(6, 200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, tau_n_series, 'b-', 'LineWidth', 4); 
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',22,'interpreter','latex')
title('{(d). Labor Tax}','fontsize',24,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
ylim([0 30])

subplot(2,3,5)
ax = gca;
bar(6, 200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
bar(6, -200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
plot(x, tau_k_series, 'b-', 'LineWidth', 4); 
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',22,'interpreter','latex')
title('{(e). Capital Tax}','fontsize',24,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
ylim([-60 20])


set(gcf,'PaperPositionMode','auto')
print(gcf, '-dpng', '../Figures/f3_baseline')


%% SECOND PLOT SLIDES

figure('WindowState','maximized','Color',[1 1 1]);

% Panel (b) - Labor / Labor_ss
subplot(2,3,1)
ax = gca; 
bar(6, 200, 1, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none'); hold on
plot(x, n_series, 'b-', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
xlim([1 length(x)])
%ylabel('$n_t / n_{ss}$','fontsize',22,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
ylim([0.85 1.05])
title('{(a). Effective Labor}','fontsize',24,'interpreter','latex')


% (b) Wage Relative to SS
subplot(2,3,2)
ax = gca; 
bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, w_series, 'b-', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('','fontsize',22,'interpreter','latex')
ylim([0.98 1.05])
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
title('{(b). Wage Rate}','fontsize',24,'interpreter','latex')


% (b) Labor income
subplot(2,3,3)
ax = gca; 
bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, lab_income_series, 'b-', 'LineWidth', 4);
%plot(x, after_tax_lab_income, 'r-', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('','fontsize',22,'interpreter','latex')
% ylim([0.90 1.02])
ylim([0.90 1.9])
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
title('{(c). Labor Income}','fontsize',24,'interpreter','latex')


% Panel (a) - Capital / Output_ss
subplot(2,3,4)
ax = gca; 
bar(6, 200, 1, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none'); hold on
plot(x, k_series, 'b-', 'LineWidth', 4);
grid on
xticks(1:4:T+T_extra)
xticklabels(xtick_labels)
xtickangle(ax,0)
xlim([1 length(x)])
%ylabel('$k_t / y_{ss}$','fontsize',22,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
ylim([4.2 4.35])
title('{(d). Capital-to-Initial Output}','fontsize',24,'interpreter','latex')


subplot(2,3,5)
ax = gca; 
bar(6, 200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
bar(6, -200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
plot(x, [r_series(1:5); NaN; r_series(7:end)], 'b-', 'LineWidth', 4);
plot(x, [r_series(1:5); r_series(1); r_series(7:end)], 'b--', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',22,'interpreter','latex')
ylim([0.6 1.2])
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
title('{(e). Interest Rate}','fontsize',24,'interpreter','latex')


% subplot(2,3,6)
% bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
% plot(x, lab_tax_series, 'b-', 'LineWidth', 4);
% grid on
% xticks(1:4:T_total)
% xticklabels(xtick_labels)
% ylabel('\%','fontsize',22,'interpreter','latex')
% ylim([10 25])
% title('{(f). Labor Tax}','fontsize',25,'interpreter','latex')
% ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = 22;
% xlim([1 length(x)])

subplot(2,3,6)
ax = gca; 
bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, q_series, 'b-', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
%ylabel('\%','fontsize',22,'interpreter','latex')
ylim([0.9 1.02])
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
title('{(f). Bond Price}','fontsize',24,'interpreter','latex')
% 
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-dpng', '../Figures/f2_baseline')

return
%% SECOND PLOT OCTOBER VERSION (PRE-MANUEL FEEDBACK)

%% SECOND PLOT
wide_param = 0.003;
figure('WindowState','maximized','Color',[1 1 1]);

% Panel (b) - Labor / Labor_ss
subplot(2,3,1)
ax = gca; 
bar(6, 200, 1, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none'); hold on
plot(x, n_series, 'b-', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
xlim([1 length(x)])
%ylabel('$n_t / n_{ss}$','fontsize',22,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
ylim([0.85 1.05])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
title('{(a). Effective Labor}','fontsize',24,'interpreter','latex')


% (b) Wage Relative to SS
subplot(2,3,2)
yyaxis left
ax = gca; 
bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, w_series, 'b-', 'LineWidth', 4);
ylim([0.98 1.05])

yyaxis right
p2 = plot(x, [w_level; NaN], 'r--', 'LineWidth', 4);
ylim([0.95 1.1])
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('','fontsize',22,'interpreter','latex')
%ylim([0.98 1.05])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
title('{(b). Wage Rate}','fontsize',24,'interpreter','latex')


% (b) Labor income
subplot(2,3,3)
ax = gca; 
bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, lab_income_series, 'b-', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('','fontsize',22,'interpreter','latex')
ylim([0.90 1.02])
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
title('{(c). Labor Income}','fontsize',24,'interpreter','latex')


% Panel (a) - Capital / Output_ss
subplot(2,3,4)
ax = gca; 
bar(6, 200, 1, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none'); hold on
plot(x, k_series, 'b-', 'LineWidth', 4);
grid on
xticks(1:4:T+T_extra)
xticklabels(xtick_labels)
xtickangle(ax,0)
xlim([1 length(x)])
%ylabel('$k_t / y_{ss}$','fontsize',22,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
ylim([4.2 4.35])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
title('{(d). Capital-to-Initial Output}','fontsize',24,'interpreter','latex')


subplot(2,3,5)
ax = gca; 
yyaxis left
bar(6, 200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
bar(6, -200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
p1 = plot(x, [r_series(1:5); NaN; r_series(7:end)], 'b-', 'LineWidth', 4);
plot(x, [r_series(1:5); r_series(1); r_series(7:end)], 'b--', 'LineWidth', 4);
ylim([0.6 1.2])
ylabel('\%','fontsize',22,'interpreter','latex')

% yyaxis right
% plot(x, [int_rate; NaN], 'r--', 'LineWidth', 4);
% ylim([-4 5])

yyaxis right
p2 = plot(x, [tips_10yr; NaN], 'r--', 'LineWidth', 4);
ylim([-1 3])

% yyaxis right
% plot(x, [tips_5yr; NaN], 'r--', 'LineWidth', 4);
% ylim([-2 3])

hold off
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',22,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
title('{(e). Interest Rate}','fontsize',24,'interpreter','latex')
legend([p1 p2],'Model (left-axis) $\hspace{0.5em}$','Data (right-axis)$\hspace{0.5em}$');
set(legend,...
    'Position',[0.440118532992161 0.0110994873080163 0.160908767349009 0.04607330017689],...
    'Orientation','horizontal',...
    'Interpreter','latex',...
    'FontSize',25);
legend('boxoff');
set(gcf,'PaperPositionMode','auto')

% subplot(2,3,6)
% bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
% plot(x, lab_tax_series, 'b-', 'LineWidth', 4);
% grid on
% xticks(1:4:T_total)
% xticklabels(xtick_labels)
% ylabel('\%','fontsize',22,'interpreter','latex')
% ylim([10 25])
% title('{(f). Labor Tax}','fontsize',25,'interpreter','latex')
% ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = 22;
% xlim([1 length(x)])

subplot(2,3,6)
ax = gca; 
bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, q_series, 'b-', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
%ylabel('\%','fontsize',22,'interpreter','latex')
ylim([0.9 1.02])
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
title('{(f). Bond Price}','fontsize',24,'interpreter','latex')

set(gcf,'PaperPositionMode','auto')
print(gcf, '-dpng', '../Figures/f2_baseline')

%% SECOND PLOT IN SLIDES WITH DATA (FOR THE FED PRESENTATION)


%% SECOND PLOT WITH 4 SUBPLOTS
wide_param = 0.003;
figure('WindowState','maximized','Color',[1 1 1]);

% Panel (b) - Labor / Labor_ss
subplot(2,3,1)
ax = gca; 
bar(6, 200, 1, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none'); hold on
plot(x, n_series, 'b-', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
xlim([1 length(x)])
%ylabel('$n_t / n_{ss}$','fontsize',22,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
ylim([0.85 1.05])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
title('{(a). Effective Labor}','fontsize',24,'interpreter','latex')


% (b) Wage Relative to SS
subplot(2,3,2)
yyaxis left
ax = gca; 
bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
p1 = plot(x, w_series, 'b-', 'LineWidth', 4);
ylim([0.98 1.05])

yyaxis right
p2 = plot(x, [w_level; NaN], 'r--', 'LineWidth', 4);
ylim([0.95 1.1])
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('','fontsize',22,'interpreter','latex')
%ylim([0.98 1.05])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
title('{(b). Wage Rate}','fontsize',24,'interpreter','latex')
legend([p1 p2],'Model (left) $\hspace{0.5em}$','Data (right)$\hspace{0.5em}$');
set(legend,...
    'Location', 'northeast',...
    'Orientation','vertical',...
    'Interpreter','latex',...
    'FontSize',20);
legend('boxoff');
set(gcf,'PaperPositionMode','auto')



% (b) Labor income
subplot(2,3,3)
ax = gca; 
bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, lab_income_series, 'b-', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('','fontsize',22,'interpreter','latex')
ylim([0.90 1.02])
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
title('{(c). Labor Income}','fontsize',24,'interpreter','latex')


% Panel (a) - Capital / Output_ss
subplot(2,3,4)
ax = gca; 
bar(6, 200, 1, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none'); hold on
plot(x, k_series, 'b-', 'LineWidth', 4);
grid on
xticks(1:4:T+T_extra)
xticklabels(xtick_labels)
xtickangle(ax,0)
xlim([1 length(x)])
%ylabel('$k_t / y_{ss}$','fontsize',22,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
ylim([4.2 4.35])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
title('{(c). Capital-to-Initial Output}','fontsize',24,'interpreter','latex')


subplot(2,3,5)
ax = gca; 
yyaxis left
bar(6, 200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
bar(6, -200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
p1 = plot(x, [r_series(1:5); NaN; r_series(7:end)], 'b-', 'LineWidth', 4);
plot(x, [r_series(1:5); r_series(1); r_series(7:end)], 'b--', 'LineWidth', 4);
ylim([0.6 1.2])
ylabel('\%','fontsize',22,'interpreter','latex')

% yyaxis right
% plot(x, [int_rate; NaN], 'r--', 'LineWidth', 4);
% ylim([-4 5])

yyaxis right
p2 = plot(x, [tips_10yr; NaN], 'r--', 'LineWidth', 4);
ylim([-1 3])

% yyaxis right
% plot(x, [tips_5yr; NaN], 'r--', 'LineWidth', 4);
% ylim([-2 3])

hold off
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',22,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
title('{(d). Interest Rate}','fontsize',24,'interpreter','latex')
legend([p1 p2],'Model (left) $\hspace{0.5em}$','Data (right)$\hspace{0.5em}$');
set(legend,...
    'Location', 'northwest',...
    'Orientation','vertical',...
    'Interpreter','latex',...
    'FontSize',20);
legend('boxoff');
set(gcf,'PaperPositionMode','auto')

% subplot(2,3,6)
% bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
% plot(x, lab_tax_series, 'b-', 'LineWidth', 4);
% grid on
% xticks(1:4:T_total)
% xticklabels(xtick_labels)
% ylabel('\%','fontsize',22,'interpreter','latex')
% ylim([10 25])
% title('{(f). Labor Tax}','fontsize',25,'interpreter','latex')
% ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = 22;
% xlim([1 length(x)])

subplot(2,3,6)
ax = gca; 
bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
plot(x, q_series, 'b-', 'LineWidth', 4);
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
%ylabel('\%','fontsize',22,'interpreter','latex')
ylim([0.9 1.02])
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 19;
xlim([1 length(x)])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
title('{(f). Bond Price}','fontsize',24,'interpreter','latex')

set(gcf,'PaperPositionMode','auto')
print(gcf, '-dpng', '../Figures/f2_baseline_slides')

%% SECOND PLOT WITH 2 SUBPLOTS
wide_param = 0.003;
figure('WindowState','maximized','Color',[1 1 1]);

% Panel (b) - Labor / Labor_ss
subplot(1,2,1)
yyaxis left
ax = gca; 
bar(6, 10, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
p1 = plot(x, w_series, 'b-', 'LineWidth', 4);
ylim([0.98 1.05])

yyaxis right
p2 = plot(x, [w_level; NaN], 'r--', 'LineWidth', 4);
ylim([0.95 1.1])
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('','fontsize',23,'interpreter','latex')
%ylim([0.98 1.05])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 22;
xlim([1 length(x)])
title('{(a). Wage Rate}','fontsize',26,'interpreter','latex')
legend([p1 p2],'Model (left-axis) $\hspace{0.5em}$','Data (right-axis)$\hspace{0.5em}$');
set(legend,...
    'Location', 'northeast',...
    'Orientation','vertical',...
    'Interpreter','latex',...
    'FontSize',20);
legend('boxoff');
set(gcf,'PaperPositionMode','auto')
axis square

subplot(1,2,2)
ax = gca; 
yyaxis left
bar(6, 200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on
bar(6, -200, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
p1 = plot(x, [r_series(1:5); NaN; r_series(7:end)], 'b-', 'LineWidth', 4);
plot(x, [r_series(1:5); r_series(1); r_series(7:end)], 'b--', 'LineWidth', 4);
ylim([0.6 1.2])
ylabel('\%','fontsize',22,'interpreter','latex')

% yyaxis right
% plot(x, [int_rate; NaN], 'r--', 'LineWidth', 4);
% ylim([-4 5])

yyaxis right
p2 = plot(x, [tips_10yr; NaN], 'r--', 'LineWidth', 4);
ylim([-1 3])

% yyaxis right
% plot(x, [tips_5yr; NaN], 'r--', 'LineWidth', 4);
% ylim([-2 3])

hold off
grid on
xticks(1:4:T_total)
xticklabels(xtick_labels)
xtickangle(ax,0)
ylabel('\%','fontsize',23,'interpreter','latex')
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 22;
xlim([1 length(x)])
p = get(gca, 'Position');             % [left bottom width height]
set(gca, 'Position', [p(1)-wide_param p(2) p(3)+wide_param p(4)]);
title('{(b). Interest Rate}','fontsize',26,'interpreter','latex')
legend([p1 p2],'Model (left-axis) $\hspace{0.5em}$','Data (right-axis)$\hspace{0.5em}$');
set(legend,...
    'Location', 'northwest',...
    'Orientation','vertical',...
    'Interpreter','latex',...
    'FontSize',20);
legend('boxoff');
set(gcf,'PaperPositionMode','auto')
axis square

print(gcf, '-dpng', '../Figures/f2_baseline_slides_two')