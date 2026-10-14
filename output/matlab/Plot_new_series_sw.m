%% [SW] NEW SERIES -- objects that do not exist in the flexible-wage paper
% s_t (hours index), piw_t (nominal wage inflation), and the union labor wedge.
% Deliberately a SEPARATE figure: nothing is forced into the existing layouts.
% Reads : <repo>/output/matlab_sw/model_new_series_sw.mat
%         <repo>/output/matlab_sw/model_time_series_ratio.mat
% Writes: <repo>/output/figures_sw/f4_new_series_sw.png

clear
clc

cd("C:/Users/Leonardo Barreto/Google Drive/GoodNotes/Fifth Year/JMP_inflation/output/matlab_sw")
nw = load('model_new_series_sw.mat');
tm = load('model_time_series_ratio.mat');

s_index = nw.s_index;
piw     = nw.piw;        % row 1 is NaN: piw is defined from j=2 on
wedge   = nw.wedge;
kappa_w = nw.kappa_w;

pi_     = tm.time_series(:,2);
T_fig   = 20;                       % 2020Q2 .. 2025Q1
x       = 1:T_fig;

figure('WindowState','maximized','Color',[1 1 1]);

% (a) hours index s_t
subplot(2,2,1); hold on; box on; grid on;
plot(x, s_index(1:T_fig), 'b-', 'LineWidth', 4);
yline(1, 'k--', 'LineWidth', 2);
title('{(a). Hours index $s_t$}','fontsize',25,'interpreter','latex')
ylabel('$s_t$','fontsize',22,'interpreter','latex')
set(gca,'fontsize',18,'TickLabelInterpreter','latex')

% (b) wage inflation vs price inflation, annualised
subplot(2,2,2); hold on; box on; grid on;
plot(x, 100*((1+piw(1:T_fig)).^4 - 1), 'b-', 'LineWidth', 4);
plot(x, 100*((1+pi_(1:T_fig)).^4 - 1), 'r--', 'LineWidth', 4);
yline(0, 'k:', 'LineWidth', 1.5);
title('{(b). Wage vs price inflation}','fontsize',25,'interpreter','latex')
ylabel('Annualised \%','fontsize',22,'interpreter','latex')
legend({'$\pi^w_t$ (wage)','$\pi_t$ (price, fed)'},'interpreter','latex','fontsize',18,'Location','best')
set(gca,'fontsize',18,'TickLabelInterpreter','latex')

% (c) union labor wedge
subplot(2,2,3); hold on; box on; grid on;
plot(x, wedge(1:T_fig), 'b-', 'LineWidth', 4);
yline(0, 'k--', 'LineWidth', 2);
title('{(c). Labor wedge}','fontsize',25,'interpreter','latex')
ylabel('$s_t[\theta s_t^{1/\nu}D_t - w_tY_t]$','fontsize',18,'interpreter','latex')
set(gca,'fontsize',18,'TickLabelInterpreter','latex')

% (d) implied real wage growth, dlog w
subplot(2,2,4); hold on; box on; grid on;
w_ = tm.time_series(:,15);
dlw = [NaN; log(w_(2:T_fig)./w_(1:T_fig-1))];
plot(x, 100*4*dlw, 'b-', 'LineWidth', 4);
yline(0, 'k--', 'LineWidth', 2);
title('{(d). Real wage growth $\Delta\log w_t$}','fontsize',25,'interpreter','latex')
ylabel('Annualised \%','fontsize',22,'interpreter','latex')
set(gca,'fontsize',18,'TickLabelInterpreter','latex')

sgtitle(sprintf('Sticky-wage block, \kappa_w = %.2f', kappa_w),'fontsize',28,'interpreter','latex')
print(gcf, '-dpng', 'C:/Users/Leonardo Barreto/Google Drive/GoodNotes/Fifth Year/JMP_inflation/output/figures_sw/f4_new_series_sw')
