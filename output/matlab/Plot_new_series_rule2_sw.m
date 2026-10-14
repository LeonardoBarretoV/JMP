%% [SW26] NEW SERIES, RULE-2 BASELINE -- objects that do not exist in the
% flexible-wage paper: s_t (hours index), piw_t (nominal wage inflation), the union
% labor wedge, real wage growth, and the two new panels:
%   (e) inflation with the monetary shock switched off
%   (f) the monetary shock itself
%
% Reads : <repo>/output/matlab_rule2/model_new_series_sw.mat
%         <repo>/output/matlab_rule2/model_time_series_ratio.mat
% Writes: <repo>/output/figures_rule2/f4_new_series_rule2.png
%
% COLOUR CONVENTION (matches Plots_series_v0_rule2_sw.m, where blue = Model and
% red dashed = Data):
%   BLUE   solid      the baseline model run -- the SAME run in every panel
%   GREEN  solid      a different variable from that same run (wage inflation)
%   ORANGE dash-dot   the counterfactual with the monetary shock switched off
%   RED    dashed     DATA ONLY.  This figure plots no data, so no red appears.
%   GREY   dashed     reference lines (target, zero, unity)
%
% ON PANEL (e).  pi_endog re-runs the inflation recursion on the SOLVED baseline
% r path with epsm = 0 throughout.  It answers "what inflation would this rule
% deliver, given this real allocation, with no monetary shocks?".  It is a
% decomposition of the inflation equation holding r fixed, NOT a re-solved
% general-equilibrium counterfactual.  j = 1 is fed from the data in both series,
% so the lines coincide there by construction and separate only from j = 2 on.

clear
clc

root = "C:/Users/Leonardo Barreto/Google Drive/GoodNotes/Fifth Year/JMP_inflation";
cd(root + "/output/matlab_rule2")
nw = load('model_new_series_sw.mat');
tm = load('model_time_series_ratio.mat');

s_index  = nw.s_index;
piw      = nw.piw;          % row 1 is NaN: piw is defined from j=2 on
wedge    = nw.wedge;
kappa_w  = nw.kappa_w;
pi_solv  = nw.pi_solved;
pi_endog = nw.pi_endog;
epsm     = nw.epsm;
pi_bar   = nw.pi_bar;

pi_ = tm.time_series(:,2);          % identical to pi_solv; kept for the old mapping
T_fig = 20;                         % 2020Q2 .. 2025Q1
x = 1:T_fig;

ann = @(p) 100*((1+p).^4 - 1);      % quarterly rate -> annualised %

% ---- palette -------------------------------------------------------------
c_model   = [0.00 0.20 0.80];       % blue   : the baseline run
c_wage    = [0.00 0.50 0.25];       % green  : another variable, same run
c_noshock = [0.93 0.49 0.19];       % orange : the no-shock counterfactual
c_ref     = [0.25 0.25 0.25];       % grey   : reference lines
% c_data  = [0.85 0.10 0.10];       % red    : reserved for DATA (unused here)

figure('WindowState','maximized','Color',[1 1 1]);

% (a) hours index s_t
subplot(2,3,1); hold on; box on; grid on;
plot(x, s_index(1:T_fig), '-', 'Color', c_model, 'LineWidth', 4);
yline(1, '--', 'Color', c_ref, 'LineWidth', 2);
title('{(a). Hours index $s_t$}','fontsize',22,'interpreter','latex')
ylabel('$s_t$','fontsize',20,'interpreter','latex')
set(gca,'fontsize',16,'TickLabelInterpreter','latex')

% (b) wage inflation vs price inflation -- both from the SAME run
subplot(2,3,2); hold on; box on; grid on;
plot(x, ann(pi_(1:T_fig)), '-', 'Color', c_model, 'LineWidth', 4);
plot(x, ann(piw(1:T_fig)), '-', 'Color', c_wage,  'LineWidth', 4);
yline(0, ':', 'Color', c_ref, 'LineWidth', 1.5);
ylim([-2.5 11]);
title('{(b). Wage vs price inflation}','fontsize',22,'interpreter','latex')
ylabel('Annualised \%','fontsize',20,'interpreter','latex')
legend({'$\pi_t$ (price)','$\pi^w_t$ (wage)'},'interpreter','latex', ...
       'fontsize',15,'Location','southeast')
set(gca,'fontsize',16,'TickLabelInterpreter','latex')

% (c) union labor wedge
subplot(2,3,3); hold on; box on; grid on;
plot(x, wedge(1:T_fig), '-', 'Color', c_model, 'LineWidth', 4);
yline(0, '--', 'Color', c_ref, 'LineWidth', 2);
title('{(c). Labor wedge}','fontsize',22,'interpreter','latex')
ylabel('$s_t[\theta s_t^{1/\nu}D_t - w_tY_t]$','fontsize',16,'interpreter','latex')
set(gca,'fontsize',16,'TickLabelInterpreter','latex')

% (d) implied real wage growth, dlog w
subplot(2,3,4); hold on; box on; grid on;
w_ = tm.time_series(:,15);
dlw = [NaN; log(w_(2:T_fig)./w_(1:T_fig-1))];
plot(x, 100*4*dlw, '-', 'Color', c_model, 'LineWidth', 4);
yline(0, '--', 'Color', c_ref, 'LineWidth', 2);
title('{(d). Real wage growth $\Delta\log w_t$}','fontsize',22,'interpreter','latex')
ylabel('Annualised \%','fontsize',20,'interpreter','latex')
set(gca,'fontsize',16,'TickLabelInterpreter','latex')

% (e) [NEW] endogenous inflation with the monetary shock switched off
subplot(2,3,5); hold on; box on; grid on;
[~, pk] = max(abs(pi_solv(2:T_fig) - pi_bar)); pk = pk + 1;
shr = 100*(pi_endog(pk) - pi_bar)/(pi_solv(pk) - pi_bar);
plot(x, ann(pi_solv(1:T_fig)),  '-',  'Color', c_model,   'LineWidth', 4);
plot(x, ann(pi_endog(1:T_fig)), '-.', 'Color', c_noshock, 'LineWidth', 4);
yline(ann(pi_bar), '--', 'Color', c_ref, 'LineWidth', 2);
ylim([0 14]);
title({'(e). Inflation without the shock', ...
       [num2str(shr,'%.2f') '\% of the peak deviation']}, ...
      'fontsize',22,'interpreter','latex')
ylabel('Annualised \%','fontsize',20,'interpreter','latex')
legend({'$\pi_t$ (baseline run)','$\pi_t$ with $\varepsilon^m_t=0$','target $\bar\pi$'}, ...
       'interpreter','latex','fontsize',13,'Location','northeast')
set(gca,'fontsize',16,'TickLabelInterpreter','latex')

% (f) [NEW] the monetary shock itself
subplot(2,3,6); hold on; box on; grid on;
plot(x, 400*epsm(1:T_fig), '-', 'Color', c_model, 'LineWidth', 4);
yline(0, '--', 'Color', c_ref, 'LineWidth', 2);
title('{(f). Monetary shock $\varepsilon^m_t$}','fontsize',22,'interpreter','latex')
ylabel('Annualised p.p.','fontsize',20,'interpreter','latex')
set(gca,'fontsize',16,'TickLabelInterpreter','latex')

sgtitle(sprintf('Sticky-wage block, rule-2 baseline, $\\kappa_w = %.2f$, $\\phi_\\pi = %.2f$', ...
        kappa_w, nw.phi_pi),'fontsize',26,'interpreter','latex')

print(gcf, '-dpng', char(root + "/output/figures_rule2/f4_new_series_rule2"))
fprintf('wrote f4_new_series_rule2.png   (peak share without shocks: %.2f%%)\n', shr);
