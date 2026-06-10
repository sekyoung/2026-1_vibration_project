function result_plot()
% result_plot  5x2 summary figure for Week-4 nonlinear wing analysis.
%
%   Layout
%   -------
%   Col 1 (linear / constant-k)              | Col 2 (nonlinear, gust option 4)
%   Row 1  gust family (impulse_fg sp1)      | modified gust input u(t)
%   Row 2  tuning sweep (impulse_fg sp2)     | plunge h(t)
%   Row 3  critical response (impulse_fg sp3)| pitch theta(t)
%   Row 4  spring-1 constant k               | spring-1 tangent k1(t)
%   Row 5  equiv. bending moment, const-k    | equiv. bending moment, varying-k

here = fileparts(mfilename('fullpath'));

% ---- 1. Run source functions (suppress internal figures) -------------
gust = impulse_function_generator('plot', false);

nl = nonlinear_response(4);
close(gcf);

mout = moment_response();
close(gcf);

% ---- 2. Linear critical response (for rows 3 and 4) -----------------
mdl     = load(fullfile(here, '..', 'week3', 'wing_2dof.mat'), ...
               'A','B','C_out','D_out');
sys_lin = ss(mdl.A, mdl.B, mdl.C_out, mdl.D_out);
wg_crit = gust.wg(:, gust.crit.idx_h);
y_lin   = lsim(sys_lin, wg_crit, gust.t);   % columns: h [m], theta [rad]

% constant spring-1 stiffness (linear assumption: k1 independent of deflection)
nlp      = load(fullfile(here, 'wing_nl.mat'), 'k1');
k1_const = nlp.k1 * ones(size(gust.t));

% ---- 3. Colours (shared with source scripts) -------------------------
col_h  = [0.85 0.33 0.10];
col_th = [0.20 0.40 0.80];
col_k  = [0.30 0.60 0.30];

% ---- 4. Figure -------------------------------------------------------
figure('Name','Week-4 Results','Position',[80 60 1200 1100]);

% -- Row 1 (full width): 1-cosine gust family --------------------------
subplot(5,2,[1 2])
show = unique(round(linspace(1, numel(gust.H), 5)));
hold on
for j = show
    plot(gust.t, gust.wg(:,j), 'LineWidth',1.1, ...
         'DisplayName', sprintf('H = %.1f m', gust.H(j)));
end
plot(gust.t, gust.wg(:,gust.crit.idx_h), 'k', 'LineWidth',2, ...
     'DisplayName','critical H');
hold off
xlabel('t [s]');  ylabel('w_g [m/s]');  grid on
title('1-cosine gust family  (U_d_s/2)[1-cos(\pi s/H)]')
legend('Location','northeast')
xlim([0 max(gust.t)])

% -- Row 2, Col 1: gradient-tuning sweep -------------------------------
subplot(5,2,3)
yyaxis left
plot(gust.H, gust.resp.h_peak, '-o', 'Color',col_h, 'MarkerSize',3, ...
     'DisplayName','peak |h|');
ylabel('peak |h| [mm]')
yyaxis right
plot(gust.H, gust.resp.th_peak, '-s', 'Color',col_th, 'MarkerSize',3, ...
     'DisplayName','peak |\theta|');
ylabel('peak |\theta| [deg]')
xlabel('H [m]');  grid on
xline(gust.crit.H_h,  '--', 'Color',col_h,  'HandleVisibility','off')
xline(gust.crit.H_th, '--', 'Color',col_th, 'HandleVisibility','off')
title(sprintf('Tuning sweep:  H_h = %.2f m,  H_\\theta = %.2f m', ...
              gust.crit.H_h, gust.crit.H_th))
legend('Location','northeast')

% -- Row 2, Col 2: nonlinear h(t) --------------------------------------
subplot(5,2,4)
plot(nl.t, nl.h*1000, 'Color',col_h, 'LineWidth',1.5, 'DisplayName','h(t)')
yline(0,'Color',[.7 .7 .7],'HandleVisibility','off');  grid on
xlabel('t [s]');  ylabel('h [mm]')
title('Plunge h(t)  (nonlinear gust response)')
legend('Location','northeast')

% -- Row 3, Col 1: linear critical response (h and theta) --------------
subplot(5,2,5)
yyaxis left
plot(gust.t, y_lin(:,1)*1000, 'Color',col_h, 'LineWidth',1.5, ...
     'DisplayName','h');
ylabel('h [mm]')
yyaxis right
plot(gust.t, rad2deg(y_lin(:,2)), 'Color',col_th, 'LineWidth',1.5, ...
     'DisplayName','\theta');
ylabel('\theta [deg]')
xlabel('t [s]');  grid on
title('Linear response to critical gust')
legend('Location','northeast')

% -- Row 3, Col 2: nonlinear theta(t) ----------------------------------
subplot(5,2,6)
plot(nl.t, rad2deg(nl.theta), 'Color',col_th, 'LineWidth',1.5, ...
     'DisplayName','\theta(t)')
yline(0,'Color',[.7 .7 .7],'HandleVisibility','off');  grid on
xlabel('t [s]');  ylabel('\theta [deg]')
title('Pitch \theta(t)  (nonlinear gust response)')
legend('Location','northeast')

% -- Row 4, Col 1: constant spring stiffness (linear assumption) -------
subplot(5,2,7)
plot(gust.t, k1_const, 'Color',col_k, 'LineWidth',1.5, ...
     'DisplayName','k_1 (const)')
xlabel('t [s]');  ylabel('k_1 [N/m]');  grid on
title('Spring-1 stiffness  (linear: k_1 = const)')
legend('Location','northeast')

% -- Row 4, Col 2: nonlinear tangent stiffness k1eff(t) ----------------
subplot(5,2,8)
plot(nl.t, nl.k1eff, 'Color',col_k, 'LineWidth',1.5, ...
     'DisplayName','k_1(t)')
xlabel('t [s]');  ylabel('k_1(t) [N/m]');  grid on
title('Spring-1 tangent stiffness  k_1 + 3k_n_l_1\delta_1^2')
legend('Location','northeast')

% -- Row 5, Col 1: equivalent bending moment – constant k --------------
subplot(5,2,9)
plot(mout.t_lin, mout.M_lin, 'Color',col_h, 'LineWidth',1.5, ...
     'DisplayName','M_e_q')
yline(0,'Color',[.7 .7 .7],'HandleVisibility','off');  grid on
xlabel('t [s]');  ylabel('M_e_q [N m]')
title(sprintf('Equiv. bending moment – constant k  (L = %.2f m)', mout.L))
legend('Location','northeast')

% -- Row 5, Col 2: equivalent bending moment – varying k ---------------
subplot(5,2,10)
plot(mout.t_nl, mout.M_nl, 'Color',col_h, 'LineWidth',1.5, ...
     'DisplayName','M_e_q')
yline(0,'Color',[.7 .7 .7],'HandleVisibility','off');  grid on
xlabel('t [s]');  ylabel('M_e_q [N m]')
title('Equiv. bending moment – varying k  (cubic springs)')
legend('Location','northeast')

end
