% multisine_analysis.m  —  multisine time response + bode (2×3 figure)
%
% Row 1: time traces  — w_g(t)  |  h(t) [mm]  |  θ(t) [deg]
% Row 2: bode         — |H_h(f)|              |  |H_θ(f)|  |  phase (both)

here = fileparts(mfilename('fullpath'));
load(fullfile(here,'wing_2dof.mat'), 'A','B','C_out','D_out','fn_struct');

sys = ss(A, B, C_out, D_out);

% ---------- input signal -----------------------------------------
wg_ts = multisine_function_generator();
t     = wg_ts.Time;
wg    = wg_ts.Data;


% ---------- simulate — lsim accepts continuous sys directly ------
y_out  = lsim(sys, wg, t);
h_mm   = y_out(:,1) * 1000;
th_deg = rad2deg(y_out(:,2));

fprintf('Plunge RMS = %.3f mm\n',  rms(h_mm));
fprintf('Pitch  RMS = %.4f deg\n', rms(th_deg));

% ---------- bode curves ------------------------------------------
f  = logspace(-1, 1.7, 600);
om = 2*pi*f;
[mag, phs] = bode(sys, om);
mag_h  = squeeze(mag(1,1,:));   phs_h  = squeeze(phs(1,1,:));
mag_th = squeeze(mag(2,1,:));   phs_th = squeeze(phs(2,1,:));

fh = fn_struct(1);  ft = fn_struct(2);

% ---------- zoom to last 2 periods for time plots ----------------
Tper = 1/0.25;
i_z  = t >= t(end) - 2*Tper;
tz   = t(i_z);

col_h  = [0.85 0.33 0.1];
col_th = [0.2  0.4  0.8];
col_in = [0.2  0.4  0.8];

% ---------- 2×3 figure -------------------------------------------
figure('Name','Multisine Analysis','Position',[80 80 1200 620]);

% --- row 1: time -------------------------------------------------
subplot(2,3,1)
plot(tz, wg(i_z), 'Color',col_in, 'LineWidth',0.8)
xlabel('t [s]');  ylabel('w_g [m/s]');  title('INPUT: gust');    grid on

subplot(2,3,2)
plot(tz, h_mm(i_z),   'Color',col_h,  'LineWidth',1.2)
xlabel('t [s]');  ylabel('h [mm]');    title('TIME: plunge');    grid on

subplot(2,3,3)
plot(tz, th_deg(i_z), 'Color',col_th, 'LineWidth',1.2)
xlabel('t [s]');  ylabel('\theta [deg]'); title('TIME: pitch');  grid on

% --- row 2: bode -------------------------------------------------
subplot(2,3,4)
semilogy(f, mag_h, 'Color',col_h, 'LineWidth',1.8)
xline(fh,'--','Color',[.6 .6 .6]);  xline(ft,'--','Color',[.6 .6 .6])
set(gca,'XScale','log')
xlabel('f [Hz]');  ylabel('|H_h| [m/(m/s)]');  title('BODE: plunge mag');  grid on

subplot(2,3,5)
semilogy(f, mag_th, 'Color',col_th, 'LineWidth',1.8)
xline(fh,'--','Color',[.6 .6 .6]);  xline(ft,'--','Color',[.6 .6 .6])
set(gca,'XScale','log')
xlabel('f [Hz]');  ylabel('|H_\theta| [rad/(m/s)]');  title('BODE: pitch mag');  grid on

subplot(2,3,6)
plot(f, phs_h,  'Color',col_h,  'LineWidth',1.8, 'DisplayName','h')
hold on
plot(f, phs_th, 'Color',col_th, 'LineWidth',1.8, 'DisplayName','\theta')
xline(fh,'--','Color',[.6 .6 .6],'HandleVisibility','off')
xline(ft,'--','Color',[.6 .6 .6],'HandleVisibility','off')
set(gca,'XScale','log')
xlabel('f [Hz]');  ylabel('Phase [deg]');  title('BODE: phase');
legend('Location','southwest');  grid on
