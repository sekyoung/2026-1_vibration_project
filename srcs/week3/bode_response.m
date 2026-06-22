% bode_response.m  —  Bode plot: plunge h and pitch theta vs gust w_g

here = fileparts(mfilename('fullpath'));
load(fullfile(here,'wing_2dof.mat'), 'A','B','C_out','D_out','fn_struct');

sys = ss(A, B, C_out, D_out);
sys.OutputName = {'h','theta'};
sys.InputName  = {'w_g'};

% frequency range spanning both structural modes
f  = logspace(-1, 1.7, 600);          % 0.1 – 50 Hz
om = 2*pi*f;

[mag, phs] = bode(sys, om);
mag_h  = squeeze(mag(1,1,:));   phs_h  = squeeze(phs(1,1,:));
mag_th = squeeze(mag(2,1,:));   phs_th = squeeze(phs(2,1,:));

fh = fn_struct(1);  ft = fn_struct(2);

col_h  = [0.85 0.33 0.1];
col_th = [0.2  0.4  0.8];

figure('Name','Bode Response','Position',[100 100 760 540]);

subplot(2,1,1)
semilogy(f, mag_h,  'Color',col_h,  'LineWidth',1.8, 'DisplayName','h  (plunge)')
hold on
semilogy(f, mag_th, 'Color',col_th, 'LineWidth',1.8, 'DisplayName','\theta  (pitch)')
xline(fh,'--','Color',[.5 .5 .5],'HandleVisibility','off')
xline(ft,'--','Color',[.5 .5 .5],'HandleVisibility','off')
set(gca,'XScale','log')
ylabel('Magnitude  [m/(m/s)], [rad/(m/s)]')
title('Bode Magnitude');  legend('Location','northwest');  grid on

subplot(2,1,2)
plot(f, phs_h,  'Color',col_h,  'LineWidth',1.8, 'DisplayName','h')
hold on
plot(f, phs_th, 'Color',col_th, 'LineWidth',1.8, 'DisplayName','\theta')
xline(fh,'--','Color',[.5 .5 .5],'HandleVisibility','off')
xline(ft,'--','Color',[.5 .5 .5],'HandleVisibility','off')
set(gca,'XScale','log')
xlabel('Frequency  [Hz]');  ylabel('Phase  [deg]')
title('Bode Phase');  legend('Location','southwest');  
grid on;
exportgraphics(gcf, 'week3_bode_response.png', 'Resolution', 300)