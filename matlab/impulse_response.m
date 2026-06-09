% impulse_response.m  —  Impulse response: plunge h and pitch theta

here = fileparts(mfilename('fullpath'));
load(fullfile(here,'wing_2dof.mat'), 'A','B','C_out','D_out','fn_struct');

sys = ss(A, B, C_out, D_out);
sys.OutputName = {'h','theta'};
sys.InputName  = {'w_g'};

% time vector: long enough for transient to decay
t_end = 8;
[y, t] = impulse(sys, t_end);     % y: [Nt x 2 x 1]

h_imp  = y(:,1) * 1000;           % [mm·s]  (impulse has units output·s)
th_imp = rad2deg(y(:,2));          % [deg·s]

fh = fn_struct(1);  ft = fn_struct(2);

col_h  = [0.85 0.33 0.1];
col_th = [0.2  0.4  0.8];

figure('Name','Impulse Response','Position',[100 100 760 420]);

subplot(2,1,1)
plot(t, h_imp,  'Color',col_h,  'LineWidth',1.5)
yline(0,'Color',[.7 .7 .7])
xlabel('t  [s]');  ylabel('h  [mm\cdots]')
title(sprintf('Plunge  (f_h = %.2f Hz,  f_\\theta = %.2f Hz)', fh, ft))
grid on

subplot(2,1,2)
plot(t, th_imp, 'Color',col_th, 'LineWidth',1.5)
yline(0,'Color',[.7 .7 .7])
xlabel('t  [s]');  ylabel('\theta  [deg\cdots]')
title('Pitch')
grid on
