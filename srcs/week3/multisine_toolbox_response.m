% multisine_toolbox_response.m

here = fileparts(mfilename('fullpath'));
load(fullfile(here,'wing_2dof.mat'), 'A','B','C_out','D_out','fn_struct');

sys = ss(A, B, C_out, D_out);

% ---------- multisine design -------------------------------------
f0  = 0.25;  kk = 1:80;  fk = kk * f0;   K = length(kk);
Ak  = 1.5 / sqrt(K/2) * ones(1, K);       % flat, RMS = 1.5 m/s
phi_k = -pi * kk .* (kk-1) / K;           % Schroeder phases

fs = 400;  dt = 1/fs;  Tper = 1/f0;  nper = 8;
N  = round(fs * Tper * nper);
t  = (0:N-1)' * dt;
wg = sum(Ak .* sin(2*pi*fk.*t + phi_k), 2);

% ---------- simulate with lsim (ZOH) ----------------------------
sysd  = c2d(sys, dt, 'zoh');
y_out = lsim(sysd, wg, t);        % [N x 2]: h, theta

fprintf('Plunge RMS = %.3f mm\n',  rms(y_out(:,1))*1000);
fprintf('Pitch  RMS = %.4f deg\n', rad2deg(rms(y_out(:,2))));

% ---------- bode: continuous curve + values at excited lines -----
fc      = logspace(-1, 1.5, 600);
[Mc, ~] = bode(sys, 2*pi*fc);       % continuous analytical curve
[Mk, ~] = bode(sys, 2*pi*fk);       % at each multisine line

Hh_c = squeeze(Mc(1,1,:));          % plunge channel, continuous
Hh_k = squeeze(Mk(1,1,:))';         % plunge channel, at fk  [1xK]

% ---------- amplitudes directly from bode + design ---------------
Yamp = Hh_k .* Ak * 1000;           % output amplitude [mm] = |H|*|U|

% ---------- natural frequencies ----------------------------------
fh_nat = fn_struct(1);
ft_nat = fn_struct(2);

% ---------- 3-panel figure ---------------------------------------
figure('Name','Multisine FRF (Toolbox)','Position',[100 100 1150 380]);

subplot(1,3,1)
stem(fk, Ak, 'filled', 'MarkerSize', 3, 'Color', [0.2 0.4 0.8])
xlim([0 21]);  xlabel('f [Hz]');  ylabel('|w_g| [m/s]')
title('INPUT: multisine (flat)');  grid on

subplot(1,3,2)
semilogy(fc, Hh_c, 'Color',[0.85 0.33 0.1], 'LineWidth',1.8, ...
         'DisplayName','bode (analytical)')
hold on
scatter(fk, Hh_k, 20, [0.1 0.1 0.5], 'filled', ...
        'DisplayName','bode at f_k')
xline(fh_nat,'--','Color',[.5 .5 .5],'HandleVisibility','off')
xline(ft_nat,'--','Color',[.5 .5 .5],'HandleVisibility','off')
set(gca,'XScale','log');  xlabel('f [Hz]');  ylabel('|H_h(f)|')
title('FILTER: bode magnitude');  legend('Location','southwest');  grid on

subplot(1,3,3)
stem(fk, Yamp, 'filled', 'MarkerSize', 3, 'Color', [0.1 0.1 0.1])
xline(fh_nat,'--','Color',[.5 .5 .5])
xline(ft_nat,'--','Color',[.5 .5 .5])
xlim([0 21]);  xlabel('f [Hz]');  ylabel('|h| [mm]')
title('OUTPUT: response lines');  grid on

out_dir = fullfile(here,'..','imgs');
if ~isfolder(out_dir), mkdir(out_dir); end
saveas(gcf, fullfile(out_dir,'multisine_toolbox_response.png'));
fprintf('Saved.\n');
