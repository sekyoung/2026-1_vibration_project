function out = moment_response(L)
% moment_response  Equivalent root bending moment: impulse & multisine.
%
%   out = moment_response(L)
%       L : wing length [m]  (default 0.50)
%
%   Produces two horizontally-aligned plots.
%   Each plot shows two curves: constant k (linear) and nonlinear k (cubic).

if nargin < 1, L = 0.50; end

P     = load_params();
t_end = 8;
fs    = 2000;

% --- impulse --------------------------------------------------------
[t_i_lin, h_i_lin, th_i_lin] = linear_case(P, 1, t_end, fs);
r_imp = nonlinear_response(1);  close(gcf);

M_i_lin = equivalent_moment(h_i_lin,   th_i_lin,   P, L, false);
M_i_nl  = equivalent_moment(r_imp.h,   r_imp.theta, P, L, true);

% --- multisine ------------------------------------------------------
[t_m_lin, h_m_lin, th_m_lin] = linear_case(P, 2, t_end, fs);
r_ms = nonlinear_response(2);  close(gcf);

M_m_lin = equivalent_moment(h_m_lin, th_m_lin,  P, L, false);
M_m_nl  = equivalent_moment(r_ms.h,  r_ms.theta, P, L, true);

out = struct( ...
    'L',       L, ...
    't_i_lin', t_i_lin,  'M_i_lin', M_i_lin, ...
    't_i_nl',  r_imp.t,  'M_i_nl',  M_i_nl, ...
    't_m_lin', t_m_lin,  'M_m_lin', M_m_lin, ...
    't_m_nl',  r_ms.t,   'M_m_nl',  M_m_nl);

plot_moments(out);
end

% =====================================================================

function P = load_params()
here = fileparts(mfilename('fullpath'));
f = fullfile(here, 'wing_nl.mat');
if ~isfile(f)
    error('moment_response:model', 'wing_nl.mat not found -- run build_ss_nl.m first.');
end
P = load(f);
P.file = f;
end

% ---------------------------------------------------------------------

function [t, h, th] = linear_case(P, input_idx, t_end, fs)
% Build the same input signal as nonlinear_response/to_forcing, then lsim.
t   = (0:1/fs:t_end)';
mag = 1.0;
t0  = 0.05;

switch input_idx
    case 1                                      % impulse: short rectangular pulse
        w  = 5/fs;
        ug = mag * double(t >= t0 & t < t0 + w);

    case 2                                      % multisine (Schroeder phases)
        f0  = 0.25;  K = 80;  kk = 1:K;
        amp = 1.5 / sqrt(K/2);
        phi = -pi * kk .* (kk-1) / K;
        ug  = sum(amp .* sin(2*pi*f0*(t*kk) + phi), 2);
end

sys = ss(P.A, P.B, P.C_out, P.D_out);
y   = lsim(sys, ug, t);
h   = y(:,1);
th  = y(:,2);
end

% ---------------------------------------------------------------------

function M = equivalent_moment(h, th, P, L, nl)
d1 = h + P.r1*th;
d2 = h + P.r2*th;
if nl
    F1 = P.k1*d1 + P.knl1*d1.^3;
    F2 = P.k2*d2 + P.knl2*d2.^3;
else
    F1 = P.k1*d1;
    F2 = P.k2*d2;
end
Mb = L*(F1 + F2);
T  = P.r1*F1 + P.r2*F2;
M  = 0.5*( abs(Mb) + sqrt(Mb.^2 + T.^2) );
end

% ---------------------------------------------------------------------

function plot_moments(out)
col_lin  = [0.20 0.40 0.80];   % blue   — constant k
col_nl   = [0.85 0.33 0.10];   % orange — nonlinear k
col_diff = [0.30 0.60 0.30];   % green  — absolute difference
labels   = {'Constant k', 'Nonlinear k', '|Difference|'};

% Interpolate nonlinear (ode45 adaptive grid) onto the uniform linear grid
diff_i = abs(out.M_i_lin - interp1(out.t_i_nl, out.M_i_nl, out.t_i_lin, 'linear', 'extrap'));
diff_m = abs(out.M_m_lin - interp1(out.t_m_nl, out.M_m_nl, out.t_m_lin, 'linear', 'extrap'));

figure('Name','Equivalent bending moment','Position',[100 100 1100 420]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile
plot(out.t_i_lin, out.M_i_lin, 'Color',col_lin,  'LineWidth',1.5)
hold on
plot(out.t_i_nl,  out.M_i_nl,  'Color',col_nl,   'LineWidth',1.5)
plot(out.t_i_lin, diff_i,       'Color',col_diff, 'LineWidth',1.2, 'LineStyle','--')
yline(0,'Color',[.7 .7 .7],'HandleVisibility','off');  grid on
xlabel('t  [s]');  ylabel('M_{eq}  [N m]')
title(sprintf('Impulse response   (L = %.2f m)', out.L))
legend(labels, 'Location','northeast')

nexttile
plot(out.t_m_lin, out.M_m_lin, 'Color',col_lin,  'LineWidth',1.5)
hold on
plot(out.t_m_nl,  out.M_m_nl,  'Color',col_nl,   'LineWidth',1.5)
plot(out.t_m_lin, diff_m,       'Color',col_diff, 'LineWidth',1.2, 'LineStyle','--')
yline(0,'Color',[.7 .7 .7],'HandleVisibility','off');  grid on
xlabel('t  [s]');  ylabel('M_{eq}  [N m]')
title('Multisine response')
legend(labels, 'Location','northeast')

here = fileparts(mfilename('fullpath'));
exportgraphics(gcf, fullfile(here, 'nl_moment_response.png'), 'Resolution', 300);
fprintf('Saved: nl_moment_response.png\n');
end
