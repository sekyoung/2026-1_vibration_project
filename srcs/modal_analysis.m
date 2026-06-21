% modal_analysis.m  —  Damping ratios and natural frequencies from wing_2dof.mat
%
% Two sets of results are computed:
%   1. Structural modes  : eigenvalues of  M_mat \ Ks  (no aerodynamics)
%   2. Aeroelastic modes : eigenvalues of  A           (aero included)
%
% Eigenvalue types:
%   Complex pair  ->  wn = |lambda|,  zeta = -Re(lambda)/wn,  fn = wn/(2*pi)
%   Real negative ->  overdamped / stable pole,   reported as |lambda| [rad/s]
%   Real positive ->  divergence / unstable pole, reported as |lambda| [rad/s]

here = fileparts(mfilename('fullpath'));
load(fullfile(here, 'wing_2dof.mat'), ...
     'A', 'M_mat', 'Keff', 'Ctotal', ...
     'U', 'rho', 'chord', 'span', 'Clalpha', 'ea');

%% ── Reconstruct structural matrices ─────────────────────────────────────────
qS = 0.5 * rho * U^2 * chord * span;
Ka = [0,  qS*Clalpha;
      0, -qS*Clalpha*ea];
Ca = (qS*Clalpha/U) .* [1,   0;
                         -ea, 0];

Ks       = Keff - Ka;
C_struct = Ctotal - Ca;
Minv     = inv(M_mat);

%% ── Build structural A-matrix and get all eigenvalues ───────────────────────
A_struct = [zeros(2), eye(2);
            -Minv*Ks, -Minv*C_struct];

lam_s  = eig(A_struct);
lam_ae = eig(A);

%% ── Helper: classify and extract modal parameters ───────────────────────────
[s_wn,  s_zeta,  s_fn,  s_type,  lam_s_plot]  = classify_eigs(lam_s);
[ae_wn, ae_zeta, ae_fn, ae_type, lam_ae_plot] = classify_eigs(lam_ae);

%% ── Print structural results ─────────────────────────────────────────────────
fprintf('\n=== Structural Modes (no aerodynamics) ===\n');
fprintf('  %-6s  %-14s  %-14s  %-12s  %s\n', ...
        'Mode', 'fn [Hz]', 'wn [rad/s]', 'zeta [-]', 'Type');
for i = 1:numel(s_wn)
    if isnan(s_zeta(i))
        fprintf('  %-6d  %14.4f  %14.4f  %12s  %s\n', ...
            i, s_fn(i), s_wn(i), 'N/A', s_type{i});
    else
        fprintf('  %-6d  %14.4f  %14.4f  %12.6f  %s\n', ...
            i, s_fn(i), s_wn(i), s_zeta(i), s_type{i});
    end
end

%% ── Print aeroelastic results ────────────────────────────────────────────────
fprintf('\n=== Aeroelastic Modes (U = %.1f m/s) ===\n', U);
fprintf('  %-6s  %-14s  %-14s  %-12s  %s\n', ...
        'Mode', 'fn [Hz]', 'wn [rad/s]', 'zeta [-]', 'Type');
for i = 1:numel(ae_wn)
    if isnan(ae_zeta(i))
        fprintf('  %-6d  %14.4f  %14.4f  %12s  %s\n', ...
            i, ae_fn(i), ae_wn(i), 'N/A', ae_type{i});
    else
        fprintf('  %-6d  %14.4f  %14.4f  %12.6f  %s\n', ...
            i, ae_fn(i), ae_wn(i), ae_zeta(i), ae_type{i});
    end
end
fprintf('\n');

%% ── Store results ────────────────────────────────────────────────────────────
results.structural.fn_Hz   = s_fn;
results.structural.wn_rads = s_wn;
results.structural.zeta    = s_zeta;
results.structural.type    = s_type;

results.aeroelastic.fn_Hz   = ae_fn;
results.aeroelastic.wn_rads = ae_wn;
results.aeroelastic.zeta    = ae_zeta;
results.aeroelastic.type    = ae_type;

%% ── Pole map ─────────────────────────────────────────────────────────────────
figure('Name','Pole Map','Position',[100 100 620 460]);
hold on; grid on;

% structural: plot both upper and mirrored lower half
s_conj  = [lam_s_plot;  conj(lam_s_plot(imag(lam_s_plot) ~= 0))];
ae_conj = [lam_ae_plot; conj(lam_ae_plot(imag(lam_ae_plot) ~= 0))];

scatter(real(s_conj),  imag(s_conj),  70, 'r', '^', 'LineWidth', 1.5, ...
        'DisplayName', 'Structural');
scatter(real(ae_conj), imag(ae_conj), 70, 'b', 'o', 'filled', ...
        'DisplayName', 'Aeroelastic');

xline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(0, 'k:',  'LineWidth', 0.5, 'HandleVisibility', 'off');
xlabel('Real  [1/s]');
ylabel('Imag  [rad/s]');
title(sprintf('Pole Map  (U = %.1f m/s)', U));
legend('Location', 'best');

%% ── Local function ───────────────────────────────────────────────────────────
function [wn, zeta, fn, type, lam_keep] = classify_eigs(lam)
% Separates complex (underdamped) and real (overdamped/unstable) eigenvalues.
% Returns one entry per unique mode.
    tol = max(abs(lam)) * 1e-8;

    is_complex = abs(imag(lam)) > tol;
    lam_c = lam(is_complex & imag(lam) > 0);   % upper half-plane of conjugate pairs
    lam_r = real(lam(~is_complex));             % real poles

    % sort: underdamped by ascending frequency, real by ascending value
    [~, idx] = sort(abs(lam_c));
    lam_c = lam_c(idx);
    lam_r = sort(lam_r);

    wn_c   = abs(lam_c);
    zeta_c = -real(lam_c) ./ wn_c;
    fn_c   = wn_c / (2*pi);
    type_c = cell(numel(lam_c), 1);
    for k = 1:numel(lam_c)
        if zeta_c(k) >= 0
            type_c{k} = 'underdamped (stable)';
        else
            type_c{k} = 'underdamped (unstable/flutter)';
        end
    end

    wn_r   = abs(lam_r);
    fn_r   = wn_r / (2*pi);
    zeta_r = nan(numel(lam_r), 1);
    type_r = cell(numel(lam_r), 1);
    for k = 1:numel(lam_r)
        if lam_r(k) < 0
            type_r{k} = 'real (overdamped/stable)';
        else
            type_r{k} = 'real (divergent/unstable)';
        end
    end

    wn   = [wn_c(:);   wn_r(:)];
    zeta = [zeta_c(:); zeta_r(:)];
    fn   = [fn_c(:);   fn_r(:)];
    type = [type_c;    type_r];
    lam_keep = [lam_c(:); lam_r(:)];
end
