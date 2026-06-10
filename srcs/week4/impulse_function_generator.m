function gust = impulse_function_generator(varargin)
% impulse_function_generator  Certification discrete tuned (1-cosine) gust input.
%   Generates the FAA 14 CFR 25.341(a) / EASA CS-25.341(a) discrete gust as a
%   gust-velocity signal w_g(t) [m/s] for the Week-3 2-DOF wing model, sweeps the
%   gust gradient H to find the critical response, and plots the result.
%   Rationale, formulae, parameters and references: impulse_function_generator.md
%
%   gust = impulse_function_generator('standard','scaled')   % default
%   gust = impulse_function_generator('standard','cs25')     % literal CS-25 range

cfg = default_config();
for i = 1:2:numel(varargin)
    cfg.(varargin{i}) = varargin{i+1};
end

model         = load_model(cfg.model_file);
[H, f_gust]   = gradient_sweep(cfg, model);
Uds           = design_gust_velocity(cfg, H, model.U);
[t, wg]       = build_gust_family(H, Uds, model.U, cfg.fs, cfg.settle);
resp          = evaluate_response(model, t, wg);
crit          = find_critical(H, resp);
gust          = pack_output(t, wg, H, Uds, f_gust, resp, crit, cfg);

if cfg.plot
    make_plots(gust, model);
end
end

% =====================================================================

function cfg = default_config()
cfg.standard   = 'scaled';
cfg.model_file = '';
cfg.fs         = 2000;
cfg.nH         = 25;
cfg.settle     = 6;
cfg.Fg         = 1.0;
cfg.Uref       = [];
cfg.band       = [0.5 2.0];
cfg.plot       = true;
end

% ---------------------------------------------------------------------

function model = load_model(model_file)
if isempty(model_file)
    here = fileparts(mfilename('fullpath'));
    cand = {fullfile(here,'wing_2dof.mat'), ...
            fullfile(here,'..','week3','wing_2dof.mat')};
    for c = cand
        if isfile(c{1}), model_file = c{1}; break; end
    end
    if isempty(model_file)
        error('impulse_function_generator:model', ...
              'wing_2dof.mat not found next to script or in ../week3.');
    end
end
model = load(model_file, 'A','B','C_out','D_out','U','fn_struct');
end

% ---------------------------------------------------------------------

function [H, f_gust] = gradient_sweep(cfg, model)
switch lower(cfg.standard)
    case 'cs25'
        H = logspace(log10(9.1), log10(107), cfg.nH);
    case 'scaled'
        fn   = model.fn_struct(:)';
        f_lo = cfg.band(1) * min(fn);
        f_hi = cfg.band(2) * max(fn);
        H    = logspace(log10(model.U/(2*f_hi)), log10(model.U/(2*f_lo)), cfg.nH);
    otherwise
        error('impulse_function_generator:standard', ...
              'cfg.standard must be ''scaled'' or ''cs25''.');
end
f_gust = model.U ./ (2*H);
end

% ---------------------------------------------------------------------

function Uds = design_gust_velocity(cfg, H, U)
Uref = cfg.Uref;
if isempty(Uref)
    switch lower(cfg.standard)
        case 'cs25',   Uref = 17.07;
        case 'scaled', Uref = 0.10 * U;
    end
end
Uds = Uref * cfg.Fg .* (H/107).^(1/6);
end

% ---------------------------------------------------------------------

function [t, wg] = build_gust_family(H, Uds, U, fs, settle)
Tg   = 2*H/U;
tmax = max(Tg) + settle;
t    = (0:1/fs:tmax)';
wg   = zeros(numel(t), numel(H));
for j = 1:numel(H)
    in       = t <= Tg(j);
    wg(in,j) = (Uds(j)/2) * (1 - cos(pi*U*t(in)/H(j)));
end
end

% ---------------------------------------------------------------------

function resp = evaluate_response(model, t, wg)
sys          = ss(model.A, model.B, model.C_out, model.D_out);
nH           = size(wg,2);
resp.h_peak  = zeros(nH,1);
resp.th_peak = zeros(nH,1);
for j = 1:nH
    y               = lsim(sys, wg(:,j), t);
    resp.h_peak(j)  = max(abs(y(:,1))) * 1000;
    resp.th_peak(j) = max(abs(y(:,2))) * 180/pi;
end
end

% ---------------------------------------------------------------------

function crit = find_critical(H, resp)
[crit.h_max,  ih] = max(resp.h_peak);
[crit.th_max, it] = max(resp.th_peak);
crit.H_h   = H(ih);  crit.idx_h  = ih;
crit.H_th  = H(it);  crit.idx_th = it;
end

% ---------------------------------------------------------------------

function gust = pack_output(t, wg, H, Uds, f_gust, resp, crit, cfg)
gust.t      = t;
gust.wg     = wg;
gust.H      = H(:);
gust.Uds    = Uds(:);
gust.f_gust = f_gust(:);
gust.resp   = resp;
gust.crit   = crit;
gust.cfg    = cfg;
end

% ---------------------------------------------------------------------

function make_plots(gust, model)
col_h  = [0.85 0.33 0.10];
col_th = [0.20 0.40 0.80];
sys    = ss(model.A, model.B, model.C_out, model.D_out);

figure('Name','Discrete Tuned Gust (25.341)','Position',[100 100 820 780]);

subplot(3,1,1)
show = unique(round(linspace(1, numel(gust.H), 5)));
hold on
for j = show
    plot(gust.t, gust.wg(:,j), 'LineWidth',1.1, ...
        'DisplayName',sprintf('H = %.2f m', gust.H(j)));
end
plot(gust.t, gust.wg(:,gust.crit.idx_h), 'k', 'LineWidth',2, ...
    'DisplayName','critical (h)');
xlabel('t  [s]');  ylabel('w_g  [m/s]')
title('1-cosine discrete gust family   $U(s) = (U_{ds}/2)\cdot[1-\cos(\pi s/H)]$','Interpreter','latex');
legend('Location','northeast');  grid on;  xlim([0 max(gust.t)])

subplot(3,1,2)
yyaxis left
plot(gust.H, gust.resp.h_peak, '-o', 'Color',col_h, 'MarkerSize',3);  ylabel('peak |h|  [mm]')
yyaxis right
plot(gust.H, gust.resp.th_peak, '-s', 'Color',col_th, 'MarkerSize',3); ylabel('peak |\theta|  [deg]')
xlabel('gust gradient  H  [m]')
title(sprintf('Tuning sweep -- critical  $H_h$ = %.2f m,   $H_{\\theta}$ = %.2f m', ...
    gust.crit.H_h, gust.crit.H_th),'Interpreter','latex');
xline(gust.crit.H_h, '--','Color',col_h, 'HandleVisibility','off');
xline(gust.crit.H_th,'--','Color',col_th,'HandleVisibility','off');
grid on

subplot(3,1,3)
y = lsim(sys, gust.wg(:,gust.crit.idx_h), gust.t);
yyaxis left
plot(gust.t, y(:,1)*1000, 'Color',col_h, 'LineWidth',1.5);  ylabel('h  [mm]')
yyaxis right
plot(gust.t, rad2deg(y(:,2)), 'Color',col_th, 'LineWidth',1.5); ylabel('\theta  [deg]')
xlabel('t  [s]')
title('Response to the critical gust');  grid on
end
