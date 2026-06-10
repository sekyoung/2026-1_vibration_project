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
resp.h_all   = zeros(numel(t), nH);   % plunge histories, all gusts [mm]
resp.th_all  = zeros(numel(t), nH);   % pitch  histories, all gusts [deg]
resp.h_peak  = zeros(nH,1);
resp.th_peak = zeros(nH,1);
for j = 1:nH
    y                = lsim(sys, wg(:,j), t);
    resp.h_all(:,j)  = y(:,1) * 1000;
    resp.th_all(:,j) = y(:,2) * 180/pi;
    resp.h_peak(j)   = max(abs(resp.h_all(:,j)));
    resp.th_peak(j)  = max(abs(resp.th_all(:,j)));
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

function make_plots(gust, model) %#ok<INUSD>
nH   = numel(gust.H);
show = unique(round(linspace(1, nH, 5)));   % representative gradients (same set in all panels)
pal  = lines(numel(show));                  % one colour per shown gust, shared across panels

figure('Name','Discrete Tuned Gust (25.341)','Position',[100 100 820 800]);

ax1 = subplot(3,1,1);
draw_family(ax1, gust.t, gust.wg,         show, pal, gust.H, gust.crit.idx_h, ...
    'w_g  [m/s]', '1-cosine discrete gust family   (critical in bold)', true);

ax2 = subplot(3,1,2);
draw_family(ax2, gust.t, gust.resp.h_all,  show, pal, gust.H, gust.crit.idx_h, ...
    'h  [mm]', sprintf('Plunge response   (critical H_h = %.2f m)', gust.crit.H_h), false);

ax3 = subplot(3,1,3);
draw_family(ax3, gust.t, gust.resp.th_all, show, pal, gust.H, gust.crit.idx_th, ...
    '\theta  [deg]', sprintf('Pitch response   (critical H_\\theta = %.2f m)', gust.crit.H_th), false);
xlabel(ax3, 't  [s]');
end

% ---------------------------------------------------------------------

function draw_family(ax, t, Y, show, pal, Hvec, idx_crit, ylab, ttl, do_legend)
% Same representative gusts and colours in every panel; the critical one bold black.
hold(ax, 'on');
for i = 1:numel(show)
    plot(ax, t, Y(:,show(i)), 'Color',pal(i,:), 'LineWidth',1.0, ...
        'DisplayName',sprintf('H = %.2f m', Hvec(show(i))));
end
plot(ax, t, Y(:,idx_crit), 'k', 'LineWidth',2.2, 'DisplayName','critical');
yline(ax, 0, 'Color',[.85 .85 .85], 'HandleVisibility','off');
ylabel(ax, ylab);  grid(ax, 'on');  xlim(ax, [0 max(t)]);  title(ax, ttl);
if do_legend, legend(ax, 'Location','northeast'); end
end