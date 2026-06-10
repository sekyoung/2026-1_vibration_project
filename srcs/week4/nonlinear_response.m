function out = nonlinear_response(input_idx, h0, theta0, k1)
% nonlinear_response  Time response of the 2-DOF wing with nonlinear springs.
%
%   out = nonlinear_response(input_idx, h0, theta0, k1)
%       input_idx : 1 = impulse | 2 = multisine | 3 = unit step | 4 = gust generator
%       h0        : initial plunge   [m]    (default 0)
%       theta0    : initial pitch    [rad]  (default 0)
%       k1        : spring-1 stiffness [N/m] (default: value in wing_nl.mat)
%
%   Solves   M*q'' + C*q' + Ka*q + g(q) = Fgust*u(t)
%   with     g(q) the per-spring (linear + cubic) structural restoring force,
%            u(t) built from the chosen input by to_forcing (the "input_f_modi").
%
%   Plots a 4-row figure:  input u(t) | h(t) | theta(t) | k1(t) tangent stiffness.
%   When knl1 = knl2 = 0 the restoring is linear and k1(t) is flat at k1.

if nargin < 1, input_idx = 2;  end
if nargin < 2, h0     = 0;     end
if nargin < 3, theta0 = 0;     end
if nargin < 4, k1     = [];    end           % [] -> keep wing_nl.mat value

% ---- fixed run settings -----------------------------------------
t_end = 8;        % simulation length      [s]
fs    = 2000;     % sample rate for u(t)   [Hz]
mag   = 1.0;      % magnitude (step / impulse)
t0    = 0.05;     % onset time             [s]

P = load_plant(k1);                          % matrices + spring params
u = to_forcing(input_idx, t_end, fs, mag, t0); % input_f_modi: u(t) interpolant

% ---- nonlinear equation of motion -------------------------------
g   = @(q) spring_force(q, P);
rhs = @(t,x) [ x(3:4);
               P.M \ ( P.Fgust*u(t) - P.C*x(3:4) - P.Ka*x(1:2) - g(x(1:2)) ) ];

x0   = [h0; theta0; 0; 0];
opts = odeset('RelTol',1e-7, 'AbsTol',1e-9, 'MaxStep',5/fs);
[t, x] = ode45(rhs, [0 t_end], x0, opts);

out = pack(t, x, u, input_idx, P);
plot_response(out);
end

% =====================================================================

function P = load_plant(k1)
here = fileparts(mfilename('fullpath'));
f = fullfile(here, 'wing_nl.mat');
if ~isfile(f)
    error('nonlinear_response:model', 'wing_nl.mat not found -- run build_ss_nl.m first.');
end
S = load(f);
P.M = S.M_mat;  P.C = S.C_mat;  P.Ka = S.Ka;  P.Fgust = S.Fgust;
P.r1 = S.r1;    P.r2 = S.r2;
if isempty(k1), k1 = S.k1; end
P.k1 = k1;      P.k2 = S.k2;
P.knl1 = S.knl1; P.knl2 = S.knl2;
end

% ---------------------------------------------------------------------

function f = spring_force(q, P)
d1 = q(1) + P.r1*q(2);                        % spring-1 deflection
d2 = q(1) + P.r2*q(2);                        % spring-2 deflection
F1 = P.k1*d1 + P.knl1*d1^3;
F2 = P.k2*d2 + P.knl2*d2^3;
f  = [ F1 + F2;
       P.r1*F1 + P.r2*F2 ];                   % generalised force [h; theta]
end

% ---------------------------------------------------------------------

function u = to_forcing(choice, t_end, fs, mag, t0)
% Convert any input choice into ONE uniform object: u(t) interpolant.
tg = (0:1/fs:t_end)';

switch choice
    case 1                                    % impulse: short pulse ~ Dirac, peak = mag
        w  = 5/fs;
        ug = mag * (tg >= t0 & tg < t0 + w);

    case 2                                    % multisine (Schroeder phases), as Week 3
        f0 = 0.25;  K = 80;  kk = 1:K;
        amp = 1.5/sqrt(K/2);
        phi = -pi*kk.*(kk-1)/K;
        ug  = sum(amp .* sin(2*pi*f0*(tg*kk) + phi), 2);

    case 3                                    % unit step, arbitrary magnitude
        ug = mag * (tg >= t0);

    case 4                                    % certification gust (critical column)
        g   = impulse_function_generator('fs',fs, 'plot',false);
        sel = g.crit.idx_h;
        Fg  = griddedInterpolant(g.t, g.wg(:,sel), 'linear','nearest');
        ug  = Fg(tg);

    otherwise
        error('nonlinear_response:input', 'input_idx must be 1..4.');
end

u = griddedInterpolant(tg, ug, 'linear','nearest');
end

% ---------------------------------------------------------------------

function out = pack(t, x, u, choice, P)
names   = {'impulse','multisine','step','gust-generator'};
out.t     = t;
out.h     = x(:,1);
out.theta = x(:,2);
out.u     = u(t);
out.input = names{choice};

delta1    = out.h + P.r1*out.theta;           % spring-1 deflection history
out.k1eff = P.k1 + 3*P.knl1*delta1.^2;        % instantaneous tangent stiffness
out.P     = P;
end

% ---------------------------------------------------------------------

function plot_response(out)
col_u  = [0.40 0.40 0.40];
col_h  = [0.85 0.33 0.10];
col_th = [0.20 0.40 0.80];
col_k  = [0.30 0.60 0.30];

figure('Name','Nonlinear time response','Position',[100 100 780 860]);

subplot(4,1,1)
plot(out.t, out.u, 'Color',col_u, 'LineWidth',1.4)
ylabel('u(t)  [m/s]');  grid on
title(sprintf('Input: %s', out.input))

subplot(4,1,2)
plot(out.t, out.h*1000, 'Color',col_h, 'LineWidth',1.5)
yline(0,'Color',[.7 .7 .7]);  ylabel('h  [mm]');  grid on

subplot(4,1,3)
plot(out.t, rad2deg(out.theta), 'Color',col_th, 'LineWidth',1.5)
yline(0,'Color',[.7 .7 .7]);  ylabel('\theta  [deg]');  grid on

subplot(4,1,4)
plot(out.t, out.k1eff, 'Color',col_k, 'LineWidth',1.5)
ylabel('k_1(t)  [N/m]');  xlabel('t  [s]');  grid on
title('Spring-1 tangent stiffness  k_1 + 3 k_{nl1}\delta_1^2')
end
