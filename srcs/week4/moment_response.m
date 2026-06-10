function out = moment_response(L)
% moment_response  Equivalent root bending moment under the critical gust.
%
%   out = moment_response(L)
%       L : wing length [m]  (default 0.50, explicit)
%
%   Calls impulse_function_generator (constant-k, linear model) and
%   nonlinear_response (cubic springs, varying tangent k) to obtain h(t)
%   and theta(t), then evaluates the equivalent bending moment
%
%       M_eq(t) = 0.5*( |M_b| + sqrt(M_b^2 + T^2) )
%       M_b     = L*(F1 + F2),    T = r1*F1 + r2*F2
%       Fi      = ki*di (+ knl_i*di^3),   di = h + ri*theta
%
%   Local spring-attachment effects neglected: section loads transfer
%   rigidly to the root. Derivation and assumptions: moment_response.md

if nargin < 1, L = 0.50; end                 % wing length [m], explicit

P = load_params();

[t_lin, h_lin, th_lin] = linear_case(P);
[t_nl,  h_nl,  th_nl ] = nonlinear_case();

M_lin = equivalent_moment(h_lin, th_lin, P, L, false);
M_nl  = equivalent_moment(h_nl,  th_nl,  P, L, true);

out = struct('L',L, 't_lin',t_lin, 'M_lin',M_lin, ...
                    't_nl',t_nl,   'M_nl',M_nl);
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

function [t, h, th] = linear_case(P)
g   = impulse_function_generator('model_file',P.file, 'plot',false);
wg  = g.wg(:, g.crit.idx_h);
y   = lsim(ss(P.A, P.B, P.C_out, P.D_out), wg, g.t);
t   = g.t;  h = y(:,1);  th = y(:,2);
end

% ---------------------------------------------------------------------

function [t, h, th] = nonlinear_case()
r = nonlinear_response(4);
close(gcf);                                  % drop its internal figure
t = r.t;  h = r.h;  th = r.theta;
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
Mb = L*(F1 + F2);                            % bending  (from h)
T  = P.r1*F1 + P.r2*F2;                      % torsion  (from theta)
M  = 0.5*( abs(Mb) + sqrt(Mb.^2 + T.^2) );   % max normal stress equivalent
end

% ---------------------------------------------------------------------

function plot_moments(out)
col = [0.85 0.33 0.10];

figure('Name','Equivalent bending moment','Position',[100 100 980 360]);

subplot(1,2,1)
plot(out.t_lin, out.M_lin, 'Color',col, 'LineWidth',1.5)
yline(0,'Color',[.7 .7 .7]);  grid on
xlabel('t  [s]');  ylabel('M_{eq}  [N m]')
title(sprintf('Constant k   (L = %.2f m)', out.L))

subplot(1,2,2)
plot(out.t_nl, out.M_nl, 'Color',col, 'LineWidth',1.5)
yline(0,'Color',[.7 .7 .7]);  grid on
xlabel('t  [s]');  ylabel('M_{eq}  [N m]')
title('Varying k   (cubic springs)')
end
