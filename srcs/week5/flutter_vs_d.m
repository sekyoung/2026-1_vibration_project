function flutter_vs_d()
% Flutter speed and stability surface as a function of CG-to-AC distance d.
%
% d = x_CG - x_AC  [m]:  + when CG is aft of AC (destabilising)
%                         - when CG is forward of AC  (stabilising)
%
% NACA 2412 context: AC at ~25% chord; user specifies CG also at 25% chord
% => d_nom = 0, which maps to S = m*(0 - ea) = -m*ea = -0.10 kg.m.
%
% Relation to build_2dof parameters (build_models.m):
%   S     = m * (d - ea)        static unbalance [kg.m]
%   M_mat = [m, S; S, Ip]       mass matrix
%
% Physical constraint for PD mass matrix:
%   m*Ip > S^2  =>  |d - ea| < sqrt(Ip/m) = r_gyr
%
% Outputs:
%   Figure 1  2-D line:  U_flutter vs d
%   Figure 2  3-D surf:  max Re(lambda) over (d, U) space
%             Red curve on surface = flutter boundary (Re = 0)

% ---- fixed parameters (match build_2dof / aero_terms in build_models.m) --
m  = 2.0;   Ip = 0.01;
k  = 1200;  cd = 3.0;  d_sp = 0.25;   % d_sp = spring spacing (not CG-AC d!)

rho = 1.225;  chord = 0.30;  span = 0.50;
Clalpha = 2*pi;  ea = 0.05;            % AC forward of elastic axis [m]

% Structural matrices (independent of d or U)
Ks       = [2*k,  0;   0, k*d_sp^2/2];
C_struct = [2*cd, 0;   0, cd*d_sp^2/2];

% ---- parametric range ---------------------------------------------------
% PD constraint: |d - ea| < sqrt(Ip/m)
r_gyr = sqrt(Ip / m);            % = 0.0707 m
d_min =  ea - r_gyr + 5e-4;     % m  (add small margin)
d_max =  ea + r_gyr - 5e-4;

d_vec = linspace(d_min, d_max, 200);   % CG-to-AC distance [m]
U_vec = linspace(0.5,   50,   500);    % freestream speed  [m/s]

Nd = numel(d_vec);
Nu = numel(U_vec);

max_re_surf = nan(Nu, Nd);   % max Re(eig(A)) at each (U, d) point
U_fl        = nan(1,  Nd);   % flutter speed for each d

for id = 1:Nd
    d_cga = d_vec(id);
    S     = m * (d_cga - ea);    % static unbalance
    M     = [m, S;  S, Ip];

    if det(M) <= 0; continue; end
    Minv = inv(M);

    for iu = 1:Nu
        U  = U_vec(iu);
        qS = 0.5 * rho * U^2 * chord * span;

        Ka = [0,  qS*Clalpha;   0, -qS*Clalpha*ea];
        Ca = (qS*Clalpha / U) .* [1, 0;  -ea, 0];

        A  = [zeros(2),          eye(2);
              -Minv*(Ks + Ka),  -Minv*(C_struct + Ca)];

        max_re_surf(iu, id) = max(real(eig(A)));
    end

    % flutter crossing: first U where max Re crosses 0 from below
    col = max_re_surf(:, id)';
    cx  = find(diff(sign(col)) > 0, 1);
    if ~isempty(cx)
        U_fl(id) = interp1(col(cx:cx+1), U_vec(cx:cx+1), 0);
    end
end

d_cm = d_vec * 100;   % display in cm

% print summary
fprintf('CG-to-AC distance range: [%.2f, %.2f] cm\n', d_cm(1), d_cm(end));
fprintf('d = 0 (CG at AC, NACA2412 nom.): U_flutter = %.3f m/s\n', ...
        interp1(d_cm, U_fl, 0, 'linear', NaN));
[min_Ufl, i_min] = min(U_fl);
fprintf('Lowest flutter:  U_fl = %.3f m/s at d = %.2f cm\n', min_Ufl, d_cm(i_min));

% =========================================================================
% Figure 1 — 2-D:  U_flutter vs d
% =========================================================================
f1 = figure('Name','Flutter speed vs d', 'Position',[60 80 760 470]);
plot(d_cm, U_fl, 'b-', 'LineWidth', 2.5);
hold on;  grid on;

% mark d = 0 (CG at AC)
xline(0, 'k--', 'LineWidth', 1.3, 'HandleVisibility','off');
text(0.3, min(U_fl, [], 'omitnan') * 0.93, ...
     '  d = 0  (CG at AC)', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);

% mark original model point (d from original S = 0.015)
S_orig   = 0.015;
d_orig   = S_orig/m + ea;   % = 0.0575 m = 5.75 cm
Ufl_orig = interp1(d_cm, U_fl, d_orig*100, 'linear', NaN);
if ~isnan(Ufl_orig)
    plot(d_orig*100, Ufl_orig, 'rs', 'MarkerSize', 10, 'MarkerFaceColor','r', ...
         'DisplayName', sprintf('build\\_2dof  d=%.1fcm, U_{fl}=%.1f m/s', ...
                                d_orig*100, Ufl_orig));
    legend('Location','best','FontSize',9);
end

xlabel('d = x_{CG} - x_{AC}  [cm]     (+ : CG aft of AC)', 'FontSize', 12);
ylabel('Flutter speed  U_{flutter}  [m/s]', 'FontSize', 12);
title({'Flutter speed vs CG-to-AC distance  d', ...
       'NACA 2412  --  CG at 25% chord when d = 0'}, 'FontSize', 12);

% =========================================================================
% Figure 2 — 3-D surf:  max Re(lambda) over (d, U)
% =========================================================================
[D_cm, UU] = meshgrid(d_cm, U_vec);

f2 = figure('Name','3D stability surface (d, U)', 'Position',[130 80 1000 680]);

% surface
surf(D_cm, UU, max_re_surf, 'EdgeColor','none', 'FaceAlpha', 0.88);
hold on;

% flutter boundary: Re = 0 contour on the surface (red)
contour3(D_cm, UU, max_re_surf, [0 0], 'r-', 'LineWidth', 3.5);

% transparent Re = 0 plane for reference
xlm = [min(d_cm), max(d_cm)];
ylm = [U_vec(1),  U_vec(end)];
patch([xlm(1) xlm(2) xlm(2) xlm(1)], ...
      [ylm(1) ylm(1) ylm(2) ylm(2)], ...
      [0 0 0 0], [0.6 0.6 0.6], ...
      'FaceAlpha', 0.12, 'EdgeColor','none');

% original model marker projected onto surface
if ~isnan(Ufl_orig) && d_orig*100 >= d_cm(1) && d_orig*100 <= d_cm(end)
    z_marker = interp2(D_cm, UU, max_re_surf, d_orig*100, Ufl_orig, 'linear', NaN);
    plot3(d_orig*100, Ufl_orig, 0, 'rs', 'MarkerSize', 12, ...
          'MarkerFaceColor','r', 'LineWidth', 1.5);
end

% colormap: blue (stable) -> white (neutral) -> red (unstable)
colormap(bwr_cmap(256));
c_lim = prctile(abs(max_re_surf(:)), 95);
clim([-c_lim, c_lim]);
cb = colorbar;
cb.Label.String  = 'max  Re(\lambda)  [1/s]';
cb.Label.FontSize = 11;

xlabel('d = x_{CG} - x_{AC}  [cm]', 'FontSize', 11);
ylabel('Freestream speed  U  [m/s]', 'FontSize', 11);
zlabel('max  Re(\lambda)  [1/s]',    'FontSize', 11);
title({'Stability surface over (CG-to-AC distance  d,  speed  U)', ...
       'Red curve = flutter boundary  [Re(\lambda) = 0]'}, 'FontSize', 12);
view([-42, 28]);

% ---- save ---------------------------------------------------------------
here = fileparts(mfilename('fullpath'));
saveas(f1, fullfile(here, 'flutter_vs_d_2D.png'));
saveas(f2, fullfile(here, 'flutter_vs_d_3D.png'));
fprintf('Saved: flutter_vs_d_2D.png  and  flutter_vs_d_3D.png\n');
end

% =========================================================================
% Local helper: diverging blue-white-red colormap
% =========================================================================
function cmap = bwr_cmap(n)
n2 = ceil(n / 2);
blue_to_white = [linspace(0,1,n2)', linspace(0,1,n2)', ones(n2,1)];
white_to_red  = [ones(n2,1), linspace(1,0,n2)', linspace(1,0,n2)'];
cmap = [blue_to_white; white_to_red];
cmap = cmap(1:n, :);
end
