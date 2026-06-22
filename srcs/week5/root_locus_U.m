function root_locus_U()
% Root locus of the 2-DOF aeroelastic wing as aero.U varies.
%
% Sweeps freestream speed U from near-zero to U_max, builds the linear
% A-matrix at each U (same parameters as build_2dof in build_models.m),
% tracks the four eigenvalue branches, locates the flutter speed where
% a branch crosses the imaginary axis (Re = 0, Im ≠ 0), and saves a plot.

% ---- physical parameters (must match build_2dof / aero_terms) ----------
m  = 1.5;     Ip = 0.0084;    S  = 0.0225;
k  = 900.0;  cd = 1.0;     d  = 0.22;

rho     = 1.225;   chord   = 0.30;
span    = 0.50;    Clalpha = 6;
ea      = 0.045;
U_nom   = 12.0;    % nominal speed from aero_terms [m/s]

M_mat    = [m,   S;   S,  Ip];
Ks       = [2*k, 0;   0,  k*d^2/2];
C_struct = [2*cd, 0;  0,  cd*d^2/2];

% ---- velocity sweep ------------------------------------------------
U_vec = linspace(0.5, 50, 4000);
N     = numel(U_vec);

ev_mat = zeros(4, N);   % each column: eigenvalues sorted by imag part

for i = 1:N
    U  = U_vec(i);
    qS = 0.5 * rho * U^2 * chord * span;

    Ka   = [0,  qS*Clalpha;   0, -qS*Clalpha*ea];
    Ca   = (qS*Clalpha / U) .* [1, 0;  -ea, 0];

    Minv = inv(M_mat);
    A    = [zeros(2),        eye(2);
            -Minv*(Ks+Ka),  -Minv*(C_struct+Ca)];

    ev = eig(A);
    [~, idx] = sort(imag(ev));   % ascending imag → consistent branch order
    ev_mat(:, i) = ev(idx);
end

% ---- find flutter speed: first U where max Re crosses 0 from below ----
max_re = max(real(ev_mat), [], 1);
cross  = find(diff(sign(max_re)) > 0, 1);

U_fl  = NaN;
im_fl = NaN;
if ~isempty(cross)
    % linear interpolation between the two straddle points
    U_fl = interp1(max_re(cross:cross+1), U_vec(cross:cross+1), 0);

    % find which branch is the one crossing
    for b = 1:4
        re_b = real(ev_mat(b, :));
        if re_b(cross) < 0 && re_b(cross+1) >= 0
            im_fl = interp1(U_vec(cross:cross+1), ...
                            imag(ev_mat(b, cross:cross+1)), U_fl);
            break
        end
    end

    fprintf('Flutter speed:     U_flutter = %.4f m/s\n', U_fl);
    fprintf('Flutter frequency: %.4f rad/s  (%.4f Hz)\n', ...
            abs(im_fl), abs(im_fl)/(2*pi));
else
    fprintf('No flutter detected in U = [%.1f, %.1f] m/s.\n', ...
            U_vec(1), U_vec(end));
end

% ---- plot --------------------------------------------------------------
figure('Name', 'Root Locus vs aero.U', 'Position', [80 80 950 680]);
hold on;  grid on;  ax = gca;

colors = lines(4);
branch_labels = {'Mode 1 (lower \omega)','Mode 2 (upper \omega)',...
                 'Mode 2 (conj.)','Mode 1 (conj.)'};

for b = 1:4
   
plot(real(ev_mat(b,:)), imag(ev_mat(b,:)), '-', ...
         'Color', colors(b,:), 'LineWidth', 1.8, ...
         'DisplayName', branch_labels{b});
      %{        
    plot(real(ev_mat(b,:)), imag(ev_mat(b,:)), '-', ...
         'Color', colors(b,:), 'LineWidth', 1.8);
      %}
    % circle at low-U end
    plot(real(ev_mat(b,1)), imag(ev_mat(b,1)), 'o', ...
         'Color', colors(b,:), 'MarkerFaceColor', colors(b,:), ...
         'MarkerSize', 6, 'HandleVisibility', 'off');

    % arrow mid-branch to show direction of increasing U
    mid = round(N/2);
    quiver(real(ev_mat(b, mid)), imag(ev_mat(b, mid)), ...
           real(ev_mat(b, mid+1) - ev_mat(b, mid)), ...
           imag(ev_mat(b, mid+1) - ev_mat(b, mid)), ...
           5, 'Color', colors(b,:), 'MaxHeadSize', 3, ...
           'HandleVisibility', 'off');
end

% imaginary axis
xline(0, 'k--', 'LineWidth', 1.4, 'HandleVisibility', 'off');
text(0.5, ax.YLim(1)*0.92, 'Im axis', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);

%{
% nominal operating point (U = 12 m/s)
nom_idx = round(interp1(U_vec, 1:N, U_nom));
for b = 1:4
    plot(real(ev_mat(b, nom_idx)), imag(ev_mat(b, nom_idx)), 'ks', ...
         'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility', 'off');
end
% one legend entry for the nominal marker
plot(nan, nan, 'ks', 'MarkerSize', 8, 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Nominal  U = %.0f m/s', U_nom));
%}

% flutter marker
if ~isnan(U_fl)
    plot(0,  im_fl, 'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'r', ...
         'LineWidth', 1.2, ...
         'DisplayName', sprintf('Flutter  U_{fl} = %.2f m/s', U_fl));
    plot(0, -im_fl, 'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'r', ...
         'HandleVisibility', 'off');
    text(0.4,  im_fl, sprintf('  U_{fl} = %.2f m/s', U_fl), ...
         'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');
    text(0.4, -im_fl, sprintf('  U_{fl} = %.2f m/s', U_fl), ...
         'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');
end

xlabel('Real part  \sigma  [1/s]',          'FontSize', 12);
ylabel('Imaginary part  \omega_d  [rad/s]', 'FontSize', 12);
title({'Root Locus w.r.t. aero.U  --  2-DOF aeroelastic wing'}, ...
      'FontSize', 12);
legend('Location', 'best', 'FontSize', 9);

here = fileparts(mfilename('fullpath'));
out  = fullfile(here, 'root_locus_U.png');
exportgraphics(gcf, out, 'Resolution', 300);
fprintf('Figure saved: %s\n', out);
end
