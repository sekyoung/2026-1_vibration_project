% build_ss_nl.m
% 2-DOF aeroelastic wing (plunge h, pitch theta) with EXPLICIT per-spring stiffness.
%
% New vs build_ss.m:
%   - each spring (k1,k2) and damper (c1,c2) defined separately and assembled
%     into the matrices from their moment arms r1,r2 (= +/- d/2).
%   - optional cubic stiffness per spring (knl1,knl2) for nonlinear studies.
%   - exports BOTH the q=0 linearisation (A,B,C,D) and the pieces the nonlinear
%     solver needs (M, C, Ka, Fgust, spring params) to wing_nl.mat.
%
% Original build_ss.m and wing_2dof.mat are left untouched.

% ---------- physical parameters ----------------------------------
m   = 2.0;     % section mass                        [kg]
Ip  = 0.01;    % pitch inertia about reference       [kg*m^2]
S   = 0.015;   % static unbalance = m*x_cg           [kg*m]
d   = 0.25;    % distance between the two springs    [m]

k1  = 1200.0;  % spring 1 stiffness                  [N/m]
k2  = 1200.0;  % spring 2 stiffness                  [N/m]   (set != k1 for h-theta coupling)
knl1 = 0.0;    % spring 1 cubic coefficient          [N/m^3] (set > 0 for hardening)
knl2 = 0.0;    % spring 2 cubic coefficient          [N/m^3]

c1  = 3.0;     % damper 1 coefficient                [N*s/m]
c2  = 3.0;     % damper 2 coefficient                [N*s/m]

rho = 1.225;   % air density                         [kg/m^3]
U   = 12.0;    % freestream speed                    [m/s]
chord = 0.30;  % airfoil chord                       [m]
span  = 0.50;  % section span                        [m]
Clalpha = 2*pi;% lift-curve slope                    [1/rad]
ea  = 0.05;    % aero centre forward of reference    [m]

% ---------- spring/damper geometry -------------------------------
% A point at moment arm r deflects by  delta = h + r*theta.
r1 = +d/2;     % spring/damper 1 arm from reference  [m]
r2 = -d/2;     % spring/damper 2 arm from reference  [m]

% ---------- structural matrices from per-element contributions ----
% K_hh = sum k_i ,  K_h.theta = sum k_i r_i ,  K_theta.theta = sum k_i r_i^2
M_mat = [m, S;
         S, Ip];

Ks = [ k1     + k2     ,  k1*r1     + k2*r2     ;
       k1*r1  + k2*r2  ,  k1*r1^2   + k2*r2^2   ];   % off-diagonal = (k1-k2)*d/2

C_struct = [ c1     + c2     ,  c1*r1     + c2*r2     ;
             c1*r1  + c2*r2  ,  c1*r1^2   + c2*r2^2   ];

% ---------- aerodynamic terms (unchanged) ------------------------
qS    = 0.5 * rho * U^2 * chord * span;
Ka    = [0,  qS*Clalpha;
         0, -qS*Clalpha*ea];
Ca    = (qS*Clalpha/U) .* [1,   0;
                          -ea,  0];
Fgust = (qS*Clalpha/U) .* [-1; ea];

C_mat = C_struct + Ca;          % total damping used by the solver

% ---------- q = 0 linearisation (tangent stiffness = Ks) ---------
Minv  = inv(M_mat);
A     = [zeros(2),        eye(2);
        -Minv*(Ks+Ka),   -Minv*C_mat];
B     = [zeros(2,1);
         Minv*Fgust];
C_out = [eye(2), zeros(2)];
D_out = zeros(2,1);

fn_struct = sqrt(eig(M_mat \ Ks)) / (2*pi);

% ---------- report -----------------------------------------------
fprintf('\nStructural stiffness matrix Ks =\n'); disp(Ks);
fprintf('Coupling term K_h,theta = %.4f  (0 when k1 = k2)\n', Ks(1,2));
fprintf('Structural natural frequencies: %.3f Hz, %.3f Hz\n', fn_struct(1), fn_struct(2));

% ---------- export -----------------------------------------------
save('wing_nl.mat', ...
     'A','B','C_out','D_out','fn_struct', ...          % linearisation (for linear tools)
     'M_mat','C_mat','Ka','Fgust', ...                 % nonlinear EOM operators
     'r1','r2','k1','k2','knl1','knl2', ...             % per-spring params (nonlinear restoring)
     'U');
fprintf('Saved: wing_nl.mat\n');
