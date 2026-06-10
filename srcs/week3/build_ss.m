% build_ss.m
% 2-DOF aeroelastic wing (plunge h, pitch theta)
%
% Equation of motion:
%   M*q'' + (C+Ca)*q' + (Ks+Ka)*q = Fgust * w_g
%
% State vector: x = [h; theta; hdot; thetadot]
% Input:        u = w_g  (gust velocity, m/s)
% Outputs:      y = [h; theta]
%
% Run this script once. It saves wing_2dof.mat which Simulink loads.

% ---------- physical parameters ----------------------------------
m     = 2.0;      % section mass                         [kg]
Ip    = 0.01;     % pitch inertia about reference point  [kg·m^2]
S     = 0.015;    % static unbalance = m*x_cg            [kg·m]
k     = 1200.0;   % each spring stiffness                [N/m]
cd    = 3.0;      % each damper coefficient              [N·s/m]
d     = 0.25;     % distance between springs             [m]

rho   = 1.225;    % air density                          [kg/m^3]
U     = 12.0;     % freestream speed                     [m/s]
chord = 0.30;     % airfoil chord                        [m]
span  = 0.50;     % section span                         [m]
Clalpha = 2*pi;   % lift-curve slope                     [1/rad]
ea    = 0.05;     % aero centre forward of reference     [m]

% ---------- system matrices --------------------------------------
M_mat = [m,  S;
         S,  Ip];

C_mat = [2*cd,       0;
         0,     cd*d^2/2];

Ks = [2*k,        0;
      0,      k*d^2/2];

qS = 0.5 * rho * U^2 * chord * span;

Ka = [0,       qS*Clalpha;
      0,  -qS*Clalpha*ea];

Ca = (qS*Clalpha/U) .* [1,    0;
                        -ea,   0];

Keff   = Ks + Ka;
Ctotal = C_mat + Ca;
Fgust  = (qS*Clalpha/U) .* [-1; ea];

% ---------- state-space ------------------------------------------
% x = [q; qdot],  qdotdot = -Minv*Keff*q - Minv*Ctotal*qdot + Minv*Fgust*u
Minv = inv(M_mat);

A = [zeros(2),    eye(2);
    -Minv*Keff,  -Minv*Ctotal];

B = [zeros(2,1);
     Minv*Fgust];

C_out = [eye(2), zeros(2)];  % observe [h, theta]
D_out = zeros(2,1);

% ---------- display for quick check ------------------------------
fprintf('\n=== State-space matrices ===\n');
fprintf('A =\n'); disp(A);
fprintf('B =\n'); disp(B);
fprintf('C =\n'); disp(C_out);
fprintf('D =\n'); disp(D_out);

% Natural frequencies (structural only, no aero)
Keff_struct = Ks;
eig_struct  = eig(M_mat \ Keff_struct);
fn_struct   = sqrt(eig_struct) / (2*pi);
fprintf('Structural natural frequencies: %.3f Hz,  %.3f Hz\n', fn_struct(1), fn_struct(2));

% Full aeroelastic eigenvalues
ev = eig(A);
fprintf('\nAeroelastic eigenvalues:\n'); disp(ev);

% Build ss object (requires Control System Toolbox)
if license('test','Control_Toolbox')
    sys = ss(A, B, C_out, D_out);
    sys.InputName  = {'w_g'};
    sys.OutputName = {'h','theta'};
    sys.StateName  = {'h','theta','hdot','thetadot'};
    fprintf('ss object created: sys\n');
end

% ---------- save -------------------------------------------------
save('wing_2dof.mat', 'A','B','C_out','D_out', ...
     'M_mat','Keff','Ctotal','Fgust', ...
     'U','rho','chord','span','Clalpha','ea', ...
     'fn_struct');

fprintf('\nSaved: wing_2dof.mat\n');
