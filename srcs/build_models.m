function build_models()
% build_models  Builds both Week-4 model files in this directory:
%   wing_2dof.mat -- linear 2-DOF model (Week-3 parameters, identical output)
%   wing_nl.mat   -- per-spring model with cubic stiffness + q=0 linearisation
%
% Replaces the need for week3/build_ss.m and week4/build_ss_nl.m at runtime.
%{
CF: Difference between build_models.m and the original build_ss.m and build_ss_nl.m:
build_models.m has two functions, build_2dof() and build_nl(), that construct the same physical models as build_ss.m and build_ss_nl.m, respectively. The parameters
Comparing the two files, build_models.m's build_nl() function covers all the same physics and exports identical variables to wing_nl.mat. But there is one missing piece:

Missing from build_models.m: the diagnostic fprintf statements that build_ss_nl.m prints (lines 73–75):

The full Ks matrix display (disp(Ks))
The coupling term report: Coupling term K_h,theta = ... (0 when k1 = k2)
build_models.m only prints the natural frequencies line, while build_ss_nl.m also prints the structural stiffness matrix and the off-diagonal coupling term.

Everything else — parameters, matrix assembly, state-space construction, saved variables — is identical. If those diagnostic prints matter to you, they can be added to build_nl().
%}
here = fileparts(mfilename('fullpath'));
build_2dof(here);
build_nl(here);
end

% =====================================================================

function [qS, Ka, Ca, Fgust, aero] = aero_terms()
aero.rho     = 1.225;   % air density          [kg/m^3]
aero.U       = 10.0;    % freestream speed     [m/s]
aero.chord   = 0.30;    % airfoil chord        [m]
aero.span    = 0.50;    % section span         [m]
aero.Clalpha = 6.0;    % lift-curve slope     [1/rad]
aero.ea      = 0.045;    % aero centre fwd of reference [m]

qS    = 0.5 * aero.rho * aero.U^2 * aero.chord * aero.span;
Ka    = [0,  qS*aero.Clalpha;
         0, -qS*aero.Clalpha*aero.ea];
Ca    = (qS*aero.Clalpha/aero.U) .* [1,        0;
                                     -aero.ea, 0];
Fgust = (qS*aero.Clalpha/aero.U) .* [-1; aero.ea];
end

% ---------------------------------------------------------------------

function build_2dof(here)
% Week-3 linear model, lumped springs (k each side, symmetric).
m  = 1.5;      % section mass                 [kg]
Ip = 0.0084;     % pitch inertia                [kg*m^2]
S  = 0.0225;    % static unbalance             [kg*m]
k  = 900;   % each spring stiffness        [N/m]
cd = 1.0;      % each damper coefficient      [N*s/m]
d  = 0.22;     % distance between springs     [m]

[~, Ka, Ca, Fgust, aero] = aero_terms();

M_mat = [m, S;
         S, Ip];
C_mat = [2*cd, 0;
         0,    cd*d^2/2];
Ks    = [2*k, 0;
         0,   k*d^2/2];

Keff   = Ks + Ka;
Ctotal = C_mat + Ca;

[A, B, C_out, D_out] = to_ss(M_mat, Keff, Ctotal, Fgust);
fn_struct = sqrt(eig(M_mat \ Ks)) / (2*pi);

U = aero.U;  rho = aero.rho;  chord = aero.chord;
span = aero.span;  Clalpha = aero.Clalpha;  ea = aero.ea;

save(fullfile(here, 'wing_2dof.mat'), ...
     'A','B','C_out','D_out', ...
     'M_mat','Keff','Ctotal','Fgust', ...
     'U','rho','chord','span','Clalpha','ea', ...
     'fn_struct');
fprintf('wing_2dof.mat  fn = %.3f / %.3f Hz\n', fn_struct(1), fn_struct(2));
end

% ---------------------------------------------------------------------

function build_nl(here)
% Per-spring model with cubic stiffness (cubic lives in the solver only).
m  = 2.0;      % section mass                 [kg]
Ip = 0.1;     % pitch inertia                [kg*m^2]
S  = 0.015;    % static unbalance             [kg*m]
d  = 0.25;     % distance between springs     [m]

k1   = 200.0;  % spring 1 stiffness           [N/m]
k2   = 200.0;  % spring 2 stiffness           [N/m]   !(set != k1 for h-theta coupling)
knl1 = 5e4;    % spring 1 cubic coefficient   [N/m^3] !(set > 0 for hardening)
knl2 = 5e4;    % spring 2 cubic coefficient   [N/m^3]
c1   = 3.0;    % damper 1 coefficient         [N*s/m]
c2   = 3.0;    % damper 2 coefficient         [N*s/m]

r1 = +d/2;     % spring/damper 1 moment arm   [m]
r2 = -d/2;     % spring/damper 2 moment arm   [m]

[~, Ka, Ca, Fgust, aero] = aero_terms();

M_mat = [m, S;
         S, Ip];
Ks    = [k1    + k2,     k1*r1   + k2*r2;
         k1*r1 + k2*r2,  k1*r1^2 + k2*r2^2];
C_struct = [c1    + c2,     c1*r1   + c2*r2;
            c1*r1 + c2*r2,  c1*r1^2 + c2*r2^2];
C_mat = C_struct + Ca;

[A, B, C_out, D_out] = to_ss(M_mat, Ks + Ka, C_mat, Fgust);
fn_struct = sqrt(eig(M_mat \ Ks)) / (2*pi);
U = aero.U;

save(fullfile(here, 'wing_nl.mat'), ...
     'A','B','C_out','D_out','fn_struct', ...
     'M_mat','C_mat','Ka','Fgust', ...
     'r1','r2','k1','k2','knl1','knl2', ...
     'U');
fprintf('wing_nl.mat    fn = %.3f / %.3f Hz\n', fn_struct(1), fn_struct(2));
end

% ---------------------------------------------------------------------

function [A, B, C_out, D_out] = to_ss(M, K, C, Fgust)
Minv  = inv(M);
A     = [zeros(2),  eye(2);
         -Minv*K,  -Minv*C];
B     = [zeros(2,1);
         Minv*Fgust];
C_out = [eye(2), zeros(2)];
D_out = zeros(2,1);
end
