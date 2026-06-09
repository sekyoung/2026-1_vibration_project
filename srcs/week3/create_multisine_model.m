% create_multisine_model.m
% Builds wing_multisine_model.slx — Simulink equivalent of multisine_time_response.m
% Scope shows 3 stacked axes: w_g [m/s] | h [mm] | theta [deg]

here = fileparts(mfilename('fullpath'));
load(fullfile(here,'wing_2dof.mat'), 'A','B','C_out','D_out');

% ---------- precompute multisine (mirrors multisine_time_response.m) ----
f0    = 0.25;  kk = 1:80;  K = 80;
Ak    = 1.5/sqrt(K/2) * ones(1,K);
phi_k = -pi*kk .* (kk-1) / K;           % Schroeder phases
fs = 400;  dt = 1/fs;  Tper = 1/f0;  nper = 8;
N  = round(fs * Tper * nper);
t  = (0:N-1)' * dt;
wg_sig = sum(Ak .* sin(2*pi*kk.*t + phi_k), 2);

% Store as timeseries in base workspace (From Workspace reads this)
wg_ts = timeseries(wg_sig, t, 'Name','w_g');
assignin('base','wg_ts', wg_ts);

% ---------- build model ------------------------------------------
mdl = 'wing_multisine_model';
if bdIsLoaded(mdl), close_system(mdl,0); end
new_system(mdl);
load_system(mdl);

% Source: reads wg_ts from base workspace
add_block('simulink/Sources/From Workspace', [mdl '/wg_src'], ...
    'VariableName','wg_ts', 'SampleTime',num2str(dt), ...
    'Position',[50 160 150 190]);

% Continuous state-space (A,B,C,D from wing_2dof.mat)
add_block('simulink/Continuous/State-Space', [mdl '/wing_ss'], ...
    'A',mat2str(A), 'B',mat2str(B), ...
    'C',mat2str(C_out), 'D',mat2str(D_out), ...
    'X0','[0;0;0;0]', ...
    'Position',[220 145 370 205]);

% Demux: split [h; theta] into two scalar signals
add_block('simulink/Signal Routing/Demux', [mdl '/demux'], ...
    'Outputs','2', 'Position',[410 153 430 197]);

% Gain: h [m] → h [mm]
add_block('simulink/Math Operations/Gain', [mdl '/h_mm'], ...
    'Gain','1000', 'Position',[460 132 540 162]);

% Gain: theta [rad] → theta [deg]
add_block('simulink/Math Operations/Gain', [mdl '/th_deg'], ...
    'Gain',num2str(180/pi,'%.8f'), 'Position',[460 183 540 213]);

% Scope: 3 input ports → 3 stacked axes in one window
add_block('simulink/Sinks/Scope', [mdl '/scope'], ...
    'NumInputPorts','3', 'Position',[610 155 655 205]);

% ---------- connections ------------------------------------------
% wg fans out: → wing_ss input AND → scope port 1
add_line(mdl,'wg_src/1', 'wing_ss/1','autorouting','smart');
add_line(mdl,'wg_src/1', 'scope/1',  'autorouting','smart');

% wing output → demux → gains → scope ports 2 & 3
add_line(mdl,'wing_ss/1','demux/1',  'autorouting','smart');
add_line(mdl,'demux/1',  'h_mm/1',   'autorouting','smart');
add_line(mdl,'demux/2',  'th_deg/1', 'autorouting','smart');
add_line(mdl,'h_mm/1',   'scope/2',  'autorouting','smart');
add_line(mdl,'th_deg/1', 'scope/3',  'autorouting','smart');

% ---------- solver: fixed-step ode4 at same dt as ZOH ------------
set_param(mdl,'SolverType','Fixed-step','Solver','ode4', ...
              'FixedStep',num2str(dt),'StopTime',num2str(t(end)));

% ---------- save -------------------------------------------------
out = fullfile(here,[mdl '.slx']);
save_system(mdl, out);
close_system(mdl);
fprintf('Saved: %s\n', out);
fprintf('To run: load_system(''%s''); sim(''%s'')\n', mdl, mdl);
