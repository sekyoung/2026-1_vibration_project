function wg_ts = multisine_function_generator(f0, K, rms_amp, fs, nper)
% MULTISINE_FUNCTION_GENERATOR  Generate a Schroeder-phase multisine as timeseries.
%
%   wg_ts = multisine_function_generator()          — default params
%   wg_ts = multisine_function_generator(f0,K,rms,fs,nper)
%
%   f0      fundamental frequency [Hz]   default 0.25
%   K       number of harmonics          default 80
%   rms_amp desired RMS amplitude [m/s]  default 1.5
%   fs      sample rate [Hz]             default 400
%   nper    number of periods            default 8

if nargin < 1, f0      = 0.25; end
if nargin < 2, K       = 80;   end
if nargin < 3, rms_amp = 1.5;  end
if nargin < 4, fs      = 400;  end
if nargin < 5, nper    = 8;    end

kk    = 1:K;
Ak    = rms_amp / sqrt(K/2) * ones(1, K);   % flat amplitude, RMS = rms_amp
phi_k = -pi * kk .* (kk-1) / K;             % Schroeder phases (low crest factor)

dt   = 1/fs;
N    = round(fs * (1/f0) * nper);
t    = (0:N-1)' * dt;
wg   = sum(Ak .* sin(2*pi*kk.*t + phi_k), 2);

wg_ts = timeseries(wg, t, 'Name','w_g');
wg_ts.TimeInfo.Units = 's';
end
