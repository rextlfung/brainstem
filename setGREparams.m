%% Script for defining experimental parameters to be used for a GRE scan
% For sensitivity maps

% RF/gradient delay (sec). 
% Conservative choice that should work across all GE scanners.
psd_rf_wait = 200e-6;

% Here we extend rfRingdownTime by psd_rf_wait to ensure that 
% the subsequent wait pulse (delay block) doesn't overlap with
% the 'true' RF ringdown time (54 us).
sys = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
              'maxSlew', 120, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ...
              'rfRingdownTime', 60e-6 + psd_rf_wait, ...
              'adcDeadTime', 20e-6, ...
              'adcRasterTime', 2e-6, ...
              'blockDurationRaster', 10e-6, ...
              'B0', 3.0);
sysGE = toppe.systemspecs('maxGrad', sys.maxGrad/sys.gamma*100, ...   % G/cm
    'maxSlew', sys.maxSlew/sys.gamma/10, ...           % G/cm/ms
    'maxRF', 0.25);
CRT = 20e-6; % Common raster time of Siemens: 10e-6, GE: 4e-6;

% Spatial parameters (different from EPI scans)
res_gre = [2, 2, 2]*1e-3; % resolution (m)
N_gre = [100, 100, 100]; % acquisiton tensor size
Nx_gre = N_gre(1); Ny_gre = N_gre(2); Nz_gre = N_gre(3);
fov_gre = N_gre.*res_gre; % field of view (m)

% Other acquisition params
dwell = 4e-6;                   % ADC sample time (s)
fatChemShift = 3.5e-6;          % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz
TE = 1/fatOffresFreq + 2e-4;    % fat and water in phase for both echoes
TR = 6e-3;                      % constant TR
T1 = 1500e-3;                   % approximate T1
alpha = 180/pi*acos(exp(-TR/T1));  % flip angle (degrees)

% Sequence parameters
rfDur = 0.4e-3;
rf_phase_0 = 117;               % RF spoiling initial phase (degrees)
nCyclesSpoil = 2;               % number of spoiler cycles
Tpre = 1.0e-3;                  % prephasing trapezoid duration