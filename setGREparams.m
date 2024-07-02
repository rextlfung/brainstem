%% Script for defining experimental parameters to be used for a GRE scan
% For sensitivity maps

% RF/gradient delay (sec). 
% Conservative choice that should work across all GE scanners.
psd_rf_wait = 200e-6;

% Here we extend rfRingdownTime by psd_rf_wait to ensure that 
% the subsequent wait pulse (delay block) doesn't overlap with
% the 'true' RF ringdown time (54 us).
sys = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
              'maxSlew', 150, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ...
              'rfRingdownTime', 60e-6 + psd_rf_wait, ...
              'adcDeadTime', 40e-6, ...
              'adcRasterTime', 2e-6, ...
              'blockDurationRaster', 10e-6, ...
              'B0', 3.0);

% Spatial parameters (different from EPI scans)
fov_gre = [180, 180, 10]*1e-3;      % field of view
Nx_gre = 180; Ny_gre = Nx_gre; Nz_gre = 10; % Matrix size. Make large enough to avoid foldover

% Other acquisition params
dwell = 4e-6;                   % ADC sample time (s)
fatChemShift = 3.5e-6;          % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz
TE = 1/fatOffresFreq*[1];       % fat and water in phase for both echoes
TR = 6e-3*[1];                  % constant TR
T1 = 1500e-3;                   % approximate T1
alpha = 180/pi*acos(exp(-TR/T1));  % flip angle (degrees)

% Sequence parameters
alphaPulseDuration = 0.2e-3;
rfSpoilingInc = 117;            % RF spoiling increment
nCyclesSpoil = 2;               % number of spoiler cycles
Tpre = 1.0e-3;                  % prephasing trapezoid duration