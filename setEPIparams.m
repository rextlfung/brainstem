%% Script for defining experimental parameters to be used across acq + recon
sys = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
              'maxSlew', 120, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ...
              'rfRingdownTime', 60e-6, ...
              'adcDeadTime', 20e-6, ...
              'adcRasterTime', 2e-6, ...
              'gradRasterTime', 10e-6, ...
              'blockDurationRaster', 10e-6, ...
              'B0', 3.0);
sysGE = toppe.systemspecs('maxGrad', sys.maxGrad/sys.gamma*100, ...   % G/cm
    'maxSlew', sys.maxSlew/sys.gamma/10, ...           % G/cm/ms
    'maxRF', 0.25);
CRT = 20e-6; % Common raster time of Siemens: 10e-6, GE: 4e-6;

% Basic spatial parameters
fov = [180, 180, 10]*1e-3;           % field of view
Nx = 180; Ny = Nx; Nz = 10;          % acquisition sizes
Nsegments = 4;                      % number of segments in EPI readout
Nshots = Nz*Nsegments;              % number of shots per volume

% Basic temporal parameters
Ndummyframes = 4;                   % dummy frames to reach steady stat
NframesPerLoop = lcm(40,Nshots)/Nshots; % number of temporal frames to complete one RF spoil cycle

% CAIPI sampling parameters
Ry = 1;                             % ky undersampling factor
Rz = 1;                             % kz undersampling factor
CaipiShiftZ = 2;                    % Caipi shift factor
pf_ky = 1.0;                        % partial Fourier factor along ky
etl = ceil(pf_ky*Ny/Ry);            % echo train length
NcyclesSpoil = 3;                   % number of Gx and Gz spoiler cycles

% ADC stuff
dwell = 4e-6;                       % ADC sample time (s). For GE, must be multiple of 2us.

% Decay parameters
TE = 30e-3;                         % echo time (s)
volumeTR = Nz*300e-3;               % temporal frame rate (s)
zTR = volumeTR/Nz;                  % time to acquire a "slice" or z (s)
TR = zTR/Nsegments;                 % time between excitations (s)
T1 = 1500e-3;                       % T1 (s)

% Exciting stuff
alpha = 180/pi * acos(exp(-TR/T1)); % Ernst angle (degrees)
rfDur = 8e-3;                       % RF pulse duration (s)
rfTB  = 6;                          % RF pulse time-bandwidth product
rfSpoilingInc = 117;                % RF spoiling increment (degrees)

% Fat Sat Stuff
fatChemShift = 3.5*1e-6;                        % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz