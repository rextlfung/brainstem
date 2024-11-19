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

% Type of EPI acquitsition. "multishot" or "singleshot"
mode = "singleshot";

% Basic spatial parameters
res = [1.8 1.8 1.8]*1e-3; % resolution (m)
fov = [21.6, 21.6, 7.2]*1e-2; % field of view (m)
N = round(fov./res); % acquisiton tensor size
Nx = N(1); Ny = N(2); Nz = N(3);


% Random sampling parameters
if strcmp(mode, "singleshot")
    Ry = 2; Rz = 1;
    R = [Ry Rz];                    % Acceleration/undersampling factors in each direction
    acs = [1/16 1/16];              % Central portion of ky-kz space to fully sample
    max_ky_step = round(Ny/16);     % Maximum gap in fast PE direction
    caipi_z = 3;                    % Minimum gap in slow PE direction to prevent duplicate sampling with CAIPI shift
else
    Ry = 1; Rz = 1;
end

% Number of shots per volume
if strcmp(mode, "multishot")
    Nsegments = 4;                  % number of segments in EPI readout
    Nshots = Nz*Nsegments;          % number of shots per volume
elseif strcmp(mode, "singleshot")
    Nshots = ceil(Nz/Rz/caipi_z);           % number of shots per volume
end

% Basic temporal parameters
Ndummyframes = 4;                   % dummy frames to reach steady state for calibration
NframesPerLoop = lcm(40,Nshots)/Nshots; % number of temporal frames to complete one RF spoil cycle

% ADC stuff
dwell = 4e-6;                       % ADC sample time (s). For GE, must be multiple of 2us.

% Decay parameters
TE = 30e-3;                         % echo time (s)
volumeTR = 0.9;                     % temporal frame rate (s)
TR = volumeTR/Nshots;      % time to acquire one shot (s)
T1 = 1500e-3;                       % T1 (s)

% Exciting stuff
alpha = 180/pi * acos(exp(-TR/T1)); % Ernst angle (degrees)
rfDur = 2e-3;                       % RF pulse duration (s)
rfTB  = 6;                          % RF pulse time-bandwidth product
rf_phase_0 = 117;                   % RF spoiling initial phase (degrees)
NcyclesSpoil = 2;                   % number of Gx and Gz spoiler cycles

% Fat Sat Stuff
fatChemShift = 3.5*1e-6;                        % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz