% brainstem interleaved 3D-EPI sequence in Pulseq, calibration portion
%
% This short sequence first excites the volume to steady state, then 
% acquires many readout lines without Gy and Gz blips to:
% 1. Allow the scanner to tune receiver gains
% 2. Collect data used for EPI ghost correciton
%
% This script creates the file '3DEPImultishot_cal.seq', that can be executed directly
% on Siemens MRI scanners using the Pulseq interpreter.
% The .seq file can also be converted to a .tar file that can be executed on GE
% scanners, see main.m.
%
% The experimental parameters below are chosen such that the sequence 
% can be executed identically (to us precision) on Siemens and GE systems.
% For more information about preparing a Pulseq file for execution on GE scanners,
% see the 'Pulseq on GE' manual.
%
% Performance?
% Low temporal resolution (1s) but
% High temporal resolution (1mm isotropic)
%% Path and options
seqname = '3DEPImultishot_cal';
addpath('excitation/');
caipiPythonPath = 'caipi/';

%% Define experimental parameters
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

% Basic parameters
fov = [200, 200, 5]*1e-3;           % field of view
Nx = 200; Ny = Nx; Nz = 10;          % Matrix sizes
Nsegments = 4;                      % number of segments in EPI readout

% Basic temporal parameters
Ndummyframes = 4;                   % dummy frames to reach steady state

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
volumeTR = Nz*300e-3;                  % temporal frame rate (s)
zTR = volumeTR/Nz;                  % time to acquire a "slice" or z (s)
TR = zTR/Nsegments;                 % time between excitations (s)
T1 = 1500e-3;                       % T1 (s)

% Excitation stuff
alpha = 180/pi * acos(exp(-TR/T1)); % Ernst angle (degrees)
rfDur = 8e-3;                       % RF pulse duration (s)
rfTB  = 6;                          % RF pulse time-bandwidth product
rfSpoilingInc = 117;                % RF spoiling increment (degrees)

% Fat Sat Stuff
fatChemShift = 3.5*1e-6;                        % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz

%% Get CAIPI sampling pattern (for one shot/echo train)

% create caipi.mat, and load it
pyFile = [caipiPythonPath 'skippedcaipi_sampling.py'];
pyCmd = sprintf('python3 %s %d %d %d %d %d %d', ...
    pyFile, Ny, Nz, Ry, Rz, CaipiShiftZ, 1);

% try to call Python script from Matlab
pyenv(ExecutionMode="InProcess");
[status, cmdout] = system(pyCmd, 'LD_PRELOAD', '/lib/x86_64-linux-gnu/libstdc++.so.6');
if status == 1 % if failed
    fprintf(cmdout);
    fprintf('Open a terminal and run the following python command:\n\t%s\n', pyCmd);
    input('Then press Enter to continue');
end

load caipi  

% kz and ky indeces (multiples of deltak)
kyInds = double(indices((end-etl+1):end, 2));
kzInds = double(indices((end-etl+1):end, 1));

% ky/kz encoding blip amplitude along echo train (multiples of deltak)
kyStep = diff(kyInds);
kzStep = diff(kzInds);

%% Excitation pulse
% Reusing Jon's SMS pulseq function for volume excitation
mb = 1;
slThick = fov(3);
sliceSep = fov(3)/mb;   % center-to-center separation between SMS slices (m)
[rf, gzRF, freq] = getsmspulse(alpha, slThick, rfTB, rfDur, ...
    mb, sliceSep, sysGE, sys, ...
    'doSim', false, ...    % Plot simulated SMS slice profile
    'type', 'st', ...     % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
    'ftype', 'ls');       % filter design. 'ls' = least squares

%% Fat sat pulse (TODO!!!)

%% Define readout gradients and ADC event
% The Pulseq toolbox really shines here!

deltak = 1./fov;

% Start with the blips
gyBlip = mr.makeTrapezoid('y', sys, 'Area', max(abs(kyStep))*deltak(2)); 
gzBlip = mr.makeTrapezoid('z', sys, 'Area', max(abs(kzStep))*deltak(3)); 

% Area and duration of the biggest blip
if gyBlip.area > gzBlip.area
    maxBlipArea = gyBlip.area;
    blipDuration = mr.calcDuration(gyBlip);
else
    maxBlipArea = gzBlip.area;
    blipDuration = mr.calcDuration(gzBlip);
end

% Readout trapezoid
systmp = sys;
systmp.maxGrad = deltak(1)/dwell;  % to ensure >= Nyquist sampling
gro = mr.makeTrapezoid('x', systmp, 'Area', Nx*deltak(1) + maxBlipArea);

% ADC event
Tread = mr.calcDuration(gro) - blipDuration;
if mod(round(Tread*1e6)*1e-6, dwell)
    Tread = Tread - mod(Tread, dwell) + dwell;
end
adc = mr.makeAdc(round(Tread/dwell), sys, ...
    'Duration', Tread, ...
    'Delay', blipDuration/2);

% Split blips at block boundary.
[gyBlipUp, gyBlipDown] = mr.splitGradientAt(gyBlip, mr.calcDuration(gyBlip)/2);
gyBlipUp.delay = mr.calcDuration(gro) - mr.calcDuration(gyBlip)/2;
gyBlipDown.delay = 0;

[gzBlipUp, gzBlipDown] = mr.splitGradientAt(gzBlip, mr.calcDuration(gzBlip)/2);
gzBlipUp.delay = mr.calcDuration(gro) - mr.calcDuration(gzBlip)/2;
gzBlipDown.delay = 0;

% prephasers and spoilers
gxPre = mr.makeTrapezoid('x', sys, ...
    'Area', -gro.area/2);
Tpre = mr.calcDuration(gxPre);
gyPre = mr.makeTrapezoid('y', sys, ...
    'Area', (kyInds(1)-Ny/2)*deltak(2), ... 
    'Duration', Tpre);
gzPre = mr.makeTrapezoid('z', sys, ...
    'Area', -1/2*deltak(3), ...
    'Duration', Tpre);
gxSpoil = mr.makeTrapezoid('x', sys, ...
    'Area', Nx*deltak(1)*NcyclesSpoil * (-1)^(Ny/Nsegments - 1));
gzSpoil = mr.makeTrapezoid('z', sys, ...
    'Area', Nx*deltak(1)*NcyclesSpoil * (-1)^(Ny/Nsegments - 1));

%% Calculate delay to achieve desired TE
kyIndAtTE = find(kyInds-Ny/2/Nsegments == min(abs(kyInds-Ny/2/Nsegments)));
minTE = mr.calcDuration(gzRF) - mr.calcDuration(rf)/2 - rf.delay + ...
        mr.calcDuration(gxPre) + ...
        (kyIndAtTE-0.5) * mr.calcDuration(gro);
TEdelay = floor((TE-minTE)/sys.blockDurationRaster) * sys.blockDurationRaster;

%% Calculate delay to achieve desired TR
minTR = mr.calcDuration(gzRF) + mr.calcDuration(gzPre)...
                + Ny/Nsegments*mr.calcDuration(gro) + mr.calcDuration(gzSpoil);
TRdelay = floor((TR - minTR)/sys.blockDurationRaster)*sys.blockDurationRaster;

%% Assemble sequence
seq = mr.Sequence(sys);           

kyStepMax = max(abs(kyStep));
kzStepMax = max(abs(kzStep));

% RF spoiling stuff
rf_phase = 0;
rf_inc = 0;

for frame = -Ndummyframes:0

    % Convenience booleans for turning off adc and y gradient
    isDummyFrame = frame < 0;
    isCalFrame = frame == 0;
    
    % z-loop (move to proper kz location)
    for z = 1:Nz
        gzPreTmp = mr.scaleGrad(gzPre,ceil(Nz/2) - z);
    
        % In plane loop (2D segmented EPI)
        for seg = 1:Nsegments
            % Label the first block in each segment with the segment ID (see Pulseq on GE manual)
            % !!This segment ID is NOT to be confused with the EPI segments!!
            segmentID = 3*seg - (frame <= 0) - (frame == 0);
    
            % RF spoiling
            rf.phaseOffset = rf_phase/180*pi;
            adc.phaseOffset = rf_phase/180*pi;
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
    
            % RF excitation
            seq.addBlock(rf, gzRF, mr.makeLabel('SET', 'TRID', segmentID));
    
            % TE delay
            if TE > minTE
                seq.addBlock(mr.makeDelay(TEdelay));
            end

            % Segment delay so TE is smooth along central line of k-space
            segDelay = (seg - 1)*mr.calcDuration(gro);
            if seg > 1
                seq.addBlock(mr.makeDelay(segDelay));
            end
    
            % Move to corner of k-space and sample the first line
            gyPreSeg = mr.scaleGrad(gyPre, ((seg - 1)*gyBlip.area + gyPre.area)/gyPre.area);
            if isDummyFrame
                seq.addBlock(gxPre, gyPreSeg, gzPreTmp);
                seq.addBlock(gro, ...
                             mr.scaleGrad(gyBlipUp, Nsegments*kyStep(1)/kyStepMax));
            elseif isCalFrame
                seq.addBlock(gxPre);
                seq.addBlock(gro, adc);
            else
                seq.addBlock(gxPre, gyPreSeg, gzPreTmp);
                seq.addBlock(gro, adc,...
                             mr.scaleGrad(gyBlipUp, Nsegments*kyStep(1)/kyStepMax));
            end
    
            % Zip through k-space with EPI trajectory
            for ie = 2:(etl/Nsegments - 1)
                gybd = mr.scaleGrad(gyBlipDown, Nsegments*kyStep(ie-1)/kyStepMax);
                gybu = mr.scaleGrad(gyBlipUp, Nsegments*kyStep(ie)/kyStepMax);
                gybdu = mr.addGradients({gybd, gybu}, sys);
                if isDummyFrame
                    seq.addBlock(mr.scaleGrad(gro, (-1)^(ie-1)), gybdu);
                elseif isCalFrame
                    seq.addBlock(adc, mr.scaleGrad(gro, (-1)^(ie-1)));
                else
                    seq.addBlock(adc, mr.scaleGrad(gro, (-1)^(ie-1)), gybdu);
                end
            end
    
            % Last line
            if isDummyFrame
                seq.addBlock(mr.scaleGrad(gro, (-1)^(ie)), ...
                             mr.scaleGrad(gyBlipDown, Nsegments*kyStep(ie)/kyStepMax));
            elseif isCalFrame
                seq.addBlock(adc, ...
                             mr.scaleGrad(gro, (-1)^(ie)));
            else
                seq.addBlock(adc, ...
                             mr.scaleGrad(gro, (-1)^(ie)), ...
                             mr.scaleGrad(gyBlipDown, Nsegments*kyStep(ie)/kyStepMax));
            end
    
            % spoil
            seq.addBlock(gxSpoil, gzSpoil);
    
            % Achieve desired TR
            seq.addBlock(mr.makeDelay(TRdelay - segDelay));
        end
    end
end

%% Check sequence timing
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% Output for execution
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', seqname);
seq.write(strcat(seqname, '.seq'));

%% GE stuff
seq2ge(strcat(seqname, '.seq'), sysGE, strcat(seqname, '.tar'))
system(sprintf('tar -xvf %s', strcat(seqname, '.tar')));
figure; toppe.plotseq(sysGE, 'timeRange',[0, volumeTR]);

%% Detailed check that takes some time to run
doDetailedCheck = false;
if doDetailedCheck
    % k-space trajectory calculation and plot
    [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
    figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D k-space plot
    axis('equal'); % enforce aspect ratio for the correct trajectory display
    hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
    title('full k-space trajectory (k_x x k_y)');
    
    % Optional slow step, but useful for testing during development,
    % e.g., for the real TE, TR or for staying within slewrate limits
    rep = seq.testReport;
    fprintf([rep{:}]);
end