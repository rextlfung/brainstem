% brainstem 2D-EPI sequence in Pulseq
%
% This script creates the file 'brainstem2DEPi.seq', that can be executed directly
% on Siemens MRI scanners using the Pulseq interpreter.
% The .seq file can also be converted to a .tar file that can be executed on GE
% scanners, see main.m.
%
% The experimental parameters below are chosen such that the sequence 
% can be executed identically (to us precision) on Siemens and GE systems.
% For more information about preparing a Pulseq file for execution on GE scanners,
% see the 'Pulseq on GE' manual.
%
% Performance appears to hold up to 120x120 (does better on inside scanner)
%% Paths
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
mb = 1;                             % multiband/SMS factor
Nx = 180; Ny = Nx; Nz = 1;           % Matrix sizes
fov = [180 180 3]*1e-3;             % field of view
slThick = fov(3)/Nz;                % slice thickness

% Basic temporal parameters
Nframes = 2;                        % number of temporal frames (image volumes)
Ndummyframes = 4;                   % dummy frames to reach steady state

% CAIPI sampling parameters
Ry = 1;                             % ky undersampling factor
Rz = mb;                            % kz undersampling factor
CaipiShiftZ = 1;                    % Caipi shift factor
pf_ky = 1.0;                        % partial Fourier factor along ky
etl = ceil(pf_ky*Ny/Ry);            % echo train length
nCyclesSpoil = 3;                   % number of Gx and Gz spoiler cycles

% ADC stuff
dwell = 4e-6;                       % ADC sample time (s). For GE, must be multiple of 2us.

% Decay parameters
TE = 30e-3;                         % echo time (s)
TR = 400e-3;                        % repetition time (s)
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
pyCmd = sprintf('python %s %d %d %d %d %d %d', ...
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
sliceSep = fov(3)/mb;   % center-to-center separation between SMS slices (m)
[rf, gzRF, freq] = getsmspulse(alpha, slThick, rfTB, rfDur, ...
    mb, sliceSep, sysGE, sys, ...
    'doSim', true, ...    % Plot simulated SMS slice profile
    'type', 'st', ...     % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
    'ftype', 'ls');       % filter design. 'ls' = least squares

% Doesn't work because unable to run python for some reason
% [rf, gzRF, gzr, delay] = mr.makeSLRpulse(alpha,'sliceThickness',slThick,...
%                                       'timeBwProduct',rfTB,'duration',rfDur,...
%                                       'system', sys);

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
    'Area', -Nz/2*deltak(3), ... 
    'Duration', Tpre);
gxSpoil = mr.makeTrapezoid('x', sys, ...
    'Area', -Nx*deltak(1)*nCyclesSpoil);
gzSpoil = mr.makeTrapezoid('z', sys, ...
    'Area', -Nx*deltak(1)*nCyclesSpoil);


%% Calculate delay to achieve desired TE
kyIndAtTE = find(kyInds-Ny/2 == min(abs(kyInds-Ny/2)));
minTE = mr.calcDuration(gzRF) - mr.calcDuration(rf)/2 - rf.delay + mr.calcDuration(gxPre) + ...
        (kyIndAtTE-0.5) * mr.calcDuration(gro);
TEdelay = floor((TE-minTE)/sys.blockDurationRaster) * sys.blockDurationRaster;

%% Calculate delay to achieve desired TR
minTR = mr.calcDuration(rf) + rf.delay + mr.calcDuration(gxPre)...
                + Ny*mr.calcDuration(gro) + mr.calcDuration(gxSpoil);
TRdelay = floor((TR - minTR)/sys.blockDurationRaster)*sys.blockDurationRaster;

%% Assemble sequence
seq = mr.Sequence(sys);           

nShots = Nz/mb;
kyStepMax = max(abs(kyStep));
kzStepMax = max(abs(kzStep));

% RF spoiling stuff
rf_phase = 0;
rf_inc = 0;

% Create scanloop
for frame = -Ndummyframes:Nframes

        % Convenience booleans for turning off adc and y gradient
        isDummyFrame = frame < 0;
        isCalFrame = frame == 0;

        % Label the first block in each segment with the segment ID (see Pulseq on GE manual)
        segmentID = 3 - (frame <= 0) - (frame == 0);

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

        % Move to corner of k-space
        if mb > 1
            seq.addBlock(gxPre, gyPre, mr.scaleGrad(gzPre, 1 - 2/Nz*(frame-1)));
            seq.addBlock(gro, adc, ...
                         mr.scaleGrad(gyBlipUp, kyStep(1)/kyStepMax), ...
                         mr.scaleGrad(gzBlipUp, kzStep(1)/kzStepMax));
        else
            if isDummyFrame
                seq.addBlock(gxPre, gyPre);
                seq.addBlock(gro, ...
                             mr.scaleGrad(gyBlipUp, kyStep(1)/kyStepMax));
            elseif isCalFrame
                seq.addBlock(gxPre);
                seq.addBlock(gro, adc);
            else
                seq.addBlock(gxPre, gyPre);
                seq.addBlock(gro, adc, ...
                             mr.scaleGrad(gyBlipUp, kyStep(1)/kyStepMax));
            end
        end

        % Zip through k-space with EPI trajectory
        for ie = 2:(etl-1)
            gybd = mr.scaleGrad(gyBlipDown, kyStep(ie-1)/kyStepMax);
            gybu = mr.scaleGrad(gyBlipUp, kyStep(ie)/kyStepMax);
            gybdu = mr.addGradients({gybd, gybu}, sys);
            if mb > 1
                gzbd = mr.scaleGrad(gzBlipDown, kzStep(ie-1)/kzStepMax);
                gzbu = mr.scaleGrad(gzBlipUp, kzStep(ie)/kzStepMax);
                gzbdu = mr.addGradients({gzbd, gzbu}, sys);
                seq.addBlock(adc, mr.scaleGrad(gro, (-1)^(ie-1)), gybdu, gzbdu);
            else
                if isDummyFrame
                    seq.addBlock(mr.scaleGrad(gro, (-1)^(ie-1)), gybdu);
                elseif isCalFrame
                    seq.addBlock(adc, mr.scaleGrad(gro, (-1)^(ie-1)));
                else
                    seq.addBlock(adc, mr.scaleGrad(gro, (-1)^(ie-1)), gybdu);
                end
            end
        end

        % ?
        if mb > 1
            seq.addBlock(adc, ...
                         mr.scaleGrad(gro, (-1)^(ie)), ...
                         mr.scaleGrad(gyBlipDown, kyStep(ie)/kyStepMax), ...
                         mr.scaleGrad(gzBlipDown, kzStep(ie)/kzStepMax));
        else
            if isDummyFrame
                seq.addBlock(mr.scaleGrad(gro, (-1)^(ie)), ...
                             mr.scaleGrad(gyBlipDown, kyStep(ie)/kyStepMax));
            elseif isCalFrame
                seq.addBlock(adc, ...
                             mr.scaleGrad(gro, (-1)^(ie)));
            else
                seq.addBlock(adc, ...
                             mr.scaleGrad(gro, (-1)^(ie)), ...
                             mr.scaleGrad(gyBlipDown, kyStep(ie)/kyStepMax));
            end
         end

        % spoil
        seq.addBlock(gxSpoil, gzSpoil);

        % Achieve desired TR
        seq.addBlock(mr.makeDelay(TRdelay));
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

%% Output for execution and plot
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'brainstemEPI2D');
seq.write('brainstemEPI2D.seq')       % Write to pulseq file

% seq.plot(); %'timeRange', [0 0.2]);

%% GE stuff
seq2ge('brainstemEPI2D.seq', sysGE, 'brainstemEPI2D.tar')
system('tar -xvf brainstemEPI2D.tar')
toppe.plotseq(sysGE);

%% k-space trajectory calculation and plot
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('full k-space trajectory (k_x x k_y)');

%% Optional slow step, but useful for testing during development,
% e.g., for the real TE, TR or for staying within slewrate limits
rep = seq.testReport;
fprintf([rep{:}]);