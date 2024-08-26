% brainstem interleaved 3D-EPI sequence in Pulseq, calibration portion
%
% This short sequence first excites the volume to steady state, then 
% acquires many readout lines without Gy and Gz blips to:
% 1. Allow the scanner to tune receiver gains
% 2. Collect data used for EPI ghost correciton
%
% This script creates the file '3DEPI_cal.seq', that can be executed directly
% on Siemens MRI scanners using the Pulseq interpreter.
% The .seq file can also be converted to a .tar file that can be executed on GE
% scanners, see main.m.
%
% The experimental parameters below are chosen such that the sequence 
% can be executed identically (to us precision) on Siemens and GE systems.
% For more information about preparing a Pulseq file for execution on GE scanners,
% see the 'Pulseq on GE' manual.
%
% Modified July 12th, 2024 to acquire both polarities

%% Define experimental parameters
setEPIparams;

%% Path and options
seqname = '3DEPI_cal';
addpath('excitation/');

%% Excitation pulse
[rf, gzSS, gzSSR] = mr.makeSincPulse(alpha/180*pi,...
                                     'duration',rfDur,...
                                     'sliceThickness',fov(3),...
                                     'system',sys);
gzSS = trap4ge(gzSS,CRT,sys);
gzSSR = trap4ge(gzSSR,CRT,sys);

%% Fat-sat
fatsat.flip    = 90;      % degrees
fatsat.slThick = 1e5;     % dummy value (determines slice-select gradient, but we won't use it; just needs to be large to reduce dead time before+after rf pulse)
fatsat.tbw     = 3.5;     % time-bandwidth product
fatsat.dur     = 8.0;     % pulse duration (ms)
% RF waveform in Gauss
wav = toppe.utils.rf.makeslr(fatsat.flip, fatsat.slThick, fatsat.tbw, fatsat.dur, 1e-6, toppe.systemspecs(), ...
    'type', 'ex', ...    % fatsat pulse is a 90 so is of type 'ex', not 'st' (small-tip)
    'ftype', 'min', ...
    'writeModFile', false);
% Convert from Gauss to Hz, and interpolate to sys.rfRasterTime
rfp = rf2pulseq(wav, 4e-6, sys.rfRasterTime);
% Create pulseq object
% Try to account for the fact that makeArbitraryRf scales the pulse as follows:
% signal = signal./abs(sum(signal.*opt.dwell))*flip/(2*pi);
flip_ang = fatsat.flip/180*pi;
flipAssumed = abs(sum(rfp));
rfsat = mr.makeArbitraryRf(rfp, ...
    flip_ang*abs(sum(rfp*sys.rfRasterTime))*(2*pi), ...
    'system', sys);
rfsat.signal = rfsat.signal/max(abs(rfsat.signal))*max(abs(rfp)); % ensure correct amplitude (Hz)
rfsat.freqOffset = -fatOffresFreq;  % Hz

%% Generate temporally incoherent sampling masks
% Save sampling masks for recon
omegas = zeros(Ny,Nz,Nframes);
for frame = 1:Nframes
    % Create pseudo-random 2D sampling mask. Save for recon
    rng(frame); % A different mask per frame
    omega = randsamp2d(N(2:end), R, acs);
    omegas(:,:,frame) = omega;
end

%% Define readout gradients and ADC event
% The Pulseq toolbox really shines here!

% Define k-space spacing for fully-sampled data
deltak = 1./fov;

% Find the biggest step in ky for setting blip duration
max_ky_step = 0;
for frame = 1:Nframes
    for z = 1:Nz
        max_ky_step = max([max_ky_step, max(diff(find(omegas(:,z,frame))))]);
    end
end

% Start with the blips.
gyBlip = trap4ge(mr.scaleGrad(mr.makeTrapezoid('y', sys, 'Area', max_ky_step*deltak(2)), 1/max_ky_step),CRT,sys);
gzBlip = trap4ge(mr.makeTrapezoid('z', sys, 'Area', 0*deltak(3)),CRT,sys);  % unusued rn

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
gro = trap4ge(mr.makeTrapezoid('x', systmp, 'Area', Nx*deltak(1) + maxBlipArea),CRT,systmp);

% ADC event
Tread = mr.calcDuration(gro) - blipDuration;
if mod(round(Tread/dwell), sys.adcSamplesDivisor) % Ensure Nfid is a multiple of adcSamplesDivisor
    Tread = (round(Tread/dwell) - mod(round(Tread/dwell), sys.adcSamplesDivisor))*dwell;
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
gxPre = trap4ge(mr.makeTrapezoid('x', sys, ...
    'Area', -gro.area/2),CRT,sys);
Tpre = mr.calcDuration(gxPre);
gyPre = trap4ge(mr.makeTrapezoid('y', sys, ...
    'Area', -Ny/2*deltak(2), ... 
    'Duration', Tpre),CRT,sys);
gzPre = trap4ge(mr.makeTrapezoid('z', sys, ...
    'Area', -Nz/2*deltak(3), ...
    'Duration', Tpre),CRT,sys);
gxSpoil = trap4ge(mr.makeTrapezoid('x', sys, ...
    'Area', Nx*deltak(1)*NcyclesSpoil),CRT,sys);
gzSpoil = trap4ge(mr.makeTrapezoid('z', sys, ...
    'Area', Nz*deltak(3)*NcyclesSpoil),CRT,sys);

%% Calculate delay to achieve desired TE
minTE = mr.calcDuration(rf)/2 - rf.delay...
      + mr.calcDuration(gzSSR)...
      + mr.calcDuration(gzPre)...
      + (Ny/2 - 0.5) * mr.calcDuration(gro);
TEdelay = floor((TE - minTE)/sys.blockDurationRaster) * sys.blockDurationRaster;
if TE < minTE
    warning('Minimum achievable TE exceeds prescribed TE')
end

%% Calculate delay to achieve desired TR
minTR = mr.calcDuration(rfsat) + mr.calcDuration(gzSpoil)...
      + mr.calcDuration(gzSS) + mr.calcDuration(gzSSR)...
      + mr.calcDuration(gzPre) + Ny*mr.calcDuration(gro)...
      + mr.calcDuration(gzSpoil);
if TEdelay > 0
    minTR = minTR + TEdelay
end
TRdelay = floor((TR - minTR)/sys.blockDurationRaster)*sys.blockDurationRaster;
if TR < minTR
    warning('Minimum achievable TR exceeds prescribed TR')
end

%% Assemble sequence
seq = mr.Sequence(sys);

% RF spoiling stuff
rf_phase = 0;
rf_inc = 0;

for frame = -Ndummyframes:2

    % Convenience booleans for turning off adc and y gradient
    isDummyFrame = frame < 0;
    isCalFrame = frame >= 0;

    % Label the first block in each segment with the TRID (see Pulseq on GE manual)
    if isDummyFrame
        TRID = 1;
    elseif isCalFrame
        TRID = 2;
    end
    
    % Fat-sat
    seq.addBlock(rfsat,mr.makeLabel('SET','TRID',TRID));
    seq.addBlock(gxSpoil, gzSpoil);

    % RF spoiling
    rf.phaseOffset = rf_phase/180*pi;
    adc.phaseOffset = rf_phase/180*pi;
    rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
    rf_phase = mod(rf_phase+rf_inc, 360.0);

    % "Slice" selective RF excitation
    seq.addBlock(rf, gzSS);
    seq.addBlock(gzSSR);

    % TE delay
    if TE > minTE
        seq.addBlock(mr.makeDelay(TEdelay));
    end

    if isDummyFrame
        seq.addBlock(mr.makeDelay(mr.calcDuration(gxPre)...
                                 +mr.calcDuration(gro)*Ny...
                                 +TRdelay));
    elseif isCalFrame
        % Move to corner of k-space and sample the first line
        seq.addBlock(gxPre);
        seq.addBlock(adc, gro);
   
        % Zip through k-space with EPI trajectory
        for ie = 2:(Ny - 1)
            seq.addBlock(adc, mr.scaleGrad(gro, (-1)^(ie-1)));
        end

        % Last line
        seq.addBlock(adc, mr.scaleGrad(gro, (-1)^(ie)));
    end

    % spoil
    seq.addBlock(gxSpoil, gzSpoil);

    % Achieve desired TR
    if isCalFrame && TRdelay > 0
        seq.addBlock(mr.makeDelay(TRdelay));
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

%% Plot
figure('WindowState','maximized');
toppe.plotseq(sysGE, 'timeRange', [0, (Ndummyframes + 1)*max([TR,minTR])]);

return;
%% Detailed check that takes some time to run

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