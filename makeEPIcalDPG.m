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
caipiPythonPath = 'caipi/';

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

%% Define readout gradients and ADC event
% The Pulseq toolbox really shines here!

deltak = 1./fov;

% Start with the blips
gyBlip = trap4ge(mr.makeTrapezoid('y', sys, 'Area', max(abs(kyStep))*deltak(2)),CRT,sys); 
gzBlip = trap4ge(mr.makeTrapezoid('z', sys, 'Area', max(abs(kzStep))*deltak(3)),CRT,sys); 

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
    'Area', (kyInds(1)-Ny/2)*deltak(2), ... 
    'Duration', Tpre),CRT,sys);
gzPre = trap4ge(mr.makeTrapezoid('z', sys, ...
    'Area', -Nz/2*deltak(3), ...
    'Duration', Tpre),CRT,sys);
gxSpoil = trap4ge(mr.makeTrapezoid('x', sys, ...
    'Area', Nx*deltak(1)*NcyclesSpoil * (-1)^(Ny/Nsegments - 1)),CRT,sys);
gzSpoil = trap4ge(mr.makeTrapezoid('z', sys, ...
    'Area', Nx*deltak(1)*NcyclesSpoil * (-1)^(Ny/Nsegments - 1)),CRT,sys);

%% Calculate delay to achieve desired TE
kyIndAtTE = find(kyInds-Ny/2/Nsegments == min(abs(kyInds-Ny/2/Nsegments)));
minTE = mr.calcDuration(gzSS) + mr.calcDuration(gzSSR)...
      - mr.calcDuration(rf)/2 - rf.delay...
      + mr.calcDuration(gxPre)...
      + (kyIndAtTE-0.5) * mr.calcDuration(gro);
TEdelay = floor((TE-minTE)/sys.blockDurationRaster) * sys.blockDurationRaster;

%% Calculate delay to achieve desired TR
minTR = mr.calcDuration(rfsat) + mr.calcDuration(gzSpoil)...
      + mr.calcDuration(gzSS) + mr.calcDuration(gzSSR)...
      + mr.calcDuration(gzPre) + TEdelay...
      + Ny/Nsegments*mr.calcDuration(gro)...
      + mr.calcDuration(gzSpoil);
TRdelay = floor((TR - minTR)/sys.blockDurationRaster)*sys.blockDurationRaster;

%% Assemble sequence
seq = mr.Sequence(sys);           

kyStepMax = max(abs(kyStep));
kzStepMax = max(abs(kzStep));

% RF spoiling stuff
rf_phase = 0;
rf_inc = 0;

printf('Constructing sequence with positive readout polarity\n');
for frame = -Ndummyframes:0

    % Convenience booleans for turning off adc and y gradient
    isDummyFrame = frame < 0;
    isCalFrame = frame == 0;
    
    % No kz-encoding
    for z = floor(Nz/2)
        gzPreTmp = mr.scaleGrad(gzPre,(z - floor(Nz/2))/(Nz/2));
    
        % In plane loop (2D segmented EPI)
        for seg = 1:Nsegments
            % Label the first block in each segment with the TRID (see Pulseq on GE manual)
            TRID = 3*seg - (frame <= 0) - (frame == 0);
            
            % Fat-sat
            seq.addBlock(rfsat,mr.makeLabel('SET','TRID',TRID));
            seq.addBlock(gxSpoil, gzSpoil);
    
            % RF spoiling
            rf.phaseOffset = rf_phase/180*pi;
            adc.phaseOffset = rf_phase/180*pi;
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
    
            % "Slice" selective RF excitation
            seq.addBlock(rf,gzSS);

            % Rephase
            seq.addBlock(gzSSR);
    
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

%% Now do the same thing but with readout gradients flipped
gxPre = mr.scaleGrad(gxPre, -1);
gro = mr.scaleGrad(gro, -1);
gxSpoil = mr.scaleGrad(gxSpoil, -1);

printf('Constructing sequence with negative readout polarity\n');
for frame = -Ndummyframes:0

    % Convenience booleans for turning off adc and y gradient
    isDummyFrame = frame < 0;
    isCalFrame = frame == 0;
    
    % No kz-encoding
    for z = floor(Nz/2)
        gzPreTmp = mr.scaleGrad(gzPre,(z - floor(Nz/2))/(Nz/2));
    
        % In plane loop (2D segmented EPI)
        for seg = 1:Nsegments
            % Label the first block in each segment with the TRID (see Pulseq on GE manual)
            TRID = 3*seg - (frame <= 0) - (frame == 0);
            
            % Fat-sat
            seq.addBlock(rfsat,mr.makeLabel('SET','TRID',TRID));
            seq.addBlock(gxSpoil, gzSpoil);
    
            % RF spoiling
            rf.phaseOffset = rf_phase/180*pi;
            adc.phaseOffset = rf_phase/180*pi;
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
    
            % "Slice" selective RF excitation
            seq.addBlock(rf,gzSS);

            % Rephase
            seq.addBlock(gzSSR);
    
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
figure('WindowState','maximized');
toppe.plotseq(sysGE, 'timeRange',[0, 2*zTR]);

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