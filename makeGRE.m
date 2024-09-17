% Pulseq script for acquiring GRE data
%
% Modified June 21, 2024, for acquiring GRE data for sensitivity maps

% Define params (edit setGREparams to change)
setGREparams;

% Create a new sequence object
seq = mr.Sequence(sys);           

% Fat-sat
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

% Same slab-selective pulse as 3D EPI
[rf, gzSS, gzSSR] = mr.makeSincPulse(alpha/180*pi,...
                                     'duration',rfDur,...
                                     'sliceThickness',fov_gre(3),...
                                     'system',sys);
gzSS = trap4ge(gzSS,CRT,sys);
gzSSR = trap4ge(gzSSR,CRT,sys);

% Define other gradients and ADC events
% Cut the redaout gradient into two parts for optimal spoiler timing
deltak = 1./fov_gre;
Tread = Nx_gre*dwell;

commonRasterTime = 20e-6;   

gyPre = trap4ge(mr.makeTrapezoid('y', sys, ...
    'Area', Ny_gre*deltak(2)/2, ...   % PE1 gradient, max positive amplitude
    'Duration', Tpre), ...
    commonRasterTime, sys);
gzPre = trap4ge(mr.makeTrapezoid('z', sys, ...
    'Area', Nz_gre*deltak(3)/2, ...   % PE2 gradient, max amplitude
    'Duration', Tpre), ...
    commonRasterTime, sys);

gxtmp = mr.makeTrapezoid('x', sys, ...  % readout trapezoid, temporary object
    'Amplitude', Nx_gre*deltak(1)/Tread, ...
    'FlatTime', Tread);
gxPre = trap4ge(mr.makeTrapezoid('x', sys, ...
    'Area', -gxtmp.area/2, ...
    'Duration', Tpre), ...
    commonRasterTime, sys);

adc = mr.makeAdc(Nx_gre, sys, ...
    'Duration', Tread,...
    'Delay', gxtmp.riseTime);

% extend flat time so we can split at end of ADC dead time
gxtmp2 = trap4ge(mr.makeTrapezoid('x', sys, ...  % temporary object
    'Amplitude', Nx_gre*deltak(1)/Tread, ...
    'FlatTime', Tread + adc.deadTime), ...
    commonRasterTime, sys);
[gx, ~] = mr.splitGradientAt(gxtmp2, gxtmp2.riseTime + gxtmp2.flatTime);
%gx = gxtmp;

gzSpoil = mr.makeTrapezoid('z', sys, ...
    'Area', Nx_gre*deltak(1)*nCyclesSpoil);
gxSpoil = mr.makeExtendedTrapezoidArea('x', gxtmp.amplitude, 0, gzSpoil.area, sys);
%gxSpoil = mr.makeTrapezoid('x', sys, ...
%    'Area', Nx_gre*deltak(1)*nCyclesSpoil);

%% y/z PE steps
pe1Steps = ((0:Ny_gre-1)-Ny_gre/2)/Ny_gre*2;
pe2Steps = ((0:Nz_gre-1)-Nz_gre/2)/Nz_gre*2;

%% Calculate timing
TEmin = mr.calcDuration(rf)/2 + mr.calcDuration(gzSSR)...
      + mr.calcDuration(gxPre) + adc.delay + Nx_gre/2*dwell;
delayTE = ceil((TE-TEmin)/seq.gradRasterTime)*seq.gradRasterTime;
TRmin = mr.calcDuration(rf) + mr.calcDuration(gzSSR)...
      + delayTE + mr.calcDuration(gxPre)...
      + mr.calcDuration(gx) + mr.calcDuration(gxSpoil);
delayTR = ceil((TR-TRmin)/seq.gradRasterTime)*seq.gradRasterTime;

%% Loop over phase encodes and define sequence blocks
% iZ < 0: Dummy shots to reach steady state
% iZ = 0: ADC is turned on and used for receive gain calibration on GE scanners
% iZ > 0: Image acquisition

nDummyZLoops = 4;
rf_phase = 0;
rf_inc = 0;

lastmsg = [];
for iZ = -nDummyZLoops:Nz_gre
    isDummyTR = iZ < 0;

    for ii = 1:length(lastmsg)
        fprintf('\b');
    end
    msg = sprintf('z encode %d of %d ', iZ, Nz_gre);
    fprintf(msg);
    lastmsg = msg;

    % Fat-sat
    TRID = 2 - isDummyTR;
    seq.addBlock(rfsat,mr.makeLabel('SET','TRID',TRID));
    seq.addBlock(gzSpoil);

    for iY = 1:Ny_gre
        % Turn on y and z prephasing lobes, except during dummy scans and
        % receive gain calibration (auto prescan)
        yStep = (iZ > 0) * pe1Steps(iY);
        zStep = (iZ > 0) * pe2Steps(max(1,iZ));

        % RF spoiling
        rf.phaseOffset = rf_phase/180*pi;
        adc.phaseOffset = rf_phase/180*pi;
        rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase = mod(rf_phase+rf_inc, 360.0);
        
        % Excitation
        % Mark start of segment (block group) by adding label.
        % Subsequent blocks in block group are NOT labelled.
        seq.addBlock(rf,gzSS);
        seq.addBlock(gzSSR);
        
        % Encoding
        seq.addBlock(mr.makeDelay(delayTE));
        seq.addBlock(gxPre, ...
            mr.scaleGrad(gyPre, yStep), ...
            mr.scaleGrad(gzPre, zStep));
        if isDummyTR
            seq.addBlock(gx);
        else
            seq.addBlock(gx, adc);
        end

        % rephasing/spoiling
        seq.addBlock(gxSpoil, ...
            mr.scaleGrad(gyPre, -yStep), ...
            mr.scaleGrad(gzPre, -zStep));
        seq.addBlock(mr.makeDelay(delayTR));
    end
end

%% Check sequence timing
[ok, error_report] = seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% Output for execution
seq.setDefinition('FOV', fov_gre);
seq.setDefinition('Name', 'gre');
seq.write('gre.seq');

%% Convert to .tar file for GE
seq2ge('gre.seq', sysGE, 'gre.tar');
system('tar -xvf gre.tar');

return;
%% Plot k-space trajectory
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
figure('WindowState','maximized');
plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('full k-space trajectory (k_x x k_y)');

return;
%% Optional slow step, but useful for testing during development,
% e.g., for the real TE, TR or for staying within slewrate limits
rep = seq.testReport;
fprintf([rep{:}]);
