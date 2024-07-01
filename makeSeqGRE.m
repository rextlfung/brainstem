% writeB0.m
%
% 3D GRE B0 mapping demo sequence for Pulseq on GE v1.0 User Guide
%
% Modified June 21, 2024, for acquiring GRE data for sensitivity maps

% Define params (edit setGREparams to change)
setGREparams;

% Create a new sequence object
seq = mr.Sequence(sys);           

% Create non-selective pulse
[rf] = mr.makeBlockPulse(alpha/180*pi, sys, 'Duration', alphaPulseDuration);

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
    'Area', Nz_gre/2*deltak(3), ...   % PE2 gradient, max amplitude
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
TEmin = rf.shape_dur/2 + rf.ringdownTime + mr.calcDuration(gxPre) ...
      + adc.delay + Nx_gre/2*dwell;
delayTE = ceil((TE-TEmin)/seq.gradRasterTime)*seq.gradRasterTime;
TRmin = mr.calcDuration(rf) + delayTE + mr.calcDuration(gxPre) ...
      + mr.calcDuration(gx) + mr.calcDuration(gxSpoil);
delayTR = ceil((TR-TRmin)/seq.gradRasterTime)*seq.gradRasterTime;

%% Loop over phase encodes and define sequence blocks
% iZ < 0: Dummy shots to reach steady state
% iZ = 0: ADC is turned on and used for receive gain calibration on GE scanners
% iZ > 0: Image acquisition

nDummyZLoops = 1;
rf_phase = 0;
rf_inc = 0;

lastmsg = [];
for iZ = -nDummyZLoops:Nz_gre
    isDummyTR = iZ < 0;

    for ii = 1:length(lastmsg)
        fprintf('\b');
    end
    msg = sprintf('z encode %d of %d   ', iZ, Nz_gre);
    fprintf(msg);
    lastmsg = msg;

    for iY = 1:Ny_gre
        % Turn on y and z prephasing lobes, except during dummy scans and
        % receive gain calibration (auto prescan)
        yStep = (iZ > 0) * pe1Steps(iY);
        zStep = (iZ > 0) * pe2Steps(max(1,iZ));

        for c = 1:length(TE)

            % RF spoiling
            rf.phaseOffset = rf_phase/180*pi;
            adc.phaseOffset = rf_phase/180*pi;
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
            
            % Excitation
            % Mark start of segment (block group) by adding label.
            % Subsequent blocks in block group are NOT labelled.
            seq.addBlock(rf, mr.makeLabel('SET', 'TRID', 2-isDummyTR));
            
            % Encoding
            seq.addBlock(mr.makeDelay(delayTE(c)));
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
            seq.addBlock(mr.makeDelay(delayTR(c)));
        end
    end
end
fprintf('Sequence ready\n');

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
seq.setDefinition('FOV', fov_gre);
seq.setDefinition('Name', 'b0');
seq.write('b0.seq');

%% Plot sequence
Noffset = length(TE)*Ny_gre*(nDummyZLoops+1);
% seq.plot('timerange',[Noffset Noffset+4]*TR(1), 'timedisp', 'ms');

%% Convert to .tar file for GE
toGE = true;
if toGE
    sysGE = toppe.systemspecs('maxGrad', 5, ...   % G/cm
        'maxSlew', 20, ...               % G/cm/ms
        'maxRF', 0.05, ...               % Gauss. Must be >= peak RF in sequence.
        'maxView', Ny_gre, ...               % Determines slice/view index in data file
        'rfDeadTime', 100, ...           % us
        'rfRingdownTime', 60, ...        % us
        'adcDeadTime', 40, ...           % us
        'psd_rf_wait', 148, ...          % RF/gradient delay (us)
        'psd_grd_wait', 156);            % ADC/gradient delay (us)

    seq2ge('b0.seq', sysGE, 'b0.tar');
    system('tar -xvf b0.tar');
end

%% Plot
figure('WindowState','maximized');
toppe.plotseq(sysGE);