% Pseudo-random 3D EPI sequence
% Pauses readout gradient during blips for pulseq memory efficiency
% Continuous readout gradient 
% Rex Fung
% Last modified August 26th, 2024

%% Definte experiment parameters
setEPIparams;

%% Path and options
seqname = '3DEPI_loop_rs';
addpath('excitation');

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
    'type', 'ex', ... % fatsat pulse is a 90 so is of type 'ex', not 'st' (small-tip)
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
rfsat.freqOffset = -fatOffresFreq; % Hz

%% Generate temporally incoherent sampling masks
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
gzBlip = trap4ge(mr.makeTrapezoid('z', sys, 'Area', 0*deltak(3)),CRT,sys); % unusued rn

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
gro = trap4ge(mr.makeTrapezoid('x', systmp, 'Area', Nx*deltak(1)),CRT,sys);

% ADC event
Tread = mr.calcDuration(gro);
Tread = floor(Tread/CRT/sys.adcSamplesDivisor)*CRT*sys.adcSamplesDivisor;
adc = mr.makeAdc(round(Tread/dwell), sys, ...
    'Duration', Tread, ...
    'Delay', sysGE.adcDeadTime*1e-6);

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
      + (Ny/2/Ry - 1) * (mr.calcDuration(gro) + blipDuration)...
      + 0.5*mr.calcDuration(gro);
if TE >= minTE
    TEdelay = floor((TE - minTE)/sys.blockDurationRaster) * sys.blockDurationRaster;
else
    warning(sprintf('Minimum achievable TE (%d) exceeds prescribed TE (%d)',...
                    minTE, TE))
    TEdelay = 0;
end

%% Calculate delay to achieve desired TR
minTR = mr.calcDuration(rfsat) + mr.calcDuration(gzSpoil)...
      + mr.calcDuration(gzSS) + mr.calcDuration(gzSSR)...
      + TEdelay...
      + mr.calcDuration(gzPre)...
      + (Ny/Ry - 1)*(mr.calcDuration(gro) + blipDuration)...
      + mr.calcDuration(gro)...
      + mr.calcDuration(gzSpoil);
if TR >= minTR
    TRdelay = floor((TR - minTR)/sys.blockDurationRaster)*sys.blockDurationRaster;
else
    warning(sprintf('Minimum achievable TR (%d) exceeds prescribed TR (%d)',...
                    minTR, TR))
    TRdelay = 0;
end
%% Assemble sequence
seq = mr.Sequence(sys);           

% RF spoiling trackers
rf_phase = 0;
rf_inc = 0;

for frame = 1:Nframes
    % Load in kz-ky sampling mask
    omega = omegas(:,:,frame);

    % kz encoding loop
    z_locs = find(sum(omega,1));
    for iz = 1:length(z_locs)
        gzPreTmp = mr.scaleGrad(gzPre,-(z_locs(iz) - floor(Nz/2))/(Nz/2));
    
        % Label the first block in each "unique" section with TRID (see Pulseq on GE manual)
        TRID = 1;

        % Fat-sat
        seq.addBlock(rfsat,mr.makeLabel('SET','TRID',TRID));
        seq.addBlock(gxSpoil, gzSpoil);

        % RF spoiling
        rf.phaseOffset = rf_phase/180*pi;
        adc.phaseOffset = rf_phase/180*pi;
        rf_inc = mod(rf_inc + rfSpoilingInc, 360.0);
        rf_phase = mod(rf_phase + rf_inc, 360.0);

        % Slab-selective RF excitation + rephase
        seq.addBlock(rf,gzSS);
        seq.addBlock(gzSSR);

        % TE delay
        if TE > minTE
            seq.addBlock(mr.makeDelay(TEdelay));
        end

        % ky encoding
        y_locs = find(omega(:,z_locs(iz)));

            % Move to corner of k-space and sample the first line
            gyPreTmp = mr.scaleGrad(gyPre, (gyBlip.area*(y_locs(1) - 1) + gyPre.area)/gyPre.area);
            seq.addBlock(gxPre, gyPreTmp, gzPreTmp);
    
            % Zip through k-space with EPI trajectory
            for iy = 1:(length(y_locs) - 1)
                seq.addBlock(adc, mr.scaleGrad(gro, (-1)^(iy-1)));
                seq.addBlock(mr.scaleGrad(gyBlip, y_locs(iy + 1) - y_locs(iy)));
            end

            % Last line
            seq.addBlock(adc, mr.scaleGrad(gro, (-1)^iy));
        % end ky encoding

        % spoil
        seq.addBlock(gxSpoil, gzSpoil);

        % Achieve desired TR
        if TR > minTR
            seq.addBlock(mr.makeDelay(TRdelay));
        end
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
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', seqname);
seq.write(strcat(seqname, '.seq'));

%% GE stuff
seq2ge(strcat(seqname, '.seq'), sysGE, strcat(seqname, '.tar'))
system(sprintf('tar -xvf %s', strcat(seqname, '.tar')));

%% Plot
figure('WindowState','maximized');
toppe.plotseq(sysGE, 'timeRange',[0, TR]);

%% Detailed checks that takes some time to run
return;

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
