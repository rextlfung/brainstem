% Pseudo-random 3D EPI sequence w/ z blips for CAIPI-like shifting
% Continuous readout gradient for speed.
% Create blocks without splitting blips.
% Represent readout gradients as arbitrary instead of trapezoids.
% CAIPI added to sample multiple kz locations per RF excitation
% Rex Fung

%% Definte experiment parameters
mode = 'rand_caipi';
setEPIparams;

%% Path and options
seqname = '3DEPI_loop_rs';
addpath('excitation');

%% Excitation pulse
% Target a slightly thinner slice to alleviate aliasing
[rf, gzSS, gzSSR] = mr.makeSincPulse(alpha/180*pi,...
                                     'duration',rfDur,...
                                     'sliceThickness',0.9*fov(3),...
                                     'system',sys);
gzSS = trap4ge(gzSS,CRT,sys);
gzSS.delay = rf.delay - gzSS.riseTime; % Sync up rf pulse and slice select gradient
gzSSR = trap4ge(gzSSR,CRT,sys);

%% Fat-sat
fatsat.flip    = 90;      % degrees
fatsat.slThick = 1e5;     % dummy value (determines slice-select gradient, but we won't use it; just needs to be large to reduce dead time before+after rf pulse)
fatsat.tbw     = 3;       % time-bandwidth product
fatsat.dur     = 6.0;     % pulse duration (ms)

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
% WITHOUT CAIPI shifts, which will be applied during sequence construction
omegas = zeros(Ny,Nz,NframesPerLoop);
z_shift_indices = zeros(round(Nz/Rz/caipi_z - Nz*acs(2)),NframesPerLoop);
for frame = 1:NframesPerLoop
    % Create pseudo-random 2D sampling mask. Save for recon
    rng(frame); % A different mask per frame
    [omega, nacs_indices_samp_z] = randsamp2d_caipi2(N(2:end), R, acs, max_ky_step, caipi_z);
    omegas(:,:,frame) = omega;
    z_shift_indices(:,frame) = nacs_indices_samp_z;
end

%% Define readout gradients and ADC event
% The Pulseq toolbox really shines here!

% Define k-space spacing for fully-sampled data
deltak = 1./fov;

% Start with the blips. Ensure long enough to support the largest blips
gyBlip = mr.scaleGrad(mr.makeTrapezoid('y', sys, 'area', max_ky_step*deltak(2)), 1/max_ky_step);
gyBlip = trap4ge(gyBlip,CRT,sys);
gzBlip = mr.scaleGrad(mr.makeTrapezoid('z', sys, 'area', (caipi_z - 1)*deltak(3)), 1/(caipi_z - 1));
gzBlip = trap4ge(gzBlip,CRT,sys);

% Area and duration of the biggest blip (in y or z)
if max_ky_step*deltak(2) > (caipi_z - 1)*deltak(3) % biggest blip in y
    maxBlipArea = max_ky_step*deltak(2);
    blipDuration = mr.calcDuration(gyBlip);

    % Remake the other blips to match duration
    gzBlip = mr.makeTrapezoid('z', sys, 'area', deltak(3), 'duration', blipDuration);
    gzBlip = trap4ge(gzBlip,CRT,sys);
else % biggest blip in z
    maxBlipArea = (caipi_z - 1)*deltak(3);
    blipDuration = mr.calcDuration(gzBlip);

    % Remake the other blips to match duration
    gyBlip = mr.makeTrapezoid('y', sys, 'area', deltak(2), 'duration', blipDuration);
    gyBlip = trap4ge(gyBlip,CRT,sys);
end

% Readout trapezoid
systmp = sys;
systmp.maxGrad = deltak(1)/dwell;  % to ensure Nyquist sampling
gro = trap4ge(mr.makeTrapezoid('x', systmp, 'Area', Nx*deltak(1) + maxBlipArea),CRT,sys);

% Circularly shift gro waveform to contain blips within each block
[gro1, gro2] = mr.splitGradientAt(gro, blipDuration/2);
gro2.delay = 0;
gro1.delay = gro2.shape_dur;
gro = mr.addGradients({gro2, mr.scaleGrad(gro1, -1)}, sys);
gro1.delay = 0; % This piece is necessary at the very beginning of the readout

% ADC event
Tread = mr.calcDuration(gro) - blipDuration;
Nfid = floor(Tread/dwell/4)*4;
adc = mr.makeAdc(Nfid, 'Dwell', dwell);

% Delay blips so they play after adc stops
gyBlip.delay = Tread;
gzBlip.delay = Tread;

% Prephasers (Make duration long enough to support all 3 directions)
gxPre = trap4ge(mr.makeTrapezoid('x',sys,'Area',-(Nx*deltak(1) + maxBlipArea)/2),CRT,sys);
gyPre = trap4ge(mr.makeTrapezoid('y',sys,'Area',-Ny/2*deltak(2)),CRT,sys);
gzPre = trap4ge(mr.makeTrapezoid('z',sys,'Area',-Nz/2*deltak(3)),CRT,sys);
Tpre = max([mr.calcDuration(gxPre),mr.calcDuration(gyPre),mr.calcDuration(gzPre)]);
gxPre = trap4ge(mr.makeTrapezoid('x',sys,'Area',-(Nx*deltak(1) + maxBlipArea)/2,'Duration',Tpre),CRT,sys);
gyPre = trap4ge(mr.makeTrapezoid('y',sys,'Area',-Ny/2*deltak(2),'Duration',Tpre),CRT,sys);
gzPre = trap4ge(mr.makeTrapezoid('z',sys,'Area',-Nz/2*deltak(3),'Duration',Tpre),CRT,sys);

% Spoilers (only in x and z because ?)
gxSpoil = trap4ge(mr.makeTrapezoid('x', sys, ...
    'Area', Nx*deltak(1)*NcyclesSpoil),CRT,sys);
gzSpoil = trap4ge(mr.makeTrapezoid('z', sys, ...
    'Area', Nz*deltak(3)*NcyclesSpoil),CRT,sys);

%% Calculate delay to achieve desired TE
minTE = 0.5*mr.calcDuration(rf)...
      + mr.calcDuration(gzSSR)...
      + max([mr.calcDuration(gxPre), mr.calcDuration(gyPre), mr.calcDuration(gzPre)])...
      + (2*ceil(Ny/Ry/2)/2 - 0.5) * mr.calcDuration(gro);
if TE >= minTE
    TEdelay = floor((TE - minTE)/sys.blockDurationRaster) * sys.blockDurationRaster;
else
    warning(sprintf('Minimum achievable TE (%d) exceeds prescribed TE (%d)',...
                    minTE, TE))
    TEdelay = 0;
end

%% Calculate delay to achieve desired TR
minTR = mr.calcDuration(rfsat)...
      + max([mr.calcDuration(gxSpoil), mr.calcDuration(gzSpoil)])...
      + max([mr.calcDuration(rf), mr.calcDuration(gzSS)])...
      + mr.calcDuration(gzSSR)...
      + TEdelay...
      + max([mr.calcDuration(gxPre), mr.calcDuration(gyPre), mr.calcDuration(gzPre)])...
      + 2*ceil(Ny/Ry/2) * mr.calcDuration(gro)...
      + mr.calcDuration(gzPre)...
      + max([mr.calcDuration(gxSpoil), mr.calcDuration(gzSpoil)]);
if TR >= minTR
    TRdelay = floor((TR - minTR)/sys.blockDurationRaster)*sys.blockDurationRaster;
else
    warning(sprintf('Minimum achievable TR (%d) exceeds prescribed TR (%d)',...
                    minTR, TR))
    TRdelay = 0;
end

%% Assemble sequence
seq = mr.Sequence(sys);

% log the sequence of k-space locations sampled (ky and kz)
samp_log = zeros(NframesPerLoop, ...
                 round(Nz/caipi_z/Rz)*2*ceil(Ny/Ry/2), ...
                 2);

% RF spoiling trackers
rf_count = 1;
rf_phase = rf_phase_0;

for frame = 1:NframesPerLoop
    fprintf('Writing frame %d\n', frame)
    % Load in kz-ky sampling mask (unshifted), and indices to be shifted
    omega = omegas(:,:,frame);
    nacs_indices_samp_z = z_shift_indices(:,frame);

    % Reset sample count
    samp_count = 1;

    % kz encoding loop
    % Each "z_loc" is the starting point of a partition of kz locations
    z_locs = find(sum(omega,1));
    for z = z_locs
        % Label the first block in each "unique" section with TRID (see Pulseq on GE manual)
        TRID = 1;

        % Fat-sat
        seq.addBlock(rfsat,mr.makeLabel('SET','TRID',TRID));
        seq.addBlock(gxSpoil, gzSpoil);

        % RF spoiling
        rf_phase = mod(0.5 * rf_phase_0 * rf_count^2, 360.0);
        rf.phaseOffset = rf_phase/180*pi;
        adc.phaseOffset = rf_phase/180*pi;
        rf_count = rf_count + 1;

        % Slab-selective RF excitation + rephase
        seq.addBlock(rf,gzSS);
        seq.addBlock(gzSSR);

        % TE delay
        if TE > minTE
            seq.addBlock(mr.makeDelay(TEdelay));
        end
        
        % Begin ky encoding
        y_locs = find(omega(:,z));

        % Randomly generate order of caipi shifts for this echo train
        z_shifts = zeros(1,length(y_locs));
        if ismember(z, nacs_indices_samp_z)
            for iy = 1:caipi_z:length(y_locs)
                z_shifts(iy:(iy + caipi_z - 1)) = randperm(caipi_z) - ceil(caipi_z/2);
            end
            z_shifts = z_shifts(1:length(y_locs));
        end

        % Move to corner of k-space
        gzPreTmp = mr.scaleGrad(gzPre, (z + z_shifts(1) - Nz/2 - 1)/(-Nz/2));
        gyPreTmp = mr.scaleGrad(gyPre, (y_locs(1) - Ny/2 - 1)/(-Ny/2));
        seq.addBlock(gxPre, gyPreTmp, gzPreTmp);
    
        % Zip through k-space with EPI trajectory
        seq.addBlock(gro1);
        for iy = 1:(length(y_locs) - 1)
            % Log sampling locations
            samp_log(frame,samp_count,:) = [y_locs(iy); z + z_shifts(iy)];
            samp_count = samp_count + 1;

            % Sample
            seq.addBlock(adc, mr.scaleGrad(gro, (-1)^(iy-1)),...
                mr.scaleGrad(gyBlip, y_locs(iy + 1) - y_locs(iy)),...
                mr.scaleGrad(gzBlip, z_shifts(iy + 1) - z_shifts(iy))...
                );
        end

        % Last line
        % Log sampling locations
        samp_log(frame,samp_count,:) = [y_locs(end); z + z_shifts(end)];
        samp_count = samp_count + 1;

        % Sample
        seq.addBlock(adc, mr.scaleGrad(gro2, (-1)^iy));

        % End ky encoding

        % rephase kz encoding before spoiling
        seq.addBlock(mr.scaleGrad(gzPreTmp, -1));

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

%% Save sampling log for recon
save('samp_logs/samp_log.mat','samp_log','-v7.3');

%% Output for execution
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', seqname);
seq.write(strcat(seqname, '.seq'));

%% GE stuff
seq2ge(strcat(seqname, '.seq'), sysGE, strcat(seqname, '.tar'))
system(sprintf('tar -xvf %s', strcat(seqname, '.tar')));

%% Plot sequence
figure();
toppe.plotseq(sysGE, 'timeRange',[0, 11*TR]);
fontsize(16,'points');

return;
%% Detailed checks that takes some time to run

% k-space trajectory calculation and plot
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
figure('WindowState','maximized');
plot(ktraj(2,:),ktraj(3,:),'b'); % a 2D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold; plot(ktraj_adc(2,:),ktraj_adc(3,:),'r.'); % plot the sampling points
title('full k-space trajectory (k_y x k_z)');
xlabel('k_y'); ylabel('k_z');

return;
%% Optional slow step, but useful for testing during development,
% e.g., for the real TE, TR or for staying within slewrate limits
rep = seq.testReport;
fprintf([rep{[1:9, 11:end]}]); % print report except warnings