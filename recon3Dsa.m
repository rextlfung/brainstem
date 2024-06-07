%% Reconstruction code using scan archives instead of p-files
% Requires GE Orchestra
% Rex Fung, June 6th, 2024

%% User-defined values (CHANGE THESE)
% Scanner system info
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

% Redefine parameters for convenience
fov = [200, 200, 3]*1e-3;
Nx = 200; Ny = Nx; Nz = 3;
Nsegments = 4;
Nframes = 120;
Ncoils = 32;

% Filenames and options
fn_cal = '/mnt/storage/rexfung/20240604brainstem3DEPI/scanArchives/3DEPI_cal/data.h5';
fn_loop = '/mnt/storage/rexfung/20240604brainstem3DEPI/scanArchives/3DEPI_loop_v1/data.h5';
fn_adc = sprintf('adc/P%dadc.mod', Nx);
doSENSE = false;

%% Data loading
% Load raw data from scan archives (takes some time)
mode = 1;
switch mode
    case 1
        ksp_raw_cal = toppe.utils.loadsafile(fn_cal,'acq_order',true);
        ksp_raw = toppe.utils.loadsafile(fn_loop,'acq_order',true);
    case 2
        ksp_raw_cal = read_arc(fn_cal);
        ksp_raw = read_arc(fn_loop);
end
%% Preprocessing
% Print max real and imag parts to check for reasonable magnitude
fprintf('Max real part: %d\n', max(real(ksp_raw(:))))
fprintf('Max imag part: %d\n', max(imag(ksp_raw(:))))

% Reshape and permute calibration data (a single frame w/out blips)
ksp_cal = flip(ksp_raw_cal, 1); % tv6 flips data along FID direction
[Nfid,Ncoils,N] = size(ksp_cal);
ksp_cal = ksp_cal(:,:,1:Ny*Nz); % discard trailing data
ksp_cal = reshape(ksp_cal, [Nfid, Ncoils, Ny, Nz]);
ksp_cal = permute(ksp_cal, [1 3 4 2]); % [Nfid Ny Nz Ncoils]

% Reshape and permute loop data
ksp_rampsamp = flip(ksp_raw, 1); % tv6 flips data along FID direction
[Nfid,Ncoils,N] = size(ksp_rampsamp);
ksp_rampsamp = ksp_rampsamp(:,:,1:Ny*Nz*Nframes);
ksp_rampsamp = reshape(ksp_rampsamp, Nfid, Ncoils, Ny, Nz, Nframes);
ksp_rampsamp = permute(ksp_rampsamp, [1 3 4 5 2]); % [Nfid Ny Nz Nframes Ncoils]

% Rearrange segments/shots into proper locations
ksp_cal_reordered = zeros(size(ksp_cal));
ksp_rampsamp_reordered = zeros(size(ksp_rampsamp));
for seg = 1:Nsegments
    range = (1:Ny/Nsegments) + (seg-1)*Ny/Nsegments;

    tmp1 = ksp_cal(:,range,:,:,:);
    tmp1(:,1:2:end,:,:,:) = flip(tmp1(:,1:2:end,:,:,:),1);
    ksp_cal_reordered(:,seg:Nsegments:end,:,:,:) = tmp1;

    tmp2 = ksp_rampsamp(:,range,:,:,:);
    tmp2(:,1:2:end,:,:,:) = flip(tmp2(:,1:2:end,:,:,:),1);
    ksp_rampsamp_reordered(:,seg:Nsegments:end,:,:,:) = tmp2;
end
ksp_cal = ksp_cal_reordered;
ksp_rampsamp = ksp_rampsamp_reordered;
clear tmp1 tmp2 ksp_cal_reordered ksp_rampsamp_reordered;

% Combine y and z dimensions as blips are off for cal data
ksp_cal = reshape(ksp_cal, [Nfid, Ny, Nz*Ncoils]);

% Naively recon each frame individually
if true
    imgs_naive = zeros(Nfid, Ny, Nz, Nframes);
    parfor frame = 1:Nframes
        fprintf('Naively reconstructing frame %d\n',round(frame))
        ksp_frame = squeeze(ksp_rampsamp(:,:,:,frame,:));
        img_frame_mc = toppe.utils.ift3(ksp_frame);
        imgs_naive(:,:,:,frame) = sqrt(sum(abs(img_frame_mc).^2,4));
    end
end

%% Prepare for reconstruction
% Estimate k-space center offset due to gradient delay using the first halfo
[M, I] = max(abs(ksp_cal),[],1);
I = squeeze(I);
delay = mean(I(1:Ny/2,:),'all') - Nfid/2; % use the 1st half of each shot
fprintf('Estimated offset from center of k-space (samples): %f\n', delay);

% retrieve sample locations from .mod file with adc info
[rf,gx,gy,gz,desc,paramsint16,pramsfloat,hdr] = toppe.readmod(fn_adc);
[kxo, kxe] = toppe.utils.getk(sysGE, fn_adc, Nfid, delay);

% EPI ghost correction phase offset values
oephase_data = ksp_cal;
oephase_data(:,1:2:end,:,:) = flip(oephase_data(:,1:2:end,:,:),1); % Simulate 1-shot EPI
oephase_data = hmriutils.epi.rampsampepi2cart(oephase_data, kxo, kxe, Nx, fov(1)*100, 'nufft');
oephase_data = ifftshift(ifft(fftshift(reshape(oephase_data, [Nx, Ny, Nz*Ncoils])),Nx,1));
[a, th] = hmriutils.epi.getoephase(oephase_data,true);
fprintf('Constant phase offset (radians): %f\n', a(1));
fprintf('Linear term (radians/fov): %f\n', a(2));

%% Reconstruct multicoil k-space (each frame individually)

ksp_mc = zeros(Nx, Ny, Nz, Nframes, Ncoils);
parfor frame = 1:Nframes
    fprintf('Reconstructing multicoil k-space of frame %d\n', round(frame))

    % Extract one frame and rearrange to simulate 1-shot EPI 
    ksp_rampsamp_frame = ksp_rampsamp(:,:,:,frame,:);
    ksp_rampsamp_frame(:,1:2:end,:,:) = flip(ksp_rampsamp_frame(:,1:2:end,:,:),1);
   
    % grid
    ksp_frame = hmriutils.epi.rampsampepi2cart(ksp_rampsamp_frame, kxo, kxe, Nx, fov(1)*100, 'nufft');

    % phase correct
    ksp_frame = hmriutils.epi.epiphasecorrect(ksp_frame, a);

    % Allocate
    ksp_mc(:,:,:,frame,:) = ksp_frame;
end

%% Get sensitivity maps with PISCO
if doSENSE
    addpath('pisco/');
    ksp4maps = squeeze(ksp(:,:,:,2,:));
    ksp4maps = ifftshift(ifft(fftshift(ksp4maps),Nz,3));
    smaps = zeros(Nx, Ny, Nz, Ncoils);
    for z = 1:Nz
        [smaps(:,:,z,:), eigvals] = PISCO_senseMaps_estimation(squeeze(ksp4maps(:,:,z,:)), [Nx, Ny]);
    end

    % Add a time dimension so it can be multiplied
    smaps = reshape(smaps,[Nx Ny Nz 1 Ncoils]);
end

%% IFFT to get images
imgs_mc = ifftshift(ifft(...
                     ifft(...
                      ifft(...
                       fftshift(ksp_mc)...
                       , Nx, 1)...
                      , Ny, 2)...
                     , Nz, 3)...
                    );

%% Coil combination
if doSENSE
    imgs = sum(imgs_mc .* conj(smaps), 5);
else % root sum of squares combination
    imgs = sqrt(sum(abs(imgs_mc).^2, 5));
end

%% Viz
z = 2;
figure; im(imgs(:,:,z,:),'cbar')
title(fn_loop(1:end-2));
% saveas(gcf, strcat('figs/',fn(1:end-2),'.png'));
