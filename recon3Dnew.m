%% Reconstruction code using scan archives instead of p-files
% Requires GE Orchestra
% Rex Fung, June 17th, 2024

%% Load in params
setGREparams; setEPIparams;

% Params defined at scan
Nloops = 1; % What you set toppe CV 8 to when running the scan
Nframes = Nloops*NframesPerLoop;

% Filenames and options
fn_smaps = '/mnt/storage/rexfung/20240701epi_outside/gre.h5';
fn_cal = '/mnt/storage/rexfung/20240701epi_outside/cal.h5';
fn_loop = '/mnt/storage/rexfung/20240701epi_outside/loop.h5';
fn_adc = sprintf('adc/P%dadc.mod',Nx);
showEPIphaseDiff = true;
doSENSE = true;

%% Data loading
% Load raw data from scan archives (takes some time)
mode = 2;
switch mode
    case 1 % Depends on GE Orchestra
        ksp_raw_smaps = toppe.utils.loadsafile(fn_smaps,'acq_order',true);
        ksp_raw_cal = toppe.utils.loadsafile(fn_cal,'acq_order',true);
        ksp_raw = toppe.utils.loadsafile(fn_loop,'acq_order',true);
    case 2 % Independent of GE Orchestra
        ksp_raw_smaps = read_archive(fn_smaps);
        ksp_raw_cal = read_archive(fn_cal);
        ksp_raw = read_archive(fn_loop);

        % permute the data into the same form as case 1
        ksp_raw_smaps = permute(squeeze(ksp_raw_smaps),[2,4,1,3]);
        ksp_raw_cal = permute(squeeze(ksp_raw_cal),[2,4,1,3]);
        ksp_raw = permute(squeeze(ksp_raw),[2,4,1,3]);

        % discard leading empty data
        ksp_raw_smaps = ksp_raw_smaps(:,:,(2*size(ksp_raw_smaps,3) + 1):end); % and cal data
        ksp_raw_cal = ksp_raw_cal(:,:,(size(ksp_raw_cal,3) + 1):end);
        ksp_raw = ksp_raw(:,:,(size(ksp_raw,3) + 1):end);
end

%% Preprocessing
% Print max real and imag parts to check for reasonable magnitude
fprintf('Max real part: %d\n', max(real(ksp_raw(:))))
fprintf('Max imag part: %d\n', max(imag(ksp_raw(:))))

if doSENSE
    ksp_smaps = flip(ksp_raw_smaps, 1); % tv6 flips data along FID direction
    [Nfid,Ncoils,N] = size(ksp_smaps);
    ksp_smaps = ksp_smaps(:,:,1:Ny_gre*Nz_gre); % discard trailing data
    ksp_smaps = reshape(ksp_smaps,Nx_gre,Ncoils,Ny_gre,Nz_gre);
    ksp_smaps = permute(ksp_smaps,[1 3 4 2]); % [Nx Ny Nz Ncoils]
end

% Reshape and permute calibration data (a single frame w/out blips)
ksp_cal = flip(ksp_raw_cal, 1); % tv6 flips data along FID direction
[Nfid,Ncoils,N] = size(ksp_cal);
ksp_cal = ksp_cal(:,:,1:Ny*Nz); % discard trailing data
ksp_cal = reshape(ksp_cal,Nfid,Ncoils,Ny/Nsegments,Nsegments,Nz);
ksp_cal = permute(ksp_cal,[1 3 4 5 2]); % [Nfid Ny/Nsegments Nsegments Nz Ncoils]

% Reshape and permute loop data
ksp_rampsamp = flip(ksp_raw, 1); % tv6 flips data along FID direction
[Nfid,Ncoils,N] = size(ksp_rampsamp);
ksp_rampsamp = ksp_rampsamp(:,:,1:Ny*Nz*Nframes);
ksp_rampsamp = reshape(ksp_rampsamp,Nfid,Ncoils,Ny/Nsegments,Nsegments,Nz,Nframes);
ksp_rampsamp = permute(ksp_rampsamp,[1 3 4 5 6 2]); % [Nfid Ny/Nsegments Nsegments Nz Nframes Ncoils]

%% Prepare for reconstruction
% Estimate k-space center offset due to gradient delay
cal_data = reshape(abs(ksp_cal),Nfid,Ny/Nsegments,Nsegments*Nz*Ncoils);
cal_data(:,2:2:end) = flip(cal_data(:,2:2:end),1);
[M, I] = max(cal_data,[],1);
I = squeeze(I);
delay = mean(I,'all') - Nfid/2; delay = -2.5;
fprintf('Estimated offset from center of k-space (samples): %f\n', delay);

% retrieve sample locations from .mod file with adc info
[rf,gx,gy,gz,desc,paramsint16,pramsfloat,hdr] = toppe.readmod(fn_adc);
[kxo, kxe] = toppe.utils.getk(sysGE, fn_adc, Nfid, delay);

% Extract even number of lines (in case ETL is odd)
ETL_even = size(ksp_cal,2) - mod(size(ksp_cal,2),2);
oephase_data = ksp_cal(:,1:ETL_even,:,:,:);

% EPI ghost correction phase offset values
oephase_data = hmriutils.epi.rampsampepi2cart(oephase_data, kxo, kxe, Nx, fov(1)*100, 'nufft');
oephase_data = ifftshift(ifft(fftshift(reshape(oephase_data,Nx,ETL_even,Nsegments*Nz*Ncoils)),Nx,1));
[a, th] = hmriutils.epi.getoephase(oephase_data,showEPIphaseDiff);
fprintf('Constant phase offset (radians): %f\n', a(1));
fprintf('Linear term (radians/fov): %f\n', a(2));

%% Reconstruct multicoil k-space (each frame individually)

ksp_cart = zeros(Nx,Ny/Nsegments,Nsegments,Nz,Nframes,Ncoils); 
parfor frame = 1:Nframes
    fprintf('Reconstructing multicoil k-space of frame %d\n', round(frame))

    % Extract one frame
    ksp_rampsamp_frame = ksp_rampsamp(:,:,:,:,frame,:);
   
    % grid
    ksp_frame = hmriutils.epi.rampsampepi2cart(ksp_rampsamp_frame, kxo, kxe, Nx, fov(1)*100, 'nufft');

    % phase correct
    ksp_frame = hmriutils.epi.epiphasecorrect(ksp_frame, a);

    % Allocate
    ksp_cart(:,:,:,:,frame,:) = ksp_frame;
end

%% Rearrange segments/shots into proper locations
ksp_mc = zeros(Nx,Ny,Nz,Nframes,Ncoils);
for seg = 1:Nsegments
    tmp = squeeze(ksp_cart(:,:,seg,:,:,:));
    % tmp(:,1:2:end,:,:) = flip(tmp(:,1:2:end,:,:),1);
    ksp_mc(:,seg:Nsegments:end,:,:,:) = tmp;
end

%% Get sensitivity maps with PISCO
if doSENSE
    ksp_smaps = ifftshift(ifft(fftshift(ksp_smaps),Nz_gre,3));
    smaps = zeros(Nx_gre, Ny_gre, Nz_gre, Ncoils);
    for z = 1:Nz_gre
        [smaps(:,:,z,:), eigvals] = PISCO_senseMaps_estimation(squeeze(ksp_smaps(:,:,z,:)),...
                                    [Nx_gre, Ny_gre]);
    end

    % Crop in z to match EPI FoV
    z_start = round((fov_gre(3) - fov(3))/fov_gre(3)/2*Nz_gre + 1);
    z_end = round(Nz_gre - (fov_gre(3) - fov(3))/fov_gre(3)/2*Nz_gre);
    smaps = smaps(:,:,z_start:z_end,:);

    % Interpolate to match EPI data dimensions
    smaps_new = zeros(Nx,Ny,Nz,Ncoils);
    for coil = 1:Ncoils
        smaps_new(:,:,:,coil) = imresize3(smaps(:,:,:,coil),[Nx,Ny,Nz]);
    end
    smaps = smaps_new; clear smaps_new;

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
    img_final = sum(imgs_mc .* conj(smaps), 5);
else % root sum of squares combination
    img_final = sqrt(sum(abs(imgs_mc).^2, 5));
end

%% Compute k-space by IFT3
ksp_final = toppe.utils.ift3(img_final);

%% Viz
z = ceil(Nz/2);
figure('WindowState','maximized');
im('col',Nframes,'row',Nz,reshape(permute(img_final,[1 2 4 3]),Nx,Ny,Nframes*Nz),'cbar')
title(fn_loop(1:end-3));
ylabel('z direction'); xlabel('time')

% saveas(gcf, strcat('figs/',fn(1:end-2),'.png'));


%% Make GIF (using the gif add-on)
if false
    volumeTR = 1; % seconds
    clims = [min(img_final(:)), max(img_final(:))]; % Set the same dynamic range for each frame
    imgs_movie = permute(img_final,[2 1 3 4]);
    imgs_movie = flip(imgs_movie,2);
    
    close all;
    canvas = figure('WindowState','maximized');
    sgtitle(fn_loop(1:end-2));
    im('col',Nz,imgs_movie(:,:,:,1),clims,'cbar');
    xlabel('Frame 1');
    gif('movie.gif','DelayTime',volumeTR);
    for frame = 2:Nframes
        im('col',Nz,imgs_movie(:,:,:,frame),clims,'cbar');
        xlabel(sprintf('Frame %d',round(frame)));
        gif;
    end
end