% Code for loading and preprocessing randomly undersampled 3D EPI data
% 1. Reads in scan archive files based on 
% Rex Fung, June 17th, 2024

%% Set params and options
% Load in params
setGREparams; setEPIparams;

% Params defined at scan
Nloops = 10; Nframes = Nloops*NframesPerLoop;

% Filenames
datdir = '/mnt/storage/rexfung/20240826randEPI/';
fn_gre = strcat(datdir,'gre.h5');
fn_cal = strcat(datdir,'cal.h5');
fn_loop = strcat(datdir,'loop.h5');
fn_adc = sprintf('adc/adc%d.mod',Nx);

% Options
doSENSE = true;
showEPIphaseDiff = true;

%% Load data
% Load raw data from scan archives (takes some time)
ksp_gre_raw = read_archive(fn_gre);
ksp_cal_raw = read_archive(fn_cal);
ksp_loop_raw = read_archive(fn_loop);
[Nfid, Ncoils] = size(ksp_loop_raw,2,5);

% permute the data into the same form as case 1
ksp_gre_raw = permute(squeeze(ksp_gre_raw),[2,4,1,3]);
ksp_cal_raw = permute(squeeze(ksp_cal_raw),[2,4,1,3]);
ksp_loop_raw = permute(squeeze(ksp_loop_raw),[2,4,1,3]);

% discard leading empty data
ksp_gre_raw = ksp_gre_raw(:,:,(2*size(ksp_gre_raw,3) + 1):end); % and cal data
ksp_cal_raw = ksp_cal_raw(:,:,(size(ksp_cal_raw,3) + 1):end);
ksp_loop_raw = ksp_loop_raw(:,:,(size(ksp_loop_raw,3) + 1):end);

% Print max real and imag parts to check for reasonable magnitude
fprintf('Max real part of gre data: %d\n', max(real(ksp_gre_raw(:))))
fprintf('Max imag part of gre data: %d\n', max(imag(ksp_gre_raw(:))))
fprintf('Max real part of cal data: %d\n', max(real(ksp_cal_raw(:))))
fprintf('Max imag part of cal data: %d\n', max(imag(ksp_cal_raw(:))))
fprintf('Max real part of loop data: %d\n', max(real(ksp_loop_raw(:))))
fprintf('Max imag part of loop data: %d\n', max(imag(ksp_loop_raw(:))))

%% Get sensitivity maps with PISCO
if doSENSE
    % Reshape and permute gre data
    ksp_gre = ksp_gre_raw(:,:,1:Ny_gre*Nz_gre); % discard trailing data
    ksp_gre = reshape(ksp_gre,Nx_gre,Ncoils,Ny_gre,Nz_gre);
    ksp_gre = permute(ksp_gre,[1 3 4 2]); % [Nx Ny Nz Ncoils]

    % Only use central (high SNR) region of k-space to estimate smaps
    cal_length = 32; % Length of each dimension of the calibration data
    center_x = ceil(Nx_gre/2) + ~rem(Nx_gre,2);
    center_y = ceil(Ny_gre/2) + ~rem(Ny_gre,2);
    cal_index_x = center_x + (-floor(cal_length/2):floor(cal_length/2) - ~rem(cal_length/2,2));
    cal_index_y = center_y + (-floor(cal_length/2):floor(cal_length/2) - ~rem(cal_length/2,2));

    % Compute smaps slice-by-slice
    tmp = ifftshift(ifft(fftshift(ksp_gre),Nz_gre,3));
    smaps = zeros(Nx_gre, Ny_gre, Nz_gre, Ncoils);
    for z = 1:Nz_gre
        [smaps_tmp, eigvals] = PISCO_senseMaps_estimation(...
                                        squeeze(tmp(cal_index_x,cal_index_y,z,:)),...
                                        [Nx_gre, Ny_gre]...
                                    );

        % Normalize
        smaps_tmp = smaps_tmp./sqrt(sum(abs(smaps_tmp).^2, 3));

        % Allocate
        smaps(:,:,z,:) = smaps_tmp;
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
end

%% Compute odd/even delays using calibration (blipless) data
% Reshape and permute calibration data (a single frame w/out blips)
ksp_cal = ksp_cal_raw(:,:,1:ceil(Ny/Ry)*ceil(Nz/Rz)); % discard trailing data
ksp_cal = reshape(ksp_cal,Nfid,Ncoils,ceil(Ny/Ry), ceil(Nz/Rz));
ksp_cal = permute(ksp_cal,[1 3 4 2]); % [Nfid Ny/Ry Nz/Rz Ncoils]

% Estimate k-space center offset due to gradient delay
cal_data = reshape(abs(ksp_cal),Nfid,ceil(Ny/Ry),ceil(Nz/Rz)*Ncoils);
cal_data(:,2:2:end,:) = flip(cal_data(:,2:2:end,:),1);
[M, I] = max(cal_data,[],1);
delay = 2*mean(I,'all') - Nfid;
delay = -delay; %
fprintf('Estimated offset from center of k-space (samples): %f\n', delay);

% retrieve sample locations from .mod file with adc info
[rf,gx,gy,gz,desc,paramsint16,pramsfloat,hdr] = toppe.readmod(fn_adc);
[kxo, kxe] = toppe.utils.getk(sysGE, fn_adc, Nfid, delay);

% Extract even number of lines (in case ETL is odd)
ETL_even = size(ksp_cal,2) - mod(size(ksp_cal,2),2);
oephase_data = ksp_cal(:,1:ETL_even,:,:,:);

% EPI ghost correction phase offset values
oephase_data = hmriutils.epi.rampsampepi2cart(oephase_data, kxo, kxe, Nx, fov(1)*100, 'nufft');
oephase_data = ifftshift(ifft(fftshift(reshape(oephase_data,Nx,ETL_even,ceil(Nz/Rz)*Ncoils)),Nx,1));
[a, th] = hmriutils.epi.getoephase(oephase_data,showEPIphaseDiff);
fprintf('Constant phase offset (radians): %f\n', a(1));
fprintf('Linear term (radians/fov): %f\n', a(2));

%% Grid and apply odd/even correction to loop data
% Reshape and permute loop data
ksp_loop = ksp_loop_raw(:,:,1:ceil(Ny/Ry)*ceil(Nz/Rz)*Nframes); % discard trailing empty data
ksp_loop = reshape(ksp_loop,Nfid,Ncoils,ceil(Ny/Ry),ceil(Nz/Rz),Nframes);
ksp_loop = permute(ksp_loop,[1 3 4 5 2]); % [Nfid Ny/Ry Nz/Rz Nframes Ncoils]

% Grid along kx direction via NUFFT (takes a while)
ksp_loop_cart = hmriutils.epi.rampsampepi2cart(ksp_loop, kxo, kxe, Nx, fov(1)*100, 'nufft');

% Phase correct along kx direction
ksp_loop_cart = hmriutils.epi.epiphasecorrect(ksp_loop_cart, a);

%% Create zero-filled k-space data
ksp_zf = zeros(Nx,Ny,Nz,Nframes,Ncoils);

% (Re)generate temporally incoherent sampling masks (or load it in)
omegas = false(Ny,Nz,Nframes);
for frame = 1:Nframes
    % Create pseudo-random 2D sampling mask. Save for recon
    rng(frame); % A different mask per frame
    omega = randsamp2d(N(2:end), R, acs);
    omegas(:,:,frame) = omega;
end

% Allocate data
ksp_zf(:,repmat(omegas,1,1,1,Ncoils)) = ksp_loop_cart(:,:);

%% Quick check with RSOS recon
imgs = zeros(Nx,Ny,Nz,Nframes);
for frame = 1:Nframes
    imgs(:,:,:,frame) = sqrt(sum(abs(toppe.utils.ift3(squeeze(ksp_zf(:,:,:,frame,:)))).^2,4));
end

return;

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
    img_final = sum(imgs_mc .* conj(reshape(smaps,[Nx Ny Nz 1 Ncoils])), 5);
else % root sum of squares combination
    img_final = sqrt(sum(abs(imgs_mc).^2, 5));
end

%% Compute k-space by IFT3
ksp_final = toppe.utils.ift3(img_final);

%% Compute GRE images
img_gre = sqrt(sum(abs(toppe.utils.ift3(ksp_gre)).^2,4));

%% Viz
close all;

% Final images
figure('WindowState','maximized');
im('col',Nframes,'row',Nz,reshape(permute(img_final,[1 2 4 3]),Nx,Ny,Nframes*Nz),'cbar')
title(fn_loop(1:end-3));
ylabel('z direction'); xlabel('time')

% Sensitivity maps
if doSENSE
    figure('WindowState','maximized');
    im('col',Ncoils,'row',Nz,reshape(permute(squeeze(smaps),[1 2 4 3]),Nx,Ny,Ncoils*Nz),'cbar')
    title('Sensitivity maps');
    ylabel('z direction'); xlabel('Coils');
end

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