% Code for loading and preprocessing randomly undersampled 3D EPI data
% 1. Reads in scan archive files
% Rex Fung, June 17th, 2024

%% Set params and options
% Load in params
setGREparams; setEPIparams;

% Total number of frames
Nloops = 15; % Defined as toppe cv 8 at scanner
Nframes = Nloops*NframesPerLoop;

% Filenames
datdir = '/mnt/storage/rexfung/20240912randEPIbrain3x3/';
fn_gre = strcat(datdir,'gre.h5');
fn_cal = strcat(datdir,'cal.h5');
fn_loop = strcat(datdir,'loop.h5');

% Options
doSENSE = true;
showEPIphaseDiff = true;

%% Load data
% Load raw data from scan archives (takes some time)
ksp_gre_raw = read_archive(fn_gre);
ksp_cal_raw = read_archive(fn_cal);
ksp_loop_raw = read_archive(fn_loop);
[Nfid, Ncoils] = size(ksp_cal_raw,2,5);

% permute the data into the same form as case 1
ksp_gre_raw = permute(squeeze(ksp_gre_raw),[2,4,1,3]);
ksp_cal_raw = permute(squeeze(ksp_cal_raw),[2,4,1,3]);
ksp_loop_raw = permute(squeeze(ksp_loop_raw),[2,4,1,3]);

% discard leading empty data
ksp_gre_raw = ksp_gre_raw(:,:,(size(ksp_gre_raw,3) + 1):end);
ksp_cal_raw = ksp_cal_raw(:,:,(size(ksp_cal_raw,3) + 1):end);
ksp_loop_raw = ksp_loop_raw(:,:,(size(ksp_loop_raw,3) + 1):end);

% discard calibration portion for gre data
ksp_gre_raw = ksp_gre_raw(:,:,(Ny_gre + 1):end);

% Print max real and imag parts to check for reasonable magnitude
fprintf('Max real part of gre data: %d\n', max(real(ksp_gre_raw(:))))
fprintf('Max imag part of gre data: %d\n', max(imag(ksp_gre_raw(:))))
fprintf('Max real part of cal data: %d\n', max(real(ksp_cal_raw(:))))
fprintf('Max imag part of cal data: %d\n', max(imag(ksp_cal_raw(:))))
fprintf('Max real part of loop data: %d\n', max(real(ksp_loop_raw(:))))
fprintf('Max imag part of loop data: %d\n', max(imag(ksp_loop_raw(:))))

%% Compute odd/even delays using calibration (blipless) data
% Reshape and permute calibration data (a single frame w/out blips)
ksp_cal = ksp_cal_raw(:,:,1:ceil(Ny/Ry)*ceil(Nz/Rz)); % discard trailing data
ksp_cal = reshape(ksp_cal,Nfid,Ncoils,ceil(Ny/Ry), ceil(Nz/Rz));
ksp_cal = permute(ksp_cal,[1 3 4 2]); % [Nfid Ny/Ry Nz/Rz Ncoils]

% Estimate k-space center offset due to gradient delay
cal_data = reshape(abs(ksp_cal),Nfid,ceil(Ny/Ry),ceil(Nz/Rz),Ncoils);
cal_data(:,1:2:end,:,:) = flip(cal_data(:,1:2:end,:,:),1);
[M, I] = max(cal_data,[],1);
delay = Nfid - 2*mean(I,'all'); delay = -1.5;
fprintf('Estimated offset from center of k-space (samples): %f\n', delay);

% retrieve sample locations from .mod file with adc info
fn_adc = sprintf('adc/adc%d.mod',Nfid);
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
ksp_loop = permute(ksp_loop,[1 3 4 2 5]); % [Nfid Ny/Ry Nz/Rz Ncoils Nframes]

% Grid along kx direction via NUFFT (takes a while)
tic
    ksp_loop_cart = hmriutils.epi.rampsampepi2cart(ksp_loop, kxo, kxe, Nx, fov(1)*100, 'nufft');
toc

% Phase correct along kx direction
ksp_loop_cart = hmriutils.epi.epiphasecorrect(ksp_loop_cart, a);

%% Create zero-filled k-space data
ksp_zf = zeros(Nx,Ny,Nz,Ncoils,Nframes);

% (Re)generate temporally incoherent sampling masks (or load it in)
omegas = false(Ny,Nz,NframesPerLoop);
for frame = 1:Nframes
    rng(mod(frame - 1,NframesPerLoop) + 1); % A different mask per frame
    omega = randsamp2d(N(2:end), R, acs, max_ky_step);
    omegas(:,:,frame) = omega;
end

% Allocate data
ksp_zf(:,permute(repmat(omegas,1,1,1,Ncoils),[1 2 4 3])) = ksp_loop_cart(:,:);

%% IFFT to get images
imgs_mc = ifftshift(ifft(...
                     ifft(...
                      ifft(...
                       fftshift(ksp_zf)...
                       , Nx, 1)...
                      , Ny, 2)...
                     , Nz, 3)...
                    );

%% Get sensitivity maps with either BART or PISCO
smaps_technique = 'pisco';
if doSENSE
    % Reshape and permute gre data
    ksp_gre = ksp_gre_raw(:,:,1:Ny_gre*Nz_gre); % discard trailing data
    ksp_gre = reshape(ksp_gre,Nx_gre,Ncoils,Ny_gre,Nz_gre);
    ksp_gre = permute(ksp_gre,[1 3 4 2]); % [Nx Ny Nz Ncoils]

    if strcmp(smaps_technique, 'pisco')
        fprintf('Estimating sensitivity maps from GRE data via PISCO...\n')
        % PISCO options
        tau = 3;
        threshold = 0.05;
        kernel_shape = 1;
        FFT_nullspace_C_calculation = 1;
        PowerIteration_G_nullspace_vectors = 1;
        M = 10;
        PowerIteration_flag_convergence = 1;
        PowerIteration_flag_auto = 1;
        FFT_interpolation = 1;
        interp_zp = 24;
        gauss_win_param = 100;
        verbose = 0;

        tic
            % Only use central (high SNR) region of k-space to estimate smaps
            cal_length = 32; % Length of each dimension of the calibration data
            center_x = ceil(Nx_gre/2) + ~rem(Nx_gre,2);
            center_y = ceil(Ny_gre/2) + ~rem(Ny_gre,2);
            cal_index_x = center_x + (-floor(cal_length/2):floor(cal_length/2) - ~rem(cal_length/2,2));
            cal_index_y = center_y + (-floor(cal_length/2):floor(cal_length/2) - ~rem(cal_length/2,2));
        
            % Compute smaps slice-by-slice
            tmp = ifftshift(ifft(fftshift(ksp_gre),Nz_gre,3));
            smaps = zeros(Nx_gre, Ny_gre, Nz_gre, Ncoils);
            eigmaps = zeros(Nx_gre, Ny_gre, Nz_gre, Ncoils);
            
            for z = 1:Nz_gre
                [smaps_tmp, eigvals] = PISCO_senseMaps_estimation(...
                    squeeze(tmp(cal_index_x,cal_index_y,z,:)),...
                    [Nx_gre, Ny_gre],...
                    tau,...
                    threshold,...
                    kernel_shape,...
                    FFT_nullspace_C_calculation,...
                    PowerIteration_G_nullspace_vectors,...
                    M,...
                    PowerIteration_flag_convergence,...
                    PowerIteration_flag_auto,...
                    FFT_interpolation,...
                    interp_zp,...
                    gauss_win_param,...
                    verbose...
                );
        
                % Normalize
                smaps_tmp = smaps_tmp./sqrt(sum(abs(smaps_tmp).^2, 3));

                % Mask
                smaps_tmp(repmat(eigvals>0.1, 1, 1, 32)) = 0;
        
                % Allocate
                smaps(:,:,z,:) = smaps_tmp;
                eigmaps(:,:,z) = eigvals;
            end
        toc
    elseif strcmp(smaps_technique,'bart')
        fprintf('Estimating sensitivity maps from GRE data via BART...\n')
        tic
            [calib, emaps] = bart('ecalib -r 20', ksp_gre);
            smaps = squeeze(calib(:,:,:,:,1));
        toc
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
    
    % Align x-direction of smaps
    smaps = flip(smaps,1);
end

%% Coil combination
if doSENSE
    img_final = squeeze(sum(imgs_mc .* conj(smaps), 4));
else % root sum of squares combination
    img_final = squeeze(sqrt(sum(abs(imgs_mc).^2, 4)));
end

%% Compute k-space by IFT3
ksp_final = toppe.utils.ift3(img_final);

%% Recon GRE images
img_gre_mc = toppe.utils.ift3(ksp_gre);
img_gre = sqrt(sum(abs(img_gre_mc).^2,4));

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