% Code for loading and preprocessing randomly undersampled 3D EPI data
% With CAIPI shifting
% Rex Fung, Nov 21st, 2024

%% Set params and options
% Type of EPI acquitsition. "rand" or "rand_caipi"
mode = 'rand_caipi';

% Load in params
setGREparams;
setEPIparams;

% Total number of frames
Nloops = 1; % Defined as toppe cv 8 at scanner
Nframes = Nloops*NframesPerLoop*3;

% Filenames
datdir = '/mnt/storage/rexfung/20250117ball/';
fn_gre = strcat(datdir,'gre.h5');
fn_cal = strcat(datdir,'cal.h5');
fn_loop = strcat(datdir,'loop14.h5');
fn_samp_log = strcat(datdir,'samp_logs/14.mat');

% Options
doSENSE = false; % Takes a while
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

% discard unusual zeros in the middle of the data (currently unknown why
% they exist) TODO: FIND OUT WHY
non_zero_indices = (squeeze(sum(abs(ksp_loop_raw), [1 2])) ~= 0);
ksp_loop_raw = ksp_loop_raw(:,:,non_zero_indices);

% Print max real and imag parts to check for reasonable magnitude
fprintf('Max real part of gre data: %d\n', max(real(ksp_gre_raw(:))))
fprintf('Max imag part of gre data: %d\n', max(imag(ksp_gre_raw(:))))
fprintf('Max real part of cal data: %d\n', max(real(ksp_cal_raw(:))))
fprintf('Max imag part of cal data: %d\n', max(imag(ksp_cal_raw(:))))
fprintf('Max real part of loop data: %d\n', max(real(ksp_loop_raw(:))))
fprintf('Max imag part of loop data: %d\n', max(imag(ksp_loop_raw(:))))

%% Compute odd/even delays using calibration (blipless) data
% Reshape and permute calibration data (a single frame w/out blips)
ksp_cal = ksp_cal_raw(:,:,1:2*ceil(Ny/Ry/2)*round(Nz/caipi_z/Rz)); % discard trailing data
ksp_cal = reshape(ksp_cal,Nfid,Ncoils,2*ceil(Ny/Ry/2), round(Nz/caipi_z/Rz));
ksp_cal = permute(ksp_cal,[1 3 4 2]); % [Nfid Ny/Ry Nz/Rz Ncoils]

% Estimate k-space center offset due to gradient delay
cal_data = reshape(abs(ksp_cal),Nfid,2*ceil(Ny/Ry/2),round(Nz/caipi_z/Rz),Ncoils);
cal_data(:,1:2:end,:,:) = flip(cal_data(:,1:2:end,:,:),1);
[M, I] = max(cal_data,[],1);
delay = Nfid/2 - mean(I,'all');
fprintf('Estimated offset from center of k-space (samples): %f\n', delay);

% retrieve sample locations from .mod file with adc info
fn_adc = strcat(datdir, sprintf('adc%d.mod',Nfid));
[rf,gx,gy,gz,desc,paramsint16,pramsfloat,hdr] = toppe.readmod(fn_adc);
[kxo, kxe] = toppe.utils.getk(sysGE, fn_adc, Nfid, delay);

% Extract even number of lines (in case ETL is odd)
ETL_even = size(ksp_cal,2) - mod(size(ksp_cal,2),2);
oephase_data = ksp_cal(:,1:ETL_even,:,:,:);

% EPI ghost correction phase offset values
oephase_data = hmriutils.epi.rampsampepi2cart(oephase_data, kxo, kxe, Nx, fov(1)*100, 'nufft');
oephase_data = ifftshift(ifft(fftshift(reshape(oephase_data,Nx,ETL_even,round(Nz/caipi_z/Rz)*Ncoils)),Nx,1));
[a, th] = hmriutils.epi.getoephase(oephase_data,showEPIphaseDiff);
fprintf('Constant phase offset (radians): %f\n', a(1));
fprintf('Linear term (radians/fov): %f\n', a(2));

%% Grid and apply odd/even correction to loop data
% Reshape and permute loop data
ksp_loop = ksp_loop_raw(:,:,1:2*ceil(Ny/Ry/2)*round(Nz/caipi_z/Rz)*Nframes); % discard trailing empty data
ksp_loop = reshape(ksp_loop,Nfid,Ncoils,2*ceil(Ny/Ry/2)*round(Nz/caipi_z/Rz),Nframes);
ksp_loop = permute(ksp_loop,[1 3 2 4]); % [Nfid Ny/Ry*Nz/Rz Ncoils Nframes]

% Grid along kx direction via NUFFT (takes a while)
ksp_loop_cart = zeros([Nx,size(ksp_loop,2:ndims(ksp_loop))]);
tic
    parfor frame = 1:Nframes
        fprintf('Gridding frame %d\n', round(frame));
        tmp = squeeze(ksp_loop(:,:,:,frame));
        tmp1 = hmriutils.epi.rampsampepi2cart(tmp, kxo, kxe, Nx, fov(1)*100, 'nufft');
        ksp_loop_cart(:,:,:,frame) = tmp1;
    end
toc

% Phase correct along kx direction
ksp_loop_cart = hmriutils.epi.epiphasecorrect(ksp_loop_cart, a);

%% Create zero-filled k-space data
ksp_zf = zeros(Nx,Ny,Nz,Ncoils,Nframes);
load(fn_samp_log);

% Read through log of sample locations and allocate data
for frame = 1:Nframes
    for samp_count = 1:2*ceil(Ny/Ry/2)*round(Nz/caipi_z/Rz)
        iy = samp_log(frame,samp_count,1);
        iz = samp_log(frame,samp_count,2);
        if ksp_zf(:,iy,iz,:,frame) ~= 0
            disp(frame);
            disp([iy,iz]);
            pause;
        end
        ksp_zf(:,iy,iz,:,frame) = ksp_loop_cart(:,samp_count,:,frame);
    end
end

%% IFFT to get images
imgs_mc = ifftshift(ifft(...
                     ifft(...
                      ifft(...
                       fftshift(ksp_zf)...
                       , Nx, 1)...
                      , Ny, 2)...
                     , Nz, 3)...
                    );

%% Free up memory
% clear ksp_loop_cart ksp_loop ksp_loop_raw;

%% Get sensitivity maps with either BART or PISCO
% Reshape and permute gre data
ksp_gre = ksp_gre_raw(:,:,1:Ny_gre*Nz_gre); % discard trailing data
ksp_gre = reshape(ksp_gre,Nx_gre,Ncoils,Ny_gre,Nz_gre);
ksp_gre = permute(ksp_gre,[1 3 4 2]); % [Nx Ny Nz Ncoils]

smaps_technique = 'pisco';
if doSENSE
    if strcmp(smaps_technique, 'pisco')
        fprintf('Estimating sensitivity maps from GRE data via PISCO...\n')
        % PISCO options
        tau = 3;
        threshold = 0.05;
        kernel_shape = 1;
        FFT_nullspace_C_calculation = 1;
        PowerIteration_G_nullspace_vectors = 1;
        M = 20;
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
            smaps_raw = zeros(Nx_gre, Ny_gre, Nz_gre, Ncoils);
            eigmaps = zeros(Nx_gre, Ny_gre, Nz_gre);
            
            parfor z = 1:Nz_gre
                fprintf('estimate z = %d\n', round(z));
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
        
                % Allocate
                smaps_raw(:,:,z,:) = smaps_tmp;
                eigmaps(:,:,z,:) = eigvals;
            end
        toc
    elseif strcmp(smaps_technique,'bart')
        fprintf('Estimating sensitivity maps from GRE data via BART...\n')
        tic
            [calib, emaps] = bart('ecalib -r 20', ksp_gre);
            smaps_raw = squeeze(calib(:,:,:,:,1));
        toc
    end

    % Mask
    smaps = smaps_raw;
    % smaps(repmat(eigmaps>8*threshold,1,1,1,32)) = 0;

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
    
    % Align x-direction of smaps with EPI data
    smaps = flip(smaps,1);
end

%% Coil combination
if doSENSE
    img = squeeze(sum(imgs_mc .* conj(smaps), 4));
else % root sum of squares combination
    img = squeeze(sqrt(sum(abs(imgs_mc).^2, 4)));
end

%% Compute k-space by IFT3
ksp_final = toppe.utils.ift3(img);

%% Recon GRE images
img_gre_mc = toppe.utils.ift3(ksp_gre);
img_gre = sqrt(sum(abs(img_gre_mc).^2,4));

%% Viz
close all;

% Plot a frame
frame = size(img,ndims(img));
figure('WindowState','maximized');
im('mid3',img(:,:,:,frame),'cbar')
title(sprintf('|image|, middle 3 planes of frame %d',frame));
ylabel('y'); xlabel('x')

%% Plot Sensitivity maps
if doSENSE
    figure('WindowState','maximized');
    im('col',Ncoils,'row',Nz,reshape(permute(squeeze(smaps),[1 2 4 3]),Nx,Ny,Ncoils*Nz),'cbar')
    title('Sensitivity maps');
    ylabel('z direction'); xlabel('Coils');
end

return;
%% Save for next step of recon
save(strcat(datdir,'recon/kdata.mat'),'ksp_zf','-v7.3');
save(strcat(datdir,'recon/smaps.mat'),'smaps','-v7.3');

