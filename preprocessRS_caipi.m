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
Nframes = Nloops*NframesPerLoop;

% Filenames
datdir = '/mnt/storage/rexfung/20250330ball/';
fn_gre = strcat(datdir,'gre.h5');
fn_cal = strcat(datdir,'cal.h5');
fn_loop = strcat(datdir,'loop1.h5');
fn_samp_log = strcat(datdir,'samp_logs/loop.mat');
fn_smaps = strcat(datdir,'smaps.mat');

% Options
doSENSE = true; % Takes a while
SENSE_meth = 'pisco';
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
ksp_cal = ksp_cal_raw(:,:,1:Ny*Nz); % discard trailing data
ksp_cal = reshape(ksp_cal,Nfid,Ncoils,Ny,Nz);
ksp_cal = permute(ksp_cal,[1 3 4 2]); % [Nfid Ny Nz Ncoils]

% Estimate k-space center offset due to gradient delay
cal_data = reshape(abs(ksp_cal),Nfid,Ny,Nz,Ncoils);
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
oephase_data = ifftshift(ifft(fftshift(reshape(oephase_data,Nx,ETL_even,Nz*Ncoils)),Nx,1));
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

%% Free up memory
clear ksp_loop_cart ksp_loop ksp_loop_raw;

%% Rebuild sampling mask from samp_log
omega = false(Ny, Nz, Nframes);
for f = 1:size(samp_log, 1)
    for k = 1:size(samp_log, 2)
        omega(samp_log(f,k,1), samp_log(f,k,2), f) = true;
    end
end

%% Average acquired samples across temporal dimension
% nsamps = sum(omega,3);
% ksp_mc = sum(ksp_zf,5)./permute(repmat(nsamps,[1 1 Nx Ncoils]), [3 1 2 4]);
% ksp_mc(isnan(ksp_mc)) = 0;
ksp_mc = ksp_zf;

%% IFFT to get images
imgs_mc = zeros(Nx, Ny, Nz, Ncoils, Nframes);
for frame = 1:Nframes
    imgs_mc(:,:,:,:,frame) = toppe.utils.ift3(ksp_mc(:,:,:,:,frame));
end

%% Get sensitivity maps with either BART or PISCO
% Reshape and permute gre data
ksp_gre = ksp_gre_raw(:,:,1:Ny_gre*Nz_gre); % discard trailing data
ksp_gre = reshape(ksp_gre,Nx_gre,Ncoils,Ny_gre,Nz_gre);
ksp_gre = permute(ksp_gre,[1 3 4 2]); % [Nx Ny Nz Ncoils]

if doSENSE
    if exist(fn_smaps, 'file')
        load(fn_smaps);
    else
        fprintf('Estimating sensitivity maps from GRE data via %s...\n', SENSE_meth)
    
        % Compute sensitivity maps
        tic
            smaps_raw = makeSmaps(ksp_gre, SENSE_meth);
        toc
    
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

        % Save for next time
        save(fn_smaps, 'smaps', '-v7.3');
    end
end

%% Align x-direction of smaps with EPI data (sometimes necessary)
smaps = flip(smaps,1);

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

return;
%% Viz
close all;

figure('WindowState','maximized'); tiledlayout(1,2,'TileSpacing','tight');
% Plot temporally averaged k-space
% nexttile; im('mid3', log(abs(ksp_zf(:,:,:,16))),'cbar');
% title('|k-space|, middle 3 planes, coil #16');

% Plot total sampling mask
nexttile; im('blue0', nsamps);
title(sprintf('number of samples acquired at each location (blue = 0, white = %d)', max(nsamps(:))));
xlabel('ky'); ylabel('kz');

% Plot static image
nexttile; im('mid3',img,'cbar');
title('|image|, middle 3 planes');
xlabel('x / z'); ylabel('z / y');

% Plot a frame of time series
% frame = size(img,ndims(img));
% figure('WindowState','maximized');
% im('mid3',img(:,:,:,frame),'cbar')
% title(sprintf('|image|, middle 3 planes of frame %d',frame));
% ylabel('y'); xlabel('x')

return;

%% Plot Sensitivity maps
if doSENSE
    figure('WindowState','maximized');
    im('col',Ncoils,'row',Nz,reshape(permute(squeeze(smaps),[1 2 4 3]),Nx,Ny,Ncoils*Nz),'cbar')
    title('Sensitivity maps');
    ylabel('z direction'); xlabel('Coils');
end
%% Save for next step of recon
save(strcat(datdir,'recon/kdata.mat'),'ksp_mc','-v7.3');
save(strcat(datdir,'recon/smaps.mat'),'smaps','-v7.3');

