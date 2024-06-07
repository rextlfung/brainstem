%%
Nx = 200; Nsegments = 4;
fn = 'Pball200x200x3mm3D.7';
fn_adc = sprintf('adc/P%dadc.mod', Nx);
doSENSE = false;

% Redefine some parameters for convenience
Nz = 3; Ny = Nx; Nframes = 120; Ncoils = 32;
ETL = Ny; % echo train length (Ny)

% Debugging raw data
if true
    ksp_raw = hmriutils.io.ge.loadpfile(fn, [], [], [], 'acqOrder', true); 
    ksp_raw = flip(ksp_raw, 1); % tv6 flips data along FID direction
    [Nfid,Ncoils,N] = size(ksp_raw);
    ksp_raw = ksp_raw(:,:,1:ETL*Nz*Nframes);
    ksp_raw = reshape(ksp_raw, Nfid, Ncoils, ETL, Nz, Nframes);
    ksp_raw = permute(ksp_raw, [1 3 4 5 2]); % [Nfid ETL Nz Nframes Nc]
    % size(ksp_raw); figure; im(log(abs(ksp_raw(:,:,ceil(Nz/2),1,:))))

    % naive recon
    ksp = zeros(size(ksp_raw));
    for seg = 1:Nsegments
        range = (1:Ny/Nsegments) + (seg-1)*Ny/Nsegments;
        tmp = ksp_raw(:,range,:,:,:);
        tmp(:,1:2:end,:,:,:) = flip(tmp(:,1:2:end,:,:,:),1);
        ksp(:,seg:Nsegments:end,:,:,:) = tmp;
    end

    % naively recon each frame
    imgs_naive = zeros(Nfid, ETL, Nz, Nframes);
    for frame = 1:Nframes
        ksp_frame = squeeze(ksp(:,:,:,frame,:));
        img_frame_mc = toppe.utils.ift3(ksp_frame);
        imgs_naive(:,:,:,frame) = sqrt(sum(abs(img_frame_mc).^2,4));
    end
end

%% Load in 2D time series data frame by frame.
ksp = zeros(Nx, Ny, Nz, Nframes, Ncoils);
for frame = 1:Nframes
    isCalShot = frame == 1;

    % load a frame
    ksp_frame = hmriutils.epi.loadframeraw_ge(fn, ETL, Nz, frame, true);

    % Print max real and imag parts to ensure no clipping
    fprintf('Max real part: %d\n', max(real(ksp_frame(:))))
    fprintf('Max imag part: %d\n', max(imag(ksp_frame(:))))

    % rearrange segments into proper locations
    tmp = zeros(size(ksp_frame));
    for seg = 1:Nsegments
        range = (1:Ny/Nsegments) + (seg-1)*Ny/Nsegments;
        tmp2 = ksp_frame(:,range,:,:);
        tmp2(:,1:2:end,:,:) = flip(tmp2(:,1:2:end,:,:),1);
        tmp(:,seg:Nsegments:end,:,:) = tmp2;
    end
    tmp(:,1:2:end,:,:) = flip(tmp(:,1:2:end,:,:),1); % make it look like 1-shot EPI
    ksp_frame = tmp;

    % odd/even echo k-space sampling locations (ramp sampling)
    [rf,gx,gy,gz,desc,paramsint16,pramsfloat,hdr] = toppe.readmod(fn_adc); % commented out because doesn't match
    Nfid = size(ksp_frame,1); % hdr.rfres;

    % estimate and apply odd/even k-space delay (samples)
    if isCalShot
        ksp_cal_stacked = reshape(ksp_frame, [Nfid, Nz*ETL*Ncoils]);
        ksp_cal_stacked(:,1:2:end) = flip(ksp_cal_stacked(:,1:2:end),1);
        [M, I] = max(abs(ksp_cal_stacked),[],1);
        delay = mean(I) - 1 - Nfid/2;
        delay = 0;
    end
    [kxo, kxe] = toppe.utils.getk(sysGE, fn_adc, Nfid, delay);
    
    % grid
    ksp_frame = hmriutils.epi.rampsampepi2cart(ksp_frame, kxo, kxe, Nx, fov(1)*100, 'nufft');
    ksp_frame = squeeze(ksp_frame); % Discard slice dimensions since it's one slice

    % phase correct
    if isCalShot
       cal_data = ifftshift(ifft(fftshift(reshape(ksp_frame, [Nx, Ny*Nz, Ncoils])),Nx,1));
       [a, th] = hmriutils.epi.getoephase(cal_data,false);
    end
    ksp_frame = hmriutils.epi.epiphasecorrect(ksp_frame, a);

    % Allocate
    ksp(:,:,:,frame,:) = ksp_frame;
end

% Print EPI phase offset info
fprintf('Constant phase offset (radians): %f\n', a(1));
fprintf('Linear term (radians/fov): %f\n', a(2));

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
imgs_multicoil = ifftshift(...
    ifft(...
         ifft(...
           ifft(...
             fftshift(ksp), Nx, 1)...
           , Ny, 2)...
         , Nz, 3)...
         );

%% Coil combination
if doSENSE
    imgs = sum(imgs_multicoil .* conj(smaps), 5);
else % root sum of squares combination
    imgs = sqrt(sum(abs(imgs_multicoil).^2, 5));
end

%% Viz
figure; im(imgs,'cbar')
title(fn(1:end-2));
xlabel('space (z)');
ylabel(sprintf('time (frames at %ds volume TR)', volumeTR));
% saveas(gcf, strcat('figs/',fn(1:end-2),'.png'));
