%%
Nx = 200; Nsegments = 4; scanner = 'inside';
fn = 'Pball200x200x3mm.7';
fn_adc = sprintf('adc/Pball%dadc_%dsegs.mod', Nx, Nsegments);
doSENSE = false;

% Redefine some parameters for convenience
Ny = Nx; Nframes = 1 + 2; Ncoils = 32;
ETL = Ny; % echo train length (Ny)
Np = 1; % number of partitions

%% Debugging raw data
close all;
if false
    ksp_raw = hmriutils.io.ge.loadpfile(fn, [], [], [], 'acqOrder', true); 
    ksp_raw = flip(ksp_raw, 1); % tv6 flips data along FID direction
    [Nfid,Ncoils,N] = size(ksp_raw);
    ksp_raw = ksp_raw(:,:,1:ETL*Np*Nframes);
    ksp_raw = reshape(ksp_raw, Nfid, Ncoils, ETL, Np*Nframes);
    ksp_raw = permute(ksp_raw, [1 3 4 2]); % [Nfid ETL Np Nc]
    size(ksp_raw);
    figure; im(log(abs(ksp_raw(:,:,2,:))))

    % naive recon
    ksp = zeros([size(ksp_raw), Nsegments]);
    for seg = 1:Nsegments
        range = (1:Ny/Nsegments) + (seg-1)*Ny/Nsegments;
        tmp = ksp_raw(:,range,:,:);
        tmp(:,1:2:end,:,:) = flip(tmp(:,1:2:end,:,:),1);
        ksp(:,seg:Nsegments:end,:,:,seg) = tmp;
    end
    figure; im(log(abs(ksp(:,:,2,2,:))))
    
    imgs = ifftshift(ifft2(fftshift(ksp)));
    imgs = sum(imgs,5);
    img = sqrt(sum(imgs.^2,4));
    figure; im(img(:,:,2,:))
end

%% Load in 2D time series data frame by frame.
ksp = zeros(Nx, Ny, Nframes, Ncoils);
for frame = 1:Nframes
    isCalShot = frame == 1;

    % load a frame
    ksp_frame = hmriutils.epi.loadframeraw_ge(fn, ETL, Np, frame, true);

    % Print max real and imag parts to ensure no clipping
    fprintf('Max real part: %d\n', max(real(ksp_frame(:))))
    fprintf('Max imag part: %d\n', max(imag(ksp_frame(:))))

    % rearrange segments into proper locations
    tmp = zeros(size(ksp_frame));
    for seg = 1:Nsegments
        range = (1:Ny/Nsegments) + (seg-1)*Ny/Nsegments;
        tmp2 = ksp_frame(:,range,:);
        tmp2(:,1:2:end,:) = flip(tmp2(:,1:2:end,:),1);
        tmp(:,seg:Nsegments:end,:) = tmp2;
    end
    tmp(:,1:2:end,:) = flip(tmp(:,1:2:end,:),1); % make it look like 1-shot EPI
    ksp_frame = tmp;

    % odd/even echo k-space sampling locations (ramp sampling)
    [rf,gx,gy,gz,desc,paramsint16,pramsfloat,hdr] = toppe.readmod(fn_adc); % commented out because doesn't match
    Nfid = size(ksp_frame,1); % hdr.rfres;

    % estimate and apply odd/even k-space delay (samples)
    if isCalShot
        ksp_cal_stacked = reshape(ksp_frame, [Nfid, ETL*Ncoils]);
        ksp_cal_stacked(:,1:2:end) = flip(ksp_cal_stacked(:,1:2:end),1);
        [M, I] = max(abs(ksp_cal_stacked),[],1);
        delay = mean(I) - 1 - Nfid/2;
    end
    [kxo, kxe] = toppe.utils.getk(sysGE, fn_adc, Nfid, delay);
    
    % grid
    ksp_frame = hmriutils.epi.rampsampepi2cart(ksp_frame, kxo, kxe, Nx, fov(1)*100, 'spline');
    ksp_frame = squeeze(ksp_frame); % Discard slice dimensions since it's one slice

    % phase correct
    if isCalShot
       cal_data = ifftshift(ifft(fftshift(ksp_frame),Nx,1));
       [a, th] = hmriutils.epi.getoephase(cal_data,true);
    end
    ksp_frame = hmriutils.epi.epiphasecorrect(ksp_frame, a);

    % Allocate
    ksp(:,:,frame,:) = ksp_frame;
end

% Print EPI phase offset info
fprintf('Constant phase offset (radians): %f\n', a(1));
fprintf('Linear term (radians/fov): %f\n', a(2));

%% Get sensitivity maps with PISCO
if doSENSE
    addpath('pisco/');
    ksp_frame1 = squeeze(ksp(:,:,2,:));
    [smaps, eigvals] = PISCO_senseMaps_estimation(ksp_frame1, [Nx, Ny]);

    % smooth maps
    % sd = 3;
    % smaps = imgaussfilt(real(smaps),sd)...
    %             + 1j*imgaussfilt(imag(smaps),sd);

    % Add a time dimension so it can be multiplied
    smaps = reshape(smaps,[Nx Ny 1 Ncoils]);
end

%% IFFT to get images
imgs_multicoil = ifftshift(...
         ifft(...
           ifft(...
             fftshift(ksp), Nx, 1)...
           , Ny, 2)...
         );

%% Coil combination
if doSENSE
    imgs = sum(ifftshift(ifft2(fftshift(ksp))) .* conj(smaps), 4);
else % root sum of squares combination
    imgs = sqrt(sum(abs(imgs_multicoil).^2, 4));
end

%% Viz
close all;
figure; im(imgs(:,:,2:end),'cbar')
sgtitle(fn(1:end-2));
% saveas(gcf, strcat('figs/',fn(1:end-2),'.png'));
