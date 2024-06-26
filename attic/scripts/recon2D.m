%%
Nx = 120;
scanner = 'inside';
fn = sprintf('Pball%d%s.7', Nx, scanner);
fn_adc = sprintf('adc/Pball%dadc.mod', Nx);
doSENSE = false;

% Redefine some parameters for convenience
Ny = Nx; Nframes = 3; Ncoils = 32;
ETL = Ny; % echo train length (Ny)
Np = 1; % number of partitions

% Load in raw data
ksp_raw = hmriutils.io.ge.loadpfile(fn, [], [], [], 'acqOrder', true); 
ksp_raw = flip(ksp_raw, 1); % tv6 flips data along FID direction
[Nfid,Ncoils,N] = size(ksp_raw);
ksp_raw = ksp_raw(:,:,1:ETL*Np*Nframes);
ksp_raw = reshape(ksp_raw, Nfid, Ncoils, ETL, Np*Nframes);
ksp_raw = permute(ksp_raw, [1 3 4 2]); % [Nfid ETL Np Nc]
ksp_raw(:,2:2:end,:,:) = flip(ksp_raw(:,2:2:end,:,:),1); % Flip even lines

% Load in 2D time series data frame by frame.
ksp = zeros(Nx, Ny, Nframes, Ncoils);
for frame = 1:Nframes
    isCalShot = frame == 1;

    % load a frame
    ksp_raw_frame = hmriutils.epi.loadframeraw_ge(fn, ETL, Np, frame, true);

    % odd/even echo k-space sampling locations (ramp sampling)
    [rf,gx,gy,gz,desc,paramsint16,pramsfloat,hdr] = toppe.readmod(fn_adc); % commented out because doesn't match
    Nfid = size(ksp_raw_frame,1); % hdr.rfres;

    % estimate and apply odd/even k-space delay (samples)
    if strcmp(scanner,'inside')
        delay = -1.0;
    elseif strcmp(scanner,'outside')
        delay = 0.5;
    end
    [kxo, kxe] = toppe.utils.getk(sysGE, fn_adc, Nfid, delay);
    
    % grid
    ksp_frame = hmriutils.epi.rampsampepi2cart(ksp_raw_frame, kxo, kxe, Nx, fov(1)*100, 'spline');
    ksp_frame = squeeze(ksp_frame); % Discard slice dimensions since it's one slice

    % phase correct
    if isCalShot
       cal_data = ifftshift(ifft(fftshift(ksp_frame),Nfid,1));
       [a, th] = hmriutils.epi.getoephase(cal_data(:,1:end,:)); % only 4 lines needed
    end
    ksp_frame = hmriutils.epi.epiphasecorrect(ksp_frame, a);

    % Allocate
    ksp(:,:,frame,:) = ksp_frame;

    % Print max real and imag parts to ensure no clipping
    fprintf('Max real part: %d\n', max(real(ksp_raw_frame(:))))
    fprintf('Max imag part: %d\n', max(imag(ksp_raw_frame(:))))
end

% Print EPI phase offset info
fprintf('Constant phase offset (radians): %f\n', a(1));
fprintf('Linear term (radians/fov): %f\n', a(2));

%% Get sensitivity maps with PISCO
if doSENSE
    addpath('pisco/');
    smaps = zeros(Nx, Ny, Nframes, Ncoils);
    for frame = 1:Nframes
        ksp_frame = squeeze(ksp(:,:,frame,:));
        [smaps_frame, eigvals] = PISCO_senseMaps_estimation(ksp_frame, [Nx, Ny]);
        smaps(:,:,frame,:) = smaps_frame;
    end
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
    imgs = zeros(Nx, Ny, Nframes);
    for frame = 1:Nframes
        ksp_frame = squeeze(ksp(:,:,frame,:));
        smaps_frame = squeeze(smaps(:,:,frame,:));
        imgs(:,:,frame) = sum(ifftshift(ifft2(fftshift(ksp_frame))) .* conj(smaps_frame), 3);
    end
else % root sum of squares combination
    imgs = sqrt(sum(imgs_multicoil.^2, 4));
end

%% Viz
close all;
frame = 2;
figure;
subplot(121); im(abs(imgs(:,:,frame)),'cbar'); title('magnitude');
subplot(122); im(imag(imgs(:,:,frame)),'cbar'); title('phase');
sgtitle(sprintf('%s, delay = %d',fn(1:end-2),delay));
savefig(strcat('figs/',fn(1:end-2),'.fig'));
figure; im(log(abs(ksp(:,:,frame,:))))
