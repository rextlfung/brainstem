% load img_final;
% load smaps;
load mristack;
close all;
setEPIparams;
%% Extract one slice
img = ones(120,40,20);
% img = squeeze(img_final(:,:,size(img_final,3)/2,:));
% img = repmat(double(mristack(:,:,9)).',1,1,20);

N = size(img,1,2);
Ny = N(1); Nz = N(2);

%% Generate sampling masks
Nframes = size(img,3);
omegas = false(Ny, Nz, Nframes);

masktype = 'rand';
switch masktype
    case 'full'
        omegas = true(Ny, Nz, Nframes);
    case '1d'
        omegas = false(Ny, Nz, Nframes);
        Ry = 6;
        for iy = 1:Ry:Ny
            omegas(iy,:,:) = true;
        end
    case '2d'
        omegas = false(Ny, Nz, Nframes);
        Ry = 2; Rz = 3;
        for iy = 1:Ry:Ny
            for iz = 1:Rz:Nz
                omegas(iy,iz,:) = true;
            end
        end
    case 'caipi'
        omegas = false(Ny, Nz, Nframes);
        Ry = 2; Rz = 3; caipiShift = 0;
        for iy = 1:R:Ny
            for iz = (1 + caipiShift):Rz:Nz
                omegas(iy,iz,:) = true;
            end
            caipiShift = mod(caipiShift + 1, Rz);
        end
    case 'rand' % Randomly generate new sampling patterns every time frame
        Ry = 2; Rz = 3;
        R = [Ry Rz];                    % Acceleration/undersampling factors in each direction
        acs = [1/16 1/16];              % Central portion of ky-kz space to fully sample
        max_ky_step = round(Ny/16);     % Maximum gap in fast PE direction

        for t = 1:Nframes
            rng(t);
            omega = randsamp2d(N,R,acs,max_ky_step);
            omegas(:,:,t) = omega;
        end
end

psfs = ifftshift(ifft2(fftshift(omegas)));

% Plot sampling masks
figure('WindowState','maximized');
t = tiledlayout(1,3,"TileSpacing","tight");
nexttile;
im((-Ny/2):(Ny/2 - 1), (-Nz/2):(Nz/2 - 1), omegas(:,:,1)); title('Sampling mask')
xlabel('ky'); ylabel('kz');

% Plot PSFs
[Y, Z] = meshgrid((-Ny/2):(Ny/2 - 1), (-Nz/2):(Nz/2 - 1));

nexttile;
surf(Y',Z',abs(squeeze(psfs(:,:,1)))); title('Point spread function')
axis tight;
tmp = daspect; daspect([tmp(2), tmp(2), tmp(3)]);
xlabel('y (px)'); ylabel('z (px)'); zlabel('PSF (au)');

% Plot fully sampled images
% figure('WindowState','maximized');
% im(img(:,:,1),'cbar')
% title('Fully sampled image');

% Plot retrospectively undersampled images
ksp = ifftshift(fft2(fftshift(img)));
ksp_us = ksp;
ksp_us(~omegas) = 0;
img_us = ifftshift(ifft2(fftshift(ksp_us)));

nexttile;
im(img_us(:,:,1)); title(sprintf('Aliased image'));

return;

%% Plot sampling masks over time
figure('WindowState','maximized');
im(omegas(:,:,1:6)); title('Sampling mask')
xlabel('ky'); ylabel('kz');


%% Plot aliased images over time
figure('WindowState','maximized');
im(img_us(:,:,1:6)); title('Aliased image')
xlabel('ky'); ylabel('kz');

%% Make a movie of point spread functions
[Y, Z] = meshgrid((-Ny/2):(Ny/2 - 1), (-Nz/2):(Nz/2 - 1));
res_y = res(2)*1000; res_z = res(3)*1000; % resolutions in mm
Y = Y*res_y; Z = Z*res_z;

h = figure('WindowState','maximized');
surf(Y',Z',abs(squeeze(psfs(:,:,1))));
axis tight;
tmp = daspect; daspect([tmp(2), tmp(2), tmp(3)]);
xlabel('y (mm)'); ylabel('z (mm)'); zlabel('PSF (au)');
ax = gca;
ax.NextPlot = 'replaceChildren';

M(Nframes) = struct('cdata',[],'colormap',[]);

for t = 1:Nframes
    surf(Y',Z',abs(squeeze(psfs(:,:,t))));
    axis tight;
    tmp = daspect; daspect([tmp(2), tmp(2), tmp(3)]);
    xlabel('y (mm)'); ylabel('z (mm)'); zlabel('PSF (au)')
    title(sprintf('frame %d', round(t)));
    drawnow
    pause;
end

%% Plot 3 sampling masks and PSFs side by side
figure('WindowState','maximized');
Nframes = 3;
tiledlayout(2,Nframes,"TileSpacing","none");

for frame = 1:Nframes
    nexttile(frame);
    im((-Ny/2):(Ny/2 - 1), (-Nz/2):(Nz/2 - 1), omegas(:,:,frame)); title(sprintf('frame %d',round(frame)));
    if frame == 1
        xlabel('ky (px^{-1})'); ylabel('kz (px^{-1})');
    end

    nexttile(Nframes + frame);
    surf(Y',Z',abs(squeeze(psfs(:,:,frame))),'FaceColor','interp');
    axis tight;
    tmp = daspect; daspect([tmp(2), tmp(2), tmp(3)]);
    if frame == 1
        xlabel('y (px)'); ylabel('z (px)'); zlabel('PSF (au)');
    end
end