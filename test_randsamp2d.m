%% Set up
load mristack;
% load('/home/rexfung/github/data/myBrain/myBrain.mat')
setEPIparams;

%% Extract one slice
Nframes = 20;
% img = ones(100,100,Nframes);
img = repmat(double(mristack(41:80,:,9)).',1,1,Nframes);
% img = repmat(double(V(:,:,size(V,3)/2)).',1,1,Nframes);

N = size(img,1,2);
Ny = N(1); Nz = N(2);

%% Generate sampling masks
Nframes = size(img,3);
omegas = false(Ny, Nz, Nframes);

masktype = 'rand_caipi';
switch masktype
    case 'full'
        omegas = true(Ny, Nz, Nframes);
    case 'lowpass'
        Ry = sqrt(6); Rz = Ry;
        oemgas = false(Ny, Nz, Nframes);
        omegas((1:ceil(Ny/Ry)) + ceil(Ny/2) - ceil(Ny/Ry/2),...
               (1:ceil(Nz/Rz)) + ceil(Nz/2) - ceil(Nz/Rz/2),...
               :) = true;
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
        for iy = 1:Ry:Ny
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
            omega = randsamp2d(N,R,acs,max_ky_step);
            omegas(:,:,t) = omega;
        end
    case 'rand_caipi'
        Ry = 2; Rz = 1;
        R = [Ry Rz];                    % Acceleration/undersampling factors in each direction
        max_ky_step = round(Ny/16);     % Maximum gap in fast PE direction
        caipi_z = 3;                    % Minimum gap in slow PE direction to prevent duplicate sampling with CAIPI shift

        for t = 1:Nframes
            omega = randsamp2d_caipi(N,R,max_ky_step,caipi_z);
            omegas(:,:,t) = omega;
        end
end

% Compute acceleration
R = numel(omegas)/sum(omegas,'all');

% Select frame for plotting
frame = 4;

% Plot sampling masks
figure('WindowState','maximized');
t = tiledlayout(2,4,"TileSpacing","tight");

nexttile(1);
im((-Ny/2):(Ny/2 - 1), (-Nz/2):(Nz/2 - 1), omegas(:,:,frame));
title(sprintf('Sampling mask, R = %f', R));
xlabel('ky'); ylabel('kz');

% Plot PSFs
psfs = ifftshift(ifft2(fftshift(omegas)));
[Y, Z] = meshgrid((-Ny/2):(Ny/2 - 1), (-Nz/2):(Nz/2 - 1));

nexttile(2);
surf(Y',Z',abs(squeeze(psfs(:,:,frame))),'FaceColor','interp','EdgeColor','interp');
title('| Point spread function |');
axis tight;
tmp = daspect; daspect([tmp(2), tmp(2), tmp(3)]);
xlabel('y (px)'); ylabel('z (px)'); zlabel('PSF (au)');

% Plot retrospectively undersampled images
ksp = ifftshift(fft2(fftshift(img)));
ksp_us = ksp;
ksp_us(~omegas) = 0;
img_us = ifftshift(ifft2(fftshift(ksp_us)));

nexttile(3);
im(img_us(:,:,frame));
title('Aliased image');

% Plot fully sampled image for comparison
nexttile(4);
im(img(:,:,frame));
title('Fully sampled image')

% Plot normalized SVs compared to fully sampled
s_fs = svd(img(:,:,frame));
s_us = svd(img_us(:,:,frame));
s_fs = s_fs ./ s_fs(1);
s_us = s_us ./ s_us(1);

nexttile(5,[1 2]);
plot(1:min(Ny,Nz), s_fs, '-o');
hold on;
plot(1:min(Ny,Nz), s_us, '--o');
legend('Fully sampled','Under sampled');
title('Normalized singular values of displayed 2D image');

% Same thing for the entire space-time matrix
s_fs = svd(reshape(img,[Ny*Nz, Nframes]));
s_us = svd(reshape(img_us,[Ny*Nz, Nframes]));
s_fs = s_fs ./ s_fs(1);
s_us = s_us ./ s_us(1);

nexttile(7,[1 2]);
plot(s_fs, '-o');
hold on;
plot(s_us, '--o');
legend('Fully sampled','Under sampled');
title('Normalized singular values of space-time matrix');

return;

%% Plot sampling masks and aliased images over time
Nframes = 1;

figure('WindowState','maximized');
t = tiledlayout(2,1,"TileSpacing","tight");
nexttile;
im('row',1,(-Ny/2):(Ny/2 - 1), (-Nz/2):(Nz/2 - 1),omegas(:,:,1:Nframes)); title('Sampling mask')
xlabel('ky'); ylabel('kz');

nexttile;
im('row',1,(-Ny/2):(Ny/2 - 1), (-Nz/2):(Nz/2 - 1),img_us(:,:,1:Nframes)); title('Aliased image')
xlabel('ky'); ylabel('kz');

%% Make a movie of point spread functions
[Y, Z] = meshgrid((-Ny/2):(Ny/2 - 1), (-Nz/2):(Nz/2 - 1));
res_y = res(2)*1000; res_z = res(3)*1000; % resolutions in mm
Y = Y*res_y; Z = Z*res_z;

h = figure('WindowState','maximized');
surf(Y',Z',abs(squeeze(psfs(:,:,1))),'FaceColor','interp');
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
    im((-Ny/2):(Ny/2 - 1), (-Nz/2):(Nz/2 - 1), omegas(:,:,frame));
    title(sprintf('frame %d',round(frame)),'FontSize',20);
    if frame == 1
        xlabel('ky (px^{-1})','FontSize',20);
        ylabel('kz (px^{-1})','FontSize',20);
    end

    nexttile(Nframes + frame);
    surf(Y',Z',abs(squeeze(psfs(:,:,frame))),'FaceColor','interp');
    axis tight;
    tmp = daspect;
    daspect([tmp(2), tmp(2), tmp(3)]);
    if frame == 1
        xlabel('y (px)','FontSize',20);
        ylabel('z (px)','FontSize',20);
        zlabel('PSF (au)','FontSize',20);
    end
end

sgtitle()