% load img_final;
% load smaps;
% load mristack;
close all;

%% Extract one slice
img = ones(120, 40, 4);
% img = squeeze(img_final(:,:,size(img_final,3)/2,:));
% img = repmat(double(mristack(:,:,9)).',1,1,20);

N = size(img,1,2);
Ny = N(1); Nz = N(2);
Ry = 2; Rz = 3;
R = [Ry Rz];                    % Acceleration/undersampling factors in each direction
acs = [1/16 1/16];              % Central portion of ky-kz space to fully sample
max_ky_step = round(Ny/16);     % Maximum gap in fast PE direction

%% Generate a new sampling patterns every time frame
Nframes = size(img,3);
omegas = false(N(1),N(2),Nframes);

for t = 1:Nframes
    rng(t);
    omega = randsamp2d(N,R,acs,max_ky_step);
    omegas(:,:,t) = omega;
    
end

%% Plot sampling masks
figure('WindowState','maximized');
im(omegas); title('Sampling masks over time')
xlabel('ky'); ylabel('kz');

%% Plot point spread functions
figure('WindowState','maximized');
im(ifftshift(ifft2(fftshift(omegas)))); title('Point spread functions over time')
xlabel('y'); ylabel('z');

return;
%% Plot fully sampled images
figure('WindowState','maximized');
im('col',Nframes/2,'row',2,img,'cbar')
title('Fully sampled images');

%% Plot retrospectively undersampled images
ksp = ifftshift(fft2(fftshift(img)));
ksp_us = ksp;
ksp_us(~omegas) = 0;
img_us = ifftshift(ifft2(fftshift(ksp_us)));

figure('WindowState','maximized');
im('col',Nframes/2,'row',2,img_us,'cbar')
title(sprintf('Pseudo-randomly undersampled images. Rx = %d, Ry = %d',R(1),R(2)));

return;
%% Multicoil case (kinda irrelevant)
Ncoils = 32;

ksp_mc_us = squeeze(ksp_mc(:,:,size(ksp_mc,3)/2,:,:));
ksp_mc_us(~repmat(omegas,1,1,1,Ncoils)) = 0;
ksp_mc_us = permute(ksp_mc_us,[1 2 4 3]);
img_mc_us = ifftshift(ifft2(fftshift(ksp_mc_us)));
img_us = squeeze(sum(img_mc_us .* conj(smaps),3));

figure('WindowState','maximized');
im('col',Nframes/2,'row',2,img_us,'cbar')
title(sprintf('Pseudo-randomly undersampled images. Rx = %d, Ry = %d',R(1),R(2)));