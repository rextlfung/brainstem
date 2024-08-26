load img_final;
load smaps;
load mristack;
close all;

%% Extract one slice
img = squeeze(img_final(:,:,size(img_final,3)/2,:));
img = repmat(img,1,1,15);
% img = repmat(double(mristack(:,:,9)).',1,1,10);

N = size(img,1,2);
R = [3 3];
acs = [1/8 1/8];

Nframes = size(img,3);
xy_masks = false(N(1),N(2),Nframes);

%% Generate a new sampling patterns every time frame
figure('WindowState','maximized');
tiledlayout('flow');
for t = 1:Nframes
    rng(t);
    omega = randsamp2d(N,R,acs);
    xy_masks(:,:,t) = omega;

    % nexttile; im(omega); title(sprintf('Frame %d',round(t)))
    % xlabel('kx'); ylabel('ky');
end

%% Plot fully sampled images
figure('WindowState','maximized');
im('col',Nframes/2,'row',2,img,'cbar')
title('Fully sampled images');

%% Plot retrospectively undersampled images
ksp = ifftshift(fft2(fftshift(img)));
ksp_us = ksp;
ksp_us(~xy_masks) = 0;
img_us = ifftshift(ifft2(fftshift(ksp_us)));

figure('WindowState','maximized');
im('col',Nframes/2,'row',2,img_us,'cbar')
title(sprintf('Pseudo-randomly undersampled images. Rx = %d, Ry = %d',R(1),R(2)));

%% Multicoil case
Ncoils = 32;

ksp_mc_us = squeeze(ksp_mc(:,:,size(ksp_mc,3)/2,:,:));
ksp_mc_us(~repmat(xy_masks,1,1,1,Ncoils)) = 0;
ksp_mc_us = permute(ksp_mc_us,[1 2 4 3]);
img_mc_us = ifftshift(ifft2(fftshift(ksp_mc_us)));
img_us = squeeze(sum(img_mc_us .* conj(smaps),3));

figure('WindowState','maximized');
im('col',Nframes/2,'row',2,img_us,'cbar')
title(sprintf('Pseudo-randomly undersampled images. Rx = %d, Ry = %d',R(1),R(2)));