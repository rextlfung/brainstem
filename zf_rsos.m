%% ZF-RSOS recon
% Load data
setGREparams; setEPIparams;

kdata_path = "/home/rexfung/github/data/20241017fingertap/kdata2x3.mat";
smaps_path = "/home/rexfung/github/data/20241017fingertap/smaps.mat";

kdata = matfile(kdata_path); % ksp_zf

Nx = 120; Ny = 120; Nz = 40; Ncoils = 32; Nframes = 280;

%% Recon
tic;
delete(gcp('nocreate'));
parpool(4);
img = zeros(Nx,Ny,Nz,Nframes);
parfor frame = 1:Nframes
    fprintf('Reconstructing frame %d\n', round(frame));
    data = squeeze(kdata.ksp_zf(:,:,:,:,frame));
    img(:,:,:,frame) = squeeze(sqrt(sum(abs(toppe.utils.ift3(data)).^2, 4)));
end
delete(gcp('nocreate'));
toc;

%% Save
save('zf_rsos.mat','img','-v7.3');
