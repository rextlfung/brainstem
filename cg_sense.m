%% CG-SENSE recon with BART
% Load data
setGREparams; setEPIparams;

kdata_path = "/home/rexfung/github/data/20241017fingertap/kdata2x3.mat";
smaps_path = "/home/rexfung/github/data/20241017fingertap/smaps.mat";

kdata = matfile(kdata_path); % ksp_zf
load(smaps_path); % smaps
smaps = flip(smaps,1);

Nx = 120; Ny = 120; Nz = 40; Ncoils = 32; Nframes = 280;

%% Recon with CG-SENSE
tic;
delete(gcp('nocreate'));
parpool(8);
img = zeros(Nx,Ny,Nz,Nframes);
parfor frame = 1:Nframes
    fprintf('Reconstructing frame %d\n', round(frame));
    data = squeeze(kdata.ksp_zf(:,:,:,:,frame));
    img(:,:,:,frame) = bart('pics -l1 -r0.001', data, smaps);
end
delete(gcp('nocreate'));
toc;

%% Save
save('icg_1e-3.mat','img','-v7.3');
