%% CG-SENSE recon with BART
% Load data
setGREparams; setEPIparams;

kdata_path = "/home/rexfung/github/JuliaImageRecon/data/20241017fingertap/kdata2x3.mat";
smaps_path = "/home/rexfung/github/JuliaImageRecon/data/20241017fingertap/smaps.mat";

load(kdata_path); % ksp_zf
load(smaps_path); % smaps

[Nx, Ny, Nz, Ncoils, Nframes] = size(ksp_zf);

%% Recon with CG-SENSE
tic;
img = zeros(Nx,Ny,Nz,Nframes);
for frame = 1:Nframes
    disp(frame);
    data = squeeze(ksp_zf(:,:,:,:,frame));
    img(:,:,:,frame) = bart('pics -l1 -r 0.001', data, smaps);
end
toc;

%% Save file
save("output_cg.mat","img","-v7.3");