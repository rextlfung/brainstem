% This script allows reproducing the examples shown in the technical report 
%
% [1] R. A. Lobos, C.-C. Chan, J. P. Haldar.  PISCO Software Version 1.0 
% University of Southern California, Los Angeles, CA, Technical Report 
% USC-SIPI-458. 2023.
%
% The problem formulation and methods implemented by the associated
% software to this script were originally reported in
%
% [2] R. A. Lobos, C.-C. Chan, J. P. Haldar.  New Theory and Faster 
% Computations for Subspace-Based Sensitivity Map Estimation in 
% Multichannel MRI. Submitted.
%
% [3] R. A. Lobos, C.-C. Chan, J. P. Haldar.  Extended Version of "New 
% Theory and Faster Computations for Subspace-Based Sensitivity Map 
% Estimation in Multichannel MRI", 2023, arXiv:2302.13431.
% (https://arxiv.org/abs/2302.13431)
%
% The associated software and this script are available at
%
% http://mr.usc.edu/download/pisco/
%
% As described on that page, use of this software (or its derivatives) in
% your own work requires that you at least cite [1] and either [2] or [3].
%
% V1.1: Rodrigo A. Lobos (rlobos@umich.edu), 
% Chin-Cheng Chan (chinchen@usc.edu),
% Justin P. Haldar (jhaldar@usc.edu).
%
% May, 2023.  
%
% This software is Copyright �2023 The University of Southern California.
% All Rights Reserved. See the accompanying license.txt for additional
% license information.
%
% =========================================================================
clear;
close all;

%% Loading data
load T1_data

kData = double(kData);
[N1, N2, Nc] = size(kData);

figure; imagesc(mdisp(abs(fftshift(ifft2(ifftshift(kData))))));...
    axis image; axis tight; axis off; colormap gray;...
    title(['Data in the spatial domain']); caxis([0 1e-8]);

%% Selection of calibration data
cal_length = 32; % Length of each dimension of the calibration data

center_x = ceil(N1/2)+even(N1);
center_y = ceil(N2/2)+even(N2);

cal_index_x = center_x + [-floor(cal_length/2):floor(cal_length/2)-even(cal_length/2)];
cal_index_y = center_y + [-floor(cal_length/2):floor(cal_length/2)-even(cal_length/2)];

kCal = kData(cal_index_x,cal_index_y, :);

%% Nullspace-based algorithm parameters
dim_sens = [N1, N2]; % Desired dimensions for the estimated sensitivity maps

tau = 3; % Kernel radius

threshold = 0.08; % Theshold for C-matrix singular values

M = 10; % Number of iteration for Power Iteration

PowerIteration_flag_convergence = []; % If equal to 1 a convergence error 
%                                      is displayed for Power Iteration if 
%                                      the method has not converged for 
%                                      some voxels after the iterations 
%                                      indicated by the user. In this example
%                                      it corresponds to an empty array which 
%                                      indicates that the default value is 
%                                      being used, which is equal to 1.

PowerIteration_flag_auto = []; % If equal to 1 Power Iteration is run until
%                               convergence in case the number of
%                               iterations indicated by the user is too
%                               small. 
%                               In this example this variable corresponds 
%                               to an empty array which indicates that the
%                               default value is being used, which 
%                               is equal to 0.

interp_zp = []; % Amount of zero-padding to create the low-resolution grid 
%                 if FFT-interpolation is used. In this example it
%                 corresponds to an empty array which indicates that the
%                 default value is being used.

gauss_win_param = []; %     Parameter needed for the Gaussian apodizing 
%                           window used to generate the low-resolution 
%                           image in the FFT-based interpolation approach.
%                           This corresponds to the reciprocal value of
%                           the standard deviation of the Gaussian window. 
%                           In this example it corresponds to an empty 
%                           array which indicates that the default value is 
%                           being used.

%% PISCO techniques

% The following techniques are used if the corresponding binary variable is
% equal to 1

kernel_shape = 1; % An ellipsoidal shape is adopted for the calculation of 
%                   kernels (instead of a rectangular shape)

FFT_nullspace_C_calculation = 1; % FFT-based calculation of nullspace 
%                                  vectors of C by calculating C'*C directly
%                                  (instead of calculating C first)

PowerIteration_G_nullspace_vectors = 1; % A PowerIteration approach is used
%                                         to find nullspace vectors of the 
%                                         G matrices (instead of using SVD)


FFT_interpolation = 1; % Sensitivity maps are calculated on a small spatial
%                        grid and then interpolated to a grid with nominal 
%                        dimensions using an FFT-approach

verbose = 1; % If equal to 1 then PISCO information is displayed

%% PISCO estimation

[senseMaps, eigenValues] = PISCO_senseMaps_estimation(kCal,dim_sens,...
    tau, ...
    threshold, ...
    kernel_shape, ...
    FFT_nullspace_C_calculation, ...
    PowerIteration_G_nullspace_vectors, M, PowerIteration_flag_convergence, PowerIteration_flag_auto, ...
    FFT_interpolation, interp_zp, gauss_win_param, ...
    verbose);

%% Support mask created from the last eigenvalues of the G matrices 

threshold_mask = 0.05;

eig_mask = zeros(N1,N2);
eig_mask(find(eigenValues(:,:,end) < threshold_mask)) = 1;

% Optional masking step

senseMaps_masked= senseMaps.*eig_mask;

%% Estimated Sensitivity Maps 

figure; imagesc(mdisp(abs(senseMaps))); axis tight; axis image; axis off;...
colormap gray; title('Estimated sensitivity maps');

figure; imagesc(mdisp(abs(senseMaps_masked))); axis tight; axis image; axis off;...
    colormap gray; title('Masked sensitivity maps');

if PowerIteration_G_nullspace_vectors == 1 
    title_eig_values = 'Smallest eigenvalue of normalized G matrices (spatial map)';
    figure; imagesc(eigenValues); axis tight; axis image; colormap gray; colorbar; title(title_eig_values); 
else
    title_eig_values = 'Eigenvalues of normalized G matrices (spatial maps)';
    figure; imagesc(mdisp(eigenValues)); axis tight; axis image; colormap gray; colorbar; title(title_eig_values); 
end

figure; imagesc(eig_mask); axis tight; axis image; colormap gray; title('Support mask');

%% Auxiliary functions

function result = even(int)
    result = not(rem(int,2));
end

function tempForDisplay = mdisp(x)
% Function used to visualize multichannel images.
% If the input corresponds to a 3D array of dimensions N1 x N2 x Nc, the
% output corresponds to a 2D array that displays Nc images of dimension
% N1 x N2.

    [N1,N2,Nc] = size(x);
    f = factor(Nc);if numel(f)==1;f = [1,f];end
    tempForDisplay = reshape(permute(reshape(x,[N1,N2,prod(f(1:floor(numel(f)/2))),...
        prod(f(floor(numel(f)/2)+1:end))]),[1,3,2,4]),[N1*prod(f(1:floor(numel(f)/2))),...
        N2*prod(f(floor(numel(f)/2)+1:end))]); 
end