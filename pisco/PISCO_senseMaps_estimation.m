function [senseMaps, eigenValues] = PISCO_senseMaps_estimation(kCal, dim_sens, ...
    tau, threshold, kernel_shape, FFT_nullspace_C_calculation, PowerIteration_G_nullspace_vectors, ...
    M, PowerIteration_flag_convergence, PowerIteration_flag_auto, FFT_interpolation, interp_zp, gauss_win_param, verbose)

% Function that estimates coil sensitivity maps using the nullspace-based 
% algorithm and the PISCO techniques described in the technical report 
%
% [1] R. A. Lobos, C.-C. Chan, J. P. Haldar.  PISCO Software Version 1.0 
% University of Southern California, Los Angeles, CA, Technical Report 
% USC-SIPI-458. 2023.
%
% and in the papers 
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
% This software is available at
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
% The current implementation only covers sensitivity map estimation in the
% 2D case.
%
% =========================================================================
%
% Input parameters:
%    --kCal:                An N1_cal x N2_cal x Nc block of calibration data, 
%                           where N1_cal and N2_cal are the dimensions of 
%                           a rectangular block of Nyquist-sampled 
%                           k-space, and Nc is the number of channels
%                           in the array.
%
%    --dim_sens:            An 1x2 input array with the desired dimensions
%                           of the output sensitivity matrices.
%
%    --tau:                 Parameter (in Nyquist units) that determines
%                           the size of the k-space kernel. For a
%                           rectangular kernel the size corresponds to
%                           (2*tau+1) * (2*tau+1). For an ellipsoidal kernel it
%                           corresponds to the radius of the associated
%                           neighborhood.
%                           Defaults to tau=3 if not specified.
%
%    --threshold:           Specifies how small a singular value needs to 
%                           be (relative to the maximum singular value) 
%                           before its associated singular vector is
%                           considered to be in the nullspace of the 
%                           C-matrix. 
%                           Defaults to 0.05 if not specified.
% 
%    --kernel_shape:        Binary variable. If equal to 0, then a
%                           rectangular kernel is used. If 
%                           equal to 1, an ellipsoidal kernel is used
%                           Default to 1 if not specified.
% 
%    --FFT_nullspace_C_calculation:  Binary variable. If equal to 0, the 
%                           nullspace vectors of C are calculated from C'*C
%                           by calculating C first. If equal to 1, the 
%                           nullspace vectors of C are calculated from 
%                           C'*C, which is calculated directly using an 
%                           FFT-based approach. 
%                           Default to 1 if not specified.
% 
%    --PowerIteration_G_nullspace_vectors:  
%                           Binary variable. If equal to 0, the nullspace 
%                           vectors of the G matrices are calculated using 
%                           SVD. If equal to 1, the nullspace vectors of 
%                           the G matrices are calculated using a Power
%                           Iteration approach.
%                           Default to 1 if not specified.
%
% 
%    --M:                   Number of iterations used in the Power Iteration
%                           approach to calculate the nullspace vectors
%                           of the G matrices.
%                           Default to M = 10 if not specified.
%
%    --PowerIteration_flag_convergence:  
%                           Binary variable. If equal to 1 a convergence error 
%                           is displayed for Power Iteration if the method 
%                           has not converged for some voxels after the 
%                           iterations indicated by the user.
%                           Default to 1 if not specified.
%
%    --PowerIteration_flag_auto:  
%                           Binary variable. If equal to 1 Power Iteration 
%                           is run until convergence in case the number of
%                           iterations indicated by the user is too small.
%                           Default to 0 if not specified.
% 
%    --FFT_interpolation:   Binary variable. If equal to 0, no 
%                           interpolation is used. If equal to 1, the 
%                           FFT-based interpolation is used.
%
%    --interp_zp:           Amount of zero-padding to create the
%                           low-resolution grid if FFT-interpolation is
%                           used. The low-resolution grid has dimensions
%                           (N1_cal + interp_zp)x(N2_cal + interp_zp)xNc.
%                           Default to 24 if not specified.
%
%    --gauss_win_param:     Parameter needed for the Gaussian apodizing 
%                           window used to generate the low-resolution 
%                           image in the FFT-based interpolation approach.
%                           This corresponds to the reciprocal value of
%                           the standard deviation of the Gaussian window. 
%                           Default to 100 if not specified.
%
%    --verbose:             Binary variable. If equal to 1, then PISCO
%                           information is displayed. It indicates which
%                           PISCO techniques are employed and the
%                           computation times of each step in the
%                           sensitivity map estimation procedure.
%
% Output parameters:
%    --senseMaps:           An dim_sens(1) x dim_sens(2) x Nc stack corresponding to 
%                           the sensitivity maps for each channel that was
%                           present in the calibration data.
%
%
%    --eigenValues:         An dim_sens(1) x dim_sens(2)x Nc containing the value of 
%                           the eigenvalues of G(x) for each spatial location 
%                           (normalized by the kernel size). Can be used for 
%                           creating a mask describing the image support 
%                           (e.g., mask =  (eigenValues(:,:,end) < 0.08);). If 
%                           PowerIteration_G_nullspace_vectors == 1, only 
%                           the smallest eigenvalue is returned. In this 
%                           case the dimensions of eigenVal are 
%                           dim_sens(1) x dim_sens(2). If FFT_interpolation == 1, 
%                           approximations of eigenvalues are returned. 

if nargin < 3 || not(isnumeric(tau)) || not(numel(tau))
    tau = 3;
end

if nargin < 4 || not(isnumeric(threshold)) || not(numel(threshold))
    threshold = 0.05;
end

if nargin < 5 || not(isnumeric(kernel_shape)) || not(numel(kernel_shape))
    kernel_shape = 1;
end

if nargin < 6 || not(isnumeric(FFT_nullspace_C_calculation)) || not(numel(FFT_nullspace_C_calculation))
    FFT_nullspace_C_calculation = 1;
end

if nargin < 7 || not(isnumeric(PowerIteration_G_nullspace_vectors)) || not(numel(PowerIteration_G_nullspace_vectors))
    PowerIteration_G_nullspace_vectors = 1;
end

if nargin < 8 || not(isnumeric(M)) || not(numel(M))
    M = 10;
end

if nargin < 9 || not(isnumeric(PowerIteration_flag_convergence)) || not(numel(PowerIteration_flag_convergence))
    PowerIteration_flag_convergence = 0; % EDITED
end

if nargin < 10 || not(isnumeric(PowerIteration_flag_auto)) || not(numel(PowerIteration_flag_auto))
    PowerIteration_flag_auto = 0;
end

if nargin < 11 || not(isnumeric(FFT_interpolation)) || not(numel(FFT_interpolation))
    FFT_interpolation = 1;
end

if nargin < 12 || not(isnumeric(interp_zp)) || not(numel(interp_zp))
    interp_zp = 24;
end

if nargin < 13 || not(isnumeric(gauss_win_param)) || not(numel(gauss_win_param))
    gauss_win_param = 100;
end

if nargin < 14
    verbose = 1;
end

if verbose == 1
    
    if kernel_shape ==0
    kernel_shape_q = 'Rectangular';
    else
        kernel_shape_q = 'Ellipsoidal';
    end

    if FFT_nullspace_C_calculation == 0
        FFT_nullspace_C_calculation_q = 'No';
    else
        FFT_nullspace_C_calculation_q = 'Yes';
    end

    if FFT_interpolation == 0
        FFT_interpolation_q = 'No';
    else
        FFT_interpolation_q = 'Yes';
    end

    if PowerIteration_G_nullspace_vectors == 0
        PowerIteration_nullspace_vectors_q = 'No';
    else
        PowerIteration_nullspace_vectors_q = 'Yes';
    end
    
    disp('Selected PISCO techniques:')
    disp('=======================')
    disp(['Kernel shape : ' kernel_shape_q])
    disp(['FFT-based calculation of nullspace vectors of C : ' FFT_nullspace_C_calculation_q])
    disp(['FFT-based interpolation : ' FFT_interpolation_q])
    disp(['PowerIteration-based nullspace estimation for G matrices : ' PowerIteration_nullspace_vectors_q])
    disp('=======================')

end

t_null = tic;

% ==== Nullspace-based algorithm Steps (1) and (2)  ====

% Calculation of nullspace vectors of C 

t_null_vecs = tic;

U = nullspace_vectors_C_matrix(kCal, tau, threshold, kernel_shape, FFT_nullspace_C_calculation);

t_null_vecs = toc(t_null_vecs);

if verbose == 1

    if FFT_nullspace_C_calculation == 0
        aux_word = 'Calculating C first';
    else
        aux_word = 'FFT-based direct calculation of ChC';
    end
    
    disp('=======================')
    disp('PISCO computation times (secs):')
    disp('=======================')
    disp(['Time nullspace vectors of C (' aux_word ') : ' num2str(t_null_vecs)]) 
    disp('=======================')

end

% ==== Nullspace-based algorithm Step (3)  ====

% Direct computation of G matrices 

t_G_matrices = tic;

G = G_matrices(kCal, dim_sens(1), dim_sens(2), tau, U, kernel_shape, FFT_interpolation, interp_zp);

t_G_matrices = toc(t_G_matrices );

Nc = size(kCal,3);
patchSize = size(U,1)/Nc;
clear U

if verbose == 1

    disp(['Time G matrices (direct calculation): ' num2str(t_G_matrices )]) 
    disp('=======================')
    
end

% ==== Nullspace-based algorithm Step (4)  ====

% Calculation of nullspace vectors of the G matrices

t_null_G = tic;

[senseMaps, eigenValues] = nullspace_vectors_G_matrix(kCal, dim_sens(1), dim_sens(2), ...
    G, patchSize, PowerIteration_G_nullspace_vectors, M, PowerIteration_flag_convergence, PowerIteration_flag_auto,...
    FFT_interpolation, gauss_win_param, verbose);

t_null_G = toc(t_null_G);

clear G

if verbose == 1

    if PowerIteration_G_nullspace_vectors == 0
        aux_word = 'Using SVD';
    else
        aux_word = 'Using Power Iteration';
    end

    disp(['Time nullspace vector G matrices (' aux_word ') : ' num2str(t_null_G)]) 
    disp('=======================')
    
end

% ==== Nullspace-based algorithm Step (5)  ====

%  Normalization 
    
senseMaps = senseMaps.*repmat((exp(-complex(0,1)*angle(senseMaps(:,:,1)))), [1 1 Nc]); %Final maps after phase referencing w.r.t the first channel

senseMaps = senseMaps./repmat(sqrt(sum(abs(senseMaps).^2,3)), [1 1 Nc]);

if verbose == 1

    disp(['Total time: ' num2str(toc(t_null))]) 
    disp('=======================')
    
end

end

%% Extra functions

function b = vect( a )
    b = a(:);
end

function result = even(int)
    result = not(rem(int,2));
end

function result = C_matrix(x, N1, N2, Nc, tau, kernel_shape)

% Function that calculates the C matrix.
%
% Input parameters:
%    --X:                   N1 x N2 x Nc Nyquist-sampled k-space data, 
%                           where N1 and N2 are the data dimensions, and Nc
%                           is the number of channels in the array.
%
%    --N1, N2:              Dimensions of the k-space data.
%
%    --Nc:                  Number of channels of the k-space data
%
%    --tau:                 Parameter (in Nyquist units) that determines
%                           the size of the k-space kernel.For a
%                           rectangular kernel the size corresponds to
%                           (2*tau+1) * (2*tau+1). For an ellipsoial kernel it
%                           corresponds to the radius of the associated
%                           neighborhood.
%                           Defaults to tau=3 if not specified.
% 
%    --kernel_shape:        Binary variable. If equal to 0, then a
%                           rectangular kernel is used. If 
%                           equal to 1, an ellipsoidal kernel is used
%                           Default to 1 if not specified.

    x = reshape(x,N1*N2,Nc);

    [in1,in2] = meshgrid(-tau:tau,-tau:tau);
    
    if kernel_shape == 1
        i = find(in1.^2+in2.^2<=tau^2);
    else
        i = [1:numel(in1)]; 
    end     

    in1 = in1(i)';
    in2 = in2(i)';

    patchSize = numel(in1);

    result = zeros((N1-2*tau-even(N1))*(N2-2*tau-even(N2)),patchSize*Nc,'like',x);

    k = 0;
    for i = tau+1+even(N1):N1-tau
        for j = tau+1+even(N2):N2-tau
            k = k+1;
            ind = sub2ind([N1,N2],i+in1,j+in2);
            result(k,:) = vect(x(ind,:));
        end
    end

end

function PhP = ChC_FFT_convolutions(X, N1, N2, Nc, tau, pad, kernel_shape)

% Function that directly calculates the matrix C'*C using an FFT-based
% approach.
%
% Input parameters:
%    --X:                   N1 x N2 x Nc Nyquist-sampled k-space data, 
%                           where N1 and N2 are the data dimensions, and Nc
%                           is the number of channels in the array.
%
%    --N1, N2:              Dimensions of the k-space data.
%
%    --Nc:                  Number of channels of the k-space data
%
%    --tau:                 Parameter (in Nyquist units) that determines
%                           the size of the k-space kernel.For a
%                           rectangular kernel the size corresponds to
%                           (2*tau+1) * (2*tau+1). For an ellipsoial kernel it
%                           corresponds to the radius of the associated
%                           neighborhood.
%                           Defaults to tau=3 if not specified.
%
%    --pad:                 Binary variable. If equal to 1, then
%                           zero-padding is employed when calculating FFTs.
% 
% 
%    --kernel_shape:        Binary variable. If equal to 0, then a
%                           rectangular kernel is used. If 
%                           equal to 1, an ellipsoidal kernel is used
%                           Default to 1 if not specified.
%
%

[in1,in2] = meshgrid(-tau:tau,-tau:tau);
if kernel_shape == 1
    i = find(in1.^2+in2.^2<=tau^2);
else
    i = [1:numel(in1)];
end
in1 = in1(i(:));
in2 = in2(i(:));

patchSize = numel(i);

if pad
    N1n = 2^(ceil(log2(N1+2*tau))); 
    N2n = 2^(ceil(log2(N2+2*tau)));
else
    N1n = N1;
    N2n = N2;
end

inds = sub2ind([N1n,N2n], floor(N1n/2)+1-in1+in1', floor(N2n/2)+1-in2+in2');

[n2,n1] = meshgrid([-floor(N2n/2):floor(N2n/2)-even(N2n/2)]/N2n,[-floor(N1n/2):floor(N1n/2)-even(N1n/2)]/N1n);
phaseKernel = exp(complex(0,-2*pi)*(n1*(ceil(N1n/2)+tau)+n2*(ceil(N2n/2)+tau)));
cphaseKernel = exp(complex(0,-2*pi)*(n1*(ceil(N1n/2))+n2*(ceil(N2n/2))));

x = fft2(X,N1n,N2n).*phaseKernel;   

PhP = zeros(patchSize, patchSize, Nc,  Nc); 
for q = 1:Nc
    b= reshape(ifft2(conj(x(:,:,q:Nc)).*x(:,:,q).*cphaseKernel),[],Nc-q+1); 
    PhP(:,:,q:Nc,q) = reshape(b(inds,:),patchSize,patchSize,Nc-q+1);
    PhP(:,:,q,q+1:Nc) = permute(conj(PhP(:,:,q+1:Nc,q)),[2,1,4,3]);
end
PhP = reshape(permute(PhP,[1,3,2,4]), patchSize*Nc, patchSize*Nc);

end


function U = nullspace_vectors_C_matrix(kCal, tau, threshold, kernel_shape, FFT_nullspace_C_calculation)

% Function that returns the nullspace vectors of the C matrix.
%
% Input parameters:
%    --kCal:                An N1_cal x N2_cal x Nc block of calibration data, where
%                           N1_cal and N2_cal are the dimensions of a 
%                           rectangular block of Nyquist-sampled k-space, 
%                           and Nc is the number of channels in the array.
%
%    --tau:                 Parameter (in Nyquist units) that determines
%                           the size of the k-space kernel.For a
%                           rectangular kernel the size corresponds to
%                           (2*tau+1) * (2*tau+1). For an ellipsoial kernel it
%                           corresponds to the radius of the associated
%                           neighborhood.
%                           Defaults to tau=3 if not specified.
%
%    --threshold:           Specifies how small a singular value needs to 
%                           be (relative to the maximum singular value) 
%                           before its associated singular vector is
%                           considered to be in the nullspace of the 
%                           C-matrix. 
%                           Defaults to 0.05 if not specified.
% 
%    --kernel_shape:        Binary variable. If equal to 0, then a
%                           rectangular kernel is used. If 
%                           equal to 1, an ellipsoidal kernel is used
%                           Default to 1 if not specified.
% 
%    --FFT_nullspace_C_calculation:  
%                            Binary variable. If equal to 0, the 
%                            nullspace vectors of the C matrix are calculated 
%                            from C'*C by calculating C first. If equal to 1, 
%                            the nullspace vectors of the C matrix are 
%                            calculated from C'*C, which is calculated directly 
%                            using an  FFT-based approach. 
%                            Default to 1 if not specified.
%
% Output parameters:
%    --U:                    Matrix which columns corresponds to the
%                            nullspace vectors of the C matrix. 
%

if nargin < 2 || not(isnumeric(tau)) || not(numel(tau))
    tau = 3;
end

if nargin < 3 || not(isnumeric(threshold)) || not(numel(threshold))
    threshold = 0.05;
end

if nargin < 4 || not(isnumeric(kernel_shape)) || not(numel(kernel_shape))
    kernel_shape = 1;
end

if nargin < 5 || not(isnumeric(FFT_nullspace_C_calculation)) || not(numel(FFT_nullspace_C_calculation))
    FFT_nullspace_C_calculation = 1;
end

if FFT_nullspace_C_calculation == 0
    
    C = C_matrix(kCal(:), size(kCal,1), size(kCal,2), size(kCal,3), tau, kernel_shape);
       
    ChC = C'*C;
    clear C
    
else
    
    ChC = ChC_FFT_convolutions(kCal, size(kCal,1), size(kCal,2), size(kCal,3), tau, 1, kernel_shape);
      
end

[~,Sc,U] = svd(ChC,'econ');
clear ChC
sing = diag(Sc);
clear Sc
 
sing = sqrt(sing);
sing  = sing/sing(1);

Nvect = find(sing >=threshold*sing(1),1,'last');
clear sing
U = U(:, Nvect+1:end); 

end

function G = G_matrices(kCal, N1, N2, tau, U, kernel_shape, FFT_interpolation, interp_zp)

% Function that calculates the G(x) matrices directly without calculating
% H(x) first.
%
% Input parameters:
%    --kCal:                An N1_cal x N2_cal x Nc block of calibration data, where
%                           N1_cal and N2_cal are the dimensions of a 
%                           rectangular block of Nyquist-sampled k-space, 
%                           and Nc is the number of channels in the array.
%
%    --N1, N2:              The desired dimensions of the output sensitivity
%                           matrices.  
%
%    --tau:                 Parameter (in Nyquist units) that determines
%                           the size of the k-space kernel. For a
%                           rectangular kernel the size corresponds to
%                           (2*tau+1) * (2*tau+1). For an ellipsoial kernel it
%                           corresponds to the radius of the associated
%                           neighborhood.
%                           Defaults to tau=3 if not specified.

%    --U:                   Matrix which columns correspond to the
%                           nullspace vectors of the C matrix.
%
%    --kernel_shape:        Binary variable. If equal to 0, then a
%                           rectangular kernel is used. If 
%                           equal to 1, an ellipsoidal kernel is used
%                           Default to 1 if not specified.
% 
%    --FFT_interpolation:   Binary variable. If equal to 0, no interpolation
%                           is used. 
%                           If equal to 1, the FFT-based interpolation is
%                           used.
%                           Default to 1 if not specified.
%
%    --interp_zp:           Amount of zero-padding to create the
%                           low-resolution grid if FFT-interpolation is
%                           used. The low-resolution grid has dimensions
%                           (N1_acs + interp_zp)x(N2_acs + interp_zp)xNc.
%                           Default to 24 if not specified.
%
% Output parameters:
%    --G:                   An N1 x N2 x Nc x Nc array where G[i,j,:,:]
%                           corresponds to the G matrix at the (i,j) 
%                           spatial location.

if nargin < 4 || not(isnumeric(tau)) || not(numel(tau))
    tau = 3;
end

if nargin < 6 || not(isnumeric(kernel_shape)) || not(numel(kernel_shape))
    kernel_shape = 1;
end

if nargin < 7 || not(isnumeric(FFT_interpolation)) || not(numel(FFT_interpolation))
    FFT_interpolation = 1;
end

if nargin < 8 || not(isnumeric(interp_zp)) || not(numel(interp_zp))
    interp_zp = 24;
end

[N1_cal, N2_cal, Nc] = size(kCal);

[in1,in2] = meshgrid(-tau:tau,-tau:tau);

if kernel_shape == 0 

    ind = [1:numel(in1)]'; 
    
else 
    
    ind = find(in1.^2+in2.^2<=tau^2); 
    
end

in1 = in1(ind)';
in2 = in2(ind)';

patchSize = numel(in1);

in1 = in1(:);
in2 = in2(:);

eind = [patchSize:-1:1]';

G = zeros(2*(2*tau+1)* 2*(2*tau+1),Nc,Nc);
W = U*U';
clear U;
W = permute(reshape(W,patchSize,Nc,patchSize,Nc),[1,2,4,3]);

for s = 1:patchSize 
    G(sub2ind([2*(2*tau+1),2*(2*tau+1)],2*tau+1+1+in1(eind)+in1(s),2*tau+1+1+in2(eind)+in2(s)),:,:) = ...
        G(sub2ind([2*(2*tau+1),2*(2*tau+1)],2*tau+1+1+in1(eind)+in1(s),2*tau+1+1+in2(eind)+in2(s)),:,:)  + W(:,:,:,s);
end

clear W

if FFT_interpolation == 0
    
    N1_g = N1;
    N2_g = N2;
    
else
  
    if N1_cal <= N1 - interp_zp
        N1_g = N1_cal + interp_zp;
    else   
        N1_g = N1_cal;
    end

    if N2_cal <= N2 - interp_zp
        N2_g = N2_cal + interp_zp;
    else   
        N2_g = N2_cal;
    end
    
end

[n2,n1] = meshgrid([-N2_g/2:N2_g/2-1]/N2_g,[-N1_g/2:N1_g/2-1]/N1_g);
phaseKernel = exp(complex(0,-2*pi)*(n1*(N1_g-2*tau-1)+n2*(N2_g-2*tau-1)));

G = fft2(conj(reshape(G,2*(2*tau+1),2*(2*tau+1),Nc,Nc)),N1_g,N2_g).*phaseKernel; 

G = fftshift(fftshift(G,1),2);

end

function [senseMaps, eigenVal] = nullspace_vectors_G_matrix(kCal, N1, N2, ...
    G, patchSize, PowerIteration_G_nullspace_vectors, M, PowerIteration_flag_convergence, PowerIteration_flag_auto, ...
    FFT_interpolation, gauss_win_param, verbose)

% Function that calculates the nullspace vectors for each G(x) matrix. These
% vectors correspond to sensitivity maps at the x location.
%
% Input parameters:
%    --kCal:                An N1_cal x N2_cal x Nc block of calibration data, where
%                           N1_cal and N2_cal are the dimensions of a 
%                           rectangular block of Nyquist-sampled k-space, 
%                           and Nc is the number of channels in the array.
%
%    --N1, N2:              The desired dimensions of the output sensitivity
%                           matrices.  
%
%    --G:                   An N1_g x N2_g x Nc x Nc array where G[i,j,:,:]
%                           corresponds to the G matrix at the (i,j) 
%                           spatial location.
%
%    --patchSize:           Number of elements in the kernel used to calculate
%                           the nullspace vectors of the C matrix.
% 
% 
%    --PowerIteration_G_nullspace_vectors:  
%                           Binary variable. If equal to 0, the nullspace vectors  
%                           of the G matrices are calculated using SVD. If
%                           equal to 1, the nullspace vectors of the G
%                           matrices are calculated using the Power
%                           Iteration approach.
%                           Default to 1 if not specified.
% 
%    --M:                   Number of iterations used in the Power Iteration
%                           approach to calculate the nullspace vectors
%                           of the G matrices.
%                           Default to M = 10 if not specified.
%
%    --PowerIteration_flag_convergence:  
%                           Binary variable. If equal to 1 a convergence error 
%                           is displayed for Power Iteration if the method 
%                           has not converged for some voxels after the 
%                           iterations indicated by the user.
%                           Default to 1 if not specified.
%
%    --PowerIteration_flag_auto:  
%                           Binary variable. If equal to 1 Power Iteration 
%                           is run until convergence in case the number of
%                           iterations indicated by the user is too small.
%                           Default to 0 if not specified.
% 
%    --FFT_interpolation:   Binary variable. If equal to 0, no interpolation
%                           is used. If equal to 1, the FFT-based 
%                           interpolation is used.
%                           Default to 1 if not specified
%
%    --gauss_win_param:     Parameter needed for the Gaussian apodizing 
%                           window used to generate the low-resolution 
%                           image in the FFT-based interpolation approach.
%                           This corresponds to the reciprocal value of
%                           the standard deviation of the Gaussian window. 
%                           Default to 100 if not specified.  
%
%    --verbose:             Binary variable. If equal to 1, then information
%                           about the convergence of Power Iterationis is 
%                           displayed.
%
% Output parameters:
%    --senseMaps:           An N1 x N2 x Nc stack corresponding to the
%                           sensitivity maps for each channel that was
%                           present in the calibration data.
%
%    --eigenVal:            An N1 x N2 x Nc containing the value of the
%                           eigenvalues of G(x) for each location (normalized).  
%                           Can be used for creating a mask describing the 
%                           image support
%                           (e.g., mask = (eigenVal(:,:,end) < 0.08);)
%                           If PowerIteration_G_nullspace_vectors == 1, 
%                           only the smallest eigenvalue is returned. In 
%                           this case the dimensions of eigenVal are N1 x N2.
%                           If FFT_interpolation == 1, approximations of 
%                           the smallest eigenvalues are returned.
%

if nargin < 6 || not(isnumeric(PowerIteration_G_nullspace_vectors)) || not(numel(PowerIteration_G_nullspace_vectors))
    PowerIteration_G_nullspace_vectors = 1;
end

if nargin < 7 || not(isnumeric(M)) || not(numel(M))
    M = 10;
end

if nargin < 8 || not(isnumeric(PowerIteration_flag_convergence)) || not(numel(PowerIteration_flag_convergence))
    PowerIteration_flag_convergence = 1;
end

if nargin < 9 || not(isnumeric(PowerIteration_flag_auto)) || not(numel(PowerIteration_flag_auto))
    PowerIteration_flag_auto = 0;
end

if nargin < 10 || not(isnumeric(FFT_interpolation)) || not(numel(FFT_interpolation))
    FFT_interpolation = 1;
end

if nargin < 11 || not(isnumeric(gauss_win_param)) || not(numel(gauss_win_param))
    gauss_win_param = 100;
end

if nargin < 12 || not(isnumeric(verbose)) || not(numel(verbose))
    verbose = 1;
end

if PowerIteration_flag_auto == 1
    PowerIteration_flag_convergence = 0;
end

N1_g = size(G,1);
N2_g = size(G,2);
Nc = size(G,3);

senseMaps = zeros(N1_g,N2_g,Nc);

if PowerIteration_G_nullspace_vectors == 0
    
    eigenVal = zeros(N1_g,N2_g,Nc);
    
    for i = 1:N1_g
        for j = 1:N2_g 
            [~,S,Vh] = svd(squeeze(G(i,j,:,:)),'econ');
            senseMaps(i,j,:) = reshape(Vh(:,end),[1,1,Nc])*exp(-complex(0,angle(Vh(1,end))));
            eigenVal(i,j,:) = reshape(abs(diag(S)), [1 1 Nc]);
            clear S Vh
        end
    end
    
    clear G
    
    eigenVal = eigenVal/patchSize;
    
else
    
    G = G/patchSize;
    
    G_null = zeros(size(G));

    for i = 1:Nc
        G_null(:,:,i,i) = 1;
    end

    G_null = permute(G_null - G, [1 2 4 3]);

    clear G

    if PowerIteration_flag_convergence == 0 && PowerIteration_flag_auto == 0

        senseMaps = randn(size(G_null,1), size(G_null,2), size(G_null,3)) +1i*randn(size(G_null,1), size(G_null,2), size(G_null,3));
    
        senseMaps_norm = sqrt(sum(abs(senseMaps).^2, 3));
        z_out_norm = repmat(senseMaps_norm, [1 1 Nc]);
        senseMaps = senseMaps./z_out_norm;

        for m = 1:M
            senseMaps = squeeze(sum( G_null.*repmat(senseMaps, [1, 1, 1, Nc]) , 3) );

            senseMaps_norm = sqrt(sum(abs(senseMaps).^2, 3));
            z_out_norm = repmat(senseMaps_norm, [1 1 Nc]);
            senseMaps = senseMaps./z_out_norm; 
    
            if m == M
                aux1 = squeeze(sum( G_null.*repmat(senseMaps, [1, 1, 1, Nc]) , 3) );
                final_maps_norm = sqrt(sum(abs(aux1).^2, 3));
            end   
        end
        
        clear G_null 

        eigenVal = 1 - final_maps_norm;

    end

    if PowerIteration_flag_convergence == 1 || PowerIteration_flag_auto == 1

        senseMaps = randn(size(G_null,1), size(G_null,2), size(G_null,3)) +1i*randn(size(G_null,1), size(G_null,2), size(G_null,3));
        eigenVec2 = randn(size(G_null,1), size(G_null,2), size(G_null,3)) +1i*randn(size(G_null,1), size(G_null,2), size(G_null,3));
        
        senseMaps_norm = sqrt(sum(abs(senseMaps).^2, 3));
        z_out_norm = repmat(senseMaps_norm, [1 1 Nc]);
        senseMaps = senseMaps./z_out_norm;
    
        eigenVec2_norm = sqrt(sum(abs(eigenVec2).^2, 3));
        z_out_norm2 = repmat(eigenVec2_norm, [1 1 Nc]);
        eigenVec2 = eigenVec2./z_out_norm2;
    
    
        for m = 1:M
            senseMaps = squeeze(sum( G_null.*repmat(senseMaps, [1, 1, 1, Nc]) , 3) );
            eigenVec2 = squeeze(sum( G_null.*repmat(eigenVec2, [1, 1, 1, Nc]) , 3) );
            
            senseMaps_norm = sqrt(sum(abs(senseMaps).^2, 3));
            z_out_norm = repmat(senseMaps_norm, [1 1 Nc]);
            senseMaps = senseMaps./z_out_norm; 
    
            inner_prod = sum(eigenVec2.*conj(senseMaps), 3);
    
            fact_proj_mc = repmat(inner_prod, [1 1 Nc]);
    
            eigenVec2 = eigenVec2 - fact_proj_mc.*senseMaps;
   
            eigenVec2_norm = sqrt(sum(abs(eigenVec2).^2, 3));    
            z_out_norm2 = repmat(eigenVec2_norm, [1 1 Nc]);   
            eigenVec2 = eigenVec2./z_out_norm2;
    
            if m == M
                aux1 = squeeze(sum( G_null.*repmat(senseMaps, [1, 1, 1, Nc]) , 3) );
                aux2 = squeeze(sum( G_null.*repmat(eigenVec2, [1, 1, 1, Nc]) , 3) );
    
                final_maps_norm = sqrt(sum(abs(aux1).^2, 3));
                final_maps_norm2 = sqrt(sum(abs(aux2).^2, 3));
            end
        end
    
        eigen1 = final_maps_norm;
        eigen2 = final_maps_norm2;
        

        if  FFT_interpolation == 0

            eigenVal = 1 - eigen1;

            threshold_mask = 0.075;
            
            support_mask = zeros(size(eigenVal));
            support_mask(find(eigenVal < threshold_mask)) = 1;
    
            ratioEig = (eigen2./eigen1).^M;
            ratio_small = support_mask.*ratioEig;
        
            th_ratio = 0.008;
        
            ratio_small(find(ratio_small <= th_ratio)) = 0;
            ratio_small(find(ratio_small > th_ratio)) = 1;
        
            flag_convergence_PI = sum(ratio_small(:)) > 0;
            
            
            if flag_convergence_PI == 1 && PowerIteration_flag_convergence == 1
                error(['Power Iteration might have not converged for some voxels within the support after the ' int2str(M)...
                ' iterations indicated by the user. Increasing the number of iterations is recommended. You can ignore this error by setting PowerIteration_flag_convergence = 0. '...
                'The number of needed iterations for convergence can be found automatically by setting PowerIteration_flag_auto = 1. '])   
            end

            if flag_convergence_PI == 0
                if verbose == 1
                     disp(['Most likely Power Iteration has converged for all the voxels within the support after the ' int2str(M) ' iterations indicated by the user.'])
                end
            end
            
        

        if PowerIteration_flag_auto == 1 && flag_convergence_PI == 1
            if verbose == 1
                warning('off','backtrace')
                warning(['Power Iteration might have not converged for some voxels within the support after the ' int2str(M)...
                                ' iterations indicated by the user. The number of iterations for the convergence of Power Iteration will be found automatically. You can turn off this option by setting PowerIteration_flag_auto = 0. ']) 
            end
            M_auto = M+1;
            while(flag_convergence_PI == 1)
                senseMaps = squeeze(sum( G_null.*repmat(senseMaps, [1, 1, 1, Nc]) , 3) );
                eigenVec2 = squeeze(sum( G_null.*repmat(eigenVec2, [1, 1, 1, Nc]) , 3) );
                
                senseMaps_norm = sqrt(sum(abs(senseMaps).^2, 3));
                z_out_norm = repmat(senseMaps_norm, [1 1 Nc]);
                senseMaps = senseMaps./z_out_norm; 
        
                inner_prod = sum(eigenVec2.*conj(senseMaps), 3);
        
                fact_proj_mc = repmat(inner_prod, [1 1 Nc]);
        
                eigenVec2 = eigenVec2 - fact_proj_mc.*senseMaps;
       
                eigenVec2_norm = sqrt(sum(abs(eigenVec2).^2, 3));    
                z_out_norm2 = repmat(eigenVec2_norm, [1 1 Nc]);   
                eigenVec2 = eigenVec2./z_out_norm2;
        
                aux1 = squeeze(sum( G_null.*repmat(senseMaps, [1, 1, 1, Nc]) , 3) );
                aux2 = squeeze(sum( G_null.*repmat(eigenVec2, [1, 1, 1, Nc]) , 3) );
    
                final_maps_norm = sqrt(sum(abs(aux1).^2, 3));
                final_maps_norm2 = sqrt(sum(abs(aux2).^2, 3));

                eigen1 = final_maps_norm;
                eigen2 = final_maps_norm2;

                eigenVal = 1 - eigen1;
                
                support_mask = zeros(size(eigenVal));
                support_mask(find(eigenVal < threshold_mask)) = 1;
        
                ratioEig = (eigen2./eigen1).^M_auto;
                ratio_small = support_mask.*ratioEig;
            
                ratio_small(find(ratio_small <= th_ratio)) = 0;
                ratio_small(find(ratio_small > th_ratio)) = 1;
            
                flag_convergence_PI = sum(ratio_small(:)) > 0;

                M_auto = M_auto + 1;
            end
            
            if verbose == 1
                disp(['Most likely Power Iteration has converged for all the voxels within the support. ' int2str(M_auto) ' iterations were needed.'] )
            end

        end

        end

     %   clear G_null

    end

    


end

% ==== FFT-based interpolation ====

if FFT_interpolation == 1

    [N1_cal, N2_cal, ~] = size(kCal);
    
    w_sm = [0.54 - 0.46*cos(2*pi*([0:(N1_g-1)]/(N1_g-1)))].';
    w_sm2 = [0.54 - 0.46*cos(2*pi*([0:(N2_g-1)]/(N2_g-1)))].';
    w_sm = w_sm*w_sm2';
    w_sm = repmat(w_sm, [1 1 Nc]); 

    if PowerIteration_G_nullspace_vectors == 1 && (PowerIteration_flag_convergence == 1 || PowerIteration_flag_auto == 1) 

        auxVal = 1 - eigen1;

        eigenVal = abs(fftshift(fftshift(ifft2(fftshift(fft2(ifftshift(auxVal))).*w_sm(:,:,end), N1, N2), 1), 2)); 

        eigenVal = eigenVal/max(eigenVal(:));

        threshold_mask = 0.075;
        
        support_mask = zeros(size(eigenVal));
        support_mask(find(eigenVal < threshold_mask)) = 1;
    
        eigen1_us = abs(fftshift(fftshift(ifft2(fftshift(fft2(ifftshift(eigen1))).*w_sm(:,:,end), N1, N2), 1), 2)); 
    
        eigen2_us = abs(fftshift(fftshift(ifft2(fftshift(fft2(ifftshift(eigen2))).*w_sm(:,:,end), N1, N2), 1), 2)); 
    
        ratioEig = (eigen2_us./eigen1_us).^M;
        ratio_small = support_mask.*ratioEig;
    
        th_ratio = 0.008;
    
        ratio_small(find(ratio_small <= th_ratio)) = 0;
        ratio_small(find(ratio_small > th_ratio)) = 1;
    
        flag_convergence_PI = sum(ratio_small(:)) > 0;
    
        if flag_convergence_PI == 1 && PowerIteration_flag_convergence == 1
            error(['Power Iteration might have not converged for some voxels within the support after the ' int2str(M)...
                ' iterations indicated by the user. Increasing the number of iterations is recommended. You can ignore this error by setting PowerIteration_flag_convergence = 0. '...
                'The number of needed iterations for convergence can be found automatically by setting PowerIteration_flag_auto = 1. '])   
        end

        if flag_convergence_PI == 0
            if verbose == 1
                 disp(['Most likely Power Iteration has converged for all the voxels within the support after the ' int2str(M) ' iterations indicated by the user.'])
            end
        end

    if PowerIteration_flag_auto == 1 && flag_convergence_PI == 1
        if verbose == 1
            warning('off','backtrace')
            warning(['Power Iteration might have not converged for some voxels within the support after the ' int2str(M)...
                    ' iterations indicated by the user. The number of iterations for the convergence of Power Iteration will be found automatically. You can turn off this option by setting PowerIteration_flag_auto = 0.'])
        end
        M_auto = M+1;
        while(flag_convergence_PI == 1)
            senseMaps = squeeze(sum( G_null.*repmat(senseMaps, [1, 1, 1, Nc]) , 3) );
            eigenVec2 = squeeze(sum( G_null.*repmat(eigenVec2, [1, 1, 1, Nc]) , 3) );
            
            senseMaps_norm = sqrt(sum(abs(senseMaps).^2, 3));
            z_out_norm = repmat(senseMaps_norm, [1 1 Nc]);
            senseMaps = senseMaps./z_out_norm; 
    
            inner_prod = sum(eigenVec2.*conj(senseMaps), 3);
    
            fact_proj_mc = repmat(inner_prod, [1 1 Nc]);
    
            eigenVec2 = eigenVec2 - fact_proj_mc.*senseMaps;
   
            eigenVec2_norm = sqrt(sum(abs(eigenVec2).^2, 3));    
            z_out_norm2 = repmat(eigenVec2_norm, [1 1 Nc]);   
            eigenVec2 = eigenVec2./z_out_norm2;
    
            aux1 = squeeze(sum( G_null.*repmat(senseMaps, [1, 1, 1, Nc]) , 3) );
            aux2 = squeeze(sum( G_null.*repmat(eigenVec2, [1, 1, 1, Nc]) , 3) );

            final_maps_norm = sqrt(sum(abs(aux1).^2, 3));
            final_maps_norm2 = sqrt(sum(abs(aux2).^2, 3));

            eigen1 = final_maps_norm;
            eigen2 = final_maps_norm2;

            auxVal = 1 - eigen1;

            eigenVal = abs(fftshift(fftshift(ifft2(fftshift(fft2(ifftshift(auxVal))).*w_sm(:,:,end), N1, N2), 1), 2)); 

            eigenVal = eigenVal/max(eigenVal(:));
            
            support_mask = zeros(size(eigenVal));
            support_mask(find(eigenVal < threshold_mask)) = 1;

            eigen1_us = abs(fftshift(fftshift(ifft2(fftshift(fft2(ifftshift(eigen1))).*w_sm(:,:,end), N1, N2), 1), 2)); 
    
            eigen2_us = abs(fftshift(fftshift(ifft2(fftshift(fft2(ifftshift(eigen2))).*w_sm(:,:,end), N1, N2), 1), 2));     
    
            ratioEig = (eigen2_us./eigen1_us).^M_auto;
            ratio_small = support_mask.*ratioEig;
        
            ratio_small(find(ratio_small <= th_ratio)) = 0;
            ratio_small(find(ratio_small > th_ratio)) = 1;
        
            flag_convergence_PI = sum(ratio_small(:)) > 0;

            M_auto = M_auto + 1;
        end

        if verbose == 1
            disp(['Most likely Power Iteration has converged for all the voxels within the support. ' int2str(M_auto) ' iterations were needed.'] )
        end

    end

    end

    if PowerIteration_G_nullspace_vectors == 1 && PowerIteration_flag_convergence == 0 && PowerIteration_flag_auto == 0

        eigenVal = abs(fftshift(fftshift(ifft2(fftshift(fft2(ifftshift(eigenVal))).*w_sm(:,:,end), N1, N2), 1), 2)); 

        eigenVal = eigenVal/max(eigenVal(:));

    end

    if PowerIteration_G_nullspace_vectors == 0

        eigenVal = abs(fftshift(fftshift(ifft2(fftshift(fft2(ifftshift(eigenVal))).*w_sm(:,:,end), N1, N2), 1), 2)); 

        eigenVal = eigenVal/max(eigenVal(:));

    end
    
    apodizing_window = gausswin(N1_g,gauss_win_param)*gausswin(N2_g,gauss_win_param)';
    
    imLowRes_cal = zeros(N1_g,N2_g,Nc);
    imLowRes_cal(ceil(N1_g/2)+even(N1_g)+[-ceil(N1_cal/2):floor(N1_cal/2)-even(N1_cal)],ceil(N2_g/2)+even(N2_g)+[-floor(N2_cal/2):floor(N2_cal/2)-even(N2_cal)],:) = kCal;
    imLowRes_cal = fftshift(ifft2(ifftshift(imLowRes_cal.*apodizing_window)));    

    cim = sum(conj(senseMaps).*imLowRes_cal,3)./sum(abs(senseMaps).^2,3); 

    senseMaps = senseMaps.*repmat((exp(complex(0,1)*angle(cim))), [1 1 Nc]); 

    senseMaps = fftshift(fftshift(ifft2(fftshift(fft2(ifftshift(senseMaps))).*w_sm, N1, N2), 1), 2); 

    

end

end


