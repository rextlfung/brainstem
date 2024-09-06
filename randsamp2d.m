%% randsamp2d
%
% Generate a 2D pseudo random sampling pattern specified by a fully sampled
% "ACS" region. Outside the ACS region, randomly sample with higher
% probability towards the center, drawn from a Gaussian distribution.
%
% Input arguments:
% N = [Ny Nz] = sampling grid dimensions. Positive vector of integers >= 1.
% R = [Ry Rz] = acceleration factors. Positive vector of integers >= 1.
% acs = [acs_y acs_z] = portions of k-space to be fully sampled as ACS
% region. Vector of numbers between 0 and 1.
% max_ky_step = the largest tolerable gap between subsequently sampled fast
% PE (ky) locations.
%
% Output arguments:
% omega = sampling mask. Boolean matrix of size [Ny Nz].
%
% Design constraints:
% 1. Even number of samples about the center of k-space along the fast PE
% (ky) direction to ensure consistent TE
% 2. If the randomly selected ky sampling locations are separated by gaps
% exceeding max_ky_step, move the sampling location with the nearest
% neightbor(s) to the midpoint of the largest gap.
%
% Last modified Sep 5th, 2024. Rex Fung

function omega = randsamp2d(N, R, acs, max_ky_step)
    % Unpack and assert input arguments
    Ny = N(1); Nz = N(2);
    Ry = R(1); Rz = R(2);
    acs_y = acs(1); acs_z = acs(2);
    assert(Ny >= 1 && Nz >= 1, 'Dimensions must be >= 1');
    assert(Ry >= 1 && Rz >= 1, 'Acceleration factors must be >= 1');
    assert(0 <= acs_y && acs_y < 1 && 0 <= acs_z && acs_z < 1,...
           'ACS portions must be between [0, 1)');
    assert(Ry*acs_y <= 1 && Rz*acs_z <= 1, 'Acceleration cannot exceed ACS portion');

    % Compute number of ACS lines
    Ny_acs = 2 * ceil(acs_y * Ny / 2);
    Nz_acs = 2 * ceil(acs_z * Nz / 2);

    % Split indices into ACS and non-ACS regions
    acs_indices_y = (Ny/2 - Ny_acs/2 + 1):(Ny/2 + Ny_acs/2);
    nacs_indices_y = [1:(Ny/2 - Ny_acs/2);
                      (Ny/2 + Ny_acs/2 + 1):Ny]; % Split into halves for evenness
    acs_indices_z = (Nz/2 - Nz_acs/2 + 1):(Nz/2 + Nz_acs/2);
    nacs_indices_z = [1:(Nz/2 - Nz_acs/2), (Nz/2 + Nz_acs/2 + 1):Nz];

    % Compute the number of non-ACS locations to sample based on R
    Ny_nacs = ceil((Ny/Ry - Ny_acs)/2)*2; % Ensure even number
    Nz_nacs = ceil(Nz/Rz - Nz_acs);

    % Genearte Gaussian pdf of sampling non-ACS locations
    w_y = normpdf(nacs_indices_y,mean(1:Ny),std(1:Ny));
    w_z = normpdf(nacs_indices_z,mean(1:Nz),std(1:Nz));

    % Generate sampling locations along kz
    nacs_indices_samp_z = datasample(nacs_indices_z, Nz_nacs, 'Replace', false, 'Weights',w_z);
    indices_z = sort([acs_indices_z, nacs_indices_samp_z]);

    % Preallocate sampling grid
    omega = false(Ny, Nz);

    % Loop through kz locations and generate sampling locations along ky
    for z = indices_z
        nacs_indices_samp_y = zeros(2,Ny_nacs/2);
        for side = 1:2 % Sample the same amount on each half of k-space
            half_line = sort(datasample(nacs_indices_y(side,:), Ny_nacs/2,...
                                        'Replace', false,...
                                        'Weights',w_y(side,:)));

            % Limit spacing between consecutive ky lines
            [gap,maxdex] = max(diff(half_line)); % find biggest gap
            while gap > max_ky_step
                % Find ky location with nearest neighbor(s)
                [~,mindex] = min(conv(diff(half_line),[1 1], 'valid'));

                % Move sampling location to the midpoint of largest gap
                half_line(mindex + 1) = half_line(maxdex) + round(gap/2);

                % Resort and update
                half_line = sort(half_line);
                [gap,maxdex] = max(diff(half_line));
            end
            nacs_indices_samp_y(side,:) = half_line;
        end
        indices_y = sort([acs_indices_y, reshape(nacs_indices_samp_y, 1, Ny_nacs)]);
        omega(indices_y,z) = true;
    end
end