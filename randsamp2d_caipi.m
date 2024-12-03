%% randsamp2d_caipi
%
% Generate a 2D pseudo-random sampling pattern with higher
% probability towards the center, drawn from a Gaussian distribution.
%
% Additionally allows for CAIPI-like shifting to sample multiple kz
% locations per excitation.
% delay = 1;
% Input arguments:
% N = [Ny Nz] = sampling grid dimensions. Positive vector of integers >= 1.
% R = [Ry Rz] = acceleration factors. Positive vector of integers >= 1.
% max_ky_step = the largest tolerable gap between subsequently sampled fast
% PE (ky) locations. Positive integer >= 1.
% caipi_z = Number of kz locations to be sampled per shot. Used to
% constrain the minimum gap between adjacent sampled kz locations.
% Positive integer >= 1.
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
% 3. For kz, sample from the following subset of kz locations:
% 1:caipi_z:(Nz - caipi_z + 1), then apply CAIPI shift
%
% Last modified Nov 18th, 2024. Rex Fung

function omega = randsamp2d_caipi(N, R, max_ky_step, caipi_z)
    % Unpack and validate input arguments
    Ny = N(1); Nz = N(2);
    Ry = R(1); Rz = R(2);
    assert(Ny >= 1 && Nz >= 1, 'Dimensions must be >= 1');
    assert(Ry >= 1 && Rz >= 1, 'Acceleration factors must be >= 1');
    
    % Generate sampling indices to draw from
    indices_y = [1:Ny/2; (Ny/2 + 1):Ny]; % Split in to two halves
    indices_z = 1:caipi_z:(Nz - caipi_z + 1);

    % Compute the number of locations to sample based on R
    Ny_sample = ceil(Ny/Ry/2)*2; % Ensure even number in fast PE direction
    Nz_sample = ceil(length(indices_z)/Rz);

    % Genearte sampling weights using Gaussian pdf.
    % Mean and std are arbitrarility chosen for now.
    w_y = normpdf(indices_y, Ny/2 + 1, Ny/6);
    w_y = w_y./sum(w_y,'all');
    w_z = normpdf(indices_z, Nz/2 + 1, Nz/6);
    w_z = w_z./sum(w_z,'all');

    % Select sampling locations along kz
    indices_samp_z = sort(datasample(indices_z, Nz_sample,...
                                     'Replace', false,...
                                     'Weights', w_z));

    % Preallocate sampling grid
    omega = false(Ny, Nz);

    % Loop through kz locations and generate sampling locations along ky
    for z = indices_samp_z

        % Sample the same number of locations on each side of ky
        indices_samp_y = zeros(2, Ny_sample/2);
        for side = 1:2
            half_line = sort(datasample(indices_y(side,:), Ny_sample/2,...
                                        'Replace', false,...
                                        'Weights', w_y(side,:)));

            % Constrain maximum jump in ky
            [gap, maxdex] = max(diff(half_line)); % find biggest gap
            while gap > max_ky_step
                % Find ky location with nearest neighbor(s)
                [~, mindex] = min(conv(diff(half_line), [1 1], 'valid'));

                % Move sampling location to the midpoint of largest gap
                half_line(mindex + 1) = half_line(maxdex) + round(gap/2);

                % Resort and update
                half_line = sort(half_line);
                [gap, maxdex] = max(diff(half_line)); 
            end

            indices_samp_y(side,:) = half_line;
        end

        % Combine into one vector of sampling locations
        indices_samp_y = sort(reshape(indices_samp_y, 1, Ny_sample));

        % CAIPI-shifted sampling
        shifts = randperm(caipi_z);
        for i_shift = 1:length(shifts)
            shift = shifts(i_shift);
            omega(indices_samp_y(i_shift:caipi_z:length(indices_samp_y)), z + shift - 1) = true;
        end
    end
end