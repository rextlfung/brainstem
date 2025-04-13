function interactive4D(data)
% INTERACTIVE4DVOLUME_GALLERY_FAST Fast visualization of 4D fMRI data
%
% INPUT:
%   data - 4D tensor (X x Y x Z x Time)
%
% USAGE:
%   interactive4DVolume_gallery_fast(data)

    if ndims(data) ~= 4
        error('Input must be a 4D tensor (X x Y x Z x Time)');
    end

    % Setup figure
    fig = figure('Name', 'Fast 4D Slice Gallery', ...
                 'NumberTitle', 'off', ...
                 'WindowScrollWheelFcn', @scrollCallback, ...
                 'Color', 'k');
             
    colormap('gray');
    handles.data = data;
    handles.currentTime = 1;
    handles.maxTime = size(data, 4);
    guidata(fig, handles);

    % Prepare one axis only
    ax = axes('Parent', fig);
    ax.Position = [0 0 1 1];
    ax.XTick = [];
    ax.YTick = [];
    axis image off;

    handles.ax = ax;
    guidata(fig, handles);

    plotGallery(fig);

    function plotGallery(fig)
        handles = guidata(fig);
        volumeData = handles.data(:, :, :, handles.currentTime);
        numSlices = size(volumeData, 3);

        % Determine grid size
        nCols = ceil(sqrt(numSlices));
        nRows = ceil(numSlices / nCols);

        % Normalize slices
        volumeData = mat2gray(volumeData); % scale to [0,1]

        % Create a big empty canvas
        [sx, sy, ~] = size(volumeData);
        canvas = ones(nRows * sx, nCols * sy);

        % Paste each slice into the canvas
        for idx = 1:numSlices
            row = floor((idx-1) / nCols);
            col = mod((idx-1), nCols);
            xIdx = (row*sx + 1):(row*sx + sx);
            yIdx = (col*sy + 1):(col*sy + sy);
            canvas(xIdx, yIdx) = volumeData(:, :, idx);
        end

        % Show canvas
        imagesc(handles.ax, canvas);
        axis(handles.ax, 'image', 'off');
        title(handles.ax, sprintf('Time %d / %d', handles.currentTime, handles.maxTime), ...
            'Color', 'w', 'FontSize', 14);

        drawnow;
    end

    function scrollCallback(src, event)
        handles = guidata(src);
        handles.currentTime = handles.currentTime + sign(event.VerticalScrollCount);
        handles.currentTime = max(1, min(handles.maxTime, handles.currentTime)); % Clamp
        guidata(src, handles);
        plotGallery(src);
    end
end
