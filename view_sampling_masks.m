function view_sampling_masks(masks)
    % masks: 3D logical array (rows x cols x time)
    [rows, cols, T] = size(masks);

    % Create figure
    hFig = figure('Name', 'Sampling Mask Viewer', 'KeyPressFcn', @key_press);

    % Create subplots
    hAx1 = subplot(1, 2, 1); % current frame
    hImg1 = imagesc(hAx1, masks(:, :, 1));
    colormap(hAx1, gray);
    axis(hAx1, 'image');
    title(hAx1, 'Current Frame: 1');
    xlabel(hAx1, 'Slow PE direction');
    ylabel(hAx1, 'Fast PE direction');

    hAx2 = subplot(1, 2, 2); % cumulative count
    cumulativeCounts = sum(masks(:, :, 1), 3);
    hImg2 = imagesc(hAx2, cumulativeCounts);
    axis(hAx2, 'image');
    title(hAx2, 'Cumulative Sampling Counts');
    xlabel(hAx2, 'Slow PE direction');
    ylabel(hAx2, 'Fast PE direction');
    colorbar(hAx2);

    % Create colormap: blue for zero, hot for rest
    nColors = 256;
    baseMap = hot(nColors - 1); % or parula, jet, etc.
    customMap = [0 0 1; baseMap]; % blue at bottom
    colormap(hAx2, customMap);

    % Add slider
    hSlider = uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', T, 'Value', 1, ...
        'SliderStep', [1/(T-1), 1/(T-1)], ...
        'Position', [150 20 300 20], ...
        'Callback', @(src, event) update_frame(round(get(src, 'Value'))));

    % Frame update function
    function update_frame(frame)
        frame = max(1, min(T, frame)); % Clamp to valid range

        % Update current frame image
        set(hImg1, 'CData', masks(:, :, frame));
        title(hAx1, sprintf('Current Frame: %d', frame));

        % Update cumulative count image
        cumulativeCounts = sum(masks(:, :, 1:frame), 3);
        set(hImg2, 'CData', cumulativeCounts);
        title(hAx2, sprintf('Cumulative Sampling Counts (1:%d)', frame));
        caxis(hAx2, [0, max(cumulativeCounts(:)) + eps]);
    end

    % Key press handler
    function key_press(~, event)
        currentFrame = round(get(hSlider, 'Value'));
        switch event.Key
            case 'rightarrow'
                currentFrame = min(currentFrame + 1, T);
            case 'leftarrow'
                currentFrame = max(currentFrame - 1, 1);
        end
        set(hSlider, 'Value', currentFrame);
        update_frame(currentFrame);
    end

    % Initialize
    update_frame(1);
end