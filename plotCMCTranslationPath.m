function plotCMCTranslationPath(inputFile, startFrame, endFrame, animationSpeed)
%PLOTCMCTRANSLATIONPATH Plot CMC motion path from translation data in an Excel file.
%   
% Inputs:
%   inputFile     : Excel file containing translation data with columns
%                   D = radial-ulnar, E = dorsal-volar, F = proximal-distal
%   startFrame    : the starting frame for the plot
%   endFrame      : the ending frame for the plot
%   animationSpeed: playback speed for the animation (in seconds, default: 0.05)
%
% Outputs:
%   3D translation path plot for the CMC motion with or without animation.

    % Default animation speed if not provided
    if nargin < 4
        animationSpeed = 0.05;  % Default speed (0.05 seconds per frame)
    end
    
    % Check if the input file exists
    if nargin < 1
        error('Please provide an Excel file as input (e.g., plotCMCTranslationPath(''yourFile.xlsx''))');
    end
    
    % Read the data from Excel file (skip first row with headers)
    raw = readmatrix(inputFile, 'Range', 'A2');  % Start from the second row to skip headers

    % Columns mapping (radial-ulnar, dorsal-volar, proximal-distal)
    x = raw(:, 9);  % Radial-Ulnar (Column 9)
    y = raw(:, 10);  % Dorsal-Volar (Column 10)
    z = raw(:, 11);  % Proximal-Distal (Column 11)

    % Remove rows with missing or non-numeric data
    validRows = all(isfinite([x, y, z]), 2);  % Ensure no NaNs or Inf
    x = x(validRows);
    y = y(validRows);
    z = z(validRows);

    % Get the total number of frames (rows in data)
    totalFrames = length(x);
    
    % Validate startFrame
    if nargin < 2 || startFrame <= 0
        startFrame = 1;  % Default to frame 1 if not provided or invalid
    elseif startFrame > totalFrames
        error('startFrame exceeds the number of available frames.');
    end
    
    % Validate endFrame
    if nargin < 3 || endFrame > totalFrames
        endFrame = totalFrames;  % Default to the last frame if not provided or exceeds data length
    elseif endFrame < startFrame
        error('endFrame cannot be less than startFrame.');
    end

    % Subset the data based on start and end frames
    x = x(startFrame:endFrame);
    y = y(startFrame:endFrame);
    z = z(startFrame:endFrame);

    %% Check if animation should be shown (animationSpeed > 0)
    if animationSpeed > 0
        % Create a new figure for the animation
        figure;
        grid on;
        axis equal;
        view(3);  % Set the view to 3D
        xlabel('Radial-Ulnar (mm)');
        ylabel('Dorsal-Volar (mm)');
        zlabel('Proximal-Distal (mm)');
        title('CMC Translation Path Animation');

        % Initialize scatter for the start point (red dot)
        hStart = scatter3(x(1), y(1), z(1), 50, 'filled', 'r');  % Red start dot
        hold on;
        hPath = plot3(x, y, z, 'k-', 'LineWidth', 1.5);  % Path (black line)

        % Animation loop: update the position of the red start dot along the path
        for k = 1:length(x)
            % Update the position of the red start dot
            set(hStart, 'XData', x(k), 'YData', y(k), 'ZData', z(k));

            % Pause to create an animation effect (based on user-specified speed)
            pause(animationSpeed);  % Adjust the speed of the animation (in seconds)
            drawnow;  % Force the update of the plot
        end
    end

    %% After animation or if animation is disabled, show the static path plot
    figure;
    plot3(x, y, z, 'LineWidth', 1.5);
    grid on;
    axis equal;
    view(3);  % Set the view to 3D

    xlabel('Radial-Ulnar (mm)');
    ylabel('Dorsal-Volar (mm)');
    zlabel('Proximal-Distal (mm)');
    title(sprintf('CMC Translation Path (Frames %d to %d)', startFrame, endFrame));

    %% Optional: mark the start and end of the motion
    hold on;
    scatter3(x(1), y(1), z(1), 50, 'filled', 'r');  % Mark start
    scatter3(x(end), y(end), z(end), 50, 'filled', 'g');  % Mark end
    legend('Path', 'Start', 'End', 'Location', 'best');
end