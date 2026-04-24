function out = plotOverlappedFHA_CMC(matFile, axisLength, stepDraw, intersect, ...
    startFrame, endFrame, dt, overlapPercent, ignoreFirstN)
%PLOTOVERLAPPEDFHA_CMC
% Compute and plot FHA directions using overlapping frame windows.
%
% Each FHA is computed from:
%   T_start  -->  T_start+dt
%
% so each window contains:
%   dt + 1 frames
%
% Example:
%   dt = 9  -> uses 10 frames per FHA
%
% Overlap examples:
%   overlapPercent = 0
%       windows: 1-10, 11-20, 21-30, ...
%
%   overlapPercent = 50
%       windows: 1-10, 6-15, 11-20, ...
%
% Inputs:
%   matFile         : .mat file containing T_all
%   axisLength      : half-length of each plotted FHA line
%   stepDraw        : plot every Nth FHA line for clarity
%   intersect       : passed into screw.m
%   startFrame      : first frame to include
%   endFrame        : last frame to include
%   dt              : frame gap between first and last pose in each FHA window
%   overlapPercent  : overlap percentage between neighboring windows
%   ignoreFirstN    : number of initial FHA windows to ignore in plotting/averaging
%
% Outputs:
%   out.n_local
%   out.n_global
%   out.anchorPts
%   out.phi
%   out.t
%   out.valid
%   out.windowStart
%   out.windowEnd
%   out.windowStep
%   out.windowLength
%   out.ignoreFirstN
%   out.validUsedForAverage

    if nargin < 1, matFile = "CMC_T_matrices.mat"; end
    if nargin < 2, axisLength = 20; end
    if nargin < 3, stepDraw = 1; end
    if nargin < 4, intersect = 3; end
    if nargin < 5, startFrame = 1; end

    S = load(matFile);
    if ~isfield(S, 'T_all')
        error("File does not contain T_all.");
    end

    if exist('screw', 'file') ~= 2
        error("screw.m not found on MATLAB path.");
    end

    T_full = S.T_all;
    N_full = size(T_full, 3);

    if nargin < 6 || endFrame > N_full
        endFrame = N_full;
    end

    if nargin < 7
        dt = 9;  % default = 10-frame window
    end

    if nargin < 8
        overlapPercent = 50;
    end

    if nargin < 9
        ignoreFirstN = 1;   % default: ignore the first FHA
    end

    % Validate inputs
    startFrame   = max(1, round(startFrame));
    endFrame     = min(N_full, round(endFrame));
    dt           = round(dt);
    ignoreFirstN = max(0, round(ignoreFirstN));

    if startFrame >= endFrame
        error("startFrame must be smaller than endFrame.");
    end

    if dt < 1
        error("dt must be at least 1.");
    end

    if overlapPercent < 0 || overlapPercent >= 100
        error("overlapPercent must be in the range [0, 100).");
    end

    % Restrict to selected frame range
    T = T_full(:,:,startFrame:endFrame);
    N = size(T,3);
    frameIdx = startFrame:endFrame;

    % Window length in frames
    windowLength = dt + 1;

    if windowLength > N
        error("Selected frame range is too short for the chosen dt.");
    end

    % Step size between neighboring windows
    % Example: windowLength = 10, overlap = 50% --> step = 5
    windowStep = max(1, round(windowLength * (1 - overlapPercent/100)));

    % Starting indices of each FHA window inside the selected range
    windowStart = 1:windowStep:(N - dt);
    numAxes = numel(windowStart);

    % Path of the segment origin
    originGlobal = zeros(3,N);
    for k = 1:N
        originGlobal(:,k) = T(1:3,4,k);
    end

    % Preallocate outputs
    n_local_all  = zeros(3, numAxes);
    n_global_all = zeros(3, numAxes);
    phi_all      = zeros(1, numAxes);
    t_all        = zeros(1, numAxes);
    anchorPts    = zeros(3, numAxes);
    valid        = false(1, numAxes);
    windowEnd    = zeros(1, numAxes);

    % Compute FHA for each overlapping window
    for i = 1:numAxes
        k1 = windowStart(i);
        k2 = k1 + dt;
        windowEnd(i) = k2;

        T1 = T(:,:,k1);
        T2 = T(:,:,k2);

        % Relative motion from first pose to last pose in the window
        Trel = inv(T1) * T2;

        [n_local, ~, phi, t] = screw(Trel, intersect);

        if all(isfinite(n_local)) && norm(n_local) > 0 && isfinite(phi) && isfinite(t)
            R1 = T1(1:3,1:3);

            % Convert axis direction from local to global frame
            n_global = R1 * n_local;
            n_global = n_global / norm(n_global);

            % Keep direction sign visually consistent
            if i > 1 && valid(i-1)
                if dot(n_global, n_global_all(:,i-1)) < 0
                    n_global = -n_global;
                end
            end

            n_local_all(:,i)  = n_local / norm(n_local);
            n_global_all(:,i) = n_global;
            phi_all(i) = phi;
            t_all(i)   = t;

            % Anchor FHA line at the start-frame origin for plotting
            anchorPts(:,i) = originGlobal(:,k1);
            valid(i) = true;
        end
    end

    % Decide which FHA windows are used after ignoring the first few
    usedMask = valid;
    if ignoreFirstN > 0
        usedMask(1:min(ignoreFirstN, numAxes)) = false;
    end

    % Store outputs
    out.n_local            = n_local_all;
    out.n_global           = n_global_all;
    out.anchorPts          = anchorPts;
    out.phi                = phi_all;
    out.t                  = t_all;
    out.valid              = valid;
    out.windowStart        = frameIdx(windowStart);
    out.windowEnd          = frameIdx(windowEnd);
    out.windowStep         = windowStep;
    out.windowLength       = windowLength;
    out.ignoreFirstN       = ignoreFirstN;
    out.validUsedForAverage = usedMask;

    %% Plot 1: path + overlapped FHA directions
    figure; hold on; grid on; axis equal; view(3);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('3D Visualization of Finite Helical Axes During CMC Joint Motion');
    %{
title(sprintf(['CMC Motion with Overlapped FHA Directions\n' ...
        'Frames %d to %d | dt = %d | Window = %d frames | Overlap = %.1f%%'], ...
        startFrame, endFrame, dt, windowLength, overlapPercent));
    %}
    % Plot the motion path
    plot3(originGlobal(1,:), originGlobal(2,:), originGlobal(3,:), ...
        'b-', 'LineWidth', 2);

    % Plot FHA direction lines (skip the first ignored windows)
    firstPlotIdx = max(1, ignoreFirstN + 1);
    for i = firstPlotIdx:stepDraw:numAxes
        if ~usedMask(i)
            continue;
        end

        anchor = anchorPts(:,i);
        n = n_global_all(:,i);

        p1 = anchor - axisLength * n;
        p2 = anchor + axisLength * n;

        plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], ...
            'r-', 'LineWidth', 2);
    end

    legend('CMC path', 'FHA', 'Location', 'best');

    %% Plot 2: rotation magnitude
    figure;
    plot(out.windowStart, phi_all, 'LineWidth', 2); hold on;
    
    grid on;
    xlabel('Window start frame');
    ylabel('\phi (deg)');
    title('FHA Rotation Magnitude per Window');

    %% Console summary
    totalAttempted = numAxes;
    totalValid = sum(valid);
    totalUsed  = sum(usedMask);

    fprintf('Selected frame range: %d to %d\n', startFrame, endFrame);
    fprintf('dt: %d\n', dt);
    fprintf('Window length: %d frames\n', windowLength);
    fprintf('Overlap: %.1f%%\n', overlapPercent);
    fprintf('Window step: %d frames\n', windowStep);
    fprintf('Total FHA windows attempted: %d\n', totalAttempted);
    fprintf('Valid FHA windows: %d\n', totalValid);
    fprintf('Ignoring first %d FHA windows\n', ignoreFirstN);
    fprintf('Valid FHA windows used in average: %d\n', totalUsed);

    if any(usedMask)
        meanDir = mean(n_global_all(:,usedMask), 2);
        meanDir = meanDir / norm(meanDir);
        fprintf('Mean global FHA direction = [%.4f %.4f %.4f]\n', ...
            meanDir(1), meanDir(2), meanDir(3));
    else
        fprintf('No valid FHA directions remained after ignoring initial windows.\n');
    end
end