clear; clc;

%% =========================
% User settings
% ==========================
matFile = "CMC_T_matrices.mat";
outputExcel = "FHA_window_variability_results.xlsx";

% Window sizes in FRAMES
windowSizes = [5 10 15 20 30 40 50];

% Overlap percentages
overlapPercents = [0 20 40 60 80 99];

% Frame range to test
startFrame = 120;
endFrame   = 600;

% screw.m intersect option
intersect = 3;

% Ignore first N FHA windows
ignoreFirstN = 1;

% Minimum rotation magnitude threshold (deg)
% Windows with phi below this are ignored
minPhi = 0;   % try 1 or 2 later if you want stricter filtering

%% =========================
% Load data
% ==========================
S = load(matFile);
if ~isfield(S, 'T_all')
    error("The file does not contain T_all.");
end

if exist('screw','file') ~= 2
    error("screw.m not found on MATLAB path.");
end

T_full = S.T_all;
N_full = size(T_full,3);

startFrame = max(1, round(startFrame));
endFrame   = min(N_full, round(endFrame));

if startFrame >= endFrame
    error("startFrame must be smaller than endFrame.");
end

T = T_full(:,:,startFrame:endFrame);
N = size(T,3);

%% =========================
% Preallocate results
% ==========================
results = [];
rowCount = 0;

%% =========================
% Main parameter sweep
% ==========================
for w = 1:length(windowSizes)
    windowLength = windowSizes(w);
    dt = windowLength - 1;

    if dt >= N
        warning("Window size %d is too large for selected frame range. Skipping.", windowLength);
        continue;
    end

    for o = 1:length(overlapPercents)
        overlapPercent = overlapPercents(o);

        if overlapPercent < 0 || overlapPercent >= 100
            warning("Invalid overlap %.1f. Skipping.", overlapPercent);
            continue;
        end

        % Step between window starts
        windowStep = max(1, round(windowLength * (1 - overlapPercent/100)));

        % Window start indices inside selected range
        windowStartLocal = 1:windowStep:(N - dt);
        numAxesAttempted = numel(windowStartLocal);

        % Preallocate for this setting
        n_global_all = zeros(3, numAxesAttempted);
        phi_all      = zeros(1, numAxesAttempted);
        valid        = false(1, numAxesAttempted);

        % Compute FHA for each window
        for i = 1:numAxesAttempted
            k1 = windowStartLocal(i);
            k2 = k1 + dt;

            T1 = T(:,:,k1);
            T2 = T(:,:,k2);

            Trel = inv(T1) * T2;

            [n_local, ~, phi, ~] = screw(Trel, intersect);

            if all(isfinite(n_local)) && norm(n_local) > 0 && isfinite(phi)
                R1 = T1(1:3,1:3);

                % Convert axis direction to global frame
                n_global = R1 * n_local;
                n_global = n_global / norm(n_global);

                % Keep sign consistent for neighboring vectors
                if i > 1 && valid(i-1)
                    if dot(n_global, n_global_all(:,i-1)) < 0
                        n_global = -n_global;
                    end
                end

                n_global_all(:,i) = n_global;
                phi_all(i) = phi;
                valid(i) = true;
            end
        end

        % Ignore first N windows
        usedMask = valid;
        if ignoreFirstN > 0
            usedMask(1:min(ignoreFirstN, numAxesAttempted)) = false;
        end

        % Also ignore windows with very small rotation
        if minPhi > 0
            usedMask = usedMask & (phi_all >= minPhi);
        end

        usedIdx = find(usedMask);
        numValid = sum(valid);
        numUsed  = numel(usedIdx);

        % Need at least 2 FHA axes to compute consecutive variability
        if numUsed >= 2
            angleDiffDeg = zeros(1, numUsed - 1);

            for j = 1:(numUsed - 1)
                n1 = n_global_all(:, usedIdx(j));
                n2 = n_global_all(:, usedIdx(j+1));

                % Angle between axes, treating n and -n as same axis
                d = abs(dot(n1, n2));
                d = min(1, max(-1, d));   % clamp for numerical safety
                angleDiffDeg(j) = rad2deg(acos(d));
            end

            meanAngleDiff = mean(angleDiffDeg);
            stdAngleDiff  = std(angleDiffDeg);
            medianAngleDiff = median(angleDiffDeg);

            meanPhi = mean(phi_all(usedIdx));
            stdPhi  = std(phi_all(usedIdx));

            meanDir = mean(n_global_all(:,usedIdx), 2);
            meanDir = meanDir / norm(meanDir);
            meanNx = meanDir(1);
            meanNy = meanDir(2);
            meanNz = meanDir(3);
        else
            angleDiffDeg = [];
            meanAngleDiff = NaN;
            stdAngleDiff  = NaN;
            medianAngleDiff = NaN;
            meanPhi = NaN;
            stdPhi  = NaN;
            meanNx = NaN;
            meanNy = NaN;
            meanNz = NaN;
        end

        % Save one summary row
        rowCount = rowCount + 1;
        results(rowCount).startFrame = startFrame;
        results(rowCount).endFrame = endFrame;
        results(rowCount).windowLength = windowLength;
        results(rowCount).dt = dt;
        results(rowCount).overlapPercent = overlapPercent;
        results(rowCount).windowStep = windowStep;
        results(rowCount).attemptedFHA = numAxesAttempted;
        results(rowCount).validFHA = numValid;
        results(rowCount).usedFHA = numUsed;
        results(rowCount).ignoreFirstN = ignoreFirstN;
        results(rowCount).minPhi = minPhi;
        results(rowCount).meanAngleDiffDeg = meanAngleDiff;
        results(rowCount).medianAngleDiffDeg = medianAngleDiff;
        results(rowCount).stdAngleDiffDeg = stdAngleDiff;
        results(rowCount).meanPhiDeg = meanPhi;
        results(rowCount).stdPhiDeg = stdPhi;
        results(rowCount).meanNx = meanNx;
        results(rowCount).meanNy = meanNy;
        results(rowCount).meanNz = meanNz;
    end
end

%% =========================
% Convert to table and save
% ==========================
resultsTable = struct2table(results);

disp(resultsTable);

writetable(resultsTable, outputExcel, 'Sheet', 'Summary');

%% =========================
% Optional plots
% ==========================
figure; hold on; grid on;
colors = lines(length(overlapPercents));

for o = 1:length(overlapPercents)
    overlap = overlapPercents(o);
    idx = resultsTable.overlapPercent == overlap;

    plot(resultsTable.windowLength(idx), ...
         resultsTable.meanAngleDiffDeg(idx), ...
         '-o', 'LineWidth', 2, 'DisplayName', sprintf('%d%% overlap', overlap), ...
         'Color', colors(o,:));
end

xlabel('Window size (frames)');
ylabel('Mean consecutive FHA direction difference (deg)');
title('FHA Direction Variability vs Window Size');
legend('Location','best');

figure; hold on; grid on;
for o = 1:length(overlapPercents)
    overlap = overlapPercents(o);
    idx = resultsTable.overlapPercent == overlap;

    plot(resultsTable.windowLength(idx), ...
         resultsTable.usedFHA(idx), ...
         '-o', 'LineWidth', 2, 'DisplayName', sprintf('%d%% overlap', overlap), ...
         'Color', colors(o,:));
end

xlabel('Window size (frames)');
ylabel('Valid FHA windows used');
title('FHA Count vs Window Size');
legend('Location','best');

disp("Done. Results written to:");
disp(outputExcel);