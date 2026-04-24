function buildCMCTransformMatrix(dataFile)
%BUILDCMCTRANSFORMMATRIX Build the transformation matrices (T) for the CMC joint
% using translation and angle data from the provided Excel file.
% 
% Inputs:
%   dataFile : Excel file containing translation data (columns 3-6 for angles, 9-11 for translations)
% 
% Outputs:
%   Saves the transformation matrices (T_all) as a .mat file ("CMC_T_matrices.mat")

    % Read all data from the Excel file
    raw = readmatrix(dataFile);

    % Columns mapping (Flexion/Extension, Abduction/Adduction, Rotation)
    flexDeg = raw(:,3);   % Flexion/Extension (Column 3)
    abdDeg  = raw(:,4);   % Abduction/Adduction (Column 4)
    rotDeg  = raw(:,5);   % Rotation (Column 5)

    % Translation data (Radial-Ulnar, Dorsal-Volar, Proximal-Distal)
    x = raw(:,9);   % Radial-Ulnar (Column 9)
    y = raw(:,10);  % Dorsal-Volar (Column 10)
    z = raw(:,11);  % Proximal-Distal (Column 11)

    % Remove rows with missing or non-numeric data
    validRows = isfinite(flexDeg) & isfinite(abdDeg) & isfinite(rotDeg) & ...
                isfinite(x) & isfinite(y) & isfinite(z);

    flexDeg = flexDeg(validRows);
    abdDeg  = abdDeg(validRows);
    rotDeg  = rotDeg(validRows);

    x = x(validRows);
    y = y(validRows);
    z = z(validRows);

    N = length(flexDeg);

    % Preallocate arrays for rotation matrices (R), translation vectors (p), and homogeneous matrices (T)
    R_all = zeros(3,3,N);
    p_all = zeros(3,N);
    T_all = zeros(4,4,N);

    % Build the rotation matrices (R), translation vectors (p), and homogeneous transformation matrices (T)
    for k = 1:N
        % Convert degrees to radians for each angle
        flex = deg2rad(flexDeg(k));
        abd  = deg2rad(abdDeg(k));
        rot  = deg2rad(rotDeg(k));

        % Rotation about the x-axis (flexion/extension)
        Rx = [1 0 0;
              0 cos(flex) -sin(flex);
              0 sin(flex)  cos(flex)];

        % Rotation about the y-axis (abduction/adduction)
        Ry = [ cos(abd) 0 sin(abd);
               0        1 0;
              -sin(abd) 0 cos(abd)];

        % Rotation about the z-axis (rotation)
        Rz = [cos(rot) -sin(rot) 0;
              sin(rot)  cos(rot) 0;
              0         0        1];

        % Combine the rotations: Rx, Ry, Rz applied in the order Rz * Ry * Rx
        R = Rz * Ry * Rx;

        % Translation vector
        p = [x(k); y(k); z(k)];

        % Store the rotation matrix, translation vector, and homogeneous transformation matrix
        R_all(:,:,k) = R;
        p_all(:,k) = p;
        T_all(:,:,k) = [R p;
                        0 0 0 1];
    end

    %% Quick checks
    disp('First rotation matrix:');
    disp(R_all(:,:,1));

    disp('Orthogonality check:');
    disp(R_all(:,:,1)' * R_all(:,:,1));

    disp('Determinant of first rotation matrix:');
    disp(det(R_all(:,:,1)));

    disp('First translation vector:');
    disp(p_all(:,1));

    disp('First homogeneous transform T_all(:,:,1):');
    disp(T_all(:,:,1));

    %% Save the transformation matrices for later use
    save("CMC_T_matrices.mat", "R_all", "p_all", "T_all");

    disp("Finished building transformation matrices and saved to CMC_T_matrices.mat");
end