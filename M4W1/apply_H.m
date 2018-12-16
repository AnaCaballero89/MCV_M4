function Ioutput=apply_H (I, H) 
    % Applies given homography H to a given image I
    % Returns transformed image Ioutput
    % ___________________________________
    
    % Recompute the corners of the new image to fit transformed image
    % We transform the corners of the image in homogeneous coordinates:
    [rows, cols, channels] = size(I);
    C1 = [1 1 1]';
    C2 = [cols 1 1]';
    C3 = [1 rows 1]';
    C4 = [cols rows 1]';

    % Compute the transformation of the homogeneous corners using the homography
    HC1 = H*C1;
    HC2 = H*C2;
    HC3 = H*C3;
    HC4 = H*C4;

    % We normalise the transformed corners dividing by the last term (last
    % term will be equal to 1)
    HC1 = HC1/HC1(3);
    HC2 = HC2/HC2(3);
    HC3 = HC3/HC3(3);
    HC4 = HC4/HC4(3);
    
    % Pick bounding values and put Coordinates as integers
    x_min = round(min([HC1(1) HC2(1) HC3(1) HC4(1)]));
    y_min = round(min([HC1(2) HC2(2) HC3(2) HC4(2)]));
    x_max = round(max([HC1(1) HC2(1) HC3(1) HC4(1)]));
    y_max = round(max([HC1(2) HC2(2) HC3(2) HC4(2)]));

    % Generate meshgrid to fit resulting image
    [X,Y] = meshgrid(x_min:x_max, y_min:y_max);
    
    % Projected coordinates [x1, x2, 1]
    [Hrow,Hcol] =size(X);
    XYZ = [X(:) Y(:) ones(Hrow*Hcol,1)]';

    % Transform image
    H_inv = inv(H); % Compute inverse transformation
    Im_world= H_inv*XYZ; % Apply inverse transformation in projected space
    Im_world = Im_world/Im_world(3); % Transform to world coordinates
    X_im = reshape(Im_world(1,:), Hrow, Hcol);
    Y_im = reshape(Im_world(2,:), Hrow, Hcol);
    
    % Interpolate color values in the coordinates
    Ioutput = zeros(Hrow, Hcol, channels);
    for c=1:channels
        Ioutput(:,:,c) = interp2(double(I(:,:,c)), X_im, Y_im, '*bilinear');
    end
end
