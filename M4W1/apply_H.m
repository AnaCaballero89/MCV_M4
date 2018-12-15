function output=apply_H (I, H, size_output)

[rows, cols, channels] = size(I);

H=[1 0 1; 0 1 0; 1 0 1];


if(size_output ~= size(I))
    
% In the case that the output has different size respect to original image we recompute the corners of the new image

% We transform the corners of the image in homogeneous coordinates:
       
C1 = [1 1 1]';
C2 = [cols 1 1]';
C3 = [1 rows 1]';
C4 = [cols rows 1]';

% Compute the transformation of the homogeneous corners using the homography

HC1 = H*C1;
HC2 = H*C2;
HC3 = H*C3;
HC4 = H*C4;

% We normalise the transformed corners dividing by the last term (last term will be equal to 1)

HC1 = HC1/HC1(3);
HC2 = HC2/HC2(3);
HC3 = HC3/HC3(3);
HC4 = HC4/HC4(3);


x_min = round(min([HC1(1) HC2(1) HC3(1) HC4(1)]));
y_min = round(min([HC1(2) HC2(2) HC3(2) HC4(2)]));
x_max = round(max([HC1(1) HC2(1) HC3(1) HC4(1)]));
y_max = round(max([HC1(2) HC2(2) HC3(2) HC4(2)]));

else
    
% In the case that the output is equal than original image we use the original image size

x_min = 1;
y_min = 1;
x_max = cols;
y_max = rows;

end

% falta acabar, aplicar la homografia a la imatge i tornar la imatge transformada
end
